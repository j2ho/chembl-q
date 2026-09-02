# chembl_curator/splitter.py

"""Stage 7: Train/test split by sequence-identity clustering.

Clusters all passed targets (and optional external datasets) using MMseqs2
at a configurable seqid threshold (default 30%). Clusters are greedily
assigned to train/test sets while balancing the ratio of each data source
across both splits.

Per-entry sampling weight = 1 / log2(cluster_size + 1)

External FASTA entries must use dot-prefixed IDs:
    >{source}.{entry_id}

For example:
    >pdbbind.1a4k
    >biolip.10gs_VWW_A_1

Output format (tab-separated, with header):
    source  entry_id          compound        weight
    chembl  P12345            CHEMBL405346    0.17
    biolip  10gs_VWW_A_1      -               0.19
    pdbbind 1a4k              -               0.26

Also generates chembl_targets.tsv:
    uniprot  split  n_actives  n_decoys
"""

import logging
import math
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class TargetSplitter:
    """Split targets into train/test sets based on sequence identity clustering."""

    def __init__(
        self,
        seqid: float = 0.4,
        valid_frac: float = 1.0,
        threads: int = 4,
        log_level: str = "INFO",
    ):
        self.seqid = seqid
        self.valid_frac = valid_frac
        self.threads = threads
        # Filled in by run(): ChEMBL targets with a direct external homologue.
        self.blocked_targets: Set[str] = set()
        # Filled in by run(): test targets demoted for train proximity.
        self.demoted_targets: Set[str] = set()
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # -- I/O helpers -----------------------------------------------------------

    def _load_fasta(self, path: Path) -> Dict[str, str]:
        """Return {sequence_id: sequence} from a FASTA file."""
        seqs: Dict[str, str] = {}
        current_id: Optional[str] = None
        current_seq: List[str] = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        seqs[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif current_id:
                    current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq)
        return seqs

    # -- MMseqs2 clustering ----------------------------------------------------

    def _run_mmseqs_cluster(
        self, fasta: Path, tmpdir: Path
    ) -> Dict[str, str]:
        """Run MMseqs2 easy-cluster and return {member: representative} mapping."""
        cluster_base = tmpdir / "cluster"
        tmp = tmpdir / "tmp"
        tmp.mkdir(exist_ok=True)

        subprocess.run(
            [
                "mmseqs", "easy-cluster",
                str(fasta), str(cluster_base), str(tmp),
                "--min-seq-id", str(self.seqid),
                "-s", "7.5",
                "--threads", str(self.threads),
            ],
            check=True, capture_output=True,
        )

        member_to_rep: Dict[str, str] = {}
        cluster_tsv = Path(str(cluster_base) + "_cluster.tsv")
        with open(cluster_tsv) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    rep, member = parts
                    member_to_rep[member] = rep
        return member_to_rep

    def _find_external_homologues(
        self, chembl_fasta: Path, external_fasta: Path, tmpdir: Path
    ) -> Set[str]:
        """ChEMBL targets with a direct sequence homologue in the external sets.

        Answers the leakage question per target rather than per cluster: does
        *this* sequence have an external hit at or above self.seqid over at
        least 80% of its length. Hits shorter than that are domain-level and
        do not mean a model has seen the target.
        """
        hits = tmpdir / "external_hits.tsv"
        search_tmp = tmpdir / "search_tmp"
        search_tmp.mkdir(exist_ok=True)
        subprocess.run(
            [
                "mmseqs", "easy-search",
                str(chembl_fasta), str(external_fasta), str(hits), str(search_tmp),
                "--threads", str(self.threads),
                "-s", "7.5",
                "--format-output", "query,target,fident,alnlen,qlen",
            ],
            check=True, capture_output=True,
        )

        blocked: Set[str] = set()
        n_rows = 0
        with open(hits) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                n_rows += 1
                try:
                    fident, alnlen, qlen = float(parts[2]), int(parts[3]), int(parts[4])
                except ValueError:
                    continue
                if fident < self.seqid or alnlen / max(1, qlen) < 0.8:
                    continue
                blocked.add(parts[0].split(".", 1)[-1])

        self.logger.info(
            f"External homologue search: {n_rows} hits, "
            f"{len(blocked)} ChEMBL targets blocked from test "
            f"(>={self.seqid} identity over >=80% of the query)"
        )
        return blocked

    def _demote_test_targets_close_to_train(
        self,
        test_targets: Set[str],
        train_seq_fasta: Path,
        chembl_seqs: Dict[str, str],
        tmpdir: Path,
    ) -> Set[str]:
        """Test targets that are still too close to something in train.

        Clustering alone does not guarantee separation. MMseqs easy-cluster is
        greedy set-cover, so two targets at 45% identity land in different
        clusters whenever each attaches to a different representative; if one
        of those clusters is forced to train by an external homologue, the
        other stays in test and the pair leaks. A direct search after the split
        found 5 such pairs (all paralogues at 0.41-0.45 identity).

        Runs to a fixed point: demoting a target grows train, which can pull in
        further test targets. Train only ever grows, so this terminates.
        """
        demoted: Set[str] = set()
        for round_no in range(1, 11):
            remaining = test_targets - demoted
            if not remaining:
                break
            query = tmpdir / f"test_round{round_no}.fasta"
            with open(query, 'w') as f:
                for t in sorted(remaining):
                    if t in chembl_seqs:
                        f.write(f">{t}\n{chembl_seqs[t]}\n")

            db = tmpdir / f"train_round{round_no}.fasta"
            with open(db, 'w') as f:
                f.write(train_seq_fasta.read_text())
                for t in sorted(demoted):
                    if t in chembl_seqs:
                        f.write(f">demoted.{t}\n{chembl_seqs[t]}\n")

            hits = tmpdir / f"demote_hits{round_no}.tsv"
            search_tmp = tmpdir / f"demote_tmp{round_no}"
            search_tmp.mkdir(exist_ok=True)
            subprocess.run(
                [
                    "mmseqs", "easy-search",
                    str(query), str(db), str(hits), str(search_tmp),
                    "--threads", str(self.threads), "-s", "7.5",
                    "--max-seqs", "20000", "-e", "10000",
                    "--format-output", "query,target,fident,alnlen,qlen",
                ],
                check=True, capture_output=True,
            )

            new: Set[str] = set()
            with open(hits) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    try:
                        fident, alnlen, qlen = (
                            float(parts[2]), int(parts[3]), int(parts[4]))
                    except ValueError:
                        continue
                    if fident >= self.seqid and alnlen / max(1, qlen) >= 0.8:
                        new.add(parts[0])
            if not new:
                break
            demoted |= new
            self.logger.info(
                f"Separation round {round_no}: demoted {len(new)} test targets "
                f"within {self.seqid} identity of train ({len(demoted)} total)"
            )
        return demoted

    # -- Generic source parsing ------------------------------------------------

    @staticmethod
    def _parse_source(sid: str) -> str:
        """Extract source label from a dot-prefixed sequence ID."""
        dot = sid.find('.')
        if dot > 0:
            return sid[:dot]
        return 'other'

    # -- Greedy split ----------------------------------------------------------

    def _greedy_split(
        self, cluster_groups: Dict[str, Set[str]]
    ) -> Tuple[Set[str], Set[str]]:
        """Assign clusters to train/test.

        A cluster may go to test only if every member is a ChEMBL target *and*
        no member has a direct external homologue (see blocked_targets). The
        point of the test set is to benchmark models trained on PDBbind and
        BioLiP without retraining them, so a test target with an external
        homologue has already been seen and cannot measure generalisation.
        Everything else is forced to train, whatever that does to source ratios.

        Cluster membership alone is the wrong test. MMseqs clustering is greedy
        set-cover, so members are similar to the representative and not
        necessarily to each other: blocking on cluster composition rejected 311
        of 1,317 targets that have no external sequence above the cutoff at all.
        blocked_targets comes from a direct search instead.

        Note this enforces *sequence*-level separation only. Pocket similarity
        against the external sets is not checked here, and the two barely
        overlap: of 3,073 pocket-similar ChEMBL pairs only 107 are also caught
        by sequence identity. Test targets can still share a pocket with a
        PDBbind or BioLiP entry.
        """
        blocked = self.blocked_targets or set()
        eligible = {
            rep for rep, members in cluster_groups.items()
            if all(self._parse_source(m) == "chembl" for m in members)
            and not any(m.split(".", 1)[-1] in blocked for m in members)
        }
        forced_train = set(cluster_groups) - eligible
        n_blocked_targets = sum(
            sum(1 for m in cluster_groups[rep] if self._parse_source(m) == "chembl")
            for rep in forced_train
        )
        self.logger.info(
            f"Test-eligible clusters (ChEMBL-only): {len(eligible)}; "
            f"{len(forced_train)} clusters forced to train, holding "
            f"{n_blocked_targets} ChEMBL targets that share a cluster with "
            "PDBbind or BioLiP"
        )
        if not eligible:
            raise RuntimeError(
                "No ChEMBL-only clusters: every target clusters with an external "
                "entry, so no test set can be built. Lower --seqid or drop "
                "--external-fasta."
            )

        # Balance only over the eligible clusters; the rest are already placed.
        # Keep the full grouping around so the closing tally can count them:
        # reporting only the eligible subset makes it look as though the forced
        # clusters vanished, when they are simply already assigned to train.
        all_groups = cluster_groups
        cluster_groups = {rep: cluster_groups[rep] for rep in eligible}

        # Discover all sources present
        all_sources: Set[str] = set()
        for members in cluster_groups.values():
            for m in members:
                all_sources.add(self._parse_source(m))
        sources = sorted(all_sources)

        # Per-cluster source composition
        comp: Dict[str, Dict[str, int]] = {}
        for rep, members in cluster_groups.items():
            c: Dict[str, int] = defaultdict(int)
            for m in members:
                c[self._parse_source(m)] += 1
            comp[rep] = dict(c)

        totals = {s: sum(c.get(s, 0) for c in comp.values()) for s in sources}
        tgt_train = {s: totals[s] * (1.0 - self.valid_frac) for s in sources}
        tgt_valid = {s: totals[s] * self.valid_frac for s in sources}

        self.logger.info(
            "Source totals: "
            + ", ".join(f"{s}={totals[s]}" for s in sources)
        )

        # Sort clusters largest-first for better greedy packing
        sorted_reps = sorted(comp, key=lambda r: -sum(comp[r].values()))

        train_reps: Set[str] = set()
        valid_reps: Set[str] = set()
        cur_train = {s: 0 for s in sources}
        cur_valid = {s: 0 for s in sources}

        for rep in sorted_reps:
            c = comp[rep]
            imb_train = sum(
                abs(cur_train[s] + c.get(s, 0) - tgt_train[s])
                + abs(cur_valid[s] - tgt_valid[s])
                for s in sources
            )
            imb_valid = sum(
                abs(cur_train[s] - tgt_train[s])
                + abs(cur_valid[s] + c.get(s, 0) - tgt_valid[s])
                for s in sources
            )
            if imb_train <= imb_valid:
                train_reps.add(rep)
                for s in sources:
                    cur_train[s] += c.get(s, 0)
            else:
                valid_reps.add(rep)
                for s in sources:
                    cur_valid[s] += c.get(s, 0)

        train_reps |= forced_train
        self.logger.info(
            f"Train clusters: {len(train_reps)} "
            f"({len(forced_train)} of them forced by external overlap), "
            f"test clusters: {len(valid_reps)}"
        )
        final_train: Dict[str, int] = defaultdict(int)
        final_valid: Dict[str, int] = defaultdict(int)
        for rep, members in all_groups.items():
            bucket = final_valid if rep in valid_reps else final_train
            for m in members:
                bucket[self._parse_source(m)] += 1
        all_sources_final = sorted(set(final_train) | set(final_valid))
        self.logger.info(
            "Actual train: "
            + ", ".join(f"{s}={final_train[s]}" for s in all_sources_final)
        )
        self.logger.info(
            "Actual test:  "
            + ", ".join(f"{s}={final_valid[s]}" for s in all_sources_final)
        )
        return train_reps, valid_reps

    # -- Entry point -----------------------------------------------------------

    def run(
        self,
        data_dir: Path,
        passed_targets: Optional[List[str]] = None,
        external_fasta: Optional[List[Path]] = None,
        output_dir: Optional[Path] = None,
    ) -> Tuple[Path, Path, Path]:
        """Build train/test split files and chembl_targets.tsv.

        Args:
            data_dir: Root data directory.
            passed_targets: UniProt IDs to include. Reads passed_targets.txt if None.
            external_fasta: Optional list of external FASTA files.
                IDs must be dot-prefixed: >{source}.{entry_id}
            output_dir: Output directory (default: data_dir).

        Returns:
            (train.txt path, test.txt path, chembl_targets.tsv path)
        """
        data_dir = Path(data_dir)
        output_dir = Path(output_dir) if output_dir else data_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        if subprocess.run(["which", "mmseqs"], capture_output=True).returncode != 0:
            raise RuntimeError("mmseqs not found in PATH - install MMseqs2 first")

        if passed_targets is None:
            passed_file = data_dir / "passed_targets.txt"
            passed_targets = [
                l.strip() for l in passed_file.read_text().splitlines() if l.strip()
            ]

        # Load ChEMBL sequences
        chembl_fasta = data_dir / "sequences.fasta"
        if not chembl_fasta.exists():
            raise FileNotFoundError(
                f"sequences.fasta not found: {chembl_fasta}\n"
                "Run filter-proteins first."
            )

        all_chembl = self._load_fasta(chembl_fasta)
        passed_set = set(passed_targets)
        chembl_seqs = {k: v for k, v in all_chembl.items() if k in passed_set}
        self.logger.info(f"ChEMBL sequences: {len(chembl_seqs)}")

        # Load external sequences
        ext_seqs: Dict[str, str] = {}
        if external_fasta:
            for fp in external_fasta:
                ext_seqs.update(self._load_fasta(Path(fp)))
            self.logger.info(f"External sequences: {len(ext_seqs)}")

        with tempfile.TemporaryDirectory(prefix="split_") as tmpdir:
            td = Path(tmpdir)

            # Write combined FASTA (prefix ChEMBL IDs to disambiguate)
            combined = td / "combined.fasta"
            with open(combined, 'w') as f:
                for sid, seq in ext_seqs.items():
                    f.write(f">{sid}\n{seq}\n")
                for sid, seq in chembl_seqs.items():
                    f.write(f">chembl.{sid}\n{seq}\n")

            total = len(ext_seqs) + len(chembl_seqs)
            self.logger.info(
                f"Combined FASTA: {total} sequences. "
                f"Running MMseqs2 at seqid={self.seqid}..."
            )

            member_to_rep = self._run_mmseqs_cluster(combined, td)

            if ext_seqs:
                ext_fasta = td / "external.fasta"
                with open(ext_fasta, 'w') as f:
                    for sid, seq in ext_seqs.items():
                        f.write(f">{sid}\n{seq}\n")
                chembl_only = td / "chembl.fasta"
                with open(chembl_only, 'w') as f:
                    for sid, seq in chembl_seqs.items():
                        f.write(f">chembl.{sid}\n{seq}\n")
                self.blocked_targets = self._find_external_homologues(
                    chembl_only, ext_fasta, td
                )

            # Group members by representative
            cluster_groups: Dict[str, Set[str]] = defaultdict(set)
            for member, rep in member_to_rep.items():
                cluster_groups[rep].add(member)

            sizes = [len(v) for v in cluster_groups.values()]
            self.logger.info(
                f"Clusters: {len(cluster_groups)}, "
                f"sizes min={min(sizes)} max={max(sizes)} "
                f"median={sorted(sizes)[len(sizes)//2]}"
            )

            _, valid_reps = self._greedy_split(cluster_groups)

            # Clustering is not a separation guarantee; check it directly.
            test_targets = {
                m.split(".", 1)[1]
                for rep in valid_reps
                for m in cluster_groups[rep]
                if m.startswith("chembl.")
            }
            train_fasta = td / "train_side.fasta"
            with open(train_fasta, 'w') as f:
                for sid, seq in ext_seqs.items():
                    f.write(f">{sid}\n{seq}\n")
                for rep in cluster_groups:
                    if rep in valid_reps:
                        continue
                    for m in cluster_groups[rep]:
                        if m.startswith("chembl."):
                            t = m.split(".", 1)[1]
                            if t in chembl_seqs:
                                f.write(f">chembltrain.{t}\n{chembl_seqs[t]}\n")
            self.demoted_targets = self._demote_test_targets_close_to_train(
                test_targets, train_fasta, chembl_seqs, td
            )
            if self.demoted_targets:
                self.logger.info(
                    f"Demoted {len(self.demoted_targets)} test targets to train "
                    f"for being within {self.seqid} identity of a train sequence "
                    "despite landing in a ChEMBL-only cluster"
                )

        # Per-member sampling weights
        member_weight: Dict[str, float] = {}
        for rep, members in cluster_groups.items():
            w = round(1.0 / math.log2(len(members) + 1), 2)
            for m in members:
                member_weight[m] = w

        # Load ChEMBL actives from decoys.tsv (one active per line)
        chembl_actives: Dict[str, List[str]] = {}
        for uniprot in passed_targets:
            decoy_tsv = data_dir / uniprot / "decoys.tsv"
            if not decoy_tsv.exists():
                continue
            actives: List[str] = []
            with open(decoy_tsv) as f:
                next(f, None)  # skip header
                for line in f:
                    line = line.strip()
                    if line:
                        actives.append(line.split('\t')[0])
            if actives:
                chembl_actives[uniprot] = actives

        self.logger.info(
            f"ChEMBL targets with actives: {len(chembl_actives)}, "
            f"total active lines: {sum(len(v) for v in chembl_actives.values())}"
        )

        train_lines: List[str] = []
        test_lines: List[str] = []

        # External entries
        for sid in sorted(ext_seqs):
            rep = member_to_rep.get(sid, sid)
            w = member_weight.get(sid, 1.0)
            is_test = rep in valid_reps

            source = self._parse_source(sid)
            dot = sid.find('.')
            entry_id = sid[dot + 1:] if dot > 0 else sid
            line = f"{source}\t{entry_id}\t-\t{w:.2f}"

            (test_lines if is_test else train_lines).append(line)

        # ChEMBL entries (one line per active)
        for uniprot in sorted(chembl_actives):
            clu_key = f"chembl.{uniprot}"
            rep = member_to_rep.get(clu_key, clu_key)
            w = member_weight.get(clu_key, 1.0)
            is_test = rep in valid_reps and uniprot not in self.demoted_targets

            for active_id in chembl_actives[uniprot]:
                line = f"chembl\t{uniprot}\t{active_id}\t{w:.2f}"
                (test_lines if is_test else train_lines).append(line)

        header = "source\tentry_id\tcompound\tweight\n"

        train_path = output_dir / "train.txt"
        test_path = output_dir / "test.txt"
        train_path.write_text(header + '\n'.join(train_lines) + '\n')
        test_path.write_text(header + '\n'.join(test_lines) + '\n')

        self.logger.info(f"Train: {len(train_lines)} lines -> {train_path}")
        self.logger.info(f"Test:  {len(test_lines)} lines -> {test_path}")

        # Generate chembl_targets.tsv
        targets_path = self._write_chembl_targets(
            data_dir, output_dir, chembl_actives, member_to_rep, valid_reps
        )

        return train_path, test_path, targets_path

    def _write_chembl_targets(
        self,
        data_dir: Path,
        output_dir: Path,
        chembl_actives: Dict[str, List[str]],
        member_to_rep: Dict[str, str],
        valid_reps: Set[str],
    ) -> Path:
        """Write chembl_targets.tsv with per-target summary."""
        lines: List[str] = ["uniprot\tsplit\tn_actives\tn_decoys"]
        for uniprot in sorted(chembl_actives):
            clu_key = f"chembl.{uniprot}"
            rep = member_to_rep.get(clu_key, clu_key)
            split_label = "test" if rep in valid_reps else "train"

            n_actives = len(chembl_actives[uniprot])

            # Count decoys
            n_decoys = 0
            decoy_tsv = data_dir / uniprot / "decoys.tsv"
            if decoy_tsv.exists():
                with open(decoy_tsv) as f:
                    next(f, None)  # skip header
                    for line in f:
                        line = line.strip()
                        if line:
                            parts = line.split('\t')
                            if len(parts) >= 2 and parts[1]:
                                n_decoys += len(parts[1].split(';'))

            lines.append(f"{uniprot}\t{split_label}\t{n_actives}\t{n_decoys}")

        targets_path = output_dir / "chembl_targets.tsv"
        targets_path.write_text('\n'.join(lines) + '\n')
        self.logger.info(
            f"ChEMBL targets summary: {len(lines) - 1} targets -> {targets_path}"
        )
        return targets_path

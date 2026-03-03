# chembl_curator/splitter.py

"""Stage 7: Train/test split by sequence-identity clustering.

Clusters all passed targets (and optional external datasets such as PDBbind
or BioLip) using MMseqs2 at a configurable seqid threshold (default 30%).
Clusters are greedily assigned to train/test sets while balancing the ratio
of each data source across both splits.

Per-entry sampling weight = 1 / log2(cluster_size + 1)

Output format (one entry per active compound):
    chembl/{uniprot}         {chembl_id}  batch   {weight:.2f}
    pdbbind/{pdbid}          ligand       single  {weight:.2f}
    biolip/{entry}           {ligname}    single  {weight:.2f}

External FASTA entries must use prefixed IDs:
    >pdbbind.{pdbid}
    >biolip.{pdbid_LIGNAME_chain_num}
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
        seqid: float = 0.3,
        valid_frac: float = 0.1,
        threads: int = 4,
        log_level: str = "INFO",
    ):
        """
        Args:
            seqid: MMseqs2 clustering seqid threshold (default: 0.3 = 30%).
            valid_frac: Fraction of clusters assigned to test/validation set.
            threads: MMseqs2 thread count.
            log_level: Logging level.
        """
        self.seqid = seqid
        self.valid_frac = valid_frac
        self.threads = threads
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # ── I/O helpers ───────────────────────────────────────────────────────────

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

    # ── MMseqs2 clustering ────────────────────────────────────────────────────

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

    # ── Greedy split ─────────────────────────────────────────────────────────

    @staticmethod
    def _get_source(sid: str) -> str:
        if sid.startswith('chembl.'):
            return 'chembl'
        if sid.startswith('pdbbind.'):
            return 'pdbbind'
        if sid.startswith('biolip.'):
            return 'biolip'
        return 'other'

    def _greedy_split(
        self, cluster_groups: Dict[str, Set[str]]
    ) -> Tuple[Set[str], Set[str]]:
        """Assign clusters to train/test, balancing per-source ratios."""
        sources = ('chembl', 'pdbbind', 'biolip')

        # Per-cluster source composition
        comp: Dict[str, Dict[str, int]] = {}
        for rep, members in cluster_groups.items():
            c: Dict[str, int] = defaultdict(int)
            for m in members:
                src = self._get_source(m)
                if src in sources:
                    c[src] += 1
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

        self.logger.info(
            f"Train clusters: {len(train_reps)}, test clusters: {len(valid_reps)}"
        )
        self.logger.info(
            "Actual train: " + ", ".join(f"{s}={cur_train[s]}" for s in sources)
        )
        self.logger.info(
            "Actual test:  " + ", ".join(f"{s}={cur_valid[s]}" for s in sources)
        )
        return train_reps, valid_reps

    # ── Entry point ───────────────────────────────────────────────────────────

    def run(
        self,
        data_dir: Path,
        passed_targets: Optional[List[str]] = None,
        external_fasta: Optional[List[Path]] = None,
        output_dir: Optional[Path] = None,
    ) -> Tuple[Path, Path]:
        """Build train/test split files.

        Args:
            data_dir: Root data directory.
            passed_targets: UniProt IDs to include. Reads passed_targets.txt if None.
            external_fasta: Optional list of external FASTA files (PDBbind, BioLip, etc.).
                IDs must be prefixed with 'pdbbind.' or 'biolip.' as appropriate.
            output_dir: Output directory for train.txt and test.txt
                (default: data_dir).

        Returns:
            (train.txt path, test.txt path)
        """
        data_dir = Path(data_dir)
        output_dir = Path(output_dir) if output_dir else data_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        if subprocess.run(["which", "mmseqs"], capture_output=True).returncode != 0:
            raise RuntimeError("mmseqs not found in PATH — install MMseqs2 first")

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

            if sid.startswith('pdbbind.'):
                pdbid = sid[len('pdbbind.'):]
                line = f"pdbbind/{pdbid} ligand single {w:.2f}"
            elif sid.startswith('biolip.'):
                entry = sid[len('biolip.'):]
                parts = entry.split('_')
                ligname = parts[1] if len(parts) >= 3 else 'UNK'
                line = f"biolip/{entry} {ligname} single {w:.2f}"
            else:
                continue

            (test_lines if is_test else train_lines).append(line)

        # ChEMBL entries (one line per active)
        for uniprot in sorted(chembl_actives):
            clu_key = f"chembl.{uniprot}"
            rep = member_to_rep.get(clu_key, clu_key)
            w = member_weight.get(clu_key, 1.0)
            is_test = rep in valid_reps

            for active_id in chembl_actives[uniprot]:
                line = f"chembl/{uniprot} {active_id} batch {w:.2f}"
                (test_lines if is_test else train_lines).append(line)

        train_path = output_dir / "train.txt"
        test_path = output_dir / "test.txt"
        train_path.write_text('\n'.join(train_lines) + '\n')
        test_path.write_text('\n'.join(test_lines) + '\n')

        self.logger.info(f"Train: {len(train_lines)} lines → {train_path}")
        self.logger.info(f"Test:  {len(test_lines)} lines → {test_path}")

        return train_path, test_path

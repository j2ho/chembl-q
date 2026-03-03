# chembl_curator/receptor_similarity.py

"""Stage 5: Compute pairwise receptor similarity.

Two independent computations:

1. Sequence identity (MMseqs2 all-vs-all)
   Input:  sequences.fasta (written by protein_filter.py)
   Output: pairwise_seqid.tsv  (query, target, seqid)

2. Pocket RMSD (nuri TMalign, all-vs-all)
   Input:  best_structure.tsv + aligned/*.pdb + pocket_info.csv
   Output: pairwise_pocket_rmsd.tsv  (target_a, target_b, tm_score, pocket_rmsd, n_matched)
"""

import csv
import logging
import os
import subprocess
import tempfile
import time
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ── Multiprocessing workers (must be module-level for pickling) ───────────────

def _suppress_stderr() -> None:
    """Suppress C++ stderr in worker processes (nuri/TMalign retry messages)."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 2)
    os.close(devnull)


def _pocket_rmsd_worker(
    args: Tuple,
) -> Tuple[str, str, float, float, int]:
    """Compute pocket RMSD for one target pair.

    args: (target_a, target_b, ca_a, ca_b, pocket_idx_a, pocket_idx_b)
    Returns: (target_a, target_b, tm_score, pocket_rmsd, n_matched)
             pocket_rmsd = -1.0 when alignment is invalid or pocket too small.
    """
    import nuri._log_interface
    nuri._log_interface.set_log_level(4)  # suppress C++ warnings
    from nuri.tools.tm import TMAlign

    ta, tb, ca_a, ca_b, pidx_a, pidx_b = args

    if len(ca_a) < 5 or len(ca_b) < 5:
        return (ta, tb, -1.0, -1.0, 0)

    try:
        tma = TMAlign(ca_a, ca_b)
        xform, tm_score = tma.score()
    except (ValueError, RuntimeError):
        return (ta, tb, -1.0, -1.0, 0)

    if len(pidx_a) == 0 or len(pidx_b) == 0:
        return (ta, tb, float(tm_score), -1.0, 0)

    # Residue-level alignment from TM-align
    aligned = tma.aligned_pairs()
    pocket_a_set = set(pidx_a.tolist())
    pocket_b_set = set(pidx_b.tolist())

    pocket_pairs = [
        (qi, ti) for qi, ti in aligned
        if qi in pocket_a_set and ti in pocket_b_set
    ]
    n_matched = len(pocket_pairs)

    if n_matched < 3:
        return (ta, tb, float(tm_score), -1.0, n_matched)

    pairs = np.array(pocket_pairs)
    # Apply xform (4×4 homogeneous) to query pocket Cα, then RMSD vs template
    P_raw = ca_a[pairs[:, 0]]
    n = P_raw.shape[0]
    hom = np.hstack([P_raw, np.ones((n, 1))])
    P = (xform @ hom.T).T[:, :3]
    Q = ca_b[pairs[:, 1]]
    rmsd = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))

    return (ta, tb, float(tm_score), rmsd, n_matched)


# ── Main class ────────────────────────────────────────────────────────────────

class ReceptorSimilarity:
    """Compute pairwise receptor similarity (sequence identity and/or pocket RMSD)."""

    def __init__(self, log_level: str = "INFO"):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # ── Sequence identity (MMseqs2) ───────────────────────────────────────────

    def compute_seqid(
        self,
        data_dir: Path,
        fasta: Optional[Path] = None,
        threads: int = 4,
        output: Optional[Path] = None,
    ) -> Path:
        """Run MMseqs2 all-vs-all search to compute pairwise sequence identities.

        Args:
            data_dir: Root data directory.
            fasta: Input FASTA path (default: data_dir/sequences.fasta).
            threads: MMseqs2 thread count.
            output: Output TSV path (default: data_dir/pairwise_seqid.tsv).

        Returns:
            Path to output TSV (query, target, seqid columns).
        """
        data_dir = Path(data_dir)
        fasta = fasta or (data_dir / "sequences.fasta")
        output = output or (data_dir / "pairwise_seqid.tsv")

        if not fasta.exists():
            raise FileNotFoundError(f"sequences.fasta not found: {fasta}")

        if subprocess.run(["which", "mmseqs"], capture_output=True).returncode != 0:
            raise RuntimeError("mmseqs not found in PATH — install MMseqs2 first")

        self.logger.info(f"Running MMseqs2 all-vs-all on {fasta} ({threads} threads)")
        t0 = time.time()

        with tempfile.TemporaryDirectory(prefix="seqid_") as tmpdir:
            td = Path(tmpdir)
            tmp = td / "tmp"
            tmp.mkdir()

            # Create MMseqs2 database
            subprocess.run(
                ["mmseqs", "createdb", str(fasta), str(td / "seqdb")],
                check=True, capture_output=True,
            )

            # All-vs-all search (sensitivity 7.5, report up to all targets)
            subprocess.run(
                [
                    "mmseqs", "search",
                    str(td / "seqdb"), str(td / "seqdb"),
                    str(td / "resultdb"), str(tmp),
                    "--threads", str(threads),
                    "--alignment-mode", "3",
                    "-s", "7.5",
                    "--max-seqs", "100000",
                ],
                check=True, capture_output=True,
            )

            # Convert to TSV (query, target, %identity 0–100)
            result_tsv = td / "result.tsv"
            subprocess.run(
                [
                    "mmseqs", "convertalis",
                    str(td / "seqdb"), str(td / "seqdb"),
                    str(td / "resultdb"), str(result_tsv),
                    "--format-output", "query,target,pident",
                ],
                check=True, capture_output=True,
            )

            # Write output (convert pident 0–100 → 0–1, skip self-hits)
            n_pairs = 0
            with open(result_tsv) as fin, open(output, 'w') as fout:
                fout.write("query\ttarget\tseqid\n")
                for line in fin:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    q, t, pident = parts[0], parts[1], parts[2]
                    if q == t:
                        continue
                    fout.write(f"{q}\t{t}\t{float(pident)/100.0:.4f}\n")
                    n_pairs += 1

        self.logger.info(
            f"Wrote {n_pairs} pairs to {output} in {time.time()-t0:.1f}s"
        )
        return output

    # ── Pocket RMSD (nuri TMalign) ────────────────────────────────────────────

    def compute_pocket_rmsd(
        self,
        data_dir: Path,
        best_structure_tsv: Optional[Path] = None,
        pocket_radius: float = 10.0,
        workers: int = 4,
        output: Optional[Path] = None,
    ) -> Path:
        """Compute all-vs-all pocket RMSD using nuri TMalign.

        Args:
            data_dir: Root data directory.
            best_structure_tsv: TSV mapping uniprot → pdbid_chain
                (default: data_dir/best_structure.tsv).
            pocket_radius: Cα atoms within this radius of ligand center define the pocket.
            workers: Number of parallel worker processes.
            output: Output TSV path (default: data_dir/pairwise_pocket_rmsd.tsv).

        Returns:
            Path to output TSV.
        """
        import nuri
        import nuri._log_interface
        nuri._log_interface.set_log_level(4)

        data_dir = Path(data_dir)
        best_structure_tsv = best_structure_tsv or (data_dir / "best_structure.tsv")
        output = output or (data_dir / "pairwise_pocket_rmsd.tsv")

        if not best_structure_tsv.exists():
            raise FileNotFoundError(
                f"best_structure.tsv not found: {best_structure_tsv}\n"
                "Run filter-proteins first."
            )

        # Load best structure mapping: uniprot → pdbid_chain (e.g. "3mms_A")
        beststr: Dict[str, str] = {}
        with open(best_structure_tsv) as f:
            next(f, None)  # skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    beststr[parts[0]] = parts[1]

        self.logger.info(f"Pre-loading structures for {len(beststr)} targets")
        t0 = time.time()

        target_data: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        skipped: List[Tuple[str, str]] = []

        for uniprot, pdbid_chain in beststr.items():
            pdb_path = data_dir / uniprot / "aligned" / f"{pdbid_chain}.pdb"
            if not pdb_path.exists():
                skipped.append((uniprot, 'pdb_missing'))
                continue

            ca = self._extract_ca(pdb_path, nuri)
            if len(ca) < 5:
                skipped.append((uniprot, 'too_few_ca'))
                continue

            lig_center = self._get_lig_center(data_dir / uniprot, pdbid_chain)
            if lig_center is None:
                skipped.append((uniprot, 'no_lig_center'))
                continue

            pocket_idx = np.where(
                np.linalg.norm(ca - lig_center, axis=1) <= pocket_radius
            )[0]
            if len(pocket_idx) == 0:
                skipped.append((uniprot, 'empty_pocket'))
                continue

            target_data[uniprot] = (ca, pocket_idx)

        valid = sorted(target_data)
        n_pairs = len(valid) * (len(valid) - 1) // 2
        self.logger.info(
            f"Loaded {len(valid)} targets, {n_pairs} pairs "
            f"(skipped {len(skipped)}) in {time.time()-t0:.1f}s"
        )

        pair_args = [
            (ta, tb, target_data[ta][0], target_data[tb][0],
             target_data[ta][1], target_data[tb][1])
            for ta, tb in combinations(valid, 2)
        ]

        self.logger.info(f"Computing pocket RMSD ({workers} workers)...")
        t0 = time.time()

        with Pool(processes=workers, initializer=_suppress_stderr) as pool:
            results: List[Tuple] = list(
                pool.imap_unordered(_pocket_rmsd_worker, pair_args, chunksize=256)
            )

        self.logger.info(f"Done in {time.time()-t0:.1f}s")

        results.sort(key=lambda x: (x[0], x[1]))
        with open(output, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(
                ['target_a', 'target_b', 'tm_score', 'pocket_rmsd', 'n_matched']
            )
            for row in results:
                writer.writerow(row)

        valid_rmsds = [r[3] for r in results if r[3] >= 0]
        if valid_rmsds:
            self.logger.info(
                f"Wrote {len(results)} pairs: {len(valid_rmsds)} with valid RMSD, "
                f"median={np.median(valid_rmsds):.2f} Å"
            )
        else:
            self.logger.warning(f"Wrote {len(results)} pairs: no valid pocket RMSDs")

        return output

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _extract_ca(self, pdb_path: Path, nuri_mod) -> np.ndarray:
        """Extract Cα coordinates from a PDB file using nuri."""
        mol = list(nuri_mod.readfile('pdb', str(pdb_path), sanitize=False))[0]
        ca = []
        for i in range(mol.num_atoms()):
            a = mol[i]
            if a.name.strip() == 'CA' and a.element_symbol == 'C':
                ca.append(a.get_pos(0))
        return np.array(ca) if ca else np.empty((0, 3))

    def _get_lig_center(
        self, target_dir: Path, pdbid_chain: str
    ) -> Optional[np.ndarray]:
        """Read ligand center from pocket_info.csv for the given pdbid_chain.

        pdbid_chain is the aligned file stem, e.g. "3mms_A".
        Matches against PDB_ID + Chain columns in pocket_info.csv.
        """
        pocket_csv = target_dir / "pocket_info.csv"
        if not pocket_csv.exists():
            return None

        pdbid, chain = pdbid_chain.split('_', 1)
        pdbid_upper = pdbid.upper()

        with open(pocket_csv) as f:
            for row in csv.DictReader(f):
                if row['PDB_ID'].upper() == pdbid_upper and row['Chain'] == chain:
                    return np.array([
                        float(row['Center_X']),
                        float(row['Center_Y']),
                        float(row['Center_Z']),
                    ])
        return None

    # ── Entry point ───────────────────────────────────────────────────────────

    def run(
        self,
        data_dir: Path,
        mode: str = "both",
        seqid_threads: int = 4,
        workers: int = 4,
        pocket_radius: float = 10.0,
    ) -> dict:
        """Run receptor similarity computation.

        Args:
            data_dir: Root data directory.
            mode: "seqid", "pocket", or "both".
            seqid_threads: Threads for MMseqs2.
            workers: Worker processes for pocket RMSD.
            pocket_radius: Pocket radius in Å.

        Returns:
            dict with paths to output files.
        """
        data_dir = Path(data_dir)
        results = {}
        if mode in ("seqid", "both"):
            results['seqid_tsv'] = self.compute_seqid(
                data_dir, threads=seqid_threads
            )
        if mode in ("pocket", "both"):
            results['pocket_rmsd_tsv'] = self.compute_pocket_rmsd(
                data_dir, pocket_radius=pocket_radius, workers=workers
            )
        return results

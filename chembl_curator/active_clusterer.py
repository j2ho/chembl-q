# chembl_curator/active_clusterer.py

"""Stage 3: Butina clustering of actives per target.

Reads actives.tsv (written by curator.py) and clusters using Morgan radius=2
2048-bit fingerprints with Tanimoto distance. One representative per cluster
is selected (best pChEMBL). Output is actives_clustered.tsv per target.
"""

import logging
import multiprocessing as mp
from pathlib import Path
from typing import List, Optional, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Cluster import Butina


def _cluster_target_worker(args: Tuple) -> int:
    """Top-level worker for multiprocessing (must be picklable)."""
    target_dir, dist_thresh = args
    clusterer = ActiveClusterer(dist_thresh=dist_thresh, log_level="WARNING")
    return clusterer.cluster_target(Path(target_dir))


class ActiveClusterer:
    """Cluster actives per target using Butina algorithm on Morgan fingerprints."""

    def __init__(self, dist_thresh: float = 0.3, log_level: str = "INFO"):
        """
        Args:
            dist_thresh: Tanimoto distance threshold (default 0.3 = similarity >= 0.7)
            log_level: Logging level
        """
        self.dist_thresh = dist_thresh
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # ── I/O ──────────────────────────────────────────────────────────────────

    def _read_actives(self, target_dir: Path) -> List[Tuple[str, float, str]]:
        """Read actives from actives.tsv (preferred) or .smi files (fallback).

        Returns list of (chembl_id, pchembl, smiles).
        """
        tsv = target_dir / "actives.tsv"
        if tsv.exists():
            return self._read_actives_tsv(tsv)

        # Fallback: read individual .smi files (no affinity available)
        smi_dir = target_dir / "comps" / "smiles"
        if not smi_dir.exists():
            return []
        actives = []
        for smi_file in sorted(smi_dir.glob("*.smi")):
            smiles = smi_file.read_text().strip()
            if smiles:
                actives.append((smi_file.stem, 0.0, smiles))
        return actives

    def _read_actives_tsv(self, path: Path) -> List[Tuple[str, float, str]]:
        actives = []
        with open(path) as f:
            next(f, None)  # skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                try:
                    actives.append((parts[0], float(parts[1]), parts[2]))
                except ValueError:
                    continue
        return actives

    # ── Clustering ───────────────────────────────────────────────────────────

    def _butina_cluster(
        self, actives: List[Tuple[str, float, str]]
    ) -> List[List[Tuple[str, float, str]]]:
        """Run Butina clustering and return list of clusters.

        Each cluster is sorted by pChEMBL descending (best first).
        """
        valid = []
        fps = []
        for chembl_id, pchembl, smiles in actives:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Invalid SMILES for {chembl_id}, skipping")
                continue
            fps.append(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048))
            valid.append((chembl_id, pchembl, smiles))

        if not valid:
            return []

        n = len(fps)
        dists = []
        for i in range(1, n):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
            dists.extend(1.0 - s for s in sims)

        clusters_idx = Butina.ClusterData(dists, n, self.dist_thresh, isDistData=True)

        clusters = []
        for cluster in clusters_idx:
            members = [valid[i] for i in cluster]
            members.sort(key=lambda x: -x[1])  # highest pChEMBL first
            clusters.append(members)
        return clusters

    # ── Public API ───────────────────────────────────────────────────────────

    def cluster_target(self, target_dir: Path) -> int:
        """Cluster actives for a single target.

        Writes actives_clustered.tsv to target_dir.

        Returns:
            Number of clusters written (0 if target was skipped).
        """
        target_dir = Path(target_dir)
        actives = self._read_actives(target_dir)
        if not actives:
            return 0

        clusters = self._butina_cluster(actives)
        if not clusters:
            return 0

        out_path = target_dir / "actives_clustered.tsv"
        with open(out_path, 'w') as f:
            f.write("chembl_id\tpchembl\tsmiles\tcluster_size\n")
            for cluster in clusters:
                rep = cluster[0]
                f.write(f"{rep[0]}\t{rep[1]:.4f}\t{rep[2]}\t{len(cluster)}\n")

        return len(clusters)

    def run(
        self,
        data_dir: Path,
        passed_targets: Optional[List[str]] = None,
        workers: int = 1,
    ) -> dict:
        """Run clustering on all passed targets.

        Args:
            data_dir: Root data directory containing per-target subdirs.
            passed_targets: List of UniProt IDs. Reads from passed_targets.txt if None.
            workers: Number of parallel workers.

        Returns:
            dict with n_targets, total_actives, total_clusters.
        """
        data_dir = Path(data_dir)

        if passed_targets is None:
            passed_file = data_dir / "passed_targets.txt"
            if not passed_file.exists():
                raise FileNotFoundError(f"passed_targets.txt not found in {data_dir}")
            passed_targets = [
                l.strip() for l in passed_file.read_text().splitlines() if l.strip()
            ]

        self.logger.info(
            f"Clustering {len(passed_targets)} targets "
            f"(dist_thresh={self.dist_thresh}, workers={workers})"
        )

        target_dirs = [data_dir / t for t in passed_targets]

        if workers == 1:
            n_clusters = [self.cluster_target(td) for td in target_dirs]
        else:
            args = [(str(td), self.dist_thresh) for td in target_dirs]
            with mp.Pool(workers) as pool:
                n_clusters = pool.map(_cluster_target_worker, args)

        total_actives = sum(
            len(self._read_actives(data_dir / t)) for t in passed_targets
        )
        total_clusters = sum(n_clusters)

        self.logger.info(
            f"Done: {total_actives} actives → {total_clusters} cluster representatives"
        )
        return {
            "n_targets": len(passed_targets),
            "total_actives": total_actives,
            "total_clusters": total_clusters,
        }

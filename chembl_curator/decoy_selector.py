# chembl_curator/decoy_selector.py

"""Stage 6: Receptor-aware decoy selection.

For each active compound in each passed target, selects up to max_decoys
property-matched, chemically dissimilar decoy compounds from the global pool.
Compounds active against receptors similar to the query target (by sequence
identity and/or pocket RMSD) are excluded from the decoy pool.

Property matching windows (defaults):
    ±50 Da MW, ±2 cLogP, ±50 Å² TPSA, ±2 HBD, ±2 HBA, ±1 aromatic ring

Input:
    compound_pool.pkl        — written by compound_pool.py
    pairwise_seqid.tsv       — optional, written by receptor_similarity.py
    pairwise_pocket_rmsd.tsv — optional, written by receptor_similarity.py

Output:
    {target_dir}/decoys.tsv  — tab-separated:
        active_chembl_id  <TAB>  decoy_chembl_ids (semicolon-separated)
"""

import logging
import math
import pickle
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set

from rdkit import DataStructs


class DecoySelector:
    """Receptor-aware, property-matched decoy selector."""

    def __init__(
        self,
        max_decoys: int = 30,
        seqid_thresh: float = 0.6,
        pocket_rmsd_thresh: float = 3.0,
        exclusion_mode: str = "or",
        tanimoto_thresh: float = 0.3,
        mw_window: float = 50.0,
        clogp_window: float = 2.0,
        tpsa_window: float = 50.0,
        hbd_window: int = 2,
        hba_window: int = 2,
        arm_ring_window: int = 1,
        max_selection_count: Optional[int] = None,
        seed: int = 42,
        log_level: str = "INFO",
    ):
        """
        Args:
            max_decoys: Maximum decoys per active compound.
            seqid_thresh: Sequence identity threshold for receptor exclusion.
            pocket_rmsd_thresh: Pocket RMSD threshold (Å) for receptor exclusion.
            exclusion_mode: "or" = exclude if seqid OR pocket RMSD matches;
                            "and" = exclude only if BOTH match.
            tanimoto_thresh: Max Tanimoto similarity between active and decoy.
            mw_window, clogp_window, tpsa_window, hbd_window, hba_window,
            arm_ring_window: Property matching windows.
            max_selection_count: Max times a compound may be used as a decoy.
                Auto-computed if None.
            seed: Random seed for reproducibility.
            log_level: Logging level.
        """
        self.max_decoys = max_decoys
        self.seqid_thresh = seqid_thresh
        self.pocket_rmsd_thresh = pocket_rmsd_thresh
        self.exclusion_mode = exclusion_mode
        self.tanimoto_thresh = tanimoto_thresh
        self.mw_window = mw_window
        self.clogp_window = clogp_window
        self.tpsa_window = tpsa_window
        self.hbd_window = hbd_window
        self.hba_window = hba_window
        self.arm_ring_window = arm_ring_window
        self.max_selection_count = max_selection_count
        self.seed = seed
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # ── Similarity loading ────────────────────────────────────────────────────

    def _load_seqid_similar(self, tsv_path: Path) -> Dict[str, Set[str]]:
        similar: Dict[str, Set[str]] = defaultdict(set)
        with open(tsv_path) as f:
            next(f, None)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                q, t, seqid = parts[0], parts[1], float(parts[2])
                if q != t and seqid >= self.seqid_thresh:
                    similar[q].add(t)
                    similar[t].add(q)
        n = sum(len(v) for v in similar.values()) // 2
        self.logger.info(
            f"Seqid similarity (>={self.seqid_thresh}): {n} pairs, "
            f"{len(similar)} targets with similar peers"
        )
        return dict(similar)

    def _load_pocket_similar(self, tsv_path: Path) -> Dict[str, Set[str]]:
        similar: Dict[str, Set[str]] = defaultdict(set)
        with open(tsv_path) as f:
            next(f, None)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                ta, tb = parts[0], parts[1]
                rmsd_str = parts[3]
                try:
                    rmsd = float(rmsd_str)
                except ValueError:
                    continue
                if 0.0 <= rmsd <= self.pocket_rmsd_thresh:
                    similar[ta].add(tb)
                    similar[tb].add(ta)
        n = sum(len(v) for v in similar.values()) // 2
        self.logger.info(
            f"Pocket RMSD similarity (<={self.pocket_rmsd_thresh} Å): {n} pairs, "
            f"{len(similar)} targets with similar pockets"
        )
        return dict(similar)

    def _build_combined_similar(
        self,
        seqid_sim: Dict[str, Set[str]],
        pocket_sim: Dict[str, Set[str]],
    ) -> Dict[str, Set[str]]:
        all_targets = set(seqid_sim) | set(pocket_sim)
        combined: Dict[str, Set[str]] = {}
        for t in all_targets:
            s = seqid_sim.get(t, set())
            p = pocket_sim.get(t, set())
            combined[t] = s | p if self.exclusion_mode == "or" else s & p
        return combined

    # ── Selection ─────────────────────────────────────────────────────────────

    def _select_decoys_for_active(
        self,
        ad: dict,
        excluded: Set[str],
        pool: Dict[str, dict],
        pool_keys: List[str],
        sel_count: Dict[str, int],
        max_sel: int,
        stats: Dict[str, int],
    ) -> List[str]:
        """Select up to max_decoys for a single active compound."""
        decoys: List[str] = []
        shuffled = pool_keys.copy()
        random.shuffle(shuffled)

        for cid in shuffled:
            if len(decoys) >= self.max_decoys:
                break

            if cid in excluded:
                stats['n_excl_receptor'] += 1
                continue

            if sel_count[cid] >= max_sel:
                stats['n_excl_count'] += 1
                continue

            cd = pool[cid]

            # Property filter
            if (
                abs(cd['mw'] - ad['mw']) > self.mw_window
                or abs(cd['clogp'] - ad['clogp']) > self.clogp_window
                or abs(cd['tpsa'] - ad['tpsa']) > self.tpsa_window
                or abs(cd['n_hbd'] - ad['n_hbd']) > self.hbd_window
                or abs(cd['n_hba'] - ad['n_hba']) > self.hba_window
                or abs(cd['n_arm_ring'] - ad['n_arm_ring']) > self.arm_ring_window
            ):
                stats['n_excl_prop'] += 1
                continue

            # Chemical dissimilarity filter
            if DataStructs.TanimotoSimilarity(ad['fp'], cd['fp']) > self.tanimoto_thresh:
                stats['n_excl_sim'] += 1
                continue

            decoys.append(cid)
            sel_count[cid] += 1

        return decoys

    # ── Entry point ───────────────────────────────────────────────────────────

    def run(
        self,
        data_dir: Path,
        pool_pkl: Optional[Path] = None,
        seqid_tsv: Optional[Path] = None,
        pocket_rmsd_tsv: Optional[Path] = None,
        passed_targets: Optional[List[str]] = None,
    ) -> dict:
        """Run decoy selection for all passed targets.

        Args:
            data_dir: Root data directory.
            pool_pkl: Compound pool pickle (default: data_dir/compound_pool.pkl).
            seqid_tsv: Pairwise seqid TSV (default: data_dir/pairwise_seqid.tsv).
            pocket_rmsd_tsv: Pairwise pocket RMSD TSV
                (default: data_dir/pairwise_pocket_rmsd.tsv).
            passed_targets: List of UniProt IDs. Reads passed_targets.txt if None.

        Returns:
            dict with selection statistics.
        """
        data_dir = Path(data_dir)
        pool_pkl = pool_pkl or (data_dir / "compound_pool.pkl")
        seqid_tsv = seqid_tsv or (data_dir / "pairwise_seqid.tsv")
        pocket_rmsd_tsv = pocket_rmsd_tsv or (data_dir / "pairwise_pocket_rmsd.tsv")

        if not pool_pkl.exists():
            raise FileNotFoundError(
                f"compound_pool.pkl not found: {pool_pkl}\n"
                "Run build-pool first."
            )

        if passed_targets is None:
            passed_file = data_dir / "passed_targets.txt"
            passed_targets = [
                l.strip() for l in passed_file.read_text().splitlines() if l.strip()
            ]

        random.seed(self.seed)

        # Load compound pool
        with open(pool_pkl, 'rb') as f:
            data = pickle.load(f)
        pool: Dict[str, dict] = data['pool']
        target_actives: Dict[str, set] = data['target_actives']
        pool_keys = list(pool)
        self.logger.info(f"Pool: {len(pool)} compounds, {len(target_actives)} targets")

        # Load receptor similarity
        seqid_sim: Dict[str, Set[str]] = {}
        if seqid_tsv.exists():
            seqid_sim = self._load_seqid_similar(seqid_tsv)
        else:
            self.logger.warning(f"No seqid TSV at {seqid_tsv} — skipping seqid exclusion")

        pocket_sim: Dict[str, Set[str]] = {}
        if pocket_rmsd_tsv.exists():
            pocket_sim = self._load_pocket_similar(pocket_rmsd_tsv)
        else:
            self.logger.warning(
                f"No pocket RMSD TSV at {pocket_rmsd_tsv} — skipping pocket exclusion"
            )

        similar_targets = self._build_combined_similar(seqid_sim, pocket_sim)

        # Auto-compute max_selection_count to distribute decoys evenly
        total_actives = sum(len(target_actives.get(t, set())) for t in passed_targets)
        max_sel = (
            self.max_selection_count
            if self.max_selection_count is not None
            else math.ceil(total_actives * self.max_decoys / len(pool)) + 1
        )
        self.logger.info(
            f"total_actives={total_actives}, max_decoys={self.max_decoys}, "
            f"max_sel_count={max_sel}, exclusion_mode={self.exclusion_mode}"
        )

        sel_count: Dict[str, int] = defaultdict(int)
        stats: Dict[str, int] = defaultdict(int)

        for uniprot in passed_targets:
            actives = list(target_actives.get(uniprot, set()))
            if not actives:
                continue

            # Build exclusion set: self + all similar receptors' actives
            excluded: Set[str] = set(actives)
            for sim_target in similar_targets.get(uniprot, set()):
                excluded |= target_actives.get(sim_target, set())

            results: List[tuple] = []
            for active_id in actives:
                if active_id not in pool:
                    continue
                stats['n_actives'] += 1

                decoys = self._select_decoys_for_active(
                    pool[active_id], excluded, pool, pool_keys,
                    sel_count, max_sel, stats,
                )
                if len(decoys) < self.max_decoys:
                    stats['n_underfilled'] += 1
                stats['n_total_decoys'] += len(decoys)
                results.append((active_id, decoys))

            out_path = data_dir / uniprot / "decoys.tsv"
            with open(out_path, 'w') as f:
                f.write("active_chembl_id\tdecoy_chembl_ids\n")
                for active_id, decoys in results:
                    f.write(f"{active_id}\t{';'.join(decoys)}\n")

        self.logger.info(
            f"Done — active-decoy pairs: {stats['n_total_decoys']}, "
            f"underfilled: {stats['n_underfilled']}\n"
            f"Exclusions: receptor={stats['n_excl_receptor']}, "
            f"count={stats['n_excl_count']}, "
            f"property={stats['n_excl_prop']}, "
            f"similarity={stats['n_excl_sim']}"
        )
        return dict(stats)

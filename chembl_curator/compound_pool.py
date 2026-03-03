# chembl_curator/compound_pool.py

"""Stage 4: Build a global compound pool from all clustered actives.

Reads actives_clustered.tsv from each passed target directory, deduplicates
by ChEMBL ID (keeping best pChEMBL entry), and computes molecular properties
and Morgan fingerprints for every unique compound.

Output: compound_pool.pkl with structure:
    {
        'pool': {
            chembl_id: {
                'smiles': str,
                'n_heavy': int,
                'mw': float,
                'clogp': float,
                'tpsa': float,
                'n_arm_ring': int,
                'n_hbd': int,
                'n_hba': int,
                'fp': rdkit ExplicitBitVect,
                'target_set': set of UniProt IDs,
            }
        },
        'target_actives': {uniprot_id: set of chembl_ids}
    }
"""

import logging
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors


class CompoundPool:
    """Build and serialize the global compound pool."""

    def __init__(self, log_level: str = "INFO"):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    @staticmethod
    def _compute_properties(smiles: str):
        """Compute properties + Morgan2 fingerprint.

        Returns (n_heavy, mw, clogp, tpsa, n_arm_ring, n_hbd, n_hba, fp) or None.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return (
            mol.GetNumHeavyAtoms(),
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            Descriptors.TPSA(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048),
        )

    def build(
        self,
        data_dir: Path,
        passed_targets: Optional[List[str]] = None,
        output: Optional[Path] = None,
    ) -> Path:
        """Build compound pool from actives_clustered.tsv files.

        Args:
            data_dir: Root data directory.
            passed_targets: List of UniProt IDs. Reads passed_targets.txt if None.
            output: Output pickle path (default: data_dir/compound_pool.pkl).

        Returns:
            Path to the written pickle file.
        """
        data_dir = Path(data_dir)

        if passed_targets is None:
            passed_file = data_dir / "passed_targets.txt"
            if not passed_file.exists():
                raise FileNotFoundError(f"passed_targets.txt not found in {data_dir}")
            passed_targets = [
                l.strip() for l in passed_file.read_text().splitlines() if l.strip()
            ]

        if output is None:
            output = data_dir / "compound_pool.pkl"

        pool: Dict[str, dict] = {}
        best_pchembl: Dict[str, float] = {}
        target_actives: Dict[str, set] = {}
        skipped = 0

        self.logger.info(f"Building compound pool from {len(passed_targets)} targets")

        for i, uniprot in enumerate(passed_targets):
            tsv_path = data_dir / uniprot / "actives_clustered.tsv"
            if not tsv_path.exists():
                skipped += 1
                continue

            target_set: set = set()
            with open(tsv_path) as f:
                next(f, None)  # skip header
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    chembl_id = parts[0]
                    try:
                        pchembl = float(parts[1])
                    except ValueError:
                        pchembl = 0.0
                    smiles = parts[2]
                    target_set.add(chembl_id)

                    if chembl_id in pool:
                        pool[chembl_id]['target_set'].add(uniprot)
                        if pchembl > best_pchembl[chembl_id]:
                            best_pchembl[chembl_id] = pchembl
                            pool[chembl_id]['smiles'] = smiles
                        continue

                    props = self._compute_properties(smiles)
                    if props is None:
                        self.logger.warning(f"Invalid SMILES for {chembl_id}")
                        continue

                    n_heavy, mw, clogp, tpsa, n_arm_ring, n_hbd, n_hba, fp = props
                    pool[chembl_id] = {
                        'smiles': smiles,
                        'n_heavy': n_heavy,
                        'mw': mw,
                        'clogp': clogp,
                        'tpsa': tpsa,
                        'n_arm_ring': n_arm_ring,
                        'n_hbd': n_hbd,
                        'n_hba': n_hba,
                        'fp': fp,
                        'target_set': {uniprot},
                    }
                    best_pchembl[chembl_id] = pchembl

            target_actives[uniprot] = target_set

            if (i + 1) % 100 == 0 or (i + 1) == len(passed_targets):
                self.logger.info(
                    f"  {i+1}/{len(passed_targets)} targets, pool size: {len(pool)}"
                )

        if skipped:
            self.logger.warning(
                f"{skipped} targets had no actives_clustered.tsv — "
                "run cluster-actives first"
            )

        data = {'pool': pool, 'target_actives': target_actives}
        with open(output, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

        self.logger.info(
            f"Saved to {output}: {len(pool)} unique compounds, "
            f"{len(target_actives)} targets"
        )

        # Compound-per-target-count distribution
        multi: Dict[int, int] = defaultdict(int)
        for d in pool.values():
            multi[len(d['target_set'])] += 1
        for n_tgt in sorted(multi):
            self.logger.info(f"  Active in {n_tgt} target(s): {multi[n_tgt]} compounds")

        return output

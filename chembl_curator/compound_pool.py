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
        'target_actives': {uniprot_id: set of chembl_ids},
        'target_ref_fps': {uniprot_id: [rdkit ExplicitBitVect, ...]}
    }

`target_actives` holds CLUSTER REPRESENTATIVES and decides which actives get
decoys drawn for them. Do not widen it: there are 2.7x more actives than
representatives and the pool cannot feed decoys for all of them.

`target_ref_fps` is the separate reference set the cross-active bar compares
against, and it is deliberately wider than `target_actives`:

  * every active in actives.tsv, not just the representatives. Stage 6 used to
    compare against representatives alone, which let 94 assignments reach 1.0
    similarity to a non-representative active of their own target.
  * every cocrystal ligand of the target, from pocket_info.csv via the CCD
    asset. These have no ChEMBL id, so no id-based exclusion can see them, yet
    a solved structure puts them in the pocket. 14 assignments were the same
    molecule as a cocrystal ligand of the target they were a decoy for.
  * every weak binder from weak.tsv: activity observed, but above
    active_max_nm. A weak binder is still a binder.

Measured INACTIVES are deliberately absent. They are confirmed non-binders, so
resembling one is a point in a decoy's favour; excluding on it would invert the
rule. They stay excluded by ChEMBL id through measured.tsv.
"""

import csv
import gzip
import logging
import math
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors


class CompoundPool:
    """Build and serialize the global compound pool."""

    # Same bar check_same_pocket used to admit the target, applied again here
    # as a guard. On the build this was written for it drops nothing: the
    # widest pocket_info spread measured across 1,253 targets was exactly
    # 10.00 A. It exists so that loosening that filter upstream cannot quietly
    # feed an off-site ligand into the cross-active reference.
    MAX_POCKET_SPREAD = 10.0

    def __init__(self, log_level: str = "INFO"):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        self._ccd: Optional[Dict[str, str]] = None

    def _ccd_smiles(self) -> Dict[str, str]:
        """comp_id -> SMILES from the bundled CCD asset, loaded once.

        Built by tools/build_ccd_smiles.py. Absent means cocrystal ligands
        contribute nothing, which is a silent weakening of the bar, so say so.
        """
        if self._ccd is not None:
            return self._ccd
        asset = Path(__file__).parent / "assets" / "ccd_smiles.tsv.gz"
        table: Dict[str, str] = {}
        if asset.exists():
            with gzip.open(asset, 'rt') as handle:
                next(handle, None)
                for line in handle:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2 and parts[1]:
                        table[parts[0]] = parts[1]
            self.logger.info(f"CCD asset: {len(table)} components")
        else:
            self.logger.warning(
                f"{asset} missing - cocrystal ligands will not enter the "
                f"cross-active reference. Run tools/build_ccd_smiles.py"
            )
        self._ccd = table
        return table

    def _pocket_centre(self, target_dir: Path, rows: List[tuple]):
        """Where the pocket is, in the aligned frame.

        Anchored on the representative structure's ligand, because that is the
        pocket the release actually ships and featurises. A plain median over
        every row is not safe: one PDB entry contributing several off-site
        ligands drags the median away from the real site, and then the distance
        guard below measures from the wrong place.

        Falls back to the per-axis median when the representative is not among
        the rows.
        """
        best = ""
        best_file = target_dir.parent / "best_structure.tsv"
        if best_file.exists():
            with open(best_file) as handle:
                next(handle, None)
                for line in handle:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2 and parts[0] == target_dir.name:
                        best = parts[1].strip().lower()
                        break

        for stem, _, centre in rows:
            if best and stem == best and centre is not None:
                return centre

        # Fallback: the medoid, not a mean or a per-axis median. Both of those
        # are positions no ligand occupies, and an off-site entry pulls them
        # away from the real site; the medoid is always an actual ligand.
        centres = [c for _, _, c in rows if c is not None]
        if not centres:
            return None
        return min(
            centres,
            key=lambda a: sum(math.dist(a, b) for b in centres),
        )

    def _cocrystal_fps(self, target_dir: Path) -> List:
        """Fingerprints of the target's same-pocket cocrystal ligands.

        Every PDB entry of the target contributes, not just the representative:
        each one is a solved structure of something bound in this pocket, and
        all of them are known binders.

        Two things this has to survive. An oligomer repeating one ligand per
        chain, handled by deduplicating on component id so the copies collapse
        to a single reference entry; a component is kept when ANY of its copies
        is in the pocket, which is why the distance test runs before the dedupe
        bookkeeping. And a structure whose ligand sits at a different site,
        handled by the distance guard.

        Note that pocket_info.csv already holds ONE row per (PDB, chain): stage
        2 reduces each structure to its best-contacting ligand before writing
        it (protein_filter.py, `max(n_contacts, n_heavy)`). So a single entry
        with a buffer molecule parked elsewhere has already been resolved, and
        check_same_pocket has already rejected targets whose sites disagree by
        more than MAX_POCKET_SPREAD. The guard here catches only what those two
        would have to stop catching.
        """
        pocket_info = target_dir / "pocket_info.csv"
        if not pocket_info.exists():
            return []

        ccd = self._ccd_smiles()
        rows: List[tuple] = []
        with open(pocket_info) as handle:
            for record in csv.DictReader(handle):
                name = (record.get('Ligand_Name') or '').strip()
                if not name:
                    continue
                stem = (record.get('Aligned_File') or '').strip().lower()
                if stem.endswith('.pdb'):
                    stem = stem[:-4]
                try:
                    centre = (float(record['Center_X']), float(record['Center_Y']),
                              float(record['Center_Z']))
                except (KeyError, TypeError, ValueError):
                    centre = None
                rows.append((stem, name, centre))
        if not rows:
            return []

        centre_ref = self._pocket_centre(target_dir, rows)

        fps: List = []
        seen: set = set()
        off_site: set = set()
        for _stem, name, centre in rows:
            if name in seen:
                continue
            if centre_ref is not None and centre is not None:
                if math.dist(centre, centre_ref) > self.MAX_POCKET_SPREAD:
                    off_site.add(name)
                    continue
            smiles = ccd.get(name)
            if not smiles:
                continue
            props = self._compute_properties(smiles)
            if props is None:
                continue
            seen.add(name)
            fps.append(props[-1])

        off_site -= seen
        if off_site:
            self.logger.warning(
                f"{target_dir.name}: {len(off_site)} cocrystal ligand(s) over "
                f"{self.MAX_POCKET_SPREAD} A from the pocket centre in every "
                f"copy, excluded from the reference set: {sorted(off_site)}"
            )
        return fps

    def _reference_fps(self, target_dir: Path) -> List:
        """The cross-active reference: all actives, cocrystals, weak binders.

        Not the cluster representatives, and not measured inactives. See the
        module docstring for why each is in or out.
        """
        fps: List = []
        for name in ("actives.tsv", "weak.tsv"):
            path = target_dir / name
            if not path.exists():
                continue
            with open(path) as handle:
                header = next(handle, "")
                columns = header.rstrip('\n').split('\t')
                try:
                    smiles_at = columns.index('smiles')
                except ValueError:
                    smiles_at = 2
                for line in handle:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) <= smiles_at:
                        continue
                    props = self._compute_properties(parts[smiles_at])
                    if props is not None:
                        fps.append(props[-1])
        fps.extend(self._cocrystal_fps(target_dir))
        return fps

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
        target_ref_fps: Dict[str, list] = {}
        skipped = 0
        n_ref_fps = 0
        n_no_ref = 0

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

            reference = self._reference_fps(data_dir / uniprot)
            target_ref_fps[uniprot] = reference
            n_ref_fps += len(reference)
            if not reference:
                n_no_ref += 1

            if (i + 1) % 100 == 0 or (i + 1) == len(passed_targets):
                self.logger.info(
                    f"  {i+1}/{len(passed_targets)} targets, pool size: {len(pool)}"
                )

        if skipped:
            self.logger.warning(
                f"{skipped} targets had no actives_clustered.tsv — "
                "run cluster-actives first"
            )

        if n_no_ref:
            self.logger.warning(
                f"{n_no_ref} targets got an empty cross-active reference set; "
                f"the bar cannot fire for them"
            )

        data = {
            'pool': pool,
            'target_actives': target_actives,
            'target_ref_fps': target_ref_fps,
        }
        with open(output, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

        self.logger.info(
            f"Saved to {output}: {len(pool)} unique compounds, "
            f"{len(target_actives)} targets, "
            f"{n_ref_fps} cross-active reference fingerprints "
            f"({n_ref_fps / max(len(target_ref_fps), 1):.1f} per target)"
        )

        # Compound-per-target-count distribution
        multi: Dict[int, int] = defaultdict(int)
        for d in pool.values():
            multi[len(d['target_set'])] += 1
        for n_tgt in sorted(multi):
            self.logger.info(f"  Active in {n_tgt} target(s): {multi[n_tgt]} compounds")

        return output

# chembl_curator/protein_filter.py

import os
import shutil
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
import requests
import numpy as np
from scipy.spatial import cKDTree
from dataclasses import dataclass


def _load_excluded_ligands() -> set:
    """Load excluded ligand codes from assets/excluded_ligands.txt."""
    path = Path(__file__).parent / "assets" / "excluded_ligands.txt"
    codes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                codes.add(line)
    return codes


EXCLUDED_LIGANDS = _load_excluded_ligands()

# Matches CurationConfig.min_heavy_atoms on the ChEMBL side. Without it a
# three-atom fragment can define a pocket, and the excluded-code list only
# catches artifacts it happens to name.
MIN_LIGAND_HEAVY_ATOMS = 5

# Mean count of target-chain protein heavy atoms within BURIAL_SHELL of each
# ligand atom. A ligand lying on the surface is not a pocket, and this also
# removes the copies align_pdb keeps from other chains: it writes every HETATM
# regardless of chain, and a 53-atom FAD grazing the target with one atom used
# to pass the contact test and go on to define the pocket.
BURIAL_SHELL = 8.0
MIN_BURIAL_NEIGHBOURS = 20.0


@dataclass
class PDBInfo:
    pdb_id: str
    method: str
    resolution: str
    chains: str  # e.g. "A=1-250,B=1-250"


@dataclass
class LigandInfo:
    ligand_name: str
    center: np.ndarray  # 3D coordinates
    chain: str
    # Heavy-atom coordinates, kept so closest-approach between two ligands can
    # be measured. Centroid distance alone is contaminated by ligand size.
    coords: Optional[np.ndarray] = None
    n_heavy: int = 0
    burial: float = 0.0      # mean protein heavy atoms within BURIAL_SHELL
    n_contacts: int = 0      # protein atoms within 4.5 A


class ProteinFilter:
    """Filter protein structures based on PDB availability and binding site analysis."""

    def __init__(self, curated_dir: Path, log_level: str = "INFO",
                 max_chain_residues: int = 1500):
        """
        Args:
            curated_dir: Directory containing target subdirectories (uniprot IDs)
            log_level: Logging level
            max_chain_residues: Skip PDB structures whose target chain exceeds
                this many residues (catches ribosomes, nanodiscs, etc.). 0 = no limit.
        """
        self.curated_dir = Path(curated_dir)
        self.max_chain_residues = max_chain_residues
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

        # Setup console handler if not already configured
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            handler.setLevel(getattr(logging, log_level.upper()))
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

    def get_target_list(self) -> List[str]:
        """Get list of all target uniprot IDs from directory."""
        targets = []
        for item in self.curated_dir.iterdir():
            if item.is_dir():
                targets.append(item.name)
        return targets

    def fetch_uniprot_info(self, uniprot_id: str) -> Tuple[List[PDBInfo], Optional[str]]:
        """
        Fetch PDB structures and canonical sequence for a UniProt ID in one API call.

        Returns:
            (pdb_list, canonical_sequence) — sequence is None if unavailable
        """
        try:
            url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}'
            response = requests.get(url, headers={'Accept': 'application/json'}, timeout=120)
            response.raise_for_status()
            data = response.json()

            pdb_list = []
            if 'uniProtKBCrossReferences' in data:
                for xref in data['uniProtKBCrossReferences']:
                    if xref['database'] == 'PDB':
                        pdb_id = xref['id']
                        properties = {p['key']: p['value'] for p in xref.get('properties', [])}
                        pdb_list.append(PDBInfo(
                            pdb_id=pdb_id,
                            method=properties.get('Method', '-'),
                            resolution=properties.get('Resolution', '-'),
                            chains=properties.get('Chains', '-')
                        ))

            sequence = None
            if 'sequence' in data and 'value' in data['sequence']:
                sequence = data['sequence']['value']

            return pdb_list, sequence

        except Exception as e:
            self.logger.error(f"Error fetching UniProt info for {uniprot_id}: {e}")
            return [], None

    def fetch_pdb_list(self, uniprot_id: str) -> List[PDBInfo]:
        """Fetch PDB structures for a given uniprot ID (backward-compatible wrapper)."""
        pdb_list, _ = self.fetch_uniprot_info(uniprot_id)
        return pdb_list

    def save_pdb_list(self, target_dir: Path, pdb_list: List[PDBInfo]):
        """Save PDB list to file."""
        pdb_list_file = target_dir / "pdbid.list"
        with open(pdb_list_file, 'w') as f:
            f.write("# PDBID method resolution chains\n")
            for pdb in pdb_list:
                f.write(f"{pdb.pdb_id}  {pdb.method}  {pdb.resolution}  {pdb.chains}\n")

    def download_pdb(self, pdb_id: str, output_dir: Path) -> bool:
        """
        Download PDB file using pdb_get command (if available) or RCSB web download.

        Args:
            pdb_id: PDB ID to download
            output_dir: Directory to save PDB file

        Returns:
            True if successful, False otherwise
        """
        pdb_file = output_dir / f"{pdb_id.lower()}.pdb"

        # Method 1: Try pdb_get (fast if available on cluster)
        if shutil.which('pdb_get'):
            try:
                result = subprocess.run(
                    ['pdb_get', pdb_id],
                    cwd=output_dir,
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                # Check if PDB file was created
                if pdb_file.exists():
                    self.logger.debug(f"Downloaded {pdb_id} using pdb_get")
                    return True

                # pdb_get might create files with different naming
                pdb_files = list(output_dir.glob(f"*{pdb_id.lower()}*"))
                if pdb_files:
                    # Rename to standard format
                    pdb_files[0].rename(pdb_file)
                    self.logger.debug(f"Downloaded {pdb_id} using pdb_get")
                    return True

            except Exception as e:
                self.logger.debug(f"pdb_get failed for {pdb_id}: {e}, trying web download")

        # Method 2: Download from RCSB PDB (fallback or if pdb_get not available)
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            self.logger.debug(f"Downloading {pdb_id} from RCSB: {url}")

            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(pdb_file, 'w') as f:
                f.write(response.text)

            if pdb_file.exists() and pdb_file.stat().st_size > 0:
                self.logger.debug(f"Downloaded {pdb_id} from RCSB")
                return True

            return False

        except Exception as e:
            self.logger.error(f"Error downloading PDB {pdb_id}: {e}")
            return False

    def download_alphafold(self, uniprot_id: str, output_dir: Path) -> bool:
        """
        Download AlphaFold structure from EBI.

        Args:
            uniprot_id: UniProt accession ID
            output_dir: Directory to save AF model

        Returns:
            True if successful, False otherwise
        """
        try:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.pdb"
            output_file = output_dir / f"AF-{uniprot_id}.pdb"

            result = subprocess.run(
                ['wget', '-O', str(output_file), url],
                capture_output=True,
                text=True,
                timeout=120
            )

            if output_file.exists() and output_file.stat().st_size > 0:
                return True

            return False

        except Exception as e:
            self.logger.error(f"Error downloading AlphaFold for {uniprot_id}: {e}")
            return False

    def parse_chain_from_pdb_info(self, pdb_info: PDBInfo) -> List[str]:
        """
        Extract chain IDs from PDB info chains string.

        Args:
            pdb_info: PDBInfo object

        Returns:
            List of chain IDs (e.g., ['A', 'B'])
        """
        chains = []
        if pdb_info.chains and pdb_info.chains != '-':
            # Format: "A=1-250,B=1-250" or "A/B=1-250"
            for part in pdb_info.chains.split(','):
                chain_part = part.split('=')[0].strip()
                # Handle multiple chains like "A/B"
                for ch in chain_part.split('/'):
                    if ch.strip():
                        chains.append(ch.strip())
        return chains

    @staticmethod
    def _count_chain_residues(pdb_file: Path, chains: List[str]) -> int:
        """Count Cα atoms in *chains* as a proxy for residue count."""
        if not chains:
            return 0
        chain_set = set(chains)
        count = 0
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    if line[21] in chain_set:
                        count += 1
        return count

    def get_ligands_from_pdb(self, pdb_file: Path, target_chains: List[str]) -> List[LigandInfo]:
        """
        Extract biologically relevant ligands from PDB file.

        Args:
            pdb_file: Path to PDB file
            target_chains: List of target chain IDs

        Returns:
            List of LigandInfo objects with ligand details
        """
        ligands = {}  # Key: (ligand_name, chain), Value: list of coords
        protein_coords = {}  # Key: chain, Value: list of heavy atom coords
        protein_resid = {}   # Key: chain, Value: parallel list of (chain, resseq)
        # One point per residue for residue-level contact counting: CB, or CA
        # for glycine. Counting atoms instead would score a large residue as
        # several contacts and bias the representative-ligand choice.
        res_rep = {}         # Key: chain, Value: {(chain, resseq): coord}

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21:22].strip()
                        if chain in target_chains:
                            # Get heavy atom coordinates (not hydrogen)
                            atom_name = line[12:16].strip()
                            if not atom_name.startswith('H'):
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])

                                if chain not in protein_coords:
                                    protein_coords[chain] = []
                                    protein_resid[chain] = []
                                protein_coords[chain].append(np.array([x, y, z]))
                                protein_resid[chain].append((chain, line[22:27]))
                                if atom_name in ('CB', 'CA'):
                                    rk = (chain, line[22:27])
                                    d = res_rep.setdefault(chain, {})
                                    if atom_name == 'CB' or rk not in d:
                                        d[rk] = np.array([x, y, z])

                    elif line.startswith('HETATM'):
                        ligand_name = line[17:20].strip()
                        chain = line[21:22].strip()

                        # Skip excluded ligands
                        if ligand_name in EXCLUDED_LIGANDS:
                            continue

                        # Get heavy atom coordinates
                        atom_name = line[12:16].strip()
                        if not atom_name.startswith('H'):
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])

                            key = (ligand_name, chain)
                            if key not in ligands:
                                ligands[key] = []
                            ligands[key].append(np.array([x, y, z]))

        except Exception as e:
            self.logger.error(f"Error parsing PDB file {pdb_file}: {e}")
            return []

        # Build KD-trees for fast contact detection (within 4A)
        prot_trees = {}
        for chain, coords in protein_coords.items():
            prot_trees[chain] = cKDTree(np.array(coords))

        ligand_infos = []
        for (ligand_name, lig_chain), lig_coords in ligands.items():
            lig_arr = np.array(lig_coords)
            has_contact = False

            for target_chain in target_chains:
                if target_chain not in prot_trees:
                    continue
                # Query: any ligand atom within 4A of any protein atom?
                dists, _ = prot_trees[target_chain].query(lig_arr, distance_upper_bound=4.0)
                if np.any(np.isfinite(dists)):
                    has_contact = True
                    break

            if len(lig_arr) < MIN_LIGAND_HEAVY_ATOMS or not has_contact:
                continue

            burial, contacts = self._burial(
                lig_arr, prot_trees, target_chains, res_rep)
            if burial < MIN_BURIAL_NEIGHBOURS:
                continue

            ligand_infos.append(LigandInfo(
                ligand_name=ligand_name,
                center=np.mean(lig_arr, axis=0),
                chain=lig_chain,
                coords=lig_arr,
                n_heavy=len(lig_arr),
                burial=burial,
                n_contacts=contacts,
            ))

        return ligand_infos

    def get_chain_ligand_contacts(self, pdb_file: Path, target_chains: List[str]) -> Dict[str, List[LigandInfo]]:
        """
        Determine which chains contact which ligands.

        Args:
            pdb_file: Path to PDB file
            target_chains: List of target chain IDs

        Returns:
            Dictionary mapping chain ID to list of ligands that contact it
        """
        ligands = {}  # Key: (ligand_name, chain), Value: list of coords
        protein_coords = {}  # Key: chain, Value: list of heavy atom coords
        protein_resid = {}   # Key: chain, Value: parallel list of (chain, resseq)
        # One point per residue for residue-level contact counting: CB, or CA
        # for glycine. Counting atoms instead would score a large residue as
        # several contacts and bias the representative-ligand choice.
        res_rep = {}         # Key: chain, Value: {(chain, resseq): coord}

        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21:22].strip()
                        if chain in target_chains:
                            # Get heavy atom coordinates (not hydrogen)
                            atom_name = line[12:16].strip()
                            if not atom_name.startswith('H'):
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])

                                if chain not in protein_coords:
                                    protein_coords[chain] = []
                                    protein_resid[chain] = []
                                protein_coords[chain].append(np.array([x, y, z]))
                                protein_resid[chain].append((chain, line[22:27]))
                                if atom_name in ('CB', 'CA'):
                                    rk = (chain, line[22:27])
                                    d = res_rep.setdefault(chain, {})
                                    if atom_name == 'CB' or rk not in d:
                                        d[rk] = np.array([x, y, z])

                    elif line.startswith('HETATM'):
                        ligand_name = line[17:20].strip()
                        chain = line[21:22].strip()

                        # Skip excluded ligands
                        if ligand_name in EXCLUDED_LIGANDS:
                            continue

                        # Get heavy atom coordinates
                        atom_name = line[12:16].strip()
                        if not atom_name.startswith('H'):
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])

                            key = (ligand_name, chain)
                            if key not in ligands:
                                ligands[key] = []
                            ligands[key].append(np.array([x, y, z]))

        except Exception as e:
            self.logger.error(f"Error parsing PDB file {pdb_file}: {e}")
            return {}

        # Build KD-trees for fast contact detection (within 4A)
        prot_trees = {}
        for chain, coords in protein_coords.items():
            prot_trees[chain] = cKDTree(np.array(coords))

        chain_ligands = {chain: [] for chain in target_chains}

        for (ligand_name, lig_chain), lig_coords in ligands.items():
            lig_arr = np.array(lig_coords)
            for target_chain in target_chains:
                if target_chain not in prot_trees:
                    continue

                # Query: any ligand atom within 4A of any protein atom?
                dists, _ = prot_trees[target_chain].query(lig_arr, distance_upper_bound=4.0)
                if len(lig_arr) < MIN_LIGAND_HEAVY_ATOMS:
                    continue
                if np.any(np.isfinite(dists)):
                    burial, contacts = self._burial(
                        lig_arr, prot_trees, [target_chain], res_rep)
                    if burial < MIN_BURIAL_NEIGHBOURS:
                        continue
                    chain_ligands[target_chain].append(LigandInfo(
                        ligand_name=ligand_name,
                        center=np.mean(lig_arr, axis=0),
                        chain=lig_chain,
                        coords=lig_arr,
                        n_heavy=len(lig_arr),
                        burial=burial,
                        n_contacts=contacts,
                    ))

        # Remove chains with no ligand contacts
        chain_ligands = {chain: ligs for chain, ligs in chain_ligands.items() if ligs}

        return chain_ligands

    @staticmethod
    def _same_pocket(a: LigandInfo, b: LigandInfo, max_gap: float) -> bool:
        """Are two ligands close enough to be occupying one site?

        Uses closest atom-atom approach, not centroid separation. Centroid
        distance scales with ligand size: in 1T2F the NAD centroid sits 9.96 A
        from an oxamate whose closest atom is 2.85 A away, stacked against the
        nicotinamide. Three probe cases (NAD/OXQ, FAD/ML2, 48V/PTR) all had
        centroid distances of 9.8-10.0 A yet closest approaches of 2.9-3.8 A,
        while a genuine two-site case (P06700, NCA/XYQ) had a near-identical
        centroid distance of 9.74 A and a closest approach of 7.29 A. Centroid
        distance cannot tell those apart; closest approach can.

        Falls back to centroid distance when coordinates are unavailable.
        """
        if a.coords is None or b.coords is None:
            return bool(np.linalg.norm(a.center - b.center) <= 10.0)
        gap = cKDTree(a.coords).query(b.coords)[0].min()
        return bool(gap <= max_gap)


    @staticmethod
    def _burial(lig_arr, prot_trees, target_chains, res_rep,
                shell: float = BURIAL_SHELL,
                res_cut: float = 8.0):
        """(burial density, number of contacting residues) for one ligand.

        burial is the mean number of protein heavy atoms within `shell` of each
        ligand atom. A ligand in a cleft is surrounded; one on the surface, or
        a copy align_pdb carried over from another chain, is not.

        The residue count uses one point per residue (CB, or CA for glycine)
        rather than every heavy atom, so a tryptophan does not outweigh an
        alanine simply by having more atoms.
        """
        counts, residues = [], set()
        for chain in target_chains:
            tree = prot_trees.get(chain)
            if tree is None:
                continue
            counts.append(np.array(
                [len(x) for x in tree.query_ball_point(lig_arr, r=shell)],
                dtype=float))
            reps = res_rep.get(chain) or {}
            if reps:
                keys = list(reps)
                pts = np.stack([reps[k] for k in keys])
                near = cKDTree(lig_arr).query(pts)[0] <= res_cut
                residues.update(k for k, ok in zip(keys, near) if ok)
        if not counts:
            return 0.0, 0
        return float(np.sum(counts, axis=0).mean()), len(residues)

    def is_single_ligand_bound(self, ligands: List[LigandInfo],
                               max_gap: float = 5.0) -> bool:
        """
        Check if a structure holds a single ligand or one shared pocket.

        Args:
            ligands: List of LigandInfo objects
            max_gap: Maximum closest atom-atom approach (A) for two ligands to
                count as occupying the same pocket. 5.0 allows van der Waals
                contact plus slack; a cofactor and its substrate analogue sit
                at 2.9-3.8 A, separate sites at 7 A and beyond.

        Returns:
            True if single ligand or all ligands share a pocket, False otherwise
        """
        if len(ligands) == 0:
            return False

        if len(ligands) == 1:
            return True

        for i in range(len(ligands)):
            for j in range(i + 1, len(ligands)):
                if not self._same_pocket(ligands[i], ligands[j], max_gap):
                    return False

        return True

    @staticmethod
    def _extract_ca_chain(pdb_path: Path, chain_id: str = "") -> np.ndarray:
        """Extract Cα coordinates from a PDB file, optionally for a specific chain.

        Args:
            pdb_path: Path to PDB file
            chain_id: If non-empty, only extract Cα from this chain.
                      If empty, extract all Cα.

        Returns:
            (N, 3) numpy array of Cα coordinates.
        """
        import nuri

        mol = list(nuri.readfile('pdb', str(pdb_path), sanitize=False))[0]
        ca = []
        for sub in mol.subs:
            if chain_id and sub.props.get('chain', '') != chain_id:
                continue
            for atom in sub:
                if atom.name.strip() == 'CA' and atom.element_symbol == 'C':
                    ca.append(atom.get_pos(0))
        return np.array(ca) if ca else np.empty((0, 3))

    def align_pdb(self, query_pdb: Path, template_pdb: Path, output_pdb: Path, target_chain: str):
        """
        Align PDB structure to template using nuri TMAlign and save a single
        target protein chain and all ligands.

        Args:
            query_pdb: Query PDB file to align
            template_pdb: Template PDB file (reference, e.g. AlphaFold model)
            output_pdb: Output aligned PDB file
            target_chain: Single chain ID to save for protein atoms (all HETATMs are saved)
        """
        try:
            import nuri._log_interface
            nuri._log_interface.set_log_level(4)  # suppress C++ warnings
            from nuri.tools.tm import TMAlign

            # Extract Cα: target chain from query, all from template
            ca_query = self._extract_ca_chain(query_pdb, chain_id=target_chain)
            ca_template = self._extract_ca_chain(template_pdb)

            if len(ca_query) < 5 or len(ca_template) < 5:
                self.logger.warning(
                    f"Too few Cα atoms for alignment: query={len(ca_query)}, "
                    f"template={len(ca_template)}"
                )
                return

            tma = TMAlign(ca_query, ca_template)
            xform, _ = tma.score()       # xform is a 4x4 affine matrix
            R = xform[:3, :3]
            t = xform[:3, 3]

            # Apply transformation and save target chain (ATOM) and all ligands (HETATM)
            with open(query_pdb, 'r') as fin, open(output_pdb, 'w') as fout:
                for line in fin:
                    if line.startswith("ATOM"):
                        chain = line[21:22].strip()
                        if chain != target_chain:
                            continue
                    elif line.startswith("HETATM"):
                        pass
                    else:
                        continue

                    try:
                        coord = np.array([
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ])
                    except ValueError:
                        fout.write(line)
                        continue

                    new_coord = R @ coord + t
                    new_line = (
                        line[:30]
                        + f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}"
                        + line[54:]
                    )
                    fout.write(new_line)

        except Exception as e:
            self.logger.error(f"Error aligning PDB {query_pdb}: {e}")

    def check_same_pocket(self, ligand_centers: List[np.ndarray], distance_threshold: float = 10.0) -> bool:
        """
        Check if all ligands are in the same pocket.

        Args:
            ligand_centers: List of ligand center coordinates
            distance_threshold: Maximum distance to be considered same pocket

        Returns:
            True if all in same pocket, False otherwise
        """
        if len(ligand_centers) <= 1:
            return True

        for i in range(len(ligand_centers)):
            for j in range(i + 1, len(ligand_centers)):
                dist = np.linalg.norm(ligand_centers[i] - ligand_centers[j])
                if dist > distance_threshold:
                    return False

        return True

    def process_target(self, uniprot_id: str) -> bool:
        """
        Process a single target through the entire protein filtering pipeline.

        Steps:
        1. Fetch PDB list
        2. Download PDBs
        3. Download AlphaFold model
        4. Detect ligand-bound structures
        5. Align structures
        6. Check pocket clustering

        Args:
            uniprot_id: Target uniprot ID

        Returns:
            True if target passes all filters (single binding site), False otherwise
        """
        target_dir = self.curated_dir / uniprot_id

        self.logger.info(f"Processing target {uniprot_id}")

        # Step 1: Fetch PDB list and canonical sequence
        pdb_list, sequence = self.fetch_uniprot_info(uniprot_id)
        if not pdb_list:
            self.logger.warning(f"No PDB structures found for {uniprot_id}")
            return False

        self.save_pdb_list(target_dir, pdb_list)
        self.logger.info(f"Found {len(pdb_list)} PDB structures for {uniprot_id}")

        # Step 2 & 3: Create pdb directory and download structures
        pdb_dir = target_dir / "pdb"
        pdb_dir.mkdir(exist_ok=True)

        # Download AlphaFold model
        af_downloaded = self.download_alphafold(uniprot_id, pdb_dir)
        if not af_downloaded:
            self.logger.warning(f"No AlphaFold model available for {uniprot_id}")
            return False

        af_model = pdb_dir / f"AF-{uniprot_id}.pdb"

        # Download PDBs, skip structures with oversized target chains
        downloaded_pdbs = []
        for pdb_info in pdb_list:
            if self.download_pdb(pdb_info.pdb_id, pdb_dir):
                if self.max_chain_residues > 0:
                    pdb_file = pdb_dir / f"{pdb_info.pdb_id.lower()}.pdb"
                    chains = self.parse_chain_from_pdb_info(pdb_info)
                    n_res = self._count_chain_residues(pdb_file, chains)
                    if n_res > self.max_chain_residues:
                        self.logger.warning(
                            f"Skipping PDB {pdb_info.pdb_id} "
                            f"({n_res} residues in target chain, "
                            f"limit {self.max_chain_residues})"
                        )
                        pdb_file.unlink()
                        continue
                downloaded_pdbs.append(pdb_info)

        if not downloaded_pdbs:
            self.logger.warning(f"No PDB files downloaded for {uniprot_id}")
            return False

        self.logger.info(f"Downloaded {len(downloaded_pdbs)} PDB files")

        # Step 4: Detect chain-ligand contacts for each structure
        # Structure: list of (pdb_info, chain_ligands_dict)
        # where chain_ligands_dict: {chain_id: [LigandInfo, ...]}
        chain_bound_pdbs = []
        for pdb_info in downloaded_pdbs:
            pdb_file = pdb_dir / f"{pdb_info.pdb_id.lower()}.pdb"
            if not pdb_file.exists():
                continue

            # Get target chains
            target_chains = self.parse_chain_from_pdb_info(pdb_info)
            if not target_chains:
                continue

            # Get chain-ligand contacts
            chain_ligands = self.get_chain_ligand_contacts(pdb_file, target_chains)

            # Only keep structures where at least one chain has ligand contacts
            # and each chain with contacts has single/clustered ligands
            valid_chains = {}
            for chain, ligands in chain_ligands.items():
                if self.is_single_ligand_bound(ligands):
                    valid_chains[chain] = ligands

            if valid_chains:
                chain_bound_pdbs.append((pdb_info, valid_chains))

        if not chain_bound_pdbs:
            self.logger.info(f"No chain-ligand bound structures for {uniprot_id}")
            return False

        # Save ligand-bound PDB list
        ligand_bound_file = target_dir / "ligand_bound_pdbs.txt"
        with open(ligand_bound_file, 'w') as f:
            for pdb_info, chain_ligands in chain_bound_pdbs:
                for chain, ligands in chain_ligands.items():
                    ligand_names = ','.join([lig.ligand_name for lig in ligands])
                    f.write(f"{pdb_info.pdb_id}  Chain {chain}: {ligand_names}\n")

        total_chain_structures = sum(len(chain_ligs) for _, chain_ligs in chain_bound_pdbs)
        self.logger.info(f"Found {total_chain_structures} chain-ligand pairs across {len(chain_bound_pdbs)} PDB structures")

        # Step 5: Align structures - create separate file for each chain
        aligned_dir = target_dir / "aligned"
        aligned_dir.mkdir(exist_ok=True)

        # Structure: list of (pdb_info, chain, ligands, aligned_pdb_path)
        aligned_structures = []
        for pdb_info, chain_ligands in chain_bound_pdbs:
            pdb_file = pdb_dir / f"{pdb_info.pdb_id.lower()}.pdb"

            # Align each chain separately
            for chain, ligands in chain_ligands.items():
                aligned_pdb = aligned_dir / f"{pdb_info.pdb_id.lower()}_{chain}.pdb"
                self.align_pdb(pdb_file, af_model, aligned_pdb, chain)

                if aligned_pdb.exists():
                    aligned_structures.append((pdb_info, chain, ligands, aligned_pdb))

        if not aligned_structures:
            self.logger.warning(f"No structures aligned for {uniprot_id}")
            return False

        self.logger.info(f"Aligned {len(aligned_structures)} chain-specific structures")

        # Step 6: Check pocket clustering
        # Re-extract ligand centers from aligned structures
        all_ligand_centers = []
        pocket_info = []

        for pdb_info, chain, original_ligands, aligned_pdb in aligned_structures:
            # Extract ligands from aligned structure
            # Use only the specific chain since the file only contains that chain
            aligned_ligands = self.get_ligands_from_pdb(aligned_pdb, [chain])

            if aligned_ligands:
                # One representative per structure. Downstream code takes the
                # pocket centre and the chemical environment from this ligand,
                # so the choice matters: in 1T2F both NAD and its neighbouring
                # oxamate are present, and NAD is the one that describes the
                # site. Ranking by contacting residues rather than atom count
                # also avoids picking a large ligand that barely touches this
                # chain, which is how a 53-atom FAD with a single contact used
                # to end up defining a pocket.
                rep = max(aligned_ligands,
                          key=lambda l: (l.n_contacts, l.n_heavy))
                all_ligand_centers.append(rep.center)
                pocket_info.append({
                    'pdb_id': pdb_info.pdb_id,
                    'chain': chain,
                    'aligned_file': aligned_pdb.name,
                    'ligand_name': rep.ligand_name,
                    'center': rep.center
                })

        # Check if all ligands are in same pocket
        is_single_pocket = self.check_same_pocket(all_ligand_centers)

        # Save pocket info
        import csv
        pocket_csv = target_dir / "pocket_info.csv"
        with open(pocket_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['PDB_ID', 'Chain', 'Aligned_File', 'Ligand_Name', 'Center_X', 'Center_Y', 'Center_Z'])
            for info in pocket_info:
                writer.writerow([
                    info['pdb_id'],
                    info['chain'],
                    info['aligned_file'],
                    info['ligand_name'],
                    f"{info['center'][0]:.3f}",
                    f"{info['center'][1]:.3f}",
                    f"{info['center'][2]:.3f}"
                ])

        if not is_single_pocket:
            self.logger.info(f"Multiple binding pockets detected for {uniprot_id}")
            return False

        # Write per-target sequence FASTA (combined into sequences.fasta at the end)
        if sequence:
            seq_fasta = target_dir / "sequence.fasta"
            with open(seq_fasta, 'w') as f:
                f.write(f">{uniprot_id}\n{sequence}\n")

        self.logger.info(f"Target {uniprot_id} passed all filters (single binding site)")
        return True

    def _get_best_structure(self, target_dir: Path) -> Optional[Tuple[str, float]]:
        """
        Return (pdbid_chain, resolution) for the best (lowest-resolution) structure.

        Reads pocket_info.csv for available aligned structures and pdbid.list for resolutions.
        pdbid_chain format matches the aligned filename stem, e.g. "3mms_A".
        """
        import csv as _csv
        pocket_csv = target_dir / "pocket_info.csv"
        if not pocket_csv.exists():
            return None

        pdbid_chains = set()
        with open(pocket_csv) as f:
            for row in _csv.DictReader(f):
                pdbid_chains.add(row['Aligned_File'].replace('.pdb', ''))

        if not pdbid_chains:
            return None

        # Parse resolutions from pdbid.list
        res_map: Dict[str, float] = {}
        pdbid_list = target_dir / "pdbid.list"
        if pdbid_list.exists():
            with open(pdbid_list) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        pdbid = parts[0].upper()
                        try:
                            res_map[pdbid] = float(parts[2].rstrip('A'))
                        except ValueError:
                            res_map[pdbid] = 99.0

        # inf, not 99.0: NMR entries report resolution as "-", fall back to 99.0,
        # and would then fail a "< 99.0" test, silently leaving NMR-only targets
        # out of best_structure.tsv even though they passed every filter.
        best_chain, best_res = None, float('inf')
        for pc in sorted(pdbid_chains):
            pdbid = pc.split('_')[0].upper()
            res = res_map.get(pdbid, 99.0)
            if res < best_res:
                best_res = res
                best_chain = pc

        return (best_chain, best_res) if best_chain else None

    def run_pipeline(self, n_processes: int = 1) -> List[str]:
        """
        Run protein filtering pipeline on all targets.

        Args:
            n_processes: Number of parallel processes (default: 1 for sequential)

        Returns:
            List of target IDs that passed filters
        """
        targets = self.get_target_list()
        self.logger.info(f"Processing {len(targets)} targets")

        passed_targets = []

        if n_processes == 1:
            # Sequential processing
            for target in targets:
                try:
                    if self.process_target(target):
                        passed_targets.append(target)
                except Exception as e:
                    self.logger.error(f"Error processing {target}: {e}")
        else:
            # Parallel processing
            import multiprocessing as mp
            with mp.Pool(n_processes) as pool:
                results = pool.map(self._process_target_wrapper, targets)
                passed_targets = [t for t, passed in zip(targets, results) if passed]

        # Save passed targets
        passed_file = self.curated_dir / "passed_targets.txt"
        with open(passed_file, 'w') as f:
            for target in passed_targets:
                f.write(f"{target}\n")

        self.logger.info(f"Pipeline complete: {len(passed_targets)}/{len(targets)} targets passed")

        # Combine per-target sequence.fasta into global sequences.fasta
        # and write best_structure.tsv for pocket RMSD computation
        self._write_global_outputs(passed_targets)

        # Optionally delete filtered-out targets
        for target in targets:
            if target not in passed_targets:
                import shutil
                target_dir = self.curated_dir / target
                if target_dir.exists():
                    shutil.rmtree(target_dir)

        return passed_targets

    def _write_global_outputs(self, passed_targets: List[str]) -> None:
        """Write sequences.fasta and best_structure.tsv for all passed targets."""
        seq_fasta = self.curated_dir / "sequences.fasta"
        best_str_tsv = self.curated_dir / "best_structure.tsv"

        n_seq = 0
        n_best = 0
        with open(seq_fasta, 'w') as sf, open(best_str_tsv, 'w') as bf:
            bf.write("uniprot\tpdbid_chain\tresolution\n")
            for target in passed_targets:
                target_dir = self.curated_dir / target

                # Append sequence
                per_seq = target_dir / "sequence.fasta"
                if per_seq.exists():
                    sf.write(per_seq.read_text())
                    n_seq += 1

                # Append best structure entry
                best = self._get_best_structure(target_dir)
                if best:
                    bf.write(f"{target}\t{best[0]}\t{best[1]:.2f}\n")
                    n_best += 1

        self.logger.info(f"Wrote sequences.fasta ({n_seq} sequences) and "
                         f"best_structure.tsv ({n_best} entries)")

    def _process_target_wrapper(self, target: str) -> bool:
        """Wrapper for multiprocessing."""
        try:
            return self.process_target(target)
        except Exception as e:
            self.logger.error(f"Error processing {target}: {e}")
            return False

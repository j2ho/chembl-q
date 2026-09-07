# chembl_curator/external_pockets.py

"""Stage 7: pocket-level leakage against PDBbind and BioLiP.

Stage 8 separates the test set from the external datasets by sequence, which
is not the same guarantee. Two proteins can share a binding site while sitting
far apart in sequence space, and on the ChEMBL side the two signals barely
overlap: of 3,073 pocket-similar target pairs only 107 are also caught at
0.4 sequence identity. A model trained on PDBbind or BioLiP has seen those
pockets whether or not the sequences align.

Two steps:

1. build   Extract the pocket of every external entry, the same way stage 5
           builds a ChEMBL pocket: residues with a heavy atom within
           pocket_radius of a ligand heavy atom, kept as (residue name, CA).
           Cached in one npz so the ~39k structures are read once.

2. compare Order-free Hungarian superposition of every ChEMBL pocket against
           every external pocket. Writes the pairs that come close, and a
           per-target best match, so a split can demote test targets without
           recomputing anything.

The entry list comes from the same external FASTA stage 8 blocks against, so
the two stages cannot disagree about what "external" means.
"""

import logging
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .pocket_align import (aligned_rmsd_prepared, pocket_from_structure,
                           prepare_pocket)
from .receptor_similarity import parse_ligand_atoms, parse_protein_residues

# Residue identity is carried as a byte code so the cache stays small; the
# order here is the on-disk contract and must not be reshuffled.
RESIDUE_CODES = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",
)
_CODE_OF = {name: i for i, name in enumerate(RESIDUE_CODES)}


def resolve_external_entry(
    source: str,
    entry_id: str,
    biolip_dir: Optional[Path],
    pdbbind_dirs: Sequence[Path],
) -> Optional[Tuple[Path, Path]]:
    """Locate (receptor, ligand) for one external entry, or None if missing.

    BioLiP ids are {pdb}_{lig}_{chain}_{serial}, and the receptor is the
    single chain the site was annotated on. PDBbind ids are the bare PDB code
    and the entry may live in either the refined or the general set.
    """
    if source == "biolip":
        if biolip_dir is None:
            return None
        parts = entry_id.split("_")
        if len(parts) != 4:
            return None
        pdb_id, _lig, chain, _serial = parts
        receptor = biolip_dir / "receptor" / f"{pdb_id}{chain}.pdb"
        ligand = biolip_dir / "ligand" / f"{entry_id}.pdb"
        return (receptor, ligand) if receptor.exists() and ligand.exists() else None

    if source == "pdbbind":
        # The general set ships only the renamed ligand. Where both exist the
        # coordinates are identical (checked on the refined set), the rename
        # only makes atom names unique, so either is fine as a source of
        # heavy-atom positions.
        for root in pdbbind_dirs:
            entry = root / entry_id
            receptor = entry / f"{entry_id}_protein.pdb"
            if not receptor.exists():
                continue
            for name in (f"{entry_id}_ligand.pdb",
                         f"{entry_id}_ligand_renamed.pdb"):
                ligand = entry / name
                if ligand.exists():
                    return receptor, ligand
        return None

    return None


def _extract_worker(args):
    """Pocket of one external entry, as (key, residue codes, CA coords)."""
    key, receptor, ligand, radius, min_residues = args
    try:
        residues = parse_protein_residues(receptor)
        ligand_xyz = parse_ligand_atoms(ligand)
        pocket = pocket_from_structure(residues, ligand_xyz, radius=radius,
                                       min_residues=min_residues)
    except (OSError, ValueError):
        return key, None, None
    if pocket is None:
        return key, None, None

    codes = np.array([_CODE_OF.get(name, 255) for name, _ in pocket],
                     dtype=np.uint8)
    coords = np.stack([ca for _, ca in pocket]).astype(np.float32)
    return key, codes, coords


_EXTERNAL: Dict[str, List] = {}


def _init_compare(external):
    _EXTERNAL["pockets"] = external


def _compare_worker(args):
    """One ChEMBL pocket against every external pocket.

    Returns the hits under the reporting cutoff plus the single best match,
    so a later threshold change does not mean recomputing 50M superpositions.
    """
    target, prepared, rmsd_report, min_matched = args
    external = _EXTERNAL["pockets"]
    n_pocket = len(prepared[1])

    hits = []
    best = None
    for key, ext in external:
        try:
            n_matched, rmsd = aligned_rmsd_prepared(prepared, ext)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if n_matched < min_matched:
            continue
        n_ext = len(ext[1])
        if best is None or rmsd < best[1]:
            best = (key, rmsd, n_matched, n_ext)
        if rmsd <= rmsd_report:
            hits.append((key, rmsd, n_matched, n_ext))
    return target, n_pocket, hits, best


class ExternalPocketLeakage:
    """Build external pockets, then score ChEMBL pockets against them."""

    def __init__(self, log_level: str = "INFO"):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, log_level.upper()))

    # ── build ─────────────────────────────────────────────────────────────

    def build(
        self,
        external_fasta: Path,
        output: Path,
        biolip_dir: Optional[Path] = None,
        pdbbind_dirs: Sequence[Path] = (),
        pocket_radius: float = 8.0,
        min_residues: int = 6,
        workers: int = 8,
    ) -> Path:
        """Extract and cache the pocket of every entry in external_fasta."""
        entries = []
        with open(external_fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                sid = line[1:].strip().split()[0]
                source, _, entry_id = sid.partition(".")
                entries.append((source, entry_id))

        self.logger.info(f"External entries listed in FASTA: {len(entries)}")

        jobs = []
        unresolved = 0
        for source, entry_id in entries:
            paths = resolve_external_entry(source, entry_id, biolip_dir,
                                           pdbbind_dirs)
            if paths is None:
                unresolved += 1
                continue
            jobs.append((f"{source}.{entry_id}", paths[0], paths[1],
                         pocket_radius, min_residues))

        if unresolved:
            self.logger.warning(
                f"{unresolved} entries have no structure files on disk and are "
                "not represented in the cache; leakage against them is unchecked"
            )
        self.logger.info(f"Extracting pockets for {len(jobs)} entries "
                         f"({workers} workers)")

        t0 = time.time()
        keys: List[str] = []
        code_blocks: List[np.ndarray] = []
        coord_blocks: List[np.ndarray] = []
        too_small = 0

        with Pool(workers) as pool:
            for i, (key, codes, coords) in enumerate(
                pool.imap_unordered(_extract_worker, jobs, chunksize=64), 1
            ):
                if codes is None:
                    too_small += 1
                else:
                    keys.append(key)
                    code_blocks.append(codes)
                    coord_blocks.append(coords)
                if i % 5000 == 0:
                    self.logger.info(f"  {i}/{len(jobs)} elapsed "
                                     f"{time.time()-t0:.0f}s")

        self.logger.info(
            f"Pockets built: {len(keys)}; {too_small} entries yielded no pocket "
            f"of at least {min_residues} residues; {time.time()-t0:.0f}s"
        )

        lengths = np.array([len(c) for c in code_blocks], dtype=np.int32)
        output = Path(output)
        np.savez_compressed(
            output,
            keys=np.array(keys),
            lengths=lengths,
            codes=np.concatenate(code_blocks) if code_blocks else np.empty(0, np.uint8),
            coords=(np.concatenate(coord_blocks) if coord_blocks
                    else np.empty((0, 3), np.float32)),
            pocket_radius=np.array(pocket_radius),
        )
        self.logger.info(f"Wrote {output} "
                         f"({output.stat().st_size / 1e6:.1f} MB)")
        return output

    # ── compare ───────────────────────────────────────────────────────────

    def compare(
        self,
        data_dir: Path,
        cache: Path,
        output: Optional[Path] = None,
        best_output: Optional[Path] = None,
        pocket_radius: float = 8.0,
        rmsd_report: float = 2.0,
        min_matched: int = 15,
        workers: int = 8,
    ) -> Tuple[Path, Path]:
        """Score every ChEMBL pocket against every external pocket.

        Writes two files: the hits under rmsd_report, and the single closest
        external pocket per target. Both are shipped with the dataset. They
        are reported, not used to filter the split: cutting the test set on
        pocket similarity removes the data-rich targets, because a fold that
        has been drugged hard is also one PDBbind and BioLiP hold many
        structures of.

        The per-target best match is what a reader needs to stratify, and it
        is carried into chembl_targets.tsv. The hits file is the detail
        behind it.
        """
        from .receptor_similarity import ReceptorSimilarity, chembl_pocket

        data_dir = Path(data_dir)
        output = Path(output or data_dir / "external_pocket_hits.tsv")
        best_output = Path(best_output or data_dir / "external_pocket_best.tsv")

        external = self.load(cache)
        self.logger.info(f"External pockets loaded: {len(external)}")

        sim = ReceptorSimilarity(log_level="ERROR")
        best_structure = data_dir / "best_structure.tsv"
        if not best_structure.exists():
            raise FileNotFoundError(f"best_structure.tsv not found: {best_structure}")

        jobs = []
        skipped = 0
        with open(best_structure) as fh:
            next(fh, None)
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                uniprot, pdbid_chain = parts[0], parts[1]
                pdb_path = data_dir / uniprot / "aligned" / f"{pdbid_chain}.pdb"
                lig_name = sim._get_lig_name(data_dir / uniprot, pdbid_chain)
                if not pdb_path.exists() or lig_name is None:
                    skipped += 1
                    continue
                pocket = chembl_pocket(pdb_path, lig_name, pocket_radius)
                if pocket is None:
                    skipped += 1
                    continue
                jobs.append((uniprot, prepare_pocket(pocket),
                             rmsd_report, min_matched))

        self.logger.info(
            f"ChEMBL pockets: {len(jobs)} ({skipped} skipped). "
            f"{len(jobs) * len(external):,} superpositions to run "
            f"({workers} workers)"
        )

        t0 = time.time()
        rows: List[Tuple] = []
        best_rows: List[Tuple] = []
        with Pool(workers, initializer=_init_compare, initargs=(external,)) as pool:
            for i, (target, n_pocket, hits, best) in enumerate(
                pool.imap_unordered(_compare_worker, jobs, chunksize=1), 1
            ):
                for key, rmsd, n_matched, n_ext in hits:
                    rows.append((target, key, rmsd, n_matched, n_pocket, n_ext))
                if best is not None:
                    key, rmsd, n_matched, n_ext = best
                    best_rows.append((target, key, rmsd, n_matched, n_pocket, n_ext))
                if i % 50 == 0:
                    done = i / len(jobs)
                    elapsed = time.time() - t0
                    self.logger.info(
                        f"  {i}/{len(jobs)} ({done:.0%}) elapsed {elapsed:.0f}s, "
                        f"ETA {elapsed / done - elapsed:.0f}s, {len(rows)} hits"
                    )

        header = ("chembl_target\texternal_entry\tpocket_rmsd\tn_matched"
                  "\tn_pocket_chembl\tn_pocket_external\n")

        rows.sort(key=lambda r: (r[0], r[2]))
        with open(output, "w") as fh:
            fh.write(header)
            for r in rows:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.3f}\t{r[3]}\t{r[4]}\t{r[5]}\n")

        best_rows.sort(key=lambda r: r[2])
        with open(best_output, "w") as fh:
            fh.write(header)
            for r in best_rows:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.3f}\t{r[3]}\t{r[4]}\t{r[5]}\n")

        self.logger.info(
            f"Done in {time.time()-t0:.0f}s: {len(rows)} hits at "
            f"rmsd <= {rmsd_report} and n_matched >= {min_matched} -> {output}; "
            f"best match per target -> {best_output}"
        )
        return output, best_output

    # ── load ──────────────────────────────────────────────────────────────

    @staticmethod
    def load(cache: Path) -> List[Tuple[str, Tuple[np.ndarray, ...]]]:
        """Rehydrate the cache into prepared pockets, ready for superposition.

        Prepared rather than as (name, CA) lists: every one of these is
        compared against every ChEMBL pocket, and each worker holds its own
        copy, so the packed form is both faster and far smaller.
        """
        data = np.load(cache, allow_pickle=False)
        keys = data["keys"]
        lengths = data["lengths"]
        codes = data["codes"]
        coords = data["coords"]

        pockets = []
        start = 0
        for key, length in zip(keys, lengths):
            end = start + int(length)
            names = [RESIDUE_CODES[c] if c < len(RESIDUE_CODES) else "UNK"
                     for c in codes[start:end]]
            xyz = coords[start:end].astype(float)
            pockets.append((str(key), prepare_pocket(list(zip(names, xyz)))))
            start = end
        return pockets

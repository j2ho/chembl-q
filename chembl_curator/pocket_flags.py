# chembl_curator/pocket_flags.py

"""Whether a target's pocket was defined on a bound cofactor.

The dataset's premise is that a target's actives bind the pocket it defines.
A tightly bound prosthetic group stays put and is not what a drug displaces,
so when one is all that any structure of the target offers, nothing in the
data supports that premise for that target.

Reported, never applied. Two reasons it is not a cut:

  * the pocket is often right anyway. For cytochromes P450, the nitric oxide
    synthase oxygenase domain, aromatase and the peroxidases, inhibitors stack
    on the cofactor's distal face. All 26 such accessions in the release turn
    out to have a drug-bound structure, which puts them in the `displaced`
    class below rather than `unverifiable`.
  * the two classes mean different things to a reader, and only one of them
    will ever change.

Measured on the 1,506-target release: 87 targets have a pocket defined
exclusively on a prosthetic group, 32 of them in the test split, 11 in the
recommended 125 and 3 in the strict 52.
"""

from pathlib import Path
from typing import List, Optional

ASSET = Path(__file__).parent / "assets" / "prosthetic_groups.txt"

# Co-substrates. They turn over, and whether an inhibitor competes with one or
# binds beside it is target-specific, so they are not prosthetic groups and
# they do not count as a drug-like ligand either. Disjoint from the prosthetic
# asset by construction; a test asserts it, which is how MGD, the molybdenum
# cofactor, got moved out of here.
CO_SUBSTRATES = frozenset({
    'NAD', 'NAI', 'NAJ', 'NAP', 'NDP', 'NAX', 'NAQ', 'NHD',
    'SAM', 'SAH', 'MTA', 'SFG', '5AD',
    'ATP', 'ADP', 'AMP', 'ANP', 'ACP', 'AGS', 'GTP', 'GDP', 'GNP',
    'COA', 'ACO', 'THF', 'UDP', 'UTP', 'CTP',
})

# Solvent, buffer, cryoprotectant, detergent and lipid. Present in a structure
# without being a ligand of it, and not evidence of a druggable site.
NOT_A_LIGAND = frozenset({
    'HOH', 'SO4', 'PO4', 'GOL', 'EDO', 'ACT', 'DMS', 'TRS', 'MES', 'EPE',
    'PEG', 'PG4', 'PGE', 'PE4', 'P6G', '1PE', 'MPD', 'BU3', 'HEZ', 'PGO',
    'CL', 'NA', 'MG', 'ZN', 'CA', 'K', 'MN', 'FE', 'FE2', 'NI', 'CU', 'CU1',
    'CD', 'CO', 'IOD', 'BR', 'NO3', 'ACY', 'FMT', 'IMD', 'AZI', 'CYN', 'OXY',
    'NH4', 'BME', 'DTT', 'MLI', 'TLA', 'CIT', 'SCN', 'SIN', 'MOH',
    'NAG', 'NDG', 'BMA', 'MAN', 'FUC', 'GAL', 'BGC', 'SIA', 'XYS',
    'MYR', 'PLM', 'OLA', 'OLB', 'STE', 'CDL', 'PEE', 'LDA', 'BOG', 'LMT',
    'LMN', 'C8E', 'CXE', 'D10', 'UND', 'OCT', 'SDS', 'JEF', 'DET', 'C15',
    'Y01', 'UNX', 'UNL',
})

# A HETATM group smaller than this is a fragment or an ion, not a ligand that
# could define a druggable site.
MIN_DRUG_LIKE_ATOMS = 12


def prosthetic_codes(path: Optional[Path] = None) -> frozenset:
    """PDB component ids of tightly bound prosthetic groups."""
    path = path or ASSET
    if not path.exists():
        return frozenset()
    with open(path) as handle:
        return frozenset(line.strip() for line in handle
                         if line.strip() and not line.startswith('#'))


def pocket_ligands(target_dir: Path) -> List[str]:
    """Component ids the target's pockets were defined on, deduplicated."""
    path = target_dir / "pocket_info.csv"
    if not path.exists():
        return []
    seen = []
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split(',')
            if len(parts) >= 4 and parts[3] and parts[3] not in seen:
                seen.append(parts[3])
    return sorted(seen)


def aligned_files(target_dir: Path) -> List[Path]:
    """The aligned chain PDBs the target's pockets were taken from."""
    path = target_dir / "pocket_info.csv"
    if not path.exists():
        return []
    out, seen = [], set()
    with open(path) as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip('\n').split(',')
            if len(parts) >= 3 and parts[2] and parts[2] not in seen:
                seen.add(parts[2])
                out.append(target_dir / "aligned" / parts[2])
    return out


def has_drug_like_structure(target_dir: Path, prosthetic: frozenset) -> bool:
    """Does any structure of the target hold a plausible drug-like ligand?

    This is what separates `displaced` from `unverifiable`. A target with such
    a structure has somewhere better to put its pocket, so the pocket will
    move once the ranking prefers it; a target without one never will.
    """
    for path in aligned_files(target_dir):
        if not path.exists():
            continue
        counts = {}
        with open(path) as handle:
            for line in handle:
                if line.startswith('HETATM'):
                    key = (line[17:20].strip(), line[21], line[22:27].strip())
                    counts[key] = counts.get(key, 0) + 1
        for (code, _, _), n_atoms in counts.items():
            if n_atoms < MIN_DRUG_LIKE_ATOMS:
                continue
            if code in prosthetic or code in CO_SUBSTRATES:
                continue
            if code in NOT_A_LIGAND:
                continue
            return True
    return False


def pocket_cofactor_status(target_dir: Path, prosthetic: frozenset) -> str:
    """'', 'unverifiable' or 'displaced'.

    Empty means the pocket is defined on something other than a prosthetic
    group, which is the ordinary case.

    `unverifiable`: every structure offers only a prosthetic group, so there
    is no structural evidence of where the actives bind, and none is coming.

    `displaced`: a drug-like ligand is present in some structure but the
    cofactor won the per-structure ranking. These pockets are wrong now and
    will move in a release that rebuilds them, so a reader comparing releases
    should expect exactly these to change.
    """
    if not prosthetic:
        return ''
    ligs = pocket_ligands(target_dir)
    if not ligs or not all(code in prosthetic for code in ligs):
        return ''
    return ('displaced' if has_drug_like_structure(target_dir, prosthetic)
            else 'unverifiable')


def pass_pocket_noncofactor(status: str) -> str:
    """1 when the pocket is not defined on a cofactor alone.

    Phrased to match the other `pass_` columns, where 1 is the target clearing
    the bar. The status column says which way it failed. Blank, never 0, when
    the asset is missing: 0 would read as "checked, the pocket is fine".
    """
    return '' if status is None else int(not status)

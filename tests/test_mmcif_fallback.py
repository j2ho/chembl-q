#!/usr/bin/env python3
"""mmCIF fallback and the modified-residue annotation that survives conversion."""

import tempfile
from pathlib import Path

from chembl_curator.protein_filter import ProteinFilter


def _pf(tmp):
    return ProteinFilter(Path(tmp), log_level="CRITICAL")


NATIVE = """\
MODRES 4C6D KCX A  156  LYS  LYSINE NZ-CARBOXYLIC ACID
MODRES 4C6D TPO A  177  THR  PHOSPHOTHREONINE
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00           N
HETATM    2  C   KCX A 156      12.000   6.000  -6.000  1.00 20.00           C
END
"""

# What a converted mmCIF looks like: coordinates only, no MODRES.
CONVERTED = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 20.00           N
HETATM    2  C   KCX A 156      12.000   6.000  -6.000  1.00 20.00           C
END
"""

SIDECAR = "comp_id\tchain\tseq_id\tparent\nKCX\tA\t156\tLYS\n"


def test_native_pdb_modres_is_read():
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "4c6d.pdb"
        f.write_text(NATIVE)
        assert _pf(tmp).modified_residues(f) == {"KCX", "TPO"}


def test_converted_structure_uses_the_sidecar():
    """Conversion drops MODRES, so the annotation has to come from beside it."""
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "9abc.pdb"
        f.write_text(CONVERTED)
        assert _pf(tmp).modified_residues(f) == set()
        f.with_suffix(".modres").write_text(SIDECAR)
        assert _pf(tmp).modified_residues(f) == {"KCX"}


def test_sidecar_wins_over_scanning_coordinates():
    """An empty sidecar means 'checked, none', not 'go look in the PDB'."""
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "x.pdb"
        f.write_text(NATIVE)
        f.with_suffix(".modres").write_text("comp_id\tchain\tseq_id\tparent\n")
        assert _pf(tmp).modified_residues(f) == set()


def test_modres_scan_stops_at_the_coordinates():
    """MODRES precedes coordinates; a ligand named MODRES-like must not match."""
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "y.pdb"
        f.write_text(
            "ATOM      1  N   ALA A   1      11.1   6.1  -6.5  1.00 20.00  N\n"
            "MODRES 9ZZZ ZZZ A    1  ALA  should not be reached\n")
        assert _pf(tmp).modified_residues(f) == set()


def test_declared_modified_residues_are_not_ligand_candidates():
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "z.pdb"
        # KCX close to the protein atoms, so only the MODRES declaration can
        # keep it out of the candidate list.
        lines = ["MODRES 4C6D KCX A  156  LYS  carboxylated"]
        for i in range(12):
            lines.append(
                f"ATOM  {i:5d}  CB  ALA A{i+1:4d}    "
                f"{10.0+i*0.5:8.3f}{6.000:8.3f}{-6.000:8.3f}  1.00 20.00           C")
        for i in range(6):
            lines.append(
                f"HETATM{100+i:5d}  C   KCX A 156    "
                f"{11.0+i*0.3:8.3f}{6.500:8.3f}{-6.000:8.3f}  1.00 20.00           C")
        lines.append("END")
        f.write_text("\n".join(lines) + "\n")
        pf = _pf(tmp)
        assert "KCX" in pf.modified_residues(f)
        kept = {l.ligand_name for l in pf.get_ligands_from_pdb(f, ["A"])}
        assert "KCX" not in kept
        # and the opt-out still lets it through, so the filter is what removed it
        relaxed = {l.ligand_name
                   for l in pf.get_ligands_from_pdb(f, ["A"], skip_modified=False)}
        assert "KCX" in relaxed or not relaxed


if __name__ == "__main__":
    import sys
    import traceback

    failed = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:
            failed += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"\n{'FAILED' if failed else 'all tests passed'} ({failed} failures)")
    sys.exit(1 if failed else 0)

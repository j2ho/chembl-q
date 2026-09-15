#!/usr/bin/env python3
"""Membrane phospholipids must not define a binding pocket, real ligands must.

A cryo-EM membrane structure ships its lipid annulus as ordered HETATM, and a
phospholipid then wins the representative-ligand vote: 26 targets in
curated_v5 were scored that way, worst case a 279-residue "pocket" built from
13 copies of 3PE.

The interesting half of this list is what it leaves out. A screen for
phosphorus plus a long chain also catches prenyl diphosphates, bisphosphonate
drugs and lipid A, which are the genuine ligands of their targets. The rule
drawn here is a glycerol backbone carrying a phosphate head and at least one
acyl or alkyl chain.
"""

import pytest

from chembl_curator.protein_filter import EXCLUDED_LIGANDS

# Confirmed against the PDB chemical component dictionary.
PHOSPHOLIPIDS = [
    ("POV", "1-palmitoyl-2-oleoyl phosphatidylcholine"),
    ("PC1", "1,2-diacyl-sn-glycero-3-phosphocholine"),
    ("PCW", "1,2-dioleoyl-sn-glycero-3-phosphocholine"),
    ("PEE", "1,2-dioleoyl-sn-glycero-3-phosphoethanolamine"),
    ("3PE", "1,2-distearoyl-sn-glycerophosphoethanolamine"),
    ("CDL", "cardiolipin"),
    ("PGW", "phosphatidylglycerol"),
    ("Q3G", "phosphatidylserine"),
    ("3PH", "phosphatidic acid"),
    ("81Q", "phosphatidylinositol"),
    ("NKP", "lysophosphatidic acid"),
]

REAL_LIGANDS = [
    ("GRG", "geranylgeranyl diphosphate, a prenyltransferase substrate"),
    ("ELU", "phytyl diphosphate"),
    ("PS7", "cyclopropyl prenyl diphosphate"),
    ("WJP", "polyprenyl diphosphate"),
    ("B29", "bisphosphonate drug"),
    ("749", "bisphosphonate drug"),
    ("SRL", "phosphonate drug"),
    ("K6L", "sphingoid phosphate, not glycerol-based"),
    ("LP4", "lipid A glycolipid, a TLR4/MD-2 ligand"),
    ("LP5", "lipid A glycolipid"),
    ("PC", "free phosphocholine head group, a hapten ligand"),
]

# Steroid hormones and retinoids share the composition of a membrane lipid
# (long, carbon-rich, no phosphorus) and are exactly what their targets bind.
HORMONES = [
    ("EST", "estradiol"),
    ("DHT", "dihydrotestosterone"),
    ("HCY", "hydrocortisone"),
    ("9CR", "9-cis-retinoic acid"),
    ("RTL", "retinol"),
    ("ACD", "arachidonic acid"),
]


@pytest.mark.parametrize("code,what", PHOSPHOLIPIDS)
def test_phospholipids_are_excluded(code, what):
    assert code in EXCLUDED_LIGANDS, f"{code} ({what}) should be excluded"


@pytest.mark.parametrize("code,what", REAL_LIGANDS + HORMONES)
def test_real_ligands_are_not_excluded(code, what):
    assert code not in EXCLUDED_LIGANDS, f"{code} ({what}) must stay selectable"


def test_no_duplicate_codes():
    from pathlib import Path

    import chembl_curator

    path = Path(chembl_curator.__file__).parent / "assets" / "excluded_ligands.txt"
    codes = [l.strip() for l in path.read_text().splitlines()
             if l.strip() and not l.startswith("#")]
    dupes = {c for c in codes if codes.count(c) > 1}
    assert not dupes, dupes
    assert len(codes) == len(EXCLUDED_LIGANDS)

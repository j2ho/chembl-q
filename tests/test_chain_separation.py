#!/usr/bin/env python3
"""
Test script to verify chain-specific PDB file creation.
"""

from pathlib import Path
from chembl_curator.protein_filter import ProteinFilter

def test_chain_separation(target_id="P13738"):
    """
    Test the new chain separation logic with P13738.
    This target has chains A, B, C, D with LMU ligand at A/B interface.
    """
    curated_dir = Path("test_curated")
    curated_dir.mkdir(exist_ok=True)

    # Create target directory
    target_dir = curated_dir / target_id
    target_dir.mkdir(exist_ok=True)

    # Initialize filter
    pf = ProteinFilter(curated_dir, log_level="INFO")

    # Process target
    print(f"\n{'='*60}")
    print(f"Testing chain separation for {target_id}")
    print(f"{'='*60}\n")

    result = pf.process_target(target_id)

    print(f"\n{'='*60}")
    print(f"Result: {'PASSED' if result else 'FAILED'}")
    print(f"{'='*60}\n")

    # Check aligned directory
    aligned_dir = target_dir / "aligned"
    if aligned_dir.exists():
        print("Aligned PDB files created:")
        for pdb_file in sorted(aligned_dir.glob("*.pdb")):
            # Count chains in the file
            chains = set()
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"):
                        chain = line[21:22].strip()
                        if chain:
                            chains.add(chain)
            print(f"  - {pdb_file.name}: chains = {sorted(chains)}")

    # Check pocket_info.csv
    pocket_csv = target_dir / "pocket_info.csv"
    if pocket_csv.exists():
        print("\nPocket info:")
        with open(pocket_csv, 'r') as f:
            for line in f:
                print(f"  {line.strip()}")

    return result

if __name__ == "__main__":
    test_chain_separation()

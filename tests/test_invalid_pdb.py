#!/usr/bin/env python3
"""Test error handling for invalid PDB ID."""

import tempfile
from pathlib import Path
from chembl_curator import ProteinFilter

def test_invalid_pdb():
    """Test that invalid PDB IDs are handled gracefully."""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Use an obviously invalid PDB ID
        invalid_pdb_id = "XXXX"

        print(f"Testing error handling with invalid PDB ID: {invalid_pdb_id}")
        print(f"Temporary directory: {tmpdir}")

        # Create ProteinFilter instance
        pf = ProteinFilter(curated_dir=tmpdir, log_level="DEBUG")

        # Try to download invalid PDB
        success = pf.download_pdb(invalid_pdb_id, tmpdir)

        # Check results
        expected_file = tmpdir / f"{invalid_pdb_id.lower()}.pdb"

        if not success and not expected_file.exists():
            print(f"✓ SUCCESS: Invalid PDB correctly rejected")
            print(f"  Download returned False as expected")
            print(f"  No file was created")
            return True
        else:
            print(f"✗ FAILED: Invalid PDB was not rejected properly")
            if success:
                print(f"  Download returned True (should be False)")
            if expected_file.exists():
                print(f"  File was created (should not exist)")
            return False

if __name__ == "__main__":
    success = test_invalid_pdb()
    exit(0 if success else 1)

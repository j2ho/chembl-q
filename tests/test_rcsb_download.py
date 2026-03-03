#!/usr/bin/env python3
"""Test RCSB web download (fallback mechanism)."""

import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch
from chembl_curator import ProteinFilter

def test_rcsb_download():
    """Test downloading from RCSB (simulating pdb_get not available)."""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        test_pdb_id = "1ATP"

        print(f"Testing RCSB web download for {test_pdb_id}...")
        print(f"Temporary directory: {tmpdir}")

        # Create ProteinFilter instance
        pf = ProteinFilter(curated_dir=tmpdir, log_level="DEBUG")

        # Mock shutil.which to return None (simulating pdb_get not available)
        with patch('shutil.which', return_value=None):
            print("Simulating pdb_get NOT available...")
            success = pf.download_pdb(test_pdb_id, tmpdir)

        # Check results
        expected_file = tmpdir / f"{test_pdb_id.lower()}.pdb"

        if success and expected_file.exists():
            file_size = expected_file.stat().st_size
            print(f"✓ SUCCESS: Downloaded {test_pdb_id} from RCSB")
            print(f"  File: {expected_file}")
            print(f"  Size: {file_size} bytes")

            # Read first few lines
            with open(expected_file, 'r') as f:
                first_lines = [f.readline().strip() for _ in range(5)]

            print(f"  First lines:")
            for line in first_lines:
                print(f"    {line}")

            # Validate PDB format
            valid_pdb = any(line.startswith(('HEADER', 'ATOM', 'HETATM', 'REMARK'))
                          for line in first_lines)

            if valid_pdb:
                print(f"✓ File appears to be a valid PDB structure")
                print(f"✓ RCSB fallback download works!")
                return True
            else:
                print(f"✗ File does not appear to be a valid PDB structure")
                return False
        else:
            print(f"✗ FAILED: Could not download {test_pdb_id} from RCSB")
            return False

if __name__ == "__main__":
    success = test_rcsb_download()
    exit(0 if success else 1)

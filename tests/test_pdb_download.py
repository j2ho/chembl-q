#!/usr/bin/env python3
"""Test PDB download with fallback mechanism."""

import tempfile
import shutil
from pathlib import Path
from chembl_curator import ProteinFilter

def test_pdb_download():
    """Test downloading a PDB file using both methods."""

    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Test with a small, well-known PDB structure
        test_pdb_id = "1ATP"

        print(f"Testing PDB download for {test_pdb_id}...")
        print(f"Temporary directory: {tmpdir}")

        # Check if pdb_get is available
        pdb_get_available = shutil.which('pdb_get') is not None
        print(f"pdb_get available: {pdb_get_available}")

        # Create ProteinFilter instance
        pf = ProteinFilter(curated_dir=tmpdir, log_level="DEBUG")

        # Test download
        success = pf.download_pdb(test_pdb_id, tmpdir)

        # Check results
        expected_file = tmpdir / f"{test_pdb_id.lower()}.pdb"

        if success and expected_file.exists():
            file_size = expected_file.stat().st_size
            print(f"✓ SUCCESS: Downloaded {test_pdb_id}")
            print(f"  File: {expected_file}")
            print(f"  Size: {file_size} bytes")

            # Read first few lines to verify it's a valid PDB
            with open(expected_file, 'r') as f:
                first_lines = [f.readline().strip() for _ in range(5)]

            print(f"  First lines:")
            for line in first_lines:
                print(f"    {line}")

            # Check if it looks like a valid PDB file
            valid_pdb = any(line.startswith(('HEADER', 'ATOM', 'HETATM', 'REMARK'))
                          for line in first_lines)

            if valid_pdb:
                print(f"✓ File appears to be a valid PDB structure")
                return True
            else:
                print(f"✗ File does not appear to be a valid PDB structure")
                return False
        else:
            print(f"✗ FAILED: Could not download {test_pdb_id}")
            if not success:
                print(f"  Download method returned False")
            if not expected_file.exists():
                print(f"  Expected file not found: {expected_file}")
            return False

if __name__ == "__main__":
    success = test_pdb_download()
    exit(0 if success else 1)

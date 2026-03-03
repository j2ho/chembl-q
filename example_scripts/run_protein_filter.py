#!/usr/bin/env python3
"""
Example script to run protein filtering pipeline on curated ChEMBL data.

This script demonstrates how to use the ProteinFilter class to:
1. Fetch PDB structures from UniProt
2. Download PDB files and AlphaFold models
3. Detect ligand-bound structures
4. Align structures using TMalign
5. Filter targets with single binding sites

Usage:
    python run_protein_filter.py --curated-dir curated_data_filtered --n-processes 8
"""

import argparse
from pathlib import Path
import sys

# Add parent directory to path to import chembl_curator
sys.path.insert(0, str(Path(__file__).parent.parent))

from chembl_curator import ProteinFilter


def main():
    parser = argparse.ArgumentParser(
        description='Filter protein structures based on PDB availability and binding site analysis'
    )
    parser.add_argument(
        '--curated-dir',
        type=str,
        required=True,
        help='Directory containing curated targets (uniprot IDs)'
    )
    parser.add_argument(
        '--n-processes',
        type=int,
        default=1,
        help='Number of parallel processes (default: 1)'
    )
    parser.add_argument(
        '--log-level',
        type=str,
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )

    args = parser.parse_args()

    print(f"Starting protein filtering pipeline...")
    print(f"Curated directory: {args.curated_dir}")
    print(f"Number of processes: {args.n_processes}")
    print()

    # Create protein filter
    protein_filter = ProteinFilter(
        curated_dir=Path(args.curated_dir),
        log_level=args.log_level
    )

    # Run pipeline
    passed_targets = protein_filter.run_pipeline(n_processes=args.n_processes)

    # Print summary
    print("\n" + "="*60)
    print("Protein filtering completed!")
    print(f"Passed targets: {len(passed_targets)}")
    print(f"Results saved to: {Path(args.curated_dir) / 'passed_targets.txt'}")
    print("="*60)

    # Print some example passed targets
    if passed_targets:
        print("\nExample passed targets:")
        for target in passed_targets[:10]:
            print(f"  - {target}")
        if len(passed_targets) > 10:
            print(f"  ... and {len(passed_targets) - 10} more")


if __name__ == '__main__':
    main()

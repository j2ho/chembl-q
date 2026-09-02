#!/usr/bin/env python3
"""Test script for new BAO format, confidence score, assay type, and pChEMBL filters."""

from chembl_curator import ChEMBLCurator, CurationConfig
from pathlib import Path

def test_new_filters():
    """Test the new filtering options."""

    # Create config with the user's specified filters
    config = CurationConfig(
        # Existing filters
        activity_thresholds={'nM': 10000.0, 'uM': 10.0},
        target_types=['SINGLE PROTEIN'],
        activity_types=['IC50', 'Ki', 'Kd', 'EC50'],
        relations=['=', '<='],
        units=['nM', 'uM'],

        # New filters as requested
        bao_formats=['BAO_0000357'],      # Single protein format
        min_confidence_score=5,            # Confidence >= 5
        assay_types=['B'],                 # Binding assays only
        min_pchembl_value=5.0,             # pChEMBL >= 5.0 (10 µM active)

        # Validity filters
        exclude_invalid_data=True,
        exclude_duplicates=True
    )

    print("Configuration:")
    print(f"  BAO Format: {config.bao_formats}")
    print(f"  Min Confidence Score: {config.min_confidence_score}")
    print(f"  Assay Types: {config.assay_types}")
    print(f"  Min pChEMBL Value: {config.min_pchembl_value}")
    print()

    # Create curator with the config
    curator = ChEMBLCurator(config=config, log_level="INFO")

    # Use existing database
    db_path = Path("./chembl_data/chembl_35.db")

    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return

    print(f"Using database: {db_path}")
    print()

    # Run the curation pipeline
    print("Running curation with new filters...")
    results = curator.run_pipeline(
        database_path=db_path,
        output_dir=Path("./curated_data_filtered")
    )

    # Print results
    print("\n" + "="*60)
    print("CURATION RESULTS")
    print("="*60)
    print(f"Total activities extracted: {results.total_activities:,}")
    print(f"Activities passing filters: {results.filtered_activities:,}")
    print(f"Filter efficiency: {results.filtered_activities/results.total_activities*100:.2f}%")
    print(f"Total proteins: {results.total_proteins:,}")
    print(f"Total compounds: {results.total_compounds:,}")
    print(f"Output directory: {results.output_directory}")
    print("="*60)

if __name__ == "__main__":
    test_new_filters()

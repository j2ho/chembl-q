#!/usr/bin/env python3

"""Test script to validate the new validity filter functionality"""

import sys
import sqlite3
from pathlib import Path
from chembl_curator.config import CurationConfig

def test_validity_filters():
    """Test the new validity filtering options"""
    
    # Test different configurations
    configs = {
        "no_filters": CurationConfig(
            require_standard_flag=False,
            exclude_invalid_data=False, 
            exclude_duplicates=False
        ),
        "default_filters": CurationConfig(
            require_standard_flag=False,
            exclude_invalid_data=True,
            exclude_duplicates=True
        ),
        "strict_filters": CurationConfig(
            require_standard_flag=True,
            exclude_invalid_data=True,
            exclude_duplicates=True
        )
    }
    
    db_path = Path("chembl_data/chembl_35.db")
    if not db_path.exists():
        print("Database not found!")
        return
        
    conn = sqlite3.connect(db_path)
    
    print("Testing validity filter configurations:")
    print("=" * 50)
    
    for config_name, config in configs.items():
        # Simulate the query building logic from curator.py
        target_types = "','".join(config.target_types)
        activity_types = "','".join(config.activity_types) 
        relations = "','".join(config.relations)
        units = "','".join(config.units)
        
        validity_conditions = []
        if config.require_standard_flag:
            validity_conditions.append("a.standard_flag = 1")
        if config.exclude_invalid_data:
            validity_conditions.append("a.data_validity_comment IS NULL")
        if config.exclude_duplicates:
            validity_conditions.append("(a.potential_duplicate IS NULL OR a.potential_duplicate = 0)")
        
        validity_filter = ""
        if validity_conditions:
            validity_filter = "AND " + " AND ".join(validity_conditions)
        
        query = f"""
        SELECT COUNT(*) 
        FROM activities a
        JOIN molecule_dictionary md ON a.molregno = md.molregno
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE 
            td.target_type IN ('{target_types}')
            AND a.standard_type IN ('{activity_types}')
            AND a.standard_relation IN ('{relations}')
            AND a.standard_units IN ('{units}')
            AND a.standard_value IS NOT NULL
            {validity_filter}
        """
        
        count = conn.execute(query).fetchone()[0]
        print(f"{config_name:15}: {count:>10,} activities")
        
        # Show what filters are applied
        filters_applied = []
        if config.require_standard_flag:
            filters_applied.append("standard_flag=1")
        if config.exclude_invalid_data:
            filters_applied.append("no_validity_comments")
        if config.exclude_duplicates:
            filters_applied.append("no_duplicates")
        
        filter_desc = ", ".join(filters_applied) if filters_applied else "none"
        print(f"               Filters: {filter_desc}")
        print()
    
    conn.close()
    
    print("✅ Validity filter testing completed!")

if __name__ == "__main__":
    test_validity_filters()

# ============================================================================
# chembl_curator/utils.py  
# ============================================================================

"""Utility functions for ChEMBL curation."""

import logging
import sys
from pathlib import Path
from typing import Dict, Any


def setup_logging(level: str = "INFO") -> logging.Logger:
    """Setup logging configuration.
    
    Args:
        level: Logging level
        
    Returns:
        Configured logger
    """
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger("chembl_curator")


def check_dependencies():
    """Check if required dependencies are installed."""
    missing = []
    
    try:
        import rdkit
    except ImportError:
        missing.append("rdkit")
    
    try:
        import pandas
    except ImportError:
        missing.append("pandas")
    
    try:
        import requests
    except ImportError:
        missing.append("requests")
    
    if missing:
        raise ImportError(f"Missing dependencies: {', '.join(missing)}")


def validate_chembl_database(db_path: Path) -> bool:
    """Validate ChEMBL SQLite database structure.
    
    Args:
        db_path: Path to database file
        
    Returns:
        True if valid ChEMBL database
    """
    import sqlite3
    
    required_tables = [
        'activities',
        'molecule_dictionary', 
        'compound_structures',
        'target_dictionary',
        'assays'
    ]
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Get table names
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = {row[0] for row in cursor.fetchall()}
        
        # Check required tables exist
        missing_tables = set(required_tables) - tables
        if missing_tables:
            logging.error(f"Missing required tables: {missing_tables}")
            return False
        
        return True
        
    except Exception as e:
        logging.error(f"Database validation failed: {e}")
        return False
    
    finally:
        if 'conn' in locals():
            conn.close()

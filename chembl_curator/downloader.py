# ============================================================================
# chembl_curator/downloader.py
# ============================================================================

"""ChEMBL database downloader utilities."""

import os
import tarfile
import sqlite3
from pathlib import Path
from typing import Optional
import requests
from tqdm import tqdm
import logging


class ChEMBLDownloader:
    """Handles downloading and setting up ChEMBL SQLite database."""
    
    def __init__(self, version: str = "35"):
        self.version = version
        self.base_url = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}"
        self.logger = logging.getLogger(__name__)
    
    def download_sqlite(
        self, 
        output_dir: Optional[Path] = None,
        force_download: bool = False
    ) -> Path:
        """Download and setup ChEMBL SQLite database.
        
        Args:
            output_dir: Directory to save database
            force_download: Force re-download even if exists
            
        Returns:
            Path to the SQLite database file
        """
        if output_dir is None:
            output_dir = Path.cwd() / "chembl_data"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        db_file = output_dir / f"chembl_{self.version}.db"
        
        # Check if database already exists
        if db_file.exists() and not force_download:
            self.logger.info(f"Database already exists: {db_file}")
            return db_file
        
        # Download archive
        archive_file = self._download_archive(output_dir, force_download)
        
        # Extract and setup database
        self._setup_database(archive_file, db_file)
        
        # Cleanup archive
        if archive_file.exists():
            archive_file.unlink()
            
        return db_file
    
    def _download_archive(self, output_dir: Path, force_download: bool) -> Path:
        """Download the SQLite archive."""
        archive_file = output_dir / f"chembl_{self.version}_sqlite.tar.gz"
        
        if archive_file.exists() and not force_download:
            self.logger.info(f"Archive already exists: {archive_file}")
            return archive_file
        
        url = f"{self.base_url}/chembl_{self.version}_sqlite.tar.gz"
        
        self.logger.info(f"Downloading from: {url}")
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(archive_file, 'wb') as f:
            with tqdm(total=total_size, unit='iB', unit_scale=True, desc="Downloading") as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
        
        self.logger.info(f"Download completed: {archive_file}")
        return archive_file
    
    def _setup_database(self, archive_file: Path, db_file: Path):
        """Extract archive and setup SQLite database."""
        self.logger.info("Extracting archive...")
        
        with tarfile.open(archive_file, 'r:gz') as tar:
            tar.extractall(archive_file.parent / "temp")
        
        # Find SQL files
        temp_dir = archive_file.parent / "temp"
        sql_files = list(temp_dir.rglob("*.sql"))
        
        if not sql_files:
            raise FileNotFoundError("No SQL files found in archive")
        
        self.logger.info(f"Setting up database from {len(sql_files)} SQL files...")
        
        # Create database
        conn = sqlite3.connect(db_file)
        
        try:
            for sql_file in sorted(sql_files):
                self.logger.info(f"Processing {sql_file.name}...")
                with open(sql_file, 'r') as f:
                    conn.executescript(f.read())
            
            conn.commit()
            self.logger.info("Database setup completed")
            
        finally:
            conn.close()
        
        # Cleanup temp directory
        import shutil
        shutil.rmtree(temp_dir)

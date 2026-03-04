# chembl_curator/downloader.py

import tarfile
import sqlite3
from pathlib import Path
from typing import Optional
import requests
from tqdm import tqdm
import logging


class ChEMBLDownloader:    
    def __init__(self, version: str = "36"):
        self.version = version
        self.base_url = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version}"
        self.logger = logging.getLogger(__name__)
    
    def download_sqlite(
        self, 
        output_dir: Optional[Path] = None,
        force_download: bool = False
    ) -> Path:
        """
        Args:
            output_dir: Where to save curated data
            force_download: Force re-download and overwrite existing files
            
        Returns:
            Path to the SQLite database file (chembl_{version}.db)
        """
        if output_dir is None:
            output_dir = Path.cwd() / "chembl_data"
        
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        db_file = output_dir / f"chembl_{self.version}.db"
        
        if db_file.exists() and not force_download:
            self.logger.info(f"Database already exists. Use pre-downloaded {db_file}")
            return db_file
        
        archive_file = self._download_archive(output_dir, force_download)

        self._setup_database(archive_file, db_file)

        if archive_file.exists():
            archive_file.unlink()

        return db_file
    
    def _download_archive(self, output_dir: Path, force_download: bool) -> Path:
        archive_file = output_dir / f"chembl_{self.version}_sqlite.tar.gz"
        
        if archive_file.exists() and not force_download:
            self.logger.info(f"ChEMBLdb raw archive file already exists: {archive_file}")
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
        
        temp_dir = archive_file.parent / "temp"
        
        db_files = list(temp_dir.rglob("*.db"))
        
        if db_files:
            source_db = db_files[0]
            self.logger.info(f"Found existing SQLite database: {source_db}")
            
            source_size = source_db.stat().st_size / (1024**3)  # GB
            self.logger.info(f"ChEMBL database size: {source_size:.2f} GB")
            
            if db_file.exists():
                self.logger.info("Removing existing ChEMBL DB")
                db_file.unlink()
            
            import shutil
            self.logger.info("Copying ChEMBL database file (this may take several minutes for large databases)...")
            shutil.copy2(source_db, db_file)
            
            if not db_file.exists():
                raise RuntimeError(f"Database copy failed: {db_file} does not exist")
                
            target_size = db_file.stat().st_size / (1024**3)
            self.logger.info(f"ChEMBL database size: {target_size:.2f} GB")
            
            if abs(source_size - target_size) > 0.01:  # Allow small differences due to rounding
                raise RuntimeError(f"Database copy verification failed: source ({source_size:.2f} GB) != target ({target_size:.2f} GB)")
            
            self._verify_database_integrity(db_file)
            self.logger.info(f"Database successfully copied and verified: {db_file}")
            
        else:
            sql_files = list(temp_dir.rglob("*.sql"))
            
            if not sql_files:
                raise FileNotFoundError("No SQL files or SQLite database found in archive")
            
            self.logger.info(f"Setting up database from {len(sql_files)} SQL files...")
            
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
        
        self._cleanup_temp_directory(temp_dir)
    
    def _verify_database_integrity(self, db_file: Path):
        """Verify SQLite database integrity - just in case! / only for distribution"""
        try:
            conn = sqlite3.connect(db_file)
            cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")
            result = cursor.fetchone()
            
            if not result:
                raise RuntimeError("Database appears empty or corrupted - no tables found")
                
            self.logger.info("Verifying database integrity...")
            cursor.execute("PRAGMA quick_check")
            integrity_result = cursor.fetchone()[0]
            
            if integrity_result != "ok":
                raise RuntimeError(f"Database integrity verification failed: {integrity_result}")

            conn.close()
            self.logger.info("Database integrity verified")
            
        except sqlite3.Error as e:
            raise RuntimeError(f"Database integrity verification failed: {e}")
    
    def _cleanup_temp_directory(self, temp_dir: Path):
        try:
            if temp_dir.exists():
                import shutil
                self.logger.info(f"Cleaning up temp: {temp_dir}")
                shutil.rmtree(temp_dir)
                self.logger.info("Temp directory cleaned up successfully")
        except Exception as e:
            self.logger.warning(f"Failed to cleanup temp directory {temp_dir}: {e}")


# ============================================================================
# chembl_curator/cli.py
# ============================================================================

"""Command line interface for ChEMBL Curator."""

import click
from pathlib import Path
from .curator import ChEMBLCurator
from .config import CurationConfig


@click.command()
@click.option('--database', '-d', type=click.Path(exists=True), 
              help='Path to ChEMBL SQLite database')
@click.option('--output', '-o', default='./curated_chembl', 
              help='Output directory for curated data')
@click.option('--download', is_flag=True, 
              help='Download ChEMBL database first')
@click.option('--log-level', default='INFO', 
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
def main(database, output, download, log_level):
    """ChEMBL bioactivity data curation tool."""
    
    curator = ChEMBLCurator(log_level=log_level)
    
    db_path = None
    if database:
        db_path = Path(database)
    elif download:
        click.echo("Downloading ChEMBL database...")
        db_path = curator.download_database()
    
    if db_path is None:
        click.echo("Error: Either provide --database path or use --download flag")
        return
    
    click.echo(f"Running curation pipeline...")
    results = curator.run_pipeline(
        database_path=db_path,
        output_dir=Path(output)
    )
    
    click.echo(f"Curation completed!")
    click.echo(f"Total compounds: {results.total_compounds}")
    click.echo(f"Total proteins: {results.total_proteins}")
    click.echo(f"Output directory: {results.output_directory}")


if __name__ == '__main__':
    main()
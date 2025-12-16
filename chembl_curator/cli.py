# chembl_curator/cli.py

import click
from pathlib import Path
from .curator import ChEMBLCurator
from .config import CurationConfig
from .protein_filter import ProteinFilter


@click.group()
def cli():
    """ChEMBL Curator - A tool for curating ChEMBL data and filtering protein structures."""
    pass


@cli.command(name='curate')
@click.option('--database', '-d', type=click.Path(exists=True),
              help='Path to ChEMBL SQLite database')
@click.option('--output', '-o', default='./curated_chembl',
              help='Output directory for curated data')
@click.option('--download', is_flag=True,
              help='Download ChEMBL database first')
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Path to configuration JSON file')
@click.option('--create-config', type=click.Path(),
              help='Create example configuration file at specified path')
@click.option('--activity-types', multiple=True,
              help='Activity types to include (e.g., Ki, Kd, IC50, EC50)')
@click.option('--relations', multiple=True,
              help='Relations to include (e.g., =, <=, <)')
@click.option('--units', multiple=True,
              help='Units to include (e.g., nM, uM, pM)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
def curate(database, output, download, config, create_config, activity_types, relations, units, log_level):
    """Curate ChEMBL database and filter ligands."""
    if create_config:
        CurationConfig.create_example_config(Path(create_config))
        click.echo(f"Example configuration created at: {create_config}")
        return

    # Load config
    curation_config = None
    if config:
        curation_config = CurationConfig.from_file(Path(config))
        click.echo(f"Loaded configuration from: {config}")
    else:
        curation_config = CurationConfig()

    # Override config with cli options if wanted
    if activity_types:
        curation_config.activity_types = list(activity_types)
    if relations:
        curation_config.relations = list(relations)
    if units:
        curation_config.units = list(units)

    curator = ChEMBLCurator(config=curation_config, log_level=log_level)

    db_path = None
    if database:
        db_path = Path(database)
    elif download:
        click.echo("Downloading ChEMBL database...")
        db_path = curator.download_database()

    if db_path is None:
        click.echo("Error: Either provide --database path (downloaded already) or use --download flag")
        return

    click.echo(f"Curation of ChEMBL database in progress...")
    results = curator.run_pipeline(
        database_path=db_path,
        output_dir=Path(output)
    )

    click.echo(f"YAY! Curation completed!")
    click.echo(f"Total compounds: {results.total_compounds}")
    click.echo(f"Total proteins: {results.total_proteins}")
    click.echo(f"Output directory: {results.output_directory}")


@cli.command(name='filter-proteins')
@click.option('--curated-dir', '-d', required=True, type=click.Path(exists=True),
              help='Directory containing curated targets (uniprot IDs)')
@click.option('--n-processes', '-n', default=1, type=int,
              help='Number of parallel processes (default: 1)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
def filter_proteins(curated_dir, n_processes, log_level):
    """
    Filter protein structures based on PDB availability and binding site analysis.

    This will:
    1. Fetch PDB structures from UniProt
    2. Download PDB files and AlphaFold models
    3. Detect ligand-bound structures
    4. Align structures using TMalign
    5. Filter targets with single binding sites
    """
    click.echo(f"Starting protein filtering pipeline...")
    click.echo(f"Curated directory: {curated_dir}")
    click.echo(f"Number of processes: {n_processes}")

    protein_filter = ProteinFilter(
        curated_dir=Path(curated_dir),
        log_level=log_level
    )

    passed_targets = protein_filter.run_pipeline(n_processes=n_processes)

    click.echo(f"\nProtein filtering completed!")
    click.echo(f"Passed targets: {len(passed_targets)}")
    click.echo(f"Results saved to: {Path(curated_dir) / 'passed_targets.txt'}")


# Keep backward compatibility with old command
@click.command()
@click.option('--database', '-d', type=click.Path(exists=True), 
              help='Path to ChEMBL SQLite database')
@click.option('--output', '-o', default='./curated_chembl', 
              help='Output directory for curated data')
@click.option('--download', is_flag=True, 
              help='Download ChEMBL database first')
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Path to configuration JSON file')
@click.option('--create-config', type=click.Path(),
              help='Create example configuration file at specified path')
@click.option('--activity-types', multiple=True,
              help='Activity types to include (e.g., Ki, Kd, IC50, EC50)')
@click.option('--relations', multiple=True,
              help='Relations to include (e.g., =, <=, <)')
@click.option('--units', multiple=True,
              help='Units to include (e.g., nM, uM, pM)')
@click.option('--log-level', default='INFO', 
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))

def main(database, output, download, config, create_config, activity_types, relations, units, log_level):
    # Example config file generation for users who have trouble setting up
    if create_config:
        CurationConfig.create_example_config(Path(create_config))
        click.echo(f"Example configuration created at: {create_config}")
        return
    
    # Load config
    curation_config = None
    if config:
        curation_config = CurationConfig.from_file(Path(config))
        click.echo(f"Loaded configuration from: {config}")
    else:
        curation_config = CurationConfig()
    
    # Override config with cli options if wanted
    if activity_types:
        curation_config.activity_types = list(activity_types)
    if relations:
        curation_config.relations = list(relations)
    if units:
        curation_config.units = list(units)
    
    curator = ChEMBLCurator(config=curation_config, log_level=log_level)
    
    db_path = None
    if database:
        db_path = Path(database)
    elif download:
        click.echo("Downloading ChEMBL database...")
        db_path = curator.download_database()
    
    if db_path is None:
        click.echo("Error: Either provide --database path (downloaded already) or use --download flag")
        return
    
    click.echo(f"Curation of ChEMBL database in progress...")
    results = curator.run_pipeline(
        database_path=db_path,
        output_dir=Path(output)
    )
    
    click.echo(f"YAY! Curation completed!")
    click.echo(f"Total compounds: {results.total_compounds}")
    click.echo(f"Total proteins: {results.total_proteins}")
    click.echo(f"Output directory: {results.output_directory}")


if __name__ == '__main__':
    cli()
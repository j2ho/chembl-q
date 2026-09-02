# chembl_curator/cli.py

import functools

import click
from pathlib import Path
from .curator import ChEMBLCurator
from .config import CurationConfig
from .protein_filter import ProteinFilter
from .active_clusterer import ActiveClusterer
from .compound_pool import CompoundPool
from .receptor_similarity import ReceptorSimilarity
from .decoy_selector import DecoySelector
from .splitter import TargetSplitter
from .utils import setup_logging


def with_logging(f):
    """Attach a root log handler before the command runs.

    Every stage class builds a module logger and sets its level, but only
    ChEMBLCurator ever called setup_logging(), so the root logger had no
    handler for the others. Stages 3-7 therefore ran silently: a stage 6 run
    reported its decoy counts and dropped every exclusion statistic on the
    floor. Configuring here keeps the library free of root-logger side effects
    while still making the CLI say what it did.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        setup_logging(kwargs.get("log_level") or "INFO")
        return f(*args, **kwargs)

    return wrapper


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
@with_logging
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
@click.option('--max-chain-residues', default=1500, type=int,
              help='Skip PDB structures whose target chain exceeds this many residues (0 = no limit)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def filter_proteins(curated_dir, n_processes, max_chain_residues, log_level):
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
        log_level=log_level,
        max_chain_residues=max_chain_residues,
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


@cli.command(name='cluster-actives')
@click.option('--data-dir', '-d', required=True, type=click.Path(exists=True),
              help='Root data directory with per-target subdirs (curated_data_filtered)')
@click.option('--dist-thresh', type=float, default=0.3, show_default=True,
              help='Butina Tanimoto distance threshold (0.3 = similarity >= 0.7)')
@click.option('--workers', '-n', type=int, default=1, show_default=True,
              help='Number of parallel worker processes')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def cluster_actives(data_dir, dist_thresh, workers, log_level):
    """Stage 3: Butina cluster actives per target.

    Reads actives.tsv (or falls back to .smi files) and writes
    actives_clustered.tsv to each target directory.
    """
    clusterer = ActiveClusterer(dist_thresh=dist_thresh, log_level=log_level)
    stats = clusterer.run(Path(data_dir), workers=workers)
    click.echo(f"Done: {stats['total_actives']} actives → {stats['total_clusters']} representatives "
               f"across {stats['n_targets']} targets")


@cli.command(name='build-pool')
@click.option('--data-dir', '-d', required=True, type=click.Path(exists=True),
              help='Root data directory')
@click.option('--output', '-o', type=click.Path(),
              help='Output pickle path (default: data_dir/compound_pool.pkl)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def build_pool(data_dir, output, log_level):
    """Stage 4: Build global compound pool from clustered actives.

    Reads actives_clustered.tsv from each passed target and writes
    compound_pool.pkl with per-compound properties and fingerprints.
    """
    pool_builder = CompoundPool(log_level=log_level)
    out_path = pool_builder.build(
        Path(data_dir),
        output=Path(output) if output else None,
    )
    click.echo(f"Compound pool saved to: {out_path}")


@cli.command(name='receptor-sim')
@click.option('--data-dir', '-d', required=True, type=click.Path(exists=True),
              help='Root data directory')
@click.option('--mode', type=click.Choice(['seqid', 'pocket', 'both']), default='both',
              show_default=True, help='Which similarity metrics to compute')
@click.option('--seqid-threads', type=int, default=4, show_default=True,
              help='MMseqs2 thread count (for seqid mode)')
@click.option('--workers', '-n', type=int, default=4, show_default=True,
              help='Worker processes for pocket RMSD computation')
@click.option('--pocket-radius', type=float, default=10.0, show_default=True,
              help='Pocket radius in Å for pocket RMSD')
@click.option('--pocket-method', default='tmalign',
              type=click.Choice(['tmalign', 'hungarian']),
              help='Pocket superposition: tmalign (sequence-ordered global '
                   'alignment) or hungarian (order-free matching). Pair with '
                   '--pocket-radius 8 for hungarian.')
@click.option('--output', '-o', type=click.Path(),
              help='Pocket RMSD output path (default: '
                   'data_dir/pairwise_pocket_rmsd.tsv). Use this to keep both '
                   'methods side by side.')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def receptor_sim(data_dir, mode, seqid_threads, workers, pocket_radius,
                 pocket_method, output, log_level):
    """Stage 5: Compute pairwise receptor similarity.

    Sequence identity via MMseqs2 all-vs-all (requires mmseqs in PATH).
    Pocket RMSD by TMalign (default) or order-free Hungarian matching.

    Outputs: pairwise_seqid.tsv and/or pairwise_pocket_rmsd.tsv
    """
    sim = ReceptorSimilarity(log_level=log_level)
    results = sim.run(
        Path(data_dir),
        mode=mode,
        seqid_threads=seqid_threads,
        workers=workers,
        pocket_radius=pocket_radius,
        pocket_method=pocket_method,
        pocket_output=Path(output) if output else None,
    )
    for key, path in results.items():
        click.echo(f"  {key}: {path}")


@cli.command(name='select-decoys')
@click.option('--data-dir', '-d', required=True, type=click.Path(exists=True),
              help='Root data directory')
@click.option('--max-decoys', type=int, default=30, show_default=True,
              help='Maximum decoys per active compound')
@click.option('--seqid-thresh', type=float, default=0.6, show_default=True,
              help='Seqid threshold for receptor exclusion')
@click.option('--pocket-rmsd-thresh', type=float, default=2.0, show_default=True,
              help='Pocket RMSD threshold (A) for receptor exclusion')
@click.option('--pocket-rmsd-tsv', type=click.Path(exists=True),
              help='Pocket similarity TSV (default: data_dir/pairwise_pocket_rmsd.tsv). '
                   'Point this at a specific method run to keep provenance explicit.')
@click.option('--min-matched-residues', type=int, default=15, show_default=True,
              help='Minimum superposed residues for a pocket RMSD to count. '
                   'A close fit over a handful of residues is not evidence.')
@click.option('--exclusion-mode', type=click.Choice(['or', 'and']), default='or',
              show_default=True,
              help='"or" = exclude if seqid OR pocket matches; "and" = both must match')
@click.option('--tanimoto-thresh', type=float, default=0.3, show_default=True,
              help='Max Tanimoto similarity between active and decoy')
@click.option('--seed', type=int, default=42, show_default=True,
              help='Random seed')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def select_decoys(data_dir, max_decoys, seqid_thresh, pocket_rmsd_thresh,
                  pocket_rmsd_tsv, min_matched_residues, exclusion_mode,
                  tanimoto_thresh, seed, log_level):
    """Stage 6: Receptor-aware decoy selection.

    For each active, selects up to max_decoys property-matched,
    chemically dissimilar decoys. Actives against similar receptors
    (by seqid and/or pocket RMSD) are excluded from the decoy pool.

    Writes decoys.tsv to each target directory.
    """
    selector = DecoySelector(
        max_decoys=max_decoys,
        seqid_thresh=seqid_thresh,
        pocket_rmsd_thresh=pocket_rmsd_thresh,
        min_matched_residues=min_matched_residues,
        exclusion_mode=exclusion_mode,
        tanimoto_thresh=tanimoto_thresh,
        seed=seed,
        log_level=log_level,
    )
    stats = selector.run(
        Path(data_dir),
        pocket_rmsd_tsv=Path(pocket_rmsd_tsv) if pocket_rmsd_tsv else None,
    )
    click.echo(f"Done: {stats.get('n_total_decoys', 0)} active-decoy pairs, "
               f"{stats.get('n_underfilled', 0)} actives underfilled")


@cli.command(name='split')
@click.option('--data-dir', '-d', required=True, type=click.Path(exists=True),
              help='Root data directory')
@click.option('--seqid', type=float, default=0.4, show_default=True,
              help='MMseqs2 clustering seqid threshold')
@click.option('--valid-frac', type=float, default=1.0, show_default=True,
              help='Fraction of clusters to assign to test set')
@click.option('--threads', type=int, default=4, show_default=True,
              help='MMseqs2 thread count')
@click.option('--external-fasta', multiple=True, type=click.Path(exists=True),
              help='External FASTA file(s) with dot-prefixed IDs (>{source}.{entry_id}). '
                   'Bundled PDBbind+BioLip FASTA is used by default.')
@click.option('--no-external', is_flag=True, default=False,
              help='Skip bundled external FASTA (ChEMBL-only split)')
@click.option('--output-dir', type=click.Path(),
              help='Output directory for train.txt/test.txt (default: data-dir)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']))
@with_logging
def split(data_dir, seqid, valid_frac, threads, external_fasta, no_external,
          output_dir, log_level):
    """Stage 7: Train/test split by sequence-identity clustering.

    Clusters all passed targets (+ external datasets) at the given seqid
    threshold. Clusters are greedily assigned to train/test while balancing
    each data source.

    By default, the bundled PDBbind+BioLip FASTA is included. Use
    --no-external for a ChEMBL-only split, or --external-fasta to provide
    your own file(s).

    Requires mmseqs in PATH.
    """
    # Resolve external FASTA: explicit > bundled default > none
    if external_fasta:
        fasta_paths = [Path(p) for p in external_fasta]
    elif no_external:
        fasta_paths = None
    else:
        bundled = Path(__file__).parent / "assets" / "external_targets.fasta"
        if bundled.exists():
            fasta_paths = [bundled]
            click.echo(f"Using bundled external FASTA: {bundled}")
        else:
            fasta_paths = None
            click.echo("Bundled external FASTA not found, running ChEMBL-only split")

    splitter = TargetSplitter(seqid=seqid, valid_frac=valid_frac,
                               threads=threads, log_level=log_level)
    train_path, test_path, targets_path = splitter.run(
        Path(data_dir),
        external_fasta=fasta_paths,
        output_dir=Path(output_dir) if output_dir else None,
    )
    click.echo(f"Train:   {train_path}")
    click.echo(f"Test:    {test_path}")
    click.echo(f"Targets: {targets_path}")


if __name__ == '__main__':
    cli()
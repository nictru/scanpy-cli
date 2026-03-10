import rich_click as click
import scanpy as sc
import os
from scanpy_cli.utils import catch_errors, logger


@click.command()
@click.argument("input", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
@click.option(
    "--genome",
    type=str,
    default=None,
    help="Filter expression to genes within this genome. For hdf5 files containing data from multiple genomic references.",
)
@click.option(
    "--gex-only/--no-gex-only",
    default=True,
    help="Whether to return expression data only or all modalities. Default: True",
)
@click.option(
    "--backup-url",
    type=str,
    default=None,
    help="URL to download the file from if the local file is not found.",
)
@catch_errors
def read_10x_h5(input, output, genome, gex_only, backup_url):
    """Read 10x-Genomics-formatted hdf5 file.

    This command reads data from a 10x-Genomics-formatted hdf5 file and saves it as an AnnData object.

    INPUT is the path to the input .h5 file.
    OUTPUT is the path where the AnnData object will be saved (.h5ad format).
    """
    logger.debug(
        "Reading 10x H5 file: %s (genome=%s, gex_only=%s)", input, genome, gex_only
    )
    if backup_url:
        logger.debug("Backup URL configured: %s", backup_url)

    adata = sc.read_10x_h5(
        input, genome=genome, gex_only=gex_only, backup_url=backup_url
    )
    logger.info("Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    adata.write(output)
    logger.info("Data successfully saved to %s", output)

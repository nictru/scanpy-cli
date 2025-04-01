import rich_click as click
import scanpy as sc
import os


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
def read_10x_h5(input, output, genome, gex_only, backup_url):
    """Read 10x-Genomics-formatted hdf5 file.

    This command reads data from a 10x-Genomics-formatted hdf5 file and saves it as an AnnData object.

    INPUT is the path to the input .h5 file.
    OUTPUT is the path where the AnnData object will be saved (.h5ad format).
    """
    # Read the data
    adata = sc.read_10x_h5(
        input, genome=genome, gex_only=gex_only, backup_url=backup_url
    )

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    # Save the AnnData object
    adata.write(output)
    click.echo(f"Data successfully saved to {output}")

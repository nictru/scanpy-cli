import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--min-counts",
    type=int,
    help="Minimum number of counts required for a cell to pass filtering.",
)
@click.option(
    "--min-genes",
    type=int,
    help="Minimum number of genes expressed required for a cell to pass filtering.",
)
@click.option(
    "--max-counts",
    type=int,
    help="Maximum number of counts required for a cell to pass filtering.",
)
@click.option(
    "--max-genes",
    type=int,
    help="Maximum number of genes expressed required for a cell to pass filtering.",
)
@click.option(
    "--input-file",
    "-i",
    required=True,
    help="Input h5ad file containing AnnData object.",
)
@click.option(
    "--output-file",
    "-o",
    required=True,
    help="Output h5ad file to save the processed AnnData object.",
)
def filter_cells(min_counts, min_genes, max_counts, max_genes, input_file, output_file):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least min_counts counts or min_genes genes expressed.
    This is to filter measurement outliers, i.e. "unreliable" observations.

    Only provide one of the optional parameters min_counts, min_genes, max_counts, max_genes per call.
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's filter_cells function
        sc.pp.filter_cells(
            adata,
            min_counts=min_counts,
            min_genes=min_genes,
            max_counts=max_counts,
            max_genes=max_genes,
        )

        # Save the result
        adata.write(output_file)

        click.echo(f"Successfully filtered cells and saved to {output_file}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

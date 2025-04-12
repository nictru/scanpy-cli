import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--min-counts",
    type=int,
    help="Minimum number of counts required for a gene to pass filtering.",
)
@click.option(
    "--min-cells",
    type=int,
    help="Minimum number of cells expressed required for a gene to pass filtering.",
)
@click.option(
    "--max-counts",
    type=int,
    help="Maximum number of counts required for a gene to pass filtering.",
)
@click.option(
    "--max-cells",
    type=int,
    help="Maximum number of cells expressed required for a gene to pass filtering.",
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
def filter_genes(min_counts, min_cells, max_counts, max_cells, input_file, output_file):
    """Filter genes based on number of cells or counts.

    Keep genes that have at least min_counts counts or are expressed in at least min_cells cells
    or have at most max_counts counts or are expressed in at most max_cells cells.

    Only provide one of the optional parameters min_counts, min_cells, max_counts, max_cells per call.
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's filter_genes function
        sc.pp.filter_genes(
            adata,
            min_counts=min_counts,
            min_cells=min_cells,
            max_counts=max_counts,
            max_cells=max_cells,
        )

        # Save the result
        adata.write(output_file)

        click.echo(f"Successfully filtered genes and saved to {output_file}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

import rich_click as click
import scanpy as sc
from scanpy_cli.utils import catch_errors, logger


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
@catch_errors
def filter_cells(min_counts, min_genes, max_counts, max_genes, input_file, output_file):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least min_counts counts or min_genes genes expressed.
    This is to filter measurement outliers, i.e. "unreliable" observations.

    Only provide one of the optional parameters min_counts, min_genes, max_counts, max_genes per call.
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    n_obs_before = adata.n_obs
    logger.debug(
        "Filtering cells: min_counts=%s, min_genes=%s, max_counts=%s, max_genes=%s",
        min_counts,
        min_genes,
        max_counts,
        max_genes,
    )

    sc.pp.filter_cells(
        adata,
        min_counts=min_counts,
        min_genes=min_genes,
        max_counts=max_counts,
        max_genes=max_genes,
    )

    logger.info("Retained %d / %d cells after filtering", adata.n_obs, n_obs_before)

    adata.write(output_file)
    logger.info("Successfully filtered cells and saved to %s", output_file)

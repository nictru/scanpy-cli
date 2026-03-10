import rich_click as click
import scanpy as sc
import sys
from scanpy_cli.utils import logger


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
        adata = sc.read_h5ad(input_file)
        logger.info(
            "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
        )

        n_vars_before = adata.n_vars
        logger.debug(
            "Filtering genes: min_counts=%s, min_cells=%s, max_counts=%s, max_cells=%s",
            min_counts,
            min_cells,
            max_counts,
            max_cells,
        )

        sc.pp.filter_genes(
            adata,
            min_counts=min_counts,
            min_cells=min_cells,
            max_counts=max_counts,
            max_cells=max_cells,
        )

        logger.info(
            "Retained %d / %d genes after filtering", adata.n_vars, n_vars_before
        )

        adata.write(output_file)
        logger.info("Successfully filtered genes and saved to %s", output_file)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

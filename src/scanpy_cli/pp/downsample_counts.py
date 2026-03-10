import rich_click as click
import scanpy as sc
import numpy as np
from scanpy_cli.utils import (
    catch_errors,
    decimals_option,
    round_sparse,
    round_array,
    logger,
)


@click.command()
@decimals_option
@click.option(
    "--counts-per-cell",
    type=int,
    default=None,
    help="Target total counts per cell. Cells with more counts are downsampled; cells with fewer are left unchanged. Mutually exclusive with --total-counts.",
)
@click.option(
    "--total-counts",
    type=int,
    default=None,
    help="Target total counts across all cells. The dataset is downsampled globally. Mutually exclusive with --counts-per-cell.",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed for reproducibility (default: 0).",
)
@click.option(
    "--replace",
    is_flag=True,
    default=False,
    help="Sample with replacement (default: without replacement).",
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
def downsample_counts(
    counts_per_cell,
    total_counts,
    random_state,
    replace,
    input_file,
    output_file,
    decimals,
):
    """Downsample counts from count matrix [Satija et al., 2015].

    Downsample counts such that each cell has no more than counts_per_cell total counts,
    or downsample the entire dataset to total_counts total counts. This is useful to
    normalize sequencing depth across cells or datasets.

    Exactly one of --counts-per-cell or --total-counts must be provided.
    adata.X is modified in-place.
    """
    if counts_per_cell is None and total_counts is None:
        raise click.UsageError(
            "Exactly one of --counts-per-cell or --total-counts must be provided."
        )
    if counts_per_cell is not None and total_counts is not None:
        raise click.UsageError(
            "--counts-per-cell and --total-counts are mutually exclusive."
        )

    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )
    logger.debug(
        "downsample_counts: counts_per_cell=%s, total_counts=%s, random_state=%s, replace=%s",
        counts_per_cell,
        total_counts,
        random_state,
        replace,
    )

    sc.pp.downsample_counts(
        adata,
        counts_per_cell=counts_per_cell,
        total_counts=total_counts,
        random_state=random_state,
        replace=replace,
    )

    logger.info("Downsampling complete")

    if decimals is not None:
        if hasattr(adata.X, "data"):
            adata.X = round_sparse(adata.X, decimals)
        else:
            adata.X = round_array(np.asarray(adata.X), decimals)

    adata.write(output_file)
    logger.info("Successfully downsampled counts and saved to %s", output_file)

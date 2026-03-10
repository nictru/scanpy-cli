import rich_click as click
import scanpy as sc
import numpy as np
from scanpy_cli.utils import (
    catch_errors,
    decimals_option,
    round_array,
    round_sparse,
    logger,
)


@click.command()
@decimals_option
@click.option(
    "--target-sum",
    type=float,
    default=None,
    help="Total count after normalization. If None, normalize to the median total count (default: None).",
)
@click.option(
    "--exclude-highly-expressed",
    is_flag=True,
    default=False,
    help="Exclude highly expressed genes from the normalization. Useful to control effect of very highly expressed genes.",
)
@click.option(
    "--max-fraction",
    type=float,
    default=0.05,
    help="Fraction threshold above which a gene is considered highly expressed (default: 0.05).",
)
@click.option(
    "--key-added",
    type=str,
    default=None,
    help="Name of the field in adata.obs where normalization factors are stored. If None, not stored.",
)
@click.option(
    "--layer",
    type=str,
    default=None,
    help="Layer to normalize. If None, adata.X is normalized in-place.",
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
def normalize_total(
    target_sum,
    exclude_highly_expressed,
    max_fraction,
    key_added,
    layer,
    input_file,
    output_file,
    decimals,
):
    """Normalize counts per cell [Weinreb et al., 2017].

    Normalizes the data matrix to have equal total counts per cell (library-size normalization).
    Typically applied before log1p transformation as part of a standard preprocessing workflow.

    Results stored in adata.obs[key_added] if key_added is set (normalization factors).
    adata.X (or the specified layer) is modified in-place.
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )
    logger.debug(
        "normalize_total: target_sum=%s, exclude_highly_expressed=%s, max_fraction=%s, layer=%s",
        target_sum,
        exclude_highly_expressed,
        max_fraction,
        layer,
    )

    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        key_added=key_added,
        layer=layer,
    )

    logger.info("Normalization complete")

    if decimals is not None:
        target = adata.layers[layer] if layer else adata.X
        if hasattr(target, "data"):
            rounded = round_sparse(target, decimals)
            if layer:
                adata.layers[layer] = rounded
            else:
                adata.X = rounded
        else:
            rounded = round_array(np.asarray(target), decimals)
            if layer:
                adata.layers[layer] = rounded
            else:
                adata.X = rounded
        if key_added and key_added in adata.obs.columns:
            adata.obs[key_added] = np.round(
                adata.obs[key_added].to_numpy().astype(float), decimals
            )

    adata.write(output_file)
    logger.info("Successfully normalized counts and saved to %s", output_file)

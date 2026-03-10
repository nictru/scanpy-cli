import rich_click as click
import scanpy as sc
import numpy as np
from scanpy_cli.utils import catch_errors, decimals_option, round_array, logger


@click.command()
@decimals_option
@click.option(
    "--zero-center/--no-zero-center",
    default=True,
    help="If True, subtract the mean per gene before scaling. If False, only scale to unit variance (default: True).",
)
@click.option(
    "--max-value",
    type=float,
    default=None,
    help="Clip values to this maximum after scaling. If None, no clipping is done (default: None).",
)
@click.option(
    "--layer",
    type=str,
    default=None,
    help="Layer to scale. If None, adata.X is scaled in-place.",
)
@click.option(
    "--mask-obs",
    type=str,
    default=None,
    help="Key in adata.obs for a boolean mask. Only observations where the mask is True are used to compute mean/variance.",
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
def scale(zero_center, max_value, layer, mask_obs, input_file, output_file, decimals):
    """Scale data to unit variance and optionally zero-mean [Pedregosa et al., 2011].

    Scales each gene to unit variance. Optionally zero-centers each gene.
    This is typically the last step in a standard preprocessing workflow
    before dimensionality reduction.

    adata.X (or the specified layer) is modified in-place.
    Mean and standard deviation are stored in adata.var.
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )
    logger.debug(
        "scale: zero_center=%s, max_value=%s, layer=%s, mask_obs=%s",
        zero_center,
        max_value,
        layer,
        mask_obs,
    )

    sc.pp.scale(
        adata,
        zero_center=zero_center,
        max_value=max_value,
        layer=layer,
        mask_obs=mask_obs,
    )

    logger.info("Scaling complete")

    if decimals is not None:
        target = adata.layers[layer] if layer else adata.X
        rounded = round_array(np.asarray(target), decimals)
        if layer:
            adata.layers[layer] = rounded
        else:
            adata.X = rounded

    adata.write(output_file)
    logger.info("Successfully scaled data and saved to %s", output_file)

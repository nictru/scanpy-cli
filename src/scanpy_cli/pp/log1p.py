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
    "--base",
    type=float,
    default=None,
    help="Logarithm base to use. Natural log is used if None (default: None).",
)
@click.option(
    "--chunked",
    is_flag=True,
    default=False,
    help="Process in chunks for memory-efficient computation of large datasets.",
)
@click.option(
    "--chunk-size",
    type=int,
    default=None,
    help="Number of observations per chunk. Only used when --chunked is set.",
)
@click.option(
    "--layer",
    type=str,
    default=None,
    help="Layer to apply log1p to. If None, adata.X is transformed in-place.",
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
def log1p(base, chunked, chunk_size, layer, input_file, output_file, decimals):
    """Logarithmize the data matrix.

    Computes X = log(X + 1), where log denotes the natural logarithm unless
    a different base is specified. Typically applied after normalize_total as
    part of a standard preprocessing workflow.

    adata.X (or the specified layer) is transformed in-place.
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )
    logger.debug(
        "log1p: base=%s, chunked=%s, chunk_size=%s, layer=%s",
        base,
        chunked,
        chunk_size,
        layer,
    )

    sc.pp.log1p(
        adata,
        base=base,
        chunked=chunked or None,
        chunk_size=chunk_size,
        layer=layer,
    )

    logger.info("log1p transformation complete")

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

    adata.write(output_file)
    logger.info("Successfully applied log1p and saved to %s", output_file)

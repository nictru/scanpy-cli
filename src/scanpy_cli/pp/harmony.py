import rich_click as click
import scanpy as sc
import pickle
import numpy as np
import harmonypy
from scanpy_cli.utils import catch_errors, decimals_option, round_array, logger


@click.command()
@decimals_option
@click.option(
    "--key",
    type=str,
    required=True,
    help="Key in adata.obs of the batch annotation to harmonize over.",
)
@click.option(
    "--basis",
    type=str,
    default="X_pca",
    help="Basis in adata.obsm to harmonize (default: 'X_pca').",
)
@click.option(
    "--adjusted-basis",
    type=str,
    default=None,
    help="Key in adata.obsm to store the harmonized embedding.",
)
@click.option(
    "--theta",
    type=float,
    help="Diversity clustering penalty parameter.",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
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
@click.option(
    "--embedding-output",
    type=str,
    help="Optional path to save the harmonized embedding as a pickle file.",
)
@catch_errors
def harmony(
    key,
    basis,
    adjusted_basis,
    theta,
    random_state,
    input_file,
    output_file,
    embedding_output,
    decimals,
):
    """Run Harmony batch correction [Korsunsky et al., 2019].

    Harmony is an algorithm for integrating single-cell data from multiple experiments
    or batches. It projects cells into a shared embedding in which cells group by cell
    type rather than dataset-specific conditions.

    Results are stored in the AnnData object:
    - adata.obsm['{basis}_harmony' | adjusted_basis]: Harmonized embedding
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    out_key = adjusted_basis if adjusted_basis is not None else f"{basis}_harmony"
    logger.debug(
        "Harmony: key=%s, basis=%s, out_key=%s, theta=%s",
        key,
        basis,
        out_key,
        theta,
    )

    kwargs = {}
    if theta is not None:
        kwargs["theta"] = theta

    x = adata.obsm[basis].astype(np.float64)
    harmony_out = harmonypy.run_harmony(
        x, adata.obs, key, random_state=random_state, device="cpu", **kwargs
    )
    adata.obsm[out_key] = harmony_out.Z_corr

    logger.info(
        "Stored embedding in obsm['%s'] with shape %s",
        out_key,
        adata.obsm[out_key].shape,
    )

    if decimals is not None:
        adata.obsm[out_key] = round_array(adata.obsm[out_key], decimals)

    adata.write(output_file)
    logger.info("Successfully ran Harmony integration and saved to %s", output_file)

    if embedding_output:
        embedding = adata.obsm[out_key]
        with open(embedding_output, "wb") as f:
            pickle.dump(embedding, f)
        logger.info("Successfully saved harmonized embedding to %s", embedding_output)

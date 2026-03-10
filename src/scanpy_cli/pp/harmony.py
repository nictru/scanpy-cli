import rich_click as click
import scanpy as sc
import sys
import pickle
import numpy as np
import harmonypy


@click.command()
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
def harmony(
    key,
    basis,
    adjusted_basis,
    theta,
    random_state,
    input_file,
    output_file,
    embedding_output,
):
    """Run Harmony batch correction [Korsunsky et al., 2019].

    Harmony is an algorithm for integrating single-cell data from multiple experiments
    or batches. It projects cells into a shared embedding in which cells group by cell
    type rather than dataset-specific conditions.

    Results are stored in the AnnData object:
    - adata.obsm['{basis}_harmony' | adjusted_basis]: Harmonized embedding
    """
    try:
        adata = sc.read_h5ad(input_file)

        out_key = adjusted_basis if adjusted_basis is not None else f"{basis}_harmony"

        kwargs = {}
        if theta is not None:
            kwargs["theta"] = theta

        # Call harmonypy directly to avoid scanpy wrapper's transpose mismatch
        # with harmonypy >= 0.2.0 (Z_corr already returns (N, d), no extra .T needed)
        x = adata.obsm[basis].astype(np.float64)
        harmony_out = harmonypy.run_harmony(
            x, adata.obs, key, random_state=random_state, device="cpu", **kwargs
        )
        adata.obsm[out_key] = harmony_out.Z_corr

        adata.write(output_file)
        click.echo(f"Successfully ran Harmony integration and saved to {output_file}")

        if embedding_output:
            embedding = adata.obsm[out_key]
            with open(embedding_output, "wb") as f:
                pickle.dump(embedding, f)
            click.echo(f"Successfully saved harmonized embedding to {embedding_output}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

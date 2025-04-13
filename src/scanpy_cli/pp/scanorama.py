import rich_click as click
import scanpy as sc
import scanpy.external as sce
import sys
import pickle


@click.command()
@click.option(
    "--key",
    type=str,
    required=True,
    help="Key in adata.obs that differentiates among experiments/batches.",
)
@click.option(
    "--basis",
    type=str,
    default="X_pca",
    help="The name of the field in adata.obsm where the PCA table is stored (default: 'X_pca').",
)
@click.option(
    "--adjusted-basis",
    type=str,
    default="X_scanorama",
    help="The name of the field in adata.obsm where the integrated embeddings will be stored (default: 'X_scanorama').",
)
@click.option(
    "--knn",
    type=int,
    default=20,
    help="Number of nearest neighbors to use for matching (default: 20).",
)
@click.option(
    "--sigma",
    type=float,
    default=15,
    help="Correction smoothing parameter on Gaussian kernel (default: 15).",
)
@click.option(
    "--approx/--no-approx",
    default=True,
    help="Use approximate nearest neighbors with Python annoy (default: True).",
)
@click.option(
    "--alpha",
    type=float,
    default=0.1,
    help="Alignment score minimum cutoff (default: 0.1).",
)
@click.option(
    "--batch-size",
    type=int,
    default=5000,
    help="The batch size used in the alignment vector computation (default: 5000).",
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
    help="Optional path to save the integrated embedding as a pickle file.",
)
def scanorama(
    key,
    basis,
    adjusted_basis,
    knn,
    sigma,
    approx,
    alpha,
    batch_size,
    input_file,
    output_file,
    embedding_output,
):
    """Run Scanorama integration [Hie et al., 2019].

    Scanorama is an algorithm for integrating single-cell data from multiple experiments.
    This function should be run after performing PCA but before computing the neighbor graph.

    Results are stored in the AnnData object:
    - adata.obsm[adjusted_basis]: Scanorama embeddings such that different experiments are integrated
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's external scanorama_integrate function
        sce.pp.scanorama_integrate(
            adata,
            key=key,
            basis=basis,
            adjusted_basis=adjusted_basis,
            knn=knn,
            sigma=sigma,
            approx=approx,
            alpha=alpha,
            batch_size=batch_size,
        )

        # Save the result
        adata.write(output_file)
        click.echo(f"Successfully ran Scanorama integration and saved to {output_file}")

        # Save embedding as pickle if specified
        if embedding_output:
            embedding = adata.obsm[adjusted_basis]
            with open(embedding_output, "wb") as f:
                pickle.dump(embedding, f)
            click.echo(f"Successfully saved Scanorama embedding to {embedding_output}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

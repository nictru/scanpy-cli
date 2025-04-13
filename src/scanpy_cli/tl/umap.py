import rich_click as click
import scanpy as sc
import sys
import pickle


@click.command()
@click.option(
    "--min-dist",
    type=float,
    default=0.5,
    help="Minimum distance between points in the embedding (default: 0.5).",
)
@click.option(
    "--spread",
    type=float,
    default=1.0,
    help="Scale of the embedded points (default: 1.0).",
)
@click.option(
    "--n-components",
    type=int,
    default=2,
    help="Number of dimensions for the embedding (default: 2).",
)
@click.option(
    "--maxiter",
    type=int,
    help="Number of iterations/epochs for optimization (default: None).",
)
@click.option(
    "--alpha",
    type=float,
    default=1.0,
    help="Initial learning rate for optimization (default: 1.0).",
)
@click.option(
    "--gamma",
    type=float,
    default=1.0,
    help="Weighting applied to negative samples (default: 1.0).",
)
@click.option(
    "--negative-sample-rate",
    type=int,
    default=5,
    help="Number of negative samples per positive sample (default: 5).",
)
@click.option(
    "--init-pos",
    type=str,
    default="spectral",
    help="Initialization method ('paga', 'spectral', 'random', or obsm key) (default: 'spectral').",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
)
@click.option(
    "--a",
    type=float,
    help="More specific parameter controlling the embedding (default: None).",
)
@click.option(
    "--b",
    type=float,
    help="More specific parameter controlling the embedding (default: None).",
)
@click.option(
    "--method",
    type=click.Choice(["umap", "rapids"]),
    default="umap",
    help="UMAP implementation to use (default: 'umap').",
)
@click.option(
    "--key-added",
    type=str,
    help="If specified, the embedding is stored under this key (default: 'X_umap').",
)
@click.option(
    "--neighbors-key",
    type=str,
    default="neighbors",
    help="Key for neighbors settings (default: 'neighbors').",
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
    help="Optional path to save the UMAP embedding as a pickle file.",
)
def umap(
    min_dist,
    spread,
    n_components,
    maxiter,
    alpha,
    gamma,
    negative_sample_rate,
    init_pos,
    random_state,
    a,
    b,
    method,
    key_added,
    neighbors_key,
    input_file,
    output_file,
    embedding_output,
):
    """Embed the neighborhood graph using UMAP [McInnes et al., 2018].

    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning technique
    suitable for visualizing high-dimensional data. It optimizes the embedding to best
    reflect the topology of the data, represented using a neighborhood graph.

    This requires having run neighbors() first.

    Results are stored in the AnnData object:
    - adata.obsm['X_umap' | key_added]: UMAP coordinates
    - adata.uns['umap' | key_added]: UMAP parameters
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's umap function
        sc.tl.umap(
            adata,
            min_dist=min_dist,
            spread=spread,
            n_components=n_components,
            maxiter=maxiter,
            alpha=alpha,
            gamma=gamma,
            negative_sample_rate=negative_sample_rate,
            init_pos=init_pos,
            random_state=random_state,
            a=a,
            b=b,
            method=method,
            key_added=key_added,
            neighbors_key=neighbors_key,
        )

        # Save the result
        adata.write(output_file)
        click.echo(f"Successfully computed UMAP embedding and saved to {output_file}")

        # Save embedding as pickle if specified
        if embedding_output:
            embedding_key = key_added if key_added else "X_umap"
            embedding = adata.obsm[embedding_key]
            with open(embedding_output, "wb") as f:
                pickle.dump(embedding, f)
            click.echo(f"Successfully saved UMAP embedding to {embedding_output}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

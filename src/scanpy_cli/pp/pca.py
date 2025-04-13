import rich_click as click
import scanpy as sc
import sys
import pickle


@click.command()
@click.option(
    "--n-comps",
    type=int,
    help="Number of principal components to compute (default: 50 or 1 - min dimension).",
)
@click.option(
    "--layer",
    type=str,
    help="If provided, use this element of layers for PCA.",
)
@click.option(
    "--zero-center/--no-zero-center",
    default=True,
    help="If True, compute standard PCA from covariance matrix. If False, omit zero-centering variables.",
)
@click.option(
    "--svd-solver",
    type=click.Choice(
        [
            "auto",
            "arpack",
            "randomized",
            "full",
            "tsqr",
            "covariance_eigh",
        ]
    ),
    help="SVD solver to use (default depends on data characteristics).",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
)
@click.option(
    "--mask-var",
    type=str,
    help="Variable mask (e.g., 'highly_variable') to subset the data.",
)
@click.option(
    "--use-highly-variable/--no-use-highly-variable",
    default=None,
    help="Whether to use highly variable genes only (deprecated, use --mask-var).",
)
@click.option(
    "--dtype",
    type=str,
    default="float32",
    help="Numpy data type to which to convert the result (default: 'float32').",
)
@click.option(
    "--chunked/--no-chunked",
    default=False,
    help="If True, perform incremental PCA on segments of chunk_size.",
)
@click.option(
    "--chunk-size",
    type=int,
    help="Number of observations to include in each chunk (required if chunked=True).",
)
@click.option(
    "--key-added",
    type=str,
    help="If specified, the embedding is stored under this key (default: 'X_pca').",
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
    help="Optional path to save the PCA embedding as a pickle file.",
)
def pca(
    n_comps,
    layer,
    zero_center,
    svd_solver,
    random_state,
    mask_var,
    use_highly_variable,
    dtype,
    chunked,
    chunk_size,
    key_added,
    input_file,
    output_file,
    embedding_output,
):
    """Principal component analysis [Pedregosa et al., 2011].

    Computes PCA coordinates, loadings and variance decomposition.
    Uses the implementation of scikit-learn [Pedregosa et al., 2011].

    Results are stored in the AnnData object:
    - adata.obsm['X_pca' | key_added]: PCA representation
    - adata.varm['PCs' | key_added]: Principal components containing the loadings
    - adata.uns['pca' | key_added]['variance_ratio']: Ratio of explained variance
    - adata.uns['pca' | key_added]['variance']: Explained variance
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's pca function
        sc.pp.pca(
            adata,
            n_comps=n_comps,
            layer=layer,
            zero_center=zero_center,
            svd_solver=svd_solver,
            random_state=random_state,
            mask_var=mask_var,
            use_highly_variable=use_highly_variable,
            dtype=dtype,
            chunked=chunked,
            chunk_size=chunk_size,
            key_added=key_added,
        )

        # Save the result
        adata.write(output_file)
        click.echo(f"Successfully computed PCA and saved to {output_file}")

        # Save embedding as pickle if specified
        if embedding_output:
            embedding_key = key_added if key_added else "X_pca"
            embedding = adata.obsm[embedding_key]
            with open(embedding_output, "wb") as f:
                pickle.dump(embedding, f)
            click.echo(f"Successfully saved PCA embedding to {embedding_output}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

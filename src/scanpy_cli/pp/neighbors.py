import rich_click as click
import scanpy as sc
import sys
from scanpy_cli.utils import decimals_option, round_sparse, logger


@click.command()
@decimals_option
@click.option(
    "--n-neighbors",
    "-n",
    type=int,
    default=15,
    help="The number of neighbor cells (default: 15).",
)
@click.option(
    "--n-pcs",
    "-p",
    type=int,
    help="Use this many PCs. If n_pcs=0 use .X if use_rep is None.",
)
@click.option(
    "--use-rep",
    "-r",
    help="Use the indicated representation. 'X' or any key for .obsm is valid.",
)
@click.option(
    "--knn/--no-knn",
    default=True,
    help="If True, use a hard threshold to restrict neighbors (default: True).",
)
@click.option(
    "--method",
    type=click.Choice(["umap", "gauss"]),
    default="umap",
    help="Method for computing connectivities (default: umap).",
)
@click.option(
    "--metric",
    type=click.Choice(
        [
            # Standard metrics
            "cityblock",
            "cosine",
            "euclidean",
            "l1",
            "l2",
            "manhattan",
            # Additional metrics
            "braycurtis",
            "canberra",
            "chebyshev",
            "correlation",
            "dice",
            "hamming",
            "jaccard",
            "kulsinski",
            "mahalanobis",
            "minkowski",
            "rogerstanimoto",
            "russellrao",
            "seuclidean",
            "sokalmichener",
            "sokalsneath",
            "sqeuclidean",
            "yule",
        ]
    ),
    default="euclidean",
    help="Distance metric to use (default: euclidean).",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
)
@click.option(
    "--key-added",
    type=str,
    help="If specified, the neighbors data is stored under this key.",
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
def neighbors(
    n_neighbors,
    n_pcs,
    use_rep,
    knn,
    method,
    metric,
    random_state,
    key_added,
    input_file,
    output_file,
    decimals,
):
    """Compute a neighborhood graph of observations [McInnes et al., 2018].

    Computes the nearest neighbors, distances and connectivities of cells
    in the dataset. The neighbor search efficiency relies on UMAP [McInnes et al., 2018].

    Results are stored in the AnnData object:
    - adata.obsp['distances']: Distance matrix
    - adata.obsp['connectivities']: Weighted adjacency matrix
    - adata.uns['neighbors']: Parameters used
    """
    try:
        adata = sc.read_h5ad(input_file)
        logger.info(
            "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
        )

        logger.debug(
            "Neighbors: n_neighbors=%s, n_pcs=%s, use_rep=%s, method=%s, metric=%s",
            n_neighbors,
            n_pcs,
            use_rep,
            method,
            metric,
        )
        if key_added:
            logger.debug("Storing neighbors under key '%s'", key_added)

        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            use_rep=use_rep,
            knn=knn,
            method=method,
            metric=metric,
            random_state=random_state,
            key_added=key_added,
        )

        conn_key = f"{key_added}_connectivities" if key_added else "connectivities"
        dist_key = f"{key_added}_distances" if key_added else "distances"
        if decimals is not None:
            adata.obsp[conn_key] = round_sparse(adata.obsp[conn_key], decimals)
            adata.obsp[dist_key] = round_sparse(adata.obsp[dist_key], decimals)
        adata.write(output_file)

        logger.info("Successfully computed neighbors and saved to %s", output_file)
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

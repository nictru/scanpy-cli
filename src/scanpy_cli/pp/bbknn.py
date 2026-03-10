import rich_click as click
import scanpy as sc
import scanpy.external as sce
from scanpy_cli.utils import catch_errors, decimals_option, round_sparse, logger


@click.command()
@decimals_option
@click.option(
    "--batch-key",
    type=str,
    default="batch",
    help="Key in adata.obs discriminating between batches (default: 'batch').",
)
@click.option(
    "--use-rep",
    type=str,
    default="X_pca",
    help="The dimensionality reduction in .obsm to use for neighbour detection (default: 'X_pca').",
)
@click.option(
    "--approx/--no-approx",
    default=True,
    help="Use approximate neighbour finding (default: True).",
)
@click.option(
    "--use-annoy/--no-use-annoy",
    default=True,
    help="Use annoy for neighbour finding when approx=True (default: True).",
)
@click.option(
    "--metric",
    type=str,
    default="euclidean",
    help="Distance metric to use (default: 'euclidean').",
)
@click.option(
    "--neighbors-within-batch",
    type=int,
    default=3,
    help="How many top neighbours to report for each batch (default: 3).",
)
@click.option(
    "--n-pcs",
    type=int,
    default=50,
    help="How many dimensions to use in the analysis (default: 50).",
)
@click.option(
    "--trim",
    type=int,
    help="Trim the neighbours of each cell to these many top connectivities.",
)
@click.option(
    "--annoy-n-trees",
    type=int,
    default=10,
    help="Number of trees to construct in the annoy forest (default: 10).",
)
@click.option(
    "--pynndescent-n-neighbors",
    type=int,
    default=30,
    help="Number of neighbours to include in the approximate neighbour graph (default: 30).",
)
@click.option(
    "--pynndescent-random-state",
    type=int,
    default=0,
    help="RNG seed to use when creating the graph (default: 0).",
)
@click.option(
    "--use-faiss/--no-use-faiss",
    default=True,
    help="Use faiss package to compute nearest neighbours when approx=False and metric='euclidean' (default: True).",
)
@click.option(
    "--set-op-mix-ratio",
    type=float,
    default=1.0,
    help="UMAP connectivity computation parameter (default: 1.0).",
)
@click.option(
    "--local-connectivity",
    type=int,
    default=1,
    help="UMAP connectivity computation parameter (default: 1).",
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
def bbknn(
    batch_key,
    use_rep,
    approx,
    use_annoy,
    metric,
    neighbors_within_batch,
    n_pcs,
    trim,
    annoy_n_trees,
    pynndescent_n_neighbors,
    pynndescent_random_state,
    use_faiss,
    set_op_mix_ratio,
    local_connectivity,
    input_file,
    output_file,
    decimals,
):
    """Run BBKNN batch correction [Polański et al., 2019].

    BBKNN (Batch Balanced kNN) alters the kNN procedure to identify each cell's top neighbours
    in each batch separately instead of the entire cell pool with no accounting for batch.
    The nearest neighbours for each batch are then merged to create a final list of neighbours
    for the cell. Aligns batches in a quick and lightweight manner.

    Results are stored in the AnnData object:
    - adata.obsp['connectivities']: Sparse matrix of connectivities
    - adata.obsp['distances']: Sparse matrix of distances
    - adata.uns['neighbors']: Neighbors information
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    logger.debug(
        "BBKNN: batch_key=%s, use_rep=%s, neighbors_within_batch=%s, metric=%s, n_pcs=%s",
        batch_key,
        use_rep,
        neighbors_within_batch,
        metric,
        n_pcs,
    )

    sce.pp.bbknn(
        adata,
        batch_key=batch_key,
        use_rep=use_rep,
        approx=approx,
        use_annoy=use_annoy,
        metric=metric,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=trim,
        annoy_n_trees=annoy_n_trees,
        pynndescent_n_neighbors=pynndescent_n_neighbors,
        pynndescent_random_state=pynndescent_random_state,
        use_faiss=use_faiss,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if decimals is not None:
        adata.obsp["connectivities"] = round_sparse(
            adata.obsp["connectivities"], decimals
        )
        adata.obsp["distances"] = round_sparse(adata.obsp["distances"], decimals)
    adata.write(output_file)
    logger.info("Successfully ran BBKNN integration and saved to %s", output_file)

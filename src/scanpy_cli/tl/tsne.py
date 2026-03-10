import pickle

import rich_click as click
import scanpy as sc
from scanpy_cli.utils import catch_errors, decimals_option, round_array, logger


@click.command()
@decimals_option
@click.option(
    "--n-pcs",
    type=int,
    default=None,
    help="Number of PCs to use. If None, uses all stored PCs (default: None).",
)
@click.option(
    "--n-components",
    type=int,
    default=2,
    help="Number of tSNE dimensions (default: 2).",
)
@click.option(
    "--use-rep",
    type=str,
    default=None,
    help="obsm key to use as feature representation. Overrides --n-pcs (default: None).",
)
@click.option(
    "--perplexity",
    type=float,
    default=30,
    help="Perplexity parameter. Related to the number of nearest neighbors (default: 30).",
)
@click.option(
    "--metric",
    type=str,
    default="euclidean",
    help="Distance metric to use (default: 'euclidean').",
)
@click.option(
    "--early-exaggeration",
    type=float,
    default=12,
    help="Controls the tightness of clusters in early iterations (default: 12).",
)
@click.option(
    "--learning-rate",
    type=float,
    default=1000,
    help="Learning rate for the optimization (default: 1000).",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
)
@click.option(
    "--n-jobs",
    type=int,
    default=None,
    help="Number of parallel jobs. Uses all CPUs if None (default: None).",
)
@click.option(
    "--key-added",
    type=str,
    default=None,
    help="If specified, the embedding is stored under this key (default: 'X_tsne').",
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
    default=None,
    help="Optional path to save the tSNE embedding as a pickle file.",
)
@catch_errors
def tsne(
    n_pcs,
    n_components,
    use_rep,
    perplexity,
    metric,
    early_exaggeration,
    learning_rate,
    random_state,
    n_jobs,
    key_added,
    input_file,
    output_file,
    embedding_output,
    decimals,
):
    """t-distributed stochastic neighborhood embedding [Maaten & Hinton, 2008].

    Embeds cells into a 2D (or higher-dimensional) space for visualization,
    preserving local neighborhood structure. A popular alternative to UMAP
    for single-cell data visualization.

    Results are stored in the AnnData object:
    - adata.obsm['X_tsne' | key_added]: tSNE coordinates
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    embedding_key = key_added if key_added else "X_tsne"
    logger.debug(
        "tsne: n_pcs=%s, n_components=%s, use_rep=%s, perplexity=%s, metric=%s, key_added=%s",
        n_pcs,
        n_components,
        use_rep,
        perplexity,
        metric,
        embedding_key,
    )

    sc.tl.tsne(
        adata,
        n_pcs=n_pcs,
        n_components=n_components,
        use_rep=use_rep,
        perplexity=perplexity,
        metric=metric,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        random_state=random_state,
        n_jobs=n_jobs,
        key_added=key_added,
    )

    logger.info(
        "Stored embedding in obsm['%s'] with shape %s",
        embedding_key,
        adata.obsm[embedding_key].shape,
    )

    if decimals is not None:
        adata.obsm[embedding_key] = round_array(adata.obsm[embedding_key], decimals)

    adata.write(output_file)
    logger.info("Successfully computed tSNE embedding and saved to %s", output_file)

    if embedding_output:
        embedding = adata.obsm[embedding_key]
        with open(embedding_output, "wb") as f:
            pickle.dump(embedding, f)
        logger.info("Successfully saved tSNE embedding to %s", embedding_output)

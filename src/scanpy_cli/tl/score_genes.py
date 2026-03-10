import rich_click as click
import scanpy as sc
import numpy as np
from scanpy_cli.utils import catch_errors, decimals_option, logger


@click.command()
@decimals_option
@click.option(
    "--gene-list",
    required=True,
    type=str,
    help="Comma-separated list of gene names to score. Genes must be present in adata.var_names.",
)
@click.option(
    "--ctrl-as-ref/--no-ctrl-as-ref",
    default=True,
    help="Use control genes as the reference for scoring (default: True).",
)
@click.option(
    "--ctrl-size",
    type=int,
    default=50,
    help="Number of control genes to randomly sample per gene in gene_list (default: 50).",
)
@click.option(
    "--gene-pool",
    type=str,
    default=None,
    help="Comma-separated list of genes to draw control genes from. If None, all genes are used (default: None).",
)
@click.option(
    "--n-bins",
    type=int,
    default=25,
    help="Number of expression bins for control gene sampling (default: 25).",
)
@click.option(
    "--score-name",
    type=str,
    default="score",
    help="Key under which the scores are stored in adata.obs (default: 'score').",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
)
@click.option(
    "--layer",
    type=str,
    default=None,
    help="Layer to use for scoring. If None, adata.X is used.",
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
def score_genes(
    gene_list,
    ctrl_as_ref,
    ctrl_size,
    gene_pool,
    n_bins,
    score_name,
    random_state,
    layer,
    input_file,
    output_file,
    decimals,
):
    """Score a set of genes [Satija et al., 2015].

    Computes a per-cell score for a set of genes using a reference set of randomly
    sampled control genes with similar expression levels. The score is stored in
    adata.obs[score_name].

    Useful for scoring cell cycle genes, marker gene sets, or pathway activity.
    """
    genes = [g.strip() for g in gene_list.split(",")]
    pool = [g.strip() for g in gene_pool.split(",")] if gene_pool else None

    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )
    logger.debug(
        "score_genes: %d genes, ctrl_size=%s, n_bins=%s, score_name=%s, layer=%s",
        len(genes),
        ctrl_size,
        n_bins,
        score_name,
        layer,
    )

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        ctrl_as_ref=ctrl_as_ref,
        ctrl_size=ctrl_size,
        gene_pool=pool,
        n_bins=n_bins,
        score_name=score_name,
        random_state=random_state,
        layer=layer,
    )

    logger.info("Stored gene scores in obs['%s']", score_name)

    if decimals is not None:
        if score_name in adata.obs.columns and np.issubdtype(
            adata.obs[score_name].dtype, np.floating
        ):
            adata.obs[score_name] = np.round(adata.obs[score_name].to_numpy(), decimals)

    adata.write(output_file)
    logger.info("Successfully scored genes and saved to %s", output_file)

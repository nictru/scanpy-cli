import rich_click as click
import scanpy as sc
from scanpy_cli.utils import catch_errors, decimals_option, round_sparse, logger


@click.command()
@decimals_option
@click.option(
    "--groups",
    type=str,
    help="Key in adata.obs for grouping. If None, use the groups computed by clustering.",
)
@click.option(
    "--use-rna-velocity",
    type=bool,
    default=False,
    help="Use RNA velocity to orient edges in the abstracted graph and estimate transitions (default: False).",
)
@click.option(
    "--model",
    type=click.Choice(["v1.2", "v1.0"]),
    default="v1.2",
    help="The PAGA connectivity model to use (default: 'v1.2').",
)
@click.option(
    "--neighbors-key",
    type=str,
    default="neighbors",
    help="Key in adata.uns for neighbors (default: 'neighbors').",
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
def paga(
    groups,
    use_rna_velocity,
    model,
    neighbors_key,
    input_file,
    output_file,
    decimals,
):
    """Run PAGA (Partition-based Graph Abstraction) [Wolf et al., 2019].

    PAGA is a method for analyzing single-cell data that provides a graph-based
    abstraction of the manifold underlying the data. It can be used for:
    - Clustering
    - Trajectory inference
    - RNA velocity analysis
    - Data integration

    Results are stored in the AnnData object:
    - adata.uns['paga']: PAGA results
    - adata.uns['paga']['connectivities']: Sparse matrix with connectivities
    - adata.uns['paga']['connectivities_tree']: Sparse matrix with tree connectivities
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    logger.debug(
        "PAGA: groups=%s, model=%s, use_rna_velocity=%s, neighbors_key=%s",
        groups,
        model,
        use_rna_velocity,
        neighbors_key,
    )

    sc.tl.paga(
        adata,
        groups=groups,
        use_rna_velocity=use_rna_velocity,
        model=model,
        neighbors_key=neighbors_key,
    )

    n_nodes = adata.uns["paga"]["connectivities"].shape[0]
    logger.info("PAGA graph has %d nodes", n_nodes)

    if decimals is not None:
        paga = adata.uns["paga"]
        for key in ("connectivities", "connectivities_tree"):
            if key in paga:
                paga[key] = round_sparse(paga[key], decimals)
    adata.write(output_file)

    logger.info("Successfully ran PAGA and saved to %s", output_file)

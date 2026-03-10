import rich_click as click
import scanpy as sc
import sys
import pickle
from scanpy_cli.utils import logger


@click.command()
@click.option(
    "--resolution",
    type=float,
    default=1.0,
    help="Parameter controlling the coarseness of the clustering (default: 1.0).",
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
    default="leiden",
    help="adata.obs key under which to add the cluster labels (default: 'leiden').",
)
@click.option(
    "--use-weights/--no-weights",
    default=True,
    help="Use edge weights from the graph (default: True).",
)
@click.option(
    "--n-iterations",
    type=int,
    default=-1,
    help="Number of iterations for the Leiden algorithm (default: -1).",
)
@click.option(
    "--neighbors-key",
    type=str,
    help="Use neighbors connectivities as adjacency from this key.",
)
@click.option(
    "--obsp",
    type=str,
    help="Use .obsp[obsp] as adjacency. Can't be used with neighbors-key.",
)
@click.option(
    "--flavor",
    type=click.Choice(["leidenalg", "igraph"]),
    default="leidenalg",
    help="Which package implementation to use (default: 'leidenalg').",
)
@click.option(
    "--directed/--undirected",
    is_flag=True,
    default=None,
    help="Whether to treat the graph as directed or undirected.",
)
@click.option(
    "--restrict-to-key",
    type=str,
    help="Restrict clustering to categories within this observation key.",
)
@click.option(
    "--restrict-to-categories",
    type=str,
    help="Comma-separated list of categories to restrict the clustering to.",
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
    "--clusters-output",
    type=str,
    help="Optional path to save the cluster assignments as a pickle file.",
)
def leiden(
    resolution,
    random_state,
    key_added,
    use_weights,
    n_iterations,
    neighbors_key,
    obsp,
    flavor,
    directed,
    restrict_to_key,
    restrict_to_categories,
    input_file,
    output_file,
    clusters_output,
):
    """Cluster cells into subgroups [Traag et al., 2019].

    Cluster cells using the Leiden algorithm [Traag et al., 2019], an improved version of
    the Louvain algorithm [Blondel et al., 2008]. It was proposed for single-cell
    analysis by Levine et al. [2015].

    This requires having run neighbors() or bbknn() first.

    Results are stored in the AnnData object:
    - adata.obs[key_added]: Cluster assignments
    - adata.uns[key_added]['params']: Parameters used
    """
    try:
        adata = sc.read_h5ad(input_file)
        logger.info(
            "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
        )

        restrict_to = None
        if restrict_to_key and restrict_to_categories:
            categories_list = restrict_to_categories.split(",")
            restrict_to = (restrict_to_key, categories_list)
        elif restrict_to_key or restrict_to_categories:
            raise ValueError(
                "restrict-to-key and restrict-to-categories must be provided together"
            )

        logger.debug(
            "Leiden: resolution=%s, flavor=%s, key_added=%s, n_iterations=%s",
            resolution,
            flavor,
            key_added,
            n_iterations,
        )
        if neighbors_key:
            logger.debug("Using neighbors stored under key '%s'", neighbors_key)
        if obsp:
            logger.debug("Using adjacency from obsp['%s']", obsp)

        sc.tl.leiden(
            adata,
            resolution=resolution,
            random_state=random_state,
            key_added=key_added,
            use_weights=use_weights,
            n_iterations=n_iterations,
            neighbors_key=neighbors_key,
            obsp=obsp,
            flavor=flavor,
            directed=directed,
            restrict_to=restrict_to,
        )

        n_clusters = adata.obs[key_added].nunique()
        logger.info("Found %d clusters (key='%s')", n_clusters, key_added)

        adata.write(output_file)
        logger.info(
            "Successfully computed Leiden clustering with resolution %s and saved to %s",
            resolution,
            output_file,
        )

        if clusters_output:
            clusters = adata.obs[[key_added]]
            with open(clusters_output, "wb") as f:
                pickle.dump(clusters, f)
            logger.info("Successfully saved cluster assignments to %s", clusters_output)

    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--expected-doublet-rate",
    type=float,
    default=0.05,
    help="Where adata_sim not supplied, the estimated doublet rate for the experiment.",
)
@click.option(
    "--stdev-doublet-rate",
    type=float,
    default=0.02,
    help="Where adata_sim not supplied, uncertainty in the expected doublet rate.",
)
@click.option(
    "--sim-doublet-ratio",
    type=float,
    default=2.0,
    help="Number of doublets to simulate relative to the number of observed transcriptomes.",
)
@click.option(
    "--synthetic-doublet-umi-subsampling",
    type=float,
    default=1.0,
    help="Where adata_sim not supplied, rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMI counts from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are added and then randomly sampled at the specified rate.",
)
@click.option(
    "--knn-dist-metric",
    type=click.Choice(
        [
            "cityblock",
            "cosine",
            "euclidean",
            "l1",
            "l2",
            "manhattan",
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
    help="Distance metric used when finding nearest neighbors. For list of valid values, see the documentation for annoy (if use_approx_neighbors is True) or sklearn.neighbors.NearestNeighbors (if use_approx_neighbors is False).",
)
@click.option(
    "--normalize-variance",
    is_flag=True,
    default=True,
    help="If True, normalize the data such that each gene has a variance of 1. sklearn.decomposition.TruncatedSVD will be used for dimensionality reduction, unless mean_center is True.",
)
@click.option(
    "--log-transform",
    is_flag=True,
    default=False,
    help="Whether to use log1p() to log-transform the data prior to PCA.",
)
@click.option(
    "--mean-center",
    is_flag=True,
    default=True,
    help="If True, center the data such that each gene has a mean of 0. sklearn.decomposition.PCA will be used for dimensionality reduction.",
)
@click.option(
    "--n-prin-comps",
    type=int,
    default=30,
    help="Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.",
)
@click.option(
    "--use-approx-neighbors",
    is_flag=True,
    default=None,
    help="Use approximate nearest neighbor method (annoy) for the KNN classifier.",
)
@click.option(
    "--get-doublet-neighbor-parents",
    is_flag=True,
    default=False,
    help="If True, return (in .uns) the parent transcriptomes that generated the doublet neighbors of each observed transcriptome. This information can be used to infer the cell states that generated a given doublet state.",
)
@click.option(
    "--n-neighbors",
    type=int,
    help="Number of neighbors used to construct the KNN graph of observed transcriptomes and simulated doublets. If None, this is automatically set to np.round(0.5 * np.sqrt(n_obs)).",
)
@click.option(
    "--threshold",
    type=float,
    help="Doublet score threshold for calling a transcriptome a doublet. If None, this is set automatically by looking for the minimum between the two modes of the doublet_scores_sim_ histogram. It is best practice to check the threshold visually using the doublet_scores_sim_ histogram and/or based on co-localization of predicted doublets in a 2-D embedding.",
)
@click.option(
    "--batch-key",
    type=str,
    help="Optional obs column name discriminating between batches.",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Initial state for doublet simulation and nearest neighbors.",
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
    "--doublet-output",
    type=str,
    help="Optional path to save doublet predictions and scores as a pickle file.",
)
def scrublet(
    expected_doublet_rate,
    stdev_doublet_rate,
    sim_doublet_ratio,
    synthetic_doublet_umi_subsampling,
    knn_dist_metric,
    normalize_variance,
    log_transform,
    mean_center,
    n_prin_comps,
    use_approx_neighbors,
    get_doublet_neighbor_parents,
    n_neighbors,
    threshold,
    batch_key,
    random_state,
    input_file,
    output_file,
    doublet_output,
):
    """Predict doublets using Scrublet [Wolock et al., 2019].

    Predict cell doublets using a nearest-neighbor classifier of observed transcriptomes and simulated doublets.
    Works best if the input is a raw (unnormalized) counts matrix from a single sample or a collection of similar
    samples from the same experiment.

    Results are stored in the AnnData object:
    - adata.obs['doublet_score']: Doublet scores for each observed transcriptome
    - adata.obs['predicted_doublet']: Boolean indicating predicted doublet status
    - adata.uns['scrublet']['doublet_scores_sim']: Doublet scores for simulated doublets
    - adata.uns['scrublet']['doublet_parents']: Pairs of obs_names used to generate each simulated doublet
    - adata.uns['scrublet']['parameters']: Dictionary of Scrublet parameters
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's scrublet function
        sc.pp.scrublet(
            adata,
            expected_doublet_rate=expected_doublet_rate,
            stdev_doublet_rate=stdev_doublet_rate,
            sim_doublet_ratio=sim_doublet_ratio,
            synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling,
            knn_dist_metric=knn_dist_metric,
            normalize_variance=normalize_variance,
            log_transform=log_transform,
            mean_center=mean_center,
            n_prin_comps=n_prin_comps,
            use_approx_neighbors=use_approx_neighbors,
            get_doublet_neighbor_parents=get_doublet_neighbor_parents,
            n_neighbors=n_neighbors,
            threshold=threshold,
            batch_key=batch_key,
            random_state=random_state,
        )

        # Save the result
        adata.write(output_file)
        click.echo(
            f"Successfully ran Scrublet doublet detection and saved to {output_file}"
        )

        # Save doublet predictions and scores as pickle if specified
        if doublet_output:
            doublet_df = adata.obs[["predicted_doublet", "doublet_score"]].copy()
            doublet_df.columns = [
                "scrublet_predicted_doublet",
                "scrublet_doublet_score",
            ]
            doublet_df.to_pickle(doublet_output)
            click.echo(
                f"Successfully saved doublet predictions and scores to {doublet_output}"
            )

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

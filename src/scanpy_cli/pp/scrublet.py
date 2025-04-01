import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--expected-doublet-rate",
    type=float,
    default=0.05,
    help="Expected doublet rate for the experiment (default: 0.05).",
)
@click.option(
    "--stdev-doublet-rate",
    type=float,
    default=0.02,
    help="Uncertainty in the expected doublet rate (default: 0.02).",
)
@click.option(
    "--sim-doublet-ratio",
    type=float,
    default=2.0,
    help="Number of doublets to simulate relative to the number of observed transcriptomes (default: 2.0).",
)
@click.option(
    "--synthetic-doublet-umi-subsampling",
    type=float,
    default=1.0,
    help="Rate for sampling UMIs when creating synthetic doublets (default: 1.0).",
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
    help="Distance metric used when finding nearest neighbors (default: 'euclidean').",
)
@click.option(
    "--normalize-variance",
    is_flag=True,
    default=True,
    help="Normalize the data such that each gene has a variance of 1 (default: True).",
)
@click.option(
    "--log-transform",
    is_flag=True,
    default=False,
    help="Use log1p() to log-transform the data prior to PCA (default: False).",
)
@click.option(
    "--mean-center",
    is_flag=True,
    default=True,
    help="Center the data such that each gene has a mean of 0 (default: True).",
)
@click.option(
    "--n-prin-comps",
    type=int,
    default=30,
    help="Number of principal components used to embed the transcriptomes (default: 30).",
)
@click.option(
    "--use-approx-neighbors",
    is_flag=True,
    default=False,
    help="Use approximate nearest neighbor method (annoy) for the KNN classifier (default: False).",
)
@click.option(
    "--get-doublet-neighbor-parents",
    is_flag=True,
    default=False,
    help="Return the parent transcriptomes that generated the doublet neighbors (default: False).",
)
@click.option(
    "--n-neighbors",
    type=int,
    help="Number of neighbors used to construct the KNN graph. If None, automatically set to round(0.5 * sqrt(n_obs)).",
)
@click.option(
    "--threshold",
    type=float,
    help="Doublet score threshold for calling a transcriptome a doublet. If None, automatically determined.",
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
    help="Random seed for doublet simulation and nearest neighbors (default: 0).",
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
):
    """Run Scrublet doublet detection [Wolock et al., 2019].

    Scrublet is a method for identifying doublets in single-cell RNA-seq data.
    It works by simulating artificial doublets and comparing their expression
    profiles to those of real cells.

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
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

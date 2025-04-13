import rich_click as click
import scanpy as sc
import sys
import numpy as np
import pickle


@click.command()
@click.option(
    "--n-top-genes",
    type=int,
    help="Number of highly-variable genes to keep. Mandatory if flavor='seurat_v3'.",
)
@click.option(
    "--min-mean",
    type=float,
    default=0.0125,
    help="If n_top_genes is None, this and max_mean are the cutoffs for keeping genes. Ignored if flavor='seurat_v3'.",
)
@click.option(
    "--max-mean",
    type=float,
    default=3,
    help="If n_top_genes is None, this and min_mean are the cutoffs for keeping genes. Ignored if flavor='seurat_v3'.",
)
@click.option(
    "--min-disp",
    type=float,
    default=0.5,
    help="If n_top_genes is None, this is the cutoff for keeping genes. Ignored if flavor='seurat_v3'.",
)
@click.option(
    "--max-disp",
    type=float,
    default=np.inf,
    help="If n_top_genes is None, this is the cutoff for keeping genes. Ignored if flavor='seurat_v3'.",
)
@click.option(
    "--span",
    type=float,
    default=0.3,
    help="The fraction of the data (cells) used when estimating the variance in the loess model fit if flavor='seurat_v3'.",
)
@click.option(
    "--n-bins",
    type=int,
    default=20,
    help="Number of bins for binning the mean gene expression. Normalization is done with respect to each bin.",
)
@click.option(
    "--flavor",
    type=click.Choice(["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]),
    default="seurat",
    help="Choose the flavor for identifying highly variable genes.",
)
@click.option(
    "--batch-key",
    type=str,
    help="If specified, highly-variable genes are selected within each batch separately and merged.",
)
@click.option(
    "--layer",
    type=str,
    help="If provided, use adata.layers[layer] for expression values instead of adata.X.",
)
@click.option(
    "--subset/--no-subset",
    default=False,
    help="Inplace subset to highly-variable genes if True otherwise merely indicate highly variable genes.",
)
@click.option(
    "--inplace/--no-inplace",
    default=True,
    help="Whether to place calculated metrics in .var or return them.",
)
@click.option(
    "--check-values/--no-check-values",
    default=True,
    help="Check if counts in selected layer are integers. A Warning is returned if set to True. Only used if flavor='seurat_v3'/'seurat_v3_paper'.",
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
    "--hvg-output",
    type=str,
    help="Optional path to save the highly variable genes information as a pickle file.",
)
def highly_variable_genes(
    n_top_genes,
    min_mean,
    max_mean,
    min_disp,
    max_disp,
    span,
    n_bins,
    flavor,
    batch_key,
    layer,
    subset,
    inplace,
    check_values,
    input_file,
    output_file,
    hvg_output,
):
    """Annotate highly variable genes [Satija et al., 2015, Stuart et al., 2019, Zheng et al., 2017].

    Expects logarithmized data, except when flavor='seurat_v3'/'seurat_v3_paper', in which count data is expected.

    Results are stored in the AnnData object:
    - adata.var['highly_variable']: boolean indicator of highly-variable genes
    - adata.var['means']: means per gene
    - adata.var['dispersions']: For dispersion-based flavors, dispersions per gene
    - adata.var['dispersions_norm']: For dispersion-based flavors, normalized dispersions per gene
    - adata.var['variances']: For flavor='seurat_v3'/'seurat_v3_paper', variance per gene
    - adata.var['variances_norm']: For flavor='seurat_v3'/'seurat_v3_paper', normalized variance per gene
    - adata.var['highly_variable_rank']: For flavor='seurat_v3'/'seurat_v3_paper', rank of the gene according to normalized variance
    - adata.var['highly_variable_nbatches']: If batch_key is given, this denotes in how many batches genes are detected as HVG
    - adata.var['highly_variable_intersection']: If batch_key is given, this denotes the genes that are highly variable in all batches
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's highly_variable_genes function
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
            max_disp=max_disp,
            span=span,
            n_bins=n_bins,
            flavor=flavor,
            batch_key=batch_key,
            layer=layer,
            subset=subset,
            inplace=inplace,
            check_values=check_values,
        )

        # Save the result
        adata.write(output_file)

        # Save highly variable genes information as pickle if specified
        if hvg_output:
            hvg_info = adata.var[["highly_variable"]]
            with open(hvg_output, "wb") as f:
                pickle.dump(hvg_info, f)
            click.echo(
                f"Successfully saved highly variable genes information to {hvg_output}"
            )

        click.echo(
            f"Successfully detected highly variable genes and saved to {output_file}"
        )
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

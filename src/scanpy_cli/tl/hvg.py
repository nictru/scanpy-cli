import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--n-top-genes",
    type=int,
    help="Number of highly-variable genes to keep. If None, use min_mean and max_mean cutoffs.",
)
@click.option(
    "--min-mean",
    type=float,
    help="If n_top_genes is None, this and max_mean are the cutoffs for keeping genes.",
)
@click.option(
    "--max-mean",
    type=float,
    help="If n_top_genes is None, this and min_mean are the cutoffs for keeping genes.",
)
@click.option(
    "--min-disp",
    type=float,
    help="If n_top_genes is None, this is the cutoff for keeping genes.",
)
@click.option(
    "--batch-key",
    type=str,
    help="Key in adata.obs that differentiates among batches.",
)
@click.option(
    "--layer",
    type=str,
    help="If provided, which element of layers to use.",
)
@click.option(
    "--subset/--no-subset",
    default=False,
    help="Whether to subset the data to highly variable genes (default: False).",
)
@click.option(
    "--inplace/--no-inplace",
    default=True,
    help="Whether to modify the AnnData object inplace (default: True).",
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
def hvg(
    n_top_genes,
    min_mean,
    max_mean,
    min_disp,
    batch_key,
    layer,
    subset,
    inplace,
    input_file,
    output_file,
):
    """Detect highly variable genes [Satija et al., 2015].

    This function identifies genes that show high variability across cells.
    It can be used to reduce the dimensionality of the data by focusing
    on the most informative genes.

    Results are stored in the AnnData object:
    - adata.var['highly_variable']: Boolean indicating highly variable genes
    - adata.var['means']: Mean expression of each gene
    - adata.var['dispersions']: Dispersion of each gene
    - adata.var['dispersions_norm']: Normalized dispersion of each gene
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
            batch_key=batch_key,
            layer=layer,
            subset=subset,
            inplace=inplace,
        )

        # Save the result
        adata.write(output_file)

        click.echo(
            f"Successfully detected highly variable genes and saved to {output_file}"
        )
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

import rich_click as click
import scanpy as sc
import scanpy.external as sce
import sys


@click.command()
@click.option(
    "--key",
    type=str,
    required=True,
    help="Key in adata.obs of the batch annotation to correct over.",
)
@click.option(
    "--layer",
    type=str,
    help="If provided, which element of layers to correct.",
)
@click.option(
    "--covariates",
    type=str,
    help="Additional covariates to adjust for. Can be a comma-separated list of keys in adata.obs.",
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
def combat(
    key,
    layer,
    covariates,
    input_file,
    output_file,
):
    """Run ComBat batch correction [Johnson et al., 2007].

    ComBat is a method for adjusting for batch effects in microarray-based
    expression studies. It uses an empirical Bayes framework to adjust for
    location and scale batch effects.

    Results are stored in the AnnData object:
    - adata.layers['combat']: Batch-corrected data matrix
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Process covariates if provided
        covariates_list = None
        if covariates:
            covariates_list = covariates.split(",")

        # Call scanpy's external combat function
        sce.pp.combat(
            adata,
            key=key,
            layer=layer,
            covariates=covariates_list,
            inplace=True,  # Always True for CLI usage
        )

        # Save the result
        adata.write(output_file)

        click.echo(
            f"Successfully ran ComBat batch correction and saved to {output_file}"
        )
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

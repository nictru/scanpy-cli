import rich_click as click
import scanpy as sc
import sys
import numpy as np


@click.command()
@click.option(
    "--key",
    type=str,
    default="batch",
    help="Key to a categorical annotation from obs that will be used for batch effect removal.",
)
@click.option(
    "--covariates",
    type=str,
    help="Additional covariates besides the batch variable such as adjustment variables or biological condition. Can be a comma-separated list of keys in adata.obs.",
)
@click.option(
    "--in-layer",
    type=str,
    default="X",
    help="Layer to use as input data. If 'X', uses the main data matrix. Otherwise uses the specified layer.",
)
@click.option(
    "--out-layer",
    type=str,
    default="X",
    help="Layer to store the corrected data in. If 'X', updates the main data matrix. Otherwise stores in layers.",
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
    "--corrected-output",
    type=str,
    help="Optional path to save the batch-corrected data as a numpy file.",
)
def combat(
    key,
    covariates,
    in_layer,
    out_layer,
    input_file,
    output_file,
    corrected_output,
):
    """Run ComBat batch correction [Johnson et al., 2006, Leek et al., 2017, Pedersen, 2012].

    Corrects for batch effects by fitting linear models, gains statistical power via an EB framework
    where information is borrowed across genes. This uses the implementation combat.py [Pedersen, 2012].

    Parameters
    ----------
    key : str
        Key to a categorical annotation from obs that will be used for batch effect removal.
    covariates : str, optional
        Additional covariates besides the batch variable such as adjustment variables or biological condition.
        Can be a comma-separated list of keys in adata.obs. This parameter refers to the design matrix X in
        Equation 2.1 in Johnson et al. [2006] and to the mod argument in the original combat function in the sva R package.
        Note that not including covariates may introduce bias or lead to the removal of biological signal in unbalanced designs.
    in_layer : str
        Layer to use as input data. If 'X', uses the main data matrix. Otherwise uses the specified layer.
    out_layer : str
        Layer to store the corrected data in. If 'X', updates the main data matrix. Otherwise stores in layers.
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Store original X if needed
        original_x = None
        if in_layer != "X":
            original_x = adata.X.copy()
            adata.X = adata.layers[in_layer]

        # Process covariates if provided
        covariates_list = None
        if covariates:
            covariates_list = covariates.split(",")

        # Call scanpy's combat function
        corrected = sc.pp.combat(
            adata, key=key, covariates=covariates_list, inplace=False
        )

        if out_layer == "X":
            adata.X = corrected
        else:
            adata.layers[out_layer] = corrected
            if in_layer != "X":
                adata.X = original_x

        # Save the result
        adata.write(output_file)
        click.echo(
            f"Successfully ran ComBat batch correction and saved to {output_file}"
        )

        # Save corrected data as numpy file if specified
        if corrected_output:
            np.save(corrected_output, corrected)
            click.echo(f"Successfully saved batch-corrected data to {corrected_output}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

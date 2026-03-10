import rich_click as click
import scanpy as sc
import sys
import numpy as np
from scanpy_cli.utils import decimals_option, round_array, logger


@click.command()
@decimals_option
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
    decimals,
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
        adata = sc.read_h5ad(input_file)
        logger.info(
            "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
        )

        original_x = None
        if in_layer != "X":
            logger.debug(
                "Using layer '%s' as input (swapping into adata.X temporarily)",
                in_layer,
            )
            original_x = adata.X.copy()
            adata.X = adata.layers[in_layer]

        covariates_list = None
        if covariates:
            covariates_list = covariates.split(",")
            logger.debug("Adjusting for covariates: %s", covariates_list)

        corrected = sc.pp.combat(
            adata, key=key, covariates=covariates_list, inplace=False
        )

        if out_layer == "X":
            adata.X = corrected
            logger.debug("Stored corrected data in adata.X")
        else:
            adata.layers[out_layer] = corrected
            logger.debug("Stored corrected data in layer '%s'", out_layer)
            if in_layer != "X":
                adata.X = original_x

        if decimals is not None:
            if out_layer == "X":
                adata.X = round_array(np.asarray(adata.X), decimals)
            else:
                adata.layers[out_layer] = round_array(
                    np.asarray(adata.layers[out_layer]), decimals
                )

        adata.write(output_file)
        logger.info(
            "Successfully ran ComBat batch correction and saved to %s", output_file
        )

        if corrected_output:
            np.save(corrected_output, corrected)
            logger.info(
                "Successfully saved batch-corrected data to %s", corrected_output
            )
    except Exception as e:
        logger.error(str(e))
        sys.exit(1)

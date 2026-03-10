import rich_click as click
import scanpy as sc
import numpy as np
from scanpy_cli.utils import catch_errors, decimals_option, logger


@click.command()
@decimals_option
@click.option(
    "--qc-vars",
    type=str,
    default=None,
    help="Comma-separated list of boolean obs/var columns to use as QC variable groups (e.g. 'mt,ribo').",
)
@click.option(
    "--percent-top",
    type=str,
    default="50,100,200,500",
    help="Comma-separated integers: compute the percentage of counts in the top N genes (default: '50,100,200,500').",
)
@click.option(
    "--layer",
    type=str,
    default=None,
    help="Layer to use for QC calculation. If None, adata.X is used.",
)
@click.option(
    "--use-raw",
    is_flag=True,
    default=False,
    help="Use adata.raw.X for QC calculation.",
)
@click.option(
    "--no-log1p",
    is_flag=True,
    default=False,
    help="Skip the log1p transformation of QC metrics.",
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
def calculate_qc_metrics(
    qc_vars,
    percent_top,
    layer,
    use_raw,
    no_log1p,
    input_file,
    output_file,
    decimals,
):
    """Calculate quality control metrics [Wolf et al., 2018].

    Calculates a number of qc metrics for an AnnData object, largely based on
    calculateQCMetrics from scater [McCarthy et al., 2017]. Stores results in
    adata.obs and adata.var.

    New obs columns include: total_counts, n_genes_by_counts, pct_counts_in_top_N_genes,
    and for each qc_var: total_counts_{var}, pct_counts_{var}.
    New var columns include: total_counts, n_cells_by_counts, mean_counts, pct_dropout_by_counts.
    """
    adata = sc.read_h5ad(input_file)
    logger.info(
        "Loaded %d cells × %d genes from %s", adata.n_obs, adata.n_vars, input_file
    )

    qc_vars_list = [v.strip() for v in qc_vars.split(",")] if qc_vars else ()
    percent_top_list = [int(x.strip()) for x in percent_top.split(",")]

    logger.debug(
        "calculate_qc_metrics: qc_vars=%s, percent_top=%s, layer=%s, use_raw=%s, log1p=%s",
        qc_vars_list,
        percent_top_list,
        layer,
        use_raw,
        not no_log1p,
    )

    # percent_top values must not exceed the number of variables
    valid_percent_top = [n for n in percent_top_list if n <= adata.n_vars]
    if len(valid_percent_top) < len(percent_top_list):
        skipped = [n for n in percent_top_list if n > adata.n_vars]
        logger.warning(
            "Skipping percent_top values %s that exceed n_vars=%d",
            skipped,
            adata.n_vars,
        )
    percent_top_list = valid_percent_top or None

    obs_cols_before = set(adata.obs.columns)
    var_cols_before = set(adata.var.columns)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars_list,
        percent_top=percent_top_list,
        layer=layer,
        use_raw=use_raw,
        log1p=not no_log1p,
        inplace=True,
    )

    logger.info(
        "Added %d obs columns and %d var columns",
        len(adata.obs.columns) - len(obs_cols_before),
        len(adata.var.columns) - len(var_cols_before),
    )

    if decimals is not None:
        new_obs_cols = set(adata.obs.columns) - obs_cols_before
        new_var_cols = set(adata.var.columns) - var_cols_before
        for col in new_obs_cols:
            if np.issubdtype(adata.obs[col].dtype, np.floating):
                adata.obs[col] = np.round(adata.obs[col].to_numpy(), decimals)
        for col in new_var_cols:
            if np.issubdtype(adata.var[col].dtype, np.floating):
                adata.var[col] = np.round(adata.var[col].to_numpy(), decimals)

    adata.write(output_file)
    logger.info("Successfully calculated QC metrics and saved to %s", output_file)

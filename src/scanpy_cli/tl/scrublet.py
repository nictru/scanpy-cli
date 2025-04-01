import rich_click as click
import scanpy as sc
import scanpy.external as sce
import sys


@click.command()
@click.option(
    "--expected-doublet-rate",
    type=float,
    default=0.1,
    help="Expected doublet rate for the experiment (default: 0.1).",
)
@click.option(
    "--min-gene-count",
    type=int,
    default=2,
    help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_gene_count in more than min_cells are excluded (default: 2).",
)
@click.option(
    "--min-cells",
    type=int,
    default=3,
    help="Used for gene filtering prior to PCA. Genes expressed at fewer than min_gene_count in more than min_cells are excluded (default: 3).",
)
@click.option(
    "--min-gene-variability-pctl",
    type=float,
    default=85,
    help="Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic (default: 85).",
)
@click.option(
    "--n-prin-comps",
    type=int,
    default=30,
    help="Number of principal components to use for doublet detection (default: 30).",
)
@click.option(
    "--random-state",
    type=int,
    default=0,
    help="Random seed (default: 0).",
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
    min_gene_count,
    min_cells,
    min_gene_variability_pctl,
    n_prin_comps,
    random_state,
    input_file,
    output_file,
):
    """Run Scrublet doublet detection [Wolock et al., 2019].

    Scrublet is a method for identifying doublets in single-cell RNA-seq data.
    It works by simulating artificial doublets and comparing their expression
    profiles to those of real cells.

    Results are stored in the AnnData object:
    - adata.obs['doublet_scores']: Doublet scores for each cell
    - adata.obs['predicted_doublets']: Boolean indicating predicted doublets
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's external scrublet function
        sce.pp.scrublet(
            adata,
            expected_doublet_rate=expected_doublet_rate,
            min_gene_count=min_gene_count,
            min_cells=min_cells,
            min_gene_variability_pctl=min_gene_variability_pctl,
            n_prin_comps=n_prin_comps,
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

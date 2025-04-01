import rich_click as click
import scanpy as sc
import sys


@click.command()
@click.option(
    "--groups",
    type=str,
    help="Key in adata.obs for grouping. If None, use the groups computed by clustering.",
)
@click.option(
    "--use-rna-velocity",
    type=bool,
    default=False,
    help="Use RNA velocity to orient edges in the abstracted graph and estimate transitions (default: False).",
)
@click.option(
    "--model",
    type=click.Choice(["v1.2", "v1.0"]),
    default="v1.2",
    help="The PAGA connectivity model to use (default: 'v1.2').",
)
@click.option(
    "--neighbors-key",
    type=str,
    default="neighbors",
    help="Key in adata.uns for neighbors (default: 'neighbors').",
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
def paga(
    groups,
    use_rna_velocity,
    model,
    neighbors_key,
    input_file,
    output_file,
):
    """Run PAGA (Partition-based Graph Abstraction) [Wolf et al., 2019].

    PAGA is a method for analyzing single-cell data that provides a graph-based
    abstraction of the manifold underlying the data. It can be used for:
    - Clustering
    - Trajectory inference
    - RNA velocity analysis
    - Data integration

    Results are stored in the AnnData object:
    - adata.uns['paga']: PAGA results
    - adata.uns['paga']['connectivities']: Sparse matrix with connectivities
    - adata.uns['paga']['connectivities_tree']: Sparse matrix with tree connectivities
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's paga function
        sc.tl.paga(
            adata,
            groups=groups,
            use_rna_velocity=use_rna_velocity,
            model=model,
            neighbors_key=neighbors_key,
        )

        # Save the result
        adata.write(output_file)

        click.echo(f"Successfully ran PAGA and saved to {output_file}")
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

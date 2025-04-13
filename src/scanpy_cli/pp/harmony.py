import rich_click as click
import scanpy as sc
import scanpy.external as sce
import sys
import pickle


@click.command()
@click.option(
    "--key",
    type=str,
    required=True,
    help="Key in adata.obs of the batch annotation to harmonize over.",
)
@click.option(
    "--basis",
    type=str,
    default="X_pca",
    help="Basis in adata.obsm to harmonize (default: 'X_pca').",
)
@click.option(
    "--adjusted-basis",
    type=str,
    help="Key in adata.obsm to store the harmonized embedding.",
)
@click.option(
    "--theta",
    type=float,
    help="Diversity clustering penalty parameter.",
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
@click.option(
    "--embedding-output",
    type=str,
    help="Optional path to save the harmonized embedding as a pickle file.",
)
def harmony(
    key,
    basis,
    adjusted_basis,
    theta,
    random_state,
    input_file,
    output_file,
    embedding_output,
):
    """Run Harmony batch correction [Korsunsky et al., 2019].

    Harmony is an algorithm for integrating single-cell data from multiple experiments
    or batches. It projects cells into a shared embedding in which cells group by cell
    type rather than dataset-specific conditions.

    Results are stored in the AnnData object:
    - adata.obsm['{basis}_harmony' | adjusted_basis]: Harmonized embedding
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Call scanpy's external harmony_integrate function
        sce.pp.harmony_integrate(
            adata,
            key=key,
            basis=basis,
            adjusted_basis=adjusted_basis,
            theta=theta,
            random_state=random_state,
        )

        # Save the result
        adata.write(output_file)
        click.echo(f"Successfully ran Harmony integration and saved to {output_file}")

        # Save embedding as pickle if specified
        if embedding_output:
            embedding_key = adjusted_basis if adjusted_basis else f"{basis}_harmony"
            embedding = adata.obsm[embedding_key]
            with open(embedding_output, "wb") as f:
                pickle.dump(embedding, f)
            click.echo(f"Successfully saved harmonized embedding to {embedding_output}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

import click


@click.command()
def umap():
    """Run UMAP dimensionality reduction."""
    click.echo("Running UMAP...")

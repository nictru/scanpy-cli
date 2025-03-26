import click

@click.group()
def tl():
    """Tool commands for scanpy-cli."""
    pass

@tl.command()
def pca():
    """Run principal component analysis."""
    click.echo("Running PCA...")

@tl.command()
def umap():
    """Run UMAP dimensionality reduction."""
    click.echo("Running UMAP...")

@tl.command()
def clustering():
    """Run clustering analysis."""
    click.echo("Running clustering...") 
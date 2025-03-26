import click

@click.group()
def pl():
    """Plotting commands for scanpy-cli."""
    pass

@pl.command()
def umap():
    """Plot UMAP embedding."""
    click.echo("Plotting UMAP...")

@pl.command()
def heatmap():
    """Plot heatmap."""
    click.echo("Plotting heatmap...")

@pl.command()
def violin():
    """Plot violin plot."""
    click.echo("Plotting violin plot...") 
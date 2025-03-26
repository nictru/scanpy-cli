import click

@click.group()
def pp():
    """Preprocessing commands for scanpy-cli."""
    pass

@pp.command()
def normalize():
    """Normalize data."""
    click.echo("Normalizing data...")

@pp.command()
def filter_cells():
    """Filter cells based on various parameters."""
    click.echo("Filtering cells...")

@pp.command()
def filter_genes():
    """Filter genes based on various parameters."""
    click.echo("Filtering genes...") 
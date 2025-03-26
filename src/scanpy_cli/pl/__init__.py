import rich_click as click

from scanpy_cli.pl.umap import umap


@click.group()
def pl():
    """Plotting commands for scanpy-cli."""
    pass


pl.add_command(umap)

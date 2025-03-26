import click

from scanpy_cli.pl.umap import umap
from scanpy_cli.pl.heatmap import heatmap
from scanpy_cli.pl.violin import violin


@click.group()
def pl():
    """Plotting commands for scanpy-cli."""
    pass


pl.add_command(umap)
pl.add_command(heatmap)
pl.add_command(violin)

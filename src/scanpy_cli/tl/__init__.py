import click

from scanpy_cli.tl.pca import pca
from scanpy_cli.tl.umap import umap
from scanpy_cli.tl.clustering import clustering


@click.group()
def tl():
    """Tools commands for scanpy-cli."""
    pass


tl.add_command(pca)
tl.add_command(umap)
tl.add_command(clustering)

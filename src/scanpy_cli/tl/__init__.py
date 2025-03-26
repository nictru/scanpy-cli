import rich_click as click

from scanpy_cli.tl.umap import umap
from scanpy_cli.tl.leiden import leiden


@click.group()
def tl():
    """Tools commands for scanpy-cli."""
    pass


tl.add_command(umap)
tl.add_command(leiden)

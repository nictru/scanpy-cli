import rich_click as click

from scanpy_cli.tl.umap import umap
from scanpy_cli.tl.leiden import leiden
from scanpy_cli.tl.harmony import harmony
from scanpy_cli.tl.scrublet import scrublet
from scanpy_cli.tl.combat import combat


@click.group()
def tl():
    """Tools commands for scanpy-cli."""
    pass


tl.add_command(umap)
tl.add_command(leiden)
tl.add_command(harmony)
tl.add_command(scrublet)
tl.add_command(combat)

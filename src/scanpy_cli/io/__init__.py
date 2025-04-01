import rich_click as click
from scanpy_cli.io.read_10x_h5 import read_10x_h5


@click.group()
def io():
    """Input/Output operations for single-cell data."""
    pass


io.add_command(read_10x_h5)

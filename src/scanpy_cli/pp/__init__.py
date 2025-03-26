import rich_click as click
from scanpy_cli.pp.regress_out import regress_out
from scanpy_cli.pp.neighbors import neighbors
from scanpy_cli.pp.pca import pca


@click.group()
def pp():
    """Preprocessing commands for scanpy-cli."""
    pass


pp.add_command(regress_out)
pp.add_command(neighbors)
pp.add_command(pca)

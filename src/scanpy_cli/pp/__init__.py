import click

from scanpy_cli.pp.normalize import normalize
from scanpy_cli.pp.filter_cells import filter_cells
from scanpy_cli.pp.filter_genes import filter_genes
from scanpy_cli.pp.regress_out import regress_out
from scanpy_cli.pp.neighbors import neighbors
from scanpy_cli.pp.pca import pca


@click.group()
def pp():
    """Preprocessing commands for scanpy-cli."""
    pass


pp.add_command(normalize)
pp.add_command(filter_cells)
pp.add_command(filter_genes)
pp.add_command(regress_out)
pp.add_command(neighbors)
pp.add_command(pca)

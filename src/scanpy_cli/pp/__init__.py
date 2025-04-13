import rich_click as click
from scanpy_cli.pp.regress_out import regress_out
from scanpy_cli.pp.neighbors import neighbors
from scanpy_cli.pp.pca import pca
from scanpy_cli.pp.highly_variable_genes import highly_variable_genes
from scanpy_cli.pp.harmony import harmony
from scanpy_cli.pp.combat import combat
from scanpy_cli.pp.scrublet import scrublet
from scanpy_cli.pp.filter_genes import filter_genes
from scanpy_cli.pp.filter_cells import filter_cells
from scanpy_cli.pp.bbknn import bbknn
from scanpy_cli.pp.scanorama import scanorama


@click.group()
def pp():
    """Preprocessing commands for scanpy-cli."""
    pass


pp.add_command(scrublet)
pp.add_command(combat)
pp.add_command(harmony)
pp.add_command(highly_variable_genes)
pp.add_command(regress_out)
pp.add_command(neighbors)
pp.add_command(pca)
pp.add_command(filter_genes)
pp.add_command(filter_cells)
pp.add_command(bbknn)
pp.add_command(scanorama)

import rich_click as click

from scanpy_cli.tl.umap import umap
from scanpy_cli.tl.leiden import leiden
from scanpy_cli.tl.paga import paga
from scanpy_cli.tl.rank_genes_groups import rank_genes_groups


@click.group()
def tl():
    """Tools for processing AnnData objects."""
    pass


tl.add_command(umap)
tl.add_command(leiden)
tl.add_command(paga)
tl.add_command(rank_genes_groups)

import click
from scanpy_cli.pp import pp
from scanpy_cli.tl import tl
from scanpy_cli.pl import pl

@click.group()
def cli():
    """Scanpy command line interface for single-cell analysis."""
    pass

cli.add_command(pp)
cli.add_command(tl)
cli.add_command(pl)

def main():
    """Entry point for the scanpy-cli application."""
    cli()

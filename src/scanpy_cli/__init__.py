import rich_click as click
import importlib.metadata
from scanpy_cli.pp import pp
from scanpy_cli.tl import tl
from scanpy_cli.pl import pl
from scanpy_cli.io import io


@click.group()
@click.version_option(
    version=importlib.metadata.version("scanpy-cli"), prog_name="scanpy-cli"
)
def cli():
    """Scanpy command line interface for single-cell analysis."""
    pass


cli.add_command(pp)
cli.add_command(tl)
cli.add_command(pl)
cli.add_command(io)


def main():
    """Entry point for the scanpy-cli application."""
    cli()

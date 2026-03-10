import logging

import anndata
import rich_click as click
import importlib.metadata
from scanpy_cli.pp import pp
from scanpy_cli.tl import tl
from scanpy_cli.pl import pl
from scanpy_cli.io import io
from scanpy_cli.utils import setup_logging

anndata.settings.allow_write_nullable_strings = True


@click.group()
@click.version_option(
    version=importlib.metadata.version("scanpy-cli"), prog_name="scanpy-cli"
)
@click.option("--verbose", "-v", is_flag=True, help="Enable INFO logging.")
@click.option("--debug", is_flag=True, help="Enable DEBUG logging.")
@click.pass_context
def cli(ctx, verbose, debug):
    """Scanpy command line interface for single-cell analysis."""
    if debug:
        setup_logging(logging.DEBUG)
    elif verbose:
        setup_logging(logging.INFO)
    else:
        setup_logging(logging.WARNING)


cli.add_command(pp)
cli.add_command(tl)
cli.add_command(pl)
cli.add_command(io)


def main():
    """Entry point for the scanpy-cli application."""
    cli()

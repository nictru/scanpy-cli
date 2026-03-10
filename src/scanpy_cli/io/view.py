import os

import rich_click as click
import scanpy as sc
from rich.console import Console
from rich.table import Table
from rich import box

from scanpy_cli.utils import catch_errors, logger


def _dtype_str(arr) -> str:
    """Return a short dtype string, handling sparse datasets."""
    try:
        return str(arr.dtype)
    except AttributeError:
        return "?"


def _size_str(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} TB"


@click.command()
@click.argument("input", type=click.Path(exists=True))
@catch_errors
def view(input):
    """Display the structure of an AnnData (.h5ad) file.

    Loads INPUT in backed (read-only) mode and prints a summary of all
    stored arrays, metadata columns, and unstructured annotations without
    loading the full data into memory.

    INPUT is the path to the .h5ad file to inspect.
    """
    logger.debug("Opening %s in backed (read-only) mode", input)

    adata = sc.read_h5ad(input, backed="r")
    console = Console()

    file_size = _size_str(os.path.getsize(input))

    console.print()
    console.print(
        f"[bold cyan]AnnData[/bold cyan]  "
        f"[green]{adata.n_obs:,}[/green] obs × [green]{adata.n_vars:,}[/green] vars  "
        f"[dim]({file_size})[/dim]  [dim]{input}[/dim]"
    )
    console.print()

    # ── X matrix ────────────────────────────────────────────────────────────
    if adata.X is not None:
        x_shape = f"{adata.X.shape[0]:,} × {adata.X.shape[1]:,}"
        x_dtype = _dtype_str(adata.X)
        console.print(f"  [bold]X[/bold]          {x_shape}  [dim]{x_dtype}[/dim]")
        console.print()

    sections = [
        ("obs", adata.obs.dtypes.to_dict(), "Observation annotations"),
        ("var", adata.var.dtypes.to_dict(), "Variable annotations"),
    ]
    for section_name, dtypes_dict, title in sections:
        if dtypes_dict:
            table = Table(
                title=title,
                box=box.SIMPLE_HEAD,
                show_header=True,
                header_style="bold",
                title_style="bold magenta",
                title_justify="left",
                expand=False,
            )
            table.add_column("column", style="cyan", no_wrap=True)
            table.add_column("dtype", style="dim")
            for col, dtype in dtypes_dict.items():
                table.add_row(col, str(dtype))
            console.print(table)

    array_sections = [
        ("obsm", adata.obsm, "Observation embeddings (obsm)"),
        ("varm", adata.varm, "Variable matrices (varm)"),
        ("obsp", adata.obsp, "Observation pairwise (obsp)"),
        ("varp", adata.varp, "Variable pairwise (varp)"),
        ("layers", adata.layers, "Layers"),
    ]
    for section_name, mapping, title in array_sections:
        if len(mapping) == 0:
            continue
        table = Table(
            title=title,
            box=box.SIMPLE_HEAD,
            show_header=True,
            header_style="bold",
            title_style="bold magenta",
            title_justify="left",
            expand=False,
        )
        table.add_column("key", style="cyan", no_wrap=True)
        table.add_column("shape", style="green")
        table.add_column("dtype", style="dim")
        for key in mapping:
            arr = mapping[key]
            shape_str = " × ".join(str(d) for d in arr.shape)
            table.add_row(key, shape_str, _dtype_str(arr))
        console.print(table)

    # ── uns ─────────────────────────────────────────────────────────────────
    if adata.uns:
        table = Table(
            title="Unstructured annotations (uns)",
            box=box.SIMPLE_HEAD,
            show_header=True,
            header_style="bold",
            title_style="bold magenta",
            title_justify="left",
            expand=False,
        )
        table.add_column("key", style="cyan", no_wrap=True)
        table.add_column("type", style="dim")
        for key, val in adata.uns.items():
            table.add_row(key, type(val).__name__)
        console.print(table)

    adata.file.close()
    logger.debug("Closed %s", input)

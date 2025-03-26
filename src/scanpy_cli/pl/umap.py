import rich_click as click
import scanpy as sc
import sys
from pathlib import Path
import matplotlib.pyplot as plt


@click.command()
@click.option(
    "--input-file",
    "-i",
    required=True,
    help="Input h5ad file containing AnnData object.",
)
@click.option(
    "--output-file",
    "-o",
    required=True,
    help="Output file to save the plot. Supported formats: pdf, png, svg.",
)
@click.option(
    "--color",
    multiple=True,
    help="Keys for annotations of observations/cells or variables/genes. Can be specified multiple times.",
)
@click.option(
    "--gene-symbols",
    help="Column name in .var DataFrame that stores gene symbols.",
)
@click.option(
    "--use-raw/--no-use-raw",
    default=None,
    help="Whether to use raw attribute for coloring with gene expression.",
)
@click.option(
    "--layer",
    help="Layer to use for plotting.",
)
@click.option(
    "--dimensions",
    default="0,1",
    help="Dimensions to plot (comma-separated, e.g., '0,1'). Zero-indexed.",
)
@click.option(
    "--projection",
    type=click.Choice(["2d", "3d"]),
    default="2d",
    help="Projection type (default: '2d').",
)
@click.option(
    "--legend-loc",
    type=click.Choice(
        [
            "right margin",
            "on data",
            "none",
            "best",
            "upper right",
            "upper left",
            "lower left",
            "lower right",
            "right",
            "center left",
            "center right",
            "lower center",
            "upper center",
            "center",
        ]
    ),
    default="right margin",
    help="Location of legend (default: 'right margin').",
)
@click.option(
    "--legend-fontsize",
    help="Legend font size.",
)
@click.option(
    "--legend-fontweight",
    default="bold",
    help="Legend font weight (default: 'bold').",
)
@click.option(
    "--legend-fontoutline",
    type=float,
    help="Line width of the legend font outline in pt.",
)
@click.option(
    "--colorbar-loc",
    default="right",
    help="Location of colorbar (default: 'right').",
)
@click.option(
    "--size",
    type=float,
    help="Point size. If None, automatically computed.",
)
@click.option(
    "--color-map",
    help="Color map for continuous variables (e.g., 'viridis', 'magma').",
)
@click.option(
    "--palette",
    help="Color palette for categorical variables.",
)
@click.option(
    "--na-color",
    default="lightgray",
    help="Color for NA/missing values (default: 'lightgray').",
)
@click.option(
    "--na-in-legend/--no-na-in-legend",
    default=True,
    help="Whether to include NA in legend (default: True).",
)
@click.option(
    "--frameon/--no-frameon",
    default=True,
    help="Draw frame around scatter plot (default: True).",
)
@click.option(
    "--title",
    help="Title for the plot.",
)
@click.option(
    "--vmin",
    help="Minimum value for color scale. Can be a percentile (e.g., 'p1.5').",
)
@click.option(
    "--vmax",
    help="Maximum value for color scale. Can be a percentile (e.g., 'p98.5').",
)
@click.option(
    "--vcenter",
    help="Center value for divergent color maps.",
)
@click.option(
    "--add-outline/--no-add-outline",
    default=False,
    help="Add outline to groups of dots (default: False).",
)
@click.option(
    "--outline-color",
    default="black,white",
    help="Colors for outline (comma-separated, default: 'black,white').",
)
@click.option(
    "--outline-width",
    default="0.3,0.05",
    help="Widths for outline (comma-separated, default: '0.3,0.05').",
)
@click.option(
    "--ncols",
    type=int,
    default=4,
    help="Number of panels per row (default: 4).",
)
@click.option(
    "--wspace",
    type=float,
    help="Width space between panels.",
)
@click.option(
    "--hspace",
    type=float,
    default=0.25,
    help="Height space between panels (default: 0.25).",
)
@click.option(
    "--dpi",
    type=int,
    default=100,
    help="Figure DPI (default: 100).",
)
@click.option(
    "--figsize",
    help="Figure size in inches (comma-separated, e.g., '6,4').",
)
def umap(
    input_file,
    output_file,
    color,
    gene_symbols,
    use_raw,
    layer,
    dimensions,
    projection,
    legend_loc,
    legend_fontsize,
    legend_fontweight,
    legend_fontoutline,
    colorbar_loc,
    size,
    color_map,
    palette,
    na_color,
    na_in_legend,
    frameon,
    title,
    vmin,
    vmax,
    vcenter,
    add_outline,
    outline_color,
    outline_width,
    ncols,
    wspace,
    hspace,
    dpi,
    figsize,
):
    """Plot UMAP embeddings.

    Creates scatter plots for UMAP embeddings. This function is a wrapper around scanpy.pl.umap()
    and is designed to save the resulting plot to a file.

    The UMAP embeddings need to be computed first using the 'scanpy-cli tl umap' command.
    """
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(input_file)

        # Process dimensions
        try:
            dimensions = [int(x) for x in dimensions.split(",")]
        except ValueError:
            raise ValueError(
                "--dimensions should be comma-separated integers, e.g., '0,1'"
            )

        # Process outline colors
        if outline_color:
            try:
                outline_color = tuple(outline_color.split(","))
                if len(outline_color) != 2:
                    raise ValueError()
            except ValueError:
                raise ValueError(
                    "--outline-color should be two comma-separated colors, e.g., 'black,white'"
                )

        # Process outline width
        if outline_width:
            try:
                outline_width = tuple(float(x) for x in outline_width.split(","))
                if len(outline_width) != 2:
                    raise ValueError()
            except ValueError:
                raise ValueError(
                    "--outline-width should be two comma-separated floats, e.g., '0.3,0.05'"
                )

        # Process figsize
        if figsize:
            try:
                figsize = tuple(float(x) for x in figsize.split(","))
                if len(figsize) != 2:
                    raise ValueError()
            except ValueError:
                raise ValueError(
                    "--figsize should be two comma-separated floats, e.g., '6,4'"
                )

        # Set up figure parameters
        sc.settings.set_figure_params(dpi=dpi, figsize=figsize)

        # Call scanpy's umap plotting function
        sc.pl.umap(
            adata,
            color=color if color else None,
            gene_symbols=gene_symbols,
            use_raw=use_raw,
            layer=layer,
            dimensions=dimensions,
            projection=projection,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize,
            legend_fontweight=legend_fontweight,
            legend_fontoutline=legend_fontoutline,
            colorbar_loc=colorbar_loc,
            size=size,
            color_map=color_map,
            palette=palette,
            na_color=na_color,
            na_in_legend=na_in_legend,
            frameon=frameon,
            title=title,
            vmin=vmin,
            vmax=vmax,
            vcenter=vcenter,
            add_outline=add_outline,
            outline_color=outline_color,
            outline_width=outline_width,
            ncols=ncols,
            wspace=wspace,
            hspace=hspace,
            show=False,
            save=False,
            return_fig=True,
        )

        # Ensure output directory exists
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save the figure
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        plt.close()

        click.echo(f"Successfully created UMAP plot and saved to {output_file}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        sys.exit(1)

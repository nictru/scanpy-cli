# scanpy-cli

A command-line interface for Scanpy, a Python library for analyzing single-cell gene expression data.

## Installation

```bash
pip install scanpy-cli
```

## Usage

The scanpy-cli tool provides three main command groups:

### Preprocessing (pp)

Commands for preprocessing single-cell data:

```bash
scanpy-cli pp normalize  # Normalize data
scanpy-cli pp filter_cells  # Filter cells
scanpy-cli pp filter_genes  # Filter genes
```

### Tools (tl)

Commands for analysis tools:

```bash
scanpy-cli tl pca  # Run PCA
scanpy-cli tl umap  # Run UMAP
scanpy-cli tl clustering  # Run clustering
```

### Plotting (pl)

Commands for visualization:

```bash
scanpy-cli pl umap  # Plot UMAP
scanpy-cli pl heatmap  # Plot heatmap
scanpy-cli pl violin  # Plot violin plot
```

## Getting Help

For help on any command, use the `--help` flag:

```bash
scanpy-cli --help
scanpy-cli pp --help
scanpy-cli tl pca --help
```

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
scanpy-cli pp regress_out KEYS --input-file INPUT.h5ad --output-file OUTPUT.h5ad  # Regress out unwanted variation
```

Example of regress_out:
```bash
# Regress out cell cycle effects using S_score and G2M_score
scanpy-cli pp regress_out S_score,G2M_score -i data.h5ad -o data_regressed.h5ad

# Regress out with specified parameters
scanpy-cli pp regress_out percent_mito -l counts -j 4 -i data.h5ad -o data_regressed.h5ad

# You can use either long or short parameter names
scanpy-cli pp regress_out percent_mito --layer counts --n-jobs 4 --input-file data.h5ad --output-file data_regressed.h5ad
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
scanpy-cli pp regress_out --help
```

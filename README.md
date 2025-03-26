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
scanpy-cli pp neighbors --input-file INPUT.h5ad --output-file OUTPUT.h5ad  # Compute neighborhood graph
scanpy-cli pp pca --input-file INPUT.h5ad --output-file OUTPUT.h5ad  # Run principal component analysis
```

Example of regress_out:
```bash
# Regress out cell cycle effects using S_score and G2M_score
scanpy-cli pp regress_out S_score,G2M_score -i data.h5ad -o data_regressed.h5ad

# Regress out with specified parameters
scanpy-cli pp regress_out percent_mito -l counts -j 4 -i data.h5ad -o data_regressed.h5ad
```

Example of neighbors:
```bash
# Compute neighbors with default parameters
scanpy-cli pp neighbors -i data.h5ad -o data_neighbors.h5ad

# Compute neighbors with custom parameters
scanpy-cli pp neighbors --n-neighbors 20 --metric cosine -i data.h5ad -o data_neighbors.h5ad
```

Example of PCA:
```bash
# Run PCA with default parameters
scanpy-cli pp pca -i data.h5ad -o data_pca.h5ad

# Run PCA with custom parameters
scanpy-cli pp pca --n-comps 30 --use-highly-variable -i data.h5ad -o data_pca.h5ad
```

### Tools (tl)

Commands for analysis tools:

```bash
scanpy-cli tl umap --input-file INPUT.h5ad --output-file OUTPUT.h5ad  # Run UMAP dimensionality reduction
scanpy-cli tl leiden --input-file INPUT.h5ad --output-file OUTPUT.h5ad  # Run Leiden clustering
```

Example of UMAP:
```bash
# Run UMAP with default parameters
scanpy-cli tl umap -i data_neighbors.h5ad -o data_umap.h5ad

# Run UMAP with custom parameters
scanpy-cli tl umap --min-dist 0.3 --n-components 3 -i data_neighbors.h5ad -o data_umap.h5ad
```

Example of Leiden clustering:
```bash
# Run Leiden clustering with default parameters
scanpy-cli tl leiden -i data_neighbors.h5ad -o data_leiden.h5ad

# Run Leiden with custom resolution
scanpy-cli tl leiden --resolution 0.8 -i data_neighbors.h5ad -o data_leiden.h5ad

# Run Leiden with restricted cell types
scanpy-cli tl leiden --restrict-to-key cell_type --restrict-to-categories "T-cell,B-cell" -i data.h5ad -o data_leiden.h5ad
```

### Plotting (pl)

Commands for visualization:

```bash
scanpy-cli pl umap --input-file INPUT.h5ad --output-file OUTPUT.png  # Plot UMAP embeddings
scanpy-cli pl heatmap  # Plot heatmap
scanpy-cli pl violin  # Plot violin plot
```

Example of UMAP plotting:
```bash
# Basic UMAP plot
scanpy-cli pl umap -i data_umap.h5ad -o umap_plot.png

# UMAP colored by Leiden clusters and gene expression
scanpy-cli pl umap -i data_umap.h5ad -o umap_colored.png --color leiden --color CD4 --color CD8A

# Customized UMAP plot
scanpy-cli pl umap -i data_umap.h5ad -o umap_custom.png --color leiden --dpi 300 --figsize 8,6 --add-outline
```

## Development

### Running Tests

To run the tests, you'll need to install the package with the test dependencies:

```bash
# Install in development mode with test dependencies
pip install -e ".[testing]"

# Run the tests with pytest
pytest
```

If you're using `hatch`, you can run the tests with:

```bash
# Run all tests
hatch run test:test

# Run with coverage
hatch run test:test-cov
```

The tests use a small test dataset that's automatically generated the first time the tests are run.

## Getting Help

For help on any command, use the `--help` flag:

```bash
scanpy-cli --help
scanpy-cli pp --help
scanpy-cli tl umap --help
scanpy-cli pp neighbors --help
```

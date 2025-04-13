# scanpy-cli

A command-line interface for Scanpy, a Python library for analyzing single-cell gene expression data.

## Installation

```bash
pip install scanpy-cli
```

## Usage

The scanpy-cli tool provides three main command groups for single-cell data analysis:

### Preprocessing (pp)

Commands for preprocessing single-cell data:

- `filter_cells`: Filter cells based on counts or genes expressed
- `filter_genes`: Filter genes based on counts or cells expressing them
- `regress_out`: Regress out unwanted sources of variation
- `neighbors`: Compute neighborhood graph
- `pca`: Run principal component analysis
- `combat`: Batch effect correction using ComBat
- `harmony`: Batch effect correction using Harmony
- `scrublet`: Detect doublets in single-cell RNA-seq data
- `highly_variable_genes`: Identify highly variable genes

### Tools (tl)

Commands for analysis tools:

- `umap`: Run UMAP dimensionality reduction
- `leiden`: Run Leiden clustering
- `paga`: Run PAGA for trajectory inference
- `rank_genes_groups`: Find marker genes for clusters

### Plotting (pl)

Commands for visualization:

- `umap`: Plot UMAP embeddings

## Development

### Running Tests

To run the tests, you'll need to install the package with the test dependencies:

```bash
# Install in development mode with test dependencies
pip install -e ".[testing]"

# Run the tests with pytest
pytest
```

## Getting Help

For help on any command, use the `--help` flag:

```bash
scanpy-cli --help
scanpy-cli pp --help
scanpy-cli tl umap --help
```

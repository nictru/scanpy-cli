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

- `filter-cells`: Filter cells based on counts or genes expressed
- `filter-genes`: Filter genes based on counts or cells expressing them
- `normalize-total`: Normalize counts per cell to a target sum (library-size normalization)
- `log1p`: Logarithmize the data matrix (log(X + 1))
- `scale`: Scale data to unit variance and zero mean
- `calculate-qc-metrics`: Compute per-cell and per-gene QC metrics
- `downsample-counts`: Downsample counts to equalize sequencing depth
- `highly-variable-genes`: Identify highly variable genes
- `pca`: Run principal component analysis
- `neighbors`: Compute neighborhood graph
- `regress-out`: Regress out unwanted sources of variation
- `combat`: Batch effect correction using ComBat
- `harmony`: Batch effect correction using Harmony
- `bbknn`: Batch-balanced k-nearest neighbor graph construction
- `scanorama`: Batch effect correction using Scanorama
- `scrublet`: Detect doublets in single-cell RNA-seq data

### Tools (tl)

Commands for analysis tools:

- `tsne`: Run t-SNE dimensionality reduction
- `umap`: Run UMAP dimensionality reduction
- `leiden`: Run Leiden clustering
- `paga`: Run PAGA for trajectory inference
- `rank-genes-groups`: Find marker genes for clusters
- `score-genes`: Score a set of genes per cell

### Plotting (pl)

Commands for visualization:

- `umap`: Plot UMAP embeddings

### Input/Output (io)

Commands for reading and writing data:

- `read-10x-h5`: Read a 10x Genomics HDF5 file and save as `.h5ad`
- `view`: Display the structure of an `.h5ad` file (backed/read-only mode)

## Development

### Running Tests

To run the tests, install the package in development mode with test dependencies:

```bash
# Install in development mode with test dependencies
pip install -e ".[testing]"

# Or using uv
uv sync

# Run the tests with pytest
pytest
```

### Code Quality

The project uses [prek](https://github.com/pre-commit/pre-commit) for code quality checks:

```bash
uv run prek run --all-files
```

This runs `ruff` (linting and formatting) and `pytest`.

## Getting Help

For help on any command, use the `--help` flag:

```bash
scanpy-cli --help
scanpy-cli pp --help
scanpy-cli tl umap --help
```

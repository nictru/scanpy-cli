import scanpy as sc
import subprocess


def test_filter_genes_min_cells(raw_h5ad_path, temp_h5ad_file):
    """Test that the filter_genes command runs successfully with min_cells parameter."""
    # Get the original number of genes
    adata = sc.read_h5ad(raw_h5ad_path)
    original_n_genes = adata.n_vars

    cmd = [
        "scanpy-cli",
        "pp",
        "filter-genes",
        "--input-file",
        str(raw_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--min-cells",
        "10",  # Filter out genes expressed in fewer than 10 cells
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Filter genes command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with filtered genes
    adata_filtered = sc.read_h5ad(temp_h5ad_file)
    assert adata_filtered.n_vars < original_n_genes, "No genes were filtered out"
    assert "n_cells" in adata_filtered.var, "n_cells not found in var"


def test_filter_genes_min_counts(raw_h5ad_path, temp_h5ad_file):
    """Test that the filter_genes command runs successfully with min_counts parameter."""
    # Get the original number of genes
    adata = sc.read_h5ad(raw_h5ad_path)
    original_n_genes = adata.n_vars

    cmd = [
        "scanpy-cli",
        "pp",
        "filter-genes",
        "--input-file",
        str(raw_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--min-counts",
        "100",  # Filter out genes with fewer than 100 total counts
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Filter genes command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with filtered genes
    adata_filtered = sc.read_h5ad(temp_h5ad_file)
    assert adata_filtered.n_vars < original_n_genes, "No genes were filtered out"
    assert "n_counts" in adata_filtered.var, "n_counts not found in var"

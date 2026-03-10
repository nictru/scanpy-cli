import scanpy as sc
import subprocess


def test_filter_cells_min_genes(raw_h5ad_path, temp_h5ad_file):
    """Test that the filter_cells command runs successfully with min_genes parameter."""
    # Get the original number of cells
    adata = sc.read_h5ad(raw_h5ad_path)
    original_n_cells = adata.n_obs

    cmd = [
        "scanpy-cli",
        "pp",
        "filter-cells",
        "--input-file",
        str(raw_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--min-genes",
        "500",  # Filter out cells with fewer than 200 genes expressed
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Filter cells command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with filtered cells
    adata_filtered = sc.read_h5ad(temp_h5ad_file)
    assert adata_filtered.n_obs < original_n_cells, "No cells were filtered out"
    assert "n_genes" in adata_filtered.obs, "n_genes not found in obs"


def test_filter_cells_min_counts(raw_h5ad_path, temp_h5ad_file):
    """Test that the filter_cells command runs successfully with min_counts parameter."""
    # Get the original number of cells
    adata = sc.read_h5ad(raw_h5ad_path)
    original_n_cells = adata.n_obs

    cmd = [
        "scanpy-cli",
        "pp",
        "filter-cells",
        "--input-file",
        str(raw_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--min-counts",
        "1000",  # Filter out cells with fewer than 1000 total counts
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Filter cells command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with filtered cells
    adata_filtered = sc.read_h5ad(temp_h5ad_file)
    assert adata_filtered.n_obs < original_n_cells, "No cells were filtered out"
    assert "n_counts" in adata_filtered.obs, "n_counts not found in obs"

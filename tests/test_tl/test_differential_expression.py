import pytest
import scanpy as sc
import subprocess


@pytest.fixture
def clustered_h5ad_path(test_h5ad_path, temp_h5ad_file):
    """Create a temporary h5ad file with clustering for rank_genes_groups."""
    # Read the test data
    adata = sc.read_h5ad(test_h5ad_path)

    # Compute clustering
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


def test_rank_genes_groups_runs(clustered_h5ad_path, temp_h5ad_file):
    """Test that the rank_genes_groups command runs successfully."""
    cmd = [
        "scanpy-cli",
        "tl",
        "rank-genes-groups",
        "--input-file",
        str(clustered_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--groupby",
        "leiden",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Rank genes groups command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with rank_genes_groups results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert (
        "rank_genes_groups" in adata.uns
    ), "Rank genes groups results not found in uns"

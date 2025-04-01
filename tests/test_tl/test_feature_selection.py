import pytest
import scanpy as sc
import subprocess


@pytest.fixture
def normalized_h5ad_path(temp_h5ad_file):
    """Create a temporary h5ad file with normalized data for HVG."""
    # Create a fresh AnnData object from pbmc3k
    adata = sc.datasets.pbmc3k()

    sc.pp.log1p(adata)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


def test_highly_variable_genes_runs(normalized_h5ad_path, temp_h5ad_file):
    """Test that the highly_variable_genes command runs successfully without batch correction."""
    result = subprocess.run(
        [
            "scanpy-cli",
            "tl",
            "highly-variable-genes",
            "-i",
            normalized_h5ad_path,
            "-o",
            temp_h5ad_file,
        ],
        capture_output=True,
        text=True,
    )

    # Check that the command was successful
    assert result.returncode == 0, f"HVG command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with HVG results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "highly_variable" in adata.var, "HVG results not found in var"


def test_highly_variable_genes_with_batch_runs(batch_h5ad_path, temp_h5ad_file):
    """Test that the highly_variable_genes command runs successfully with batch correction."""
    result = subprocess.run(
        [
            "scanpy-cli",
            "tl",
            "highly-variable-genes",
            "-i",
            batch_h5ad_path,
            "-o",
            temp_h5ad_file,
            "--batch-key",
            "batch",
        ],
        capture_output=True,
        text=True,
    )

    # Check that the command was successful
    assert (
        result.returncode == 0
    ), f"HVG command with batch correction failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with HVG results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "highly_variable" in adata.var, "HVG results not found in var"
    assert (
        "highly_variable_nbatches" in adata.var
    ), "Batch-specific HVG results not found in var"
    assert (
        "highly_variable_intersection" in adata.var
    ), "Intersection HVG results not found in var"

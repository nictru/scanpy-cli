import scanpy as sc
import subprocess


def test_hvg_basic(batch_h5ad_path, temp_h5ad_file):
    """Test that the highly-variable-genes command runs successfully with basic parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "highly-variable-genes",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert (
        result.returncode == 0
    ), f"Highly variable genes command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with highly variable genes
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "highly_variable" in adata.var, "Highly variable genes not found in var"
    assert "means" in adata.var, "Gene means not found in var"
    assert "dispersions" in adata.var, "Gene dispersions not found in var"


def test_hvg_with_parameters(batch_h5ad_path, temp_h5ad_file):
    """Test that the highly-variable-genes command runs successfully with custom parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "highly-variable-genes",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--n-top-genes",
        "2000",
        "--min-mean",
        "0.0125",
        "--max-mean",
        "3",
        "--min-disp",
        "0.5",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert (
        result.returncode == 0
    ), f"Highly variable genes command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with highly variable genes
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "highly_variable" in adata.var, "Highly variable genes not found in var"
    assert "means" in adata.var, "Gene means not found in var"
    assert "dispersions" in adata.var, "Gene dispersions not found in var"


def test_hvg_with_batch(batch_h5ad_path, temp_h5ad_file):
    """Test that the highly-variable-genes command runs successfully with batch information."""
    cmd = [
        "scanpy-cli",
        "pp",
        "highly-variable-genes",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert (
        result.returncode == 0
    ), f"Highly variable genes command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with highly variable genes
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "highly_variable" in adata.var, "Highly variable genes not found in var"
    assert "means" in adata.var, "Gene means not found in var"
    assert "dispersions" in adata.var, "Gene dispersions not found in var"

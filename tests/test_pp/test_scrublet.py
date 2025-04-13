import scanpy as sc
import subprocess


def test_scrublet_basic(raw_batch_h5ad_path, temp_h5ad_file):
    """Test that the scrublet command runs successfully with basic parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "scrublet",
        "--input-file",
        str(raw_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scrublet command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with scrublet scores
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "doublet_score" in adata.obs, "Scrublet scores not found in obs"
    assert "predicted_doublet" in adata.obs, "Doublet predictions not found in obs"


def test_scrublet_with_batch(raw_batch_h5ad_path, temp_h5ad_file):
    """Test that the scrublet command runs successfully with batch information."""
    cmd = [
        "scanpy-cli",
        "pp",
        "scrublet",
        "--input-file",
        str(raw_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scrublet command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with scrublet scores
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "doublet_score" in adata.obs, "Scrublet scores not found in obs"
    assert "predicted_doublet" in adata.obs, "Doublet predictions not found in obs"

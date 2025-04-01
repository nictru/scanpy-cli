import pytest
import scanpy as sc
import subprocess
import numpy as np


@pytest.fixture
def batch_h5ad_path(test_h5ad_path, temp_h5ad_file):
    """Create a temporary h5ad file with batch information for batch correction."""
    # Read the test data
    adata = sc.read_h5ad(test_h5ad_path)

    # Set random seed for reproducibility
    np.random.seed(42)

    # Create batch column with two categories
    n_cells = adata.n_obs
    adata.obs["batch"] = np.random.choice(["batch1", "batch2"], size=n_cells)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


def test_harmony_runs(batch_h5ad_path, temp_h5ad_file):
    """Test that the harmony command runs successfully."""
    cmd = [
        "scanpy-cli",
        "tl",
        "harmony",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--basis",
        "X_pca",
        "--adjusted-basis",
        "X_harmony",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Harmony command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with Harmony results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_harmony" in adata.obsm, "Harmony results not found in obsm"


def test_combat_runs(batch_h5ad_path, temp_h5ad_file):
    """Test that the combat command runs successfully."""
    cmd = [
        "scanpy-cli",
        "tl",
        "combat",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--out-layer",
        "combat_corrected",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"ComBat command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with ComBat results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "combat_corrected" in adata.layers, "ComBat results not found in layers"

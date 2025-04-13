import scanpy as sc
import subprocess
import numpy as np


def test_combat_basic(batch_h5ad_path, temp_h5ad_file):
    """Test that the combat command runs successfully with basic parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "combat",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"ComBat command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object
    adata_corrected = sc.read_h5ad(temp_h5ad_file)
    assert adata_corrected.X is not None, "Data matrix is None after correction"


def test_combat_with_covariates(batch_h5ad_path, temp_h5ad_file):
    """Test that the combat command runs successfully with covariates."""
    # Add a test covariate
    adata = sc.read_h5ad(batch_h5ad_path)
    adata.obs["test_covariate"] = np.random.choice(
        ["A", "B", "C"], size=adata.n_obs
    )  # Binary covariate
    adata.write(batch_h5ad_path)

    cmd = [
        "scanpy-cli",
        "pp",
        "combat",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--covariates",
        "test_covariate",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"ComBat command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object
    adata_corrected = sc.read_h5ad(temp_h5ad_file)
    assert adata_corrected.X is not None, "Data matrix is None after correction"


def test_combat_with_layers(batch_h5ad_path, temp_h5ad_file):
    """Test that the combat command runs successfully with input and output layers."""
    # Add a test layer
    adata = sc.read_h5ad(batch_h5ad_path)
    adata.layers["test_layer"] = adata.X.copy()
    adata.write(batch_h5ad_path)

    cmd = [
        "scanpy-cli",
        "pp",
        "combat",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--in-layer",
        "test_layer",
        "--out-layer",
        "corrected",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"ComBat command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with the corrected layer
    adata_corrected = sc.read_h5ad(temp_h5ad_file)
    assert "corrected" in adata_corrected.layers, "Corrected layer not found"
    assert adata_corrected.layers["corrected"] is not None, "Corrected layer is None"

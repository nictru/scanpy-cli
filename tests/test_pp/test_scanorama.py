import pytest
import scanpy as sc
from pathlib import Path
import tempfile
import subprocess
import pickle
import numpy as np


@pytest.fixture
def contiguous_batch_h5ad_path(batch_h5ad_path):
    """Create a temporary h5ad file with contiguous batches."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    # Load the data
    adata = sc.read_h5ad(batch_h5ad_path)

    # Ensure batches are contiguous by sorting
    adata = adata[adata.obs["batch"].argsort()]

    # Save the prepared data
    adata.write_h5ad(tmp_path)

    yield tmp_path

    # Cleanup after test
    if tmp_path.exists():
        tmp_path.unlink()


def test_scanorama_runs(contiguous_batch_h5ad_path, temp_h5ad_file):
    """Test that the scanorama command runs successfully."""
    cmd = [
        "scanpy-cli",
        "pp",
        "scanorama",
        "--input-file",
        str(contiguous_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--basis",
        "X_pca",
        "--adjusted-basis",
        "X_scanorama",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scanorama command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with Scanorama results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_scanorama" in adata.obsm, "Scanorama results not found in obsm"


def test_scanorama_custom_parameters(contiguous_batch_h5ad_path, temp_h5ad_file):
    """Test Scanorama with custom parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "scanorama",
        "--input-file",
        str(contiguous_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--basis",
        "X_pca",
        "--adjusted-basis",
        "X_scanorama",
        "--knn",
        "30",
        "--sigma",
        "20",
        "--alpha",
        "0.2",
        "--batch-size",
        "2000",
        "--no-approx",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scanorama command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with Scanorama results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_scanorama" in adata.obsm, "Scanorama results not found in obsm"


def test_scanorama_pickle_output(contiguous_batch_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the scanorama command saves the embedding as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "scanorama_embedding.pkl"

    cmd = [
        "scanpy-cli",
        "pp",
        "scanorama",
        "--input-file",
        str(contiguous_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--basis",
        "X_pca",
        "--adjusted-basis",
        "X_scanorama",
        "--embedding-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scanorama command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pickle_embedding = pickle.load(f)

    # Check that the pickle file contains the same embedding as in the AnnData object
    assert np.array_equal(
        adata.obsm["X_scanorama"], pickle_embedding
    ), "Pickle file embedding does not match AnnData embedding"


def test_scanorama_error_handling(contiguous_batch_h5ad_path, temp_h5ad_file):
    """Test Scanorama error handling with invalid parameters."""
    # Test with invalid key
    cmd = [
        "scanpy-cli",
        "pp",
        "scanorama",
        "--input-file",
        str(contiguous_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "nonexistent_batch",
        "--basis",
        "X_pca",
        "--adjusted-basis",
        "X_scanorama",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 1, "Scanorama should fail with invalid key"
    assert "Error" in result.stderr, "Error message not found in stderr"

    # Test with invalid basis
    cmd = [
        "scanpy-cli",
        "pp",
        "scanorama",
        "--input-file",
        str(contiguous_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key",
        "batch",
        "--basis",
        "nonexistent_basis",
        "--adjusted-basis",
        "X_scanorama",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 1, "Scanorama should fail with invalid basis"
    assert "Error" in result.stderr, "Error message not found in stderr"

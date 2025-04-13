import pytest
import scanpy as sc
from pathlib import Path
import tempfile
import subprocess


@pytest.fixture
def pca_h5ad_path(batch_h5ad_path):
    """Create a temporary h5ad file with PCA computed."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    adata = sc.read_h5ad(batch_h5ad_path)
    sc.pp.pca(adata)
    adata.write_h5ad(tmp_path)

    yield tmp_path

    # Cleanup after test
    if tmp_path.exists():
        tmp_path.unlink()


def test_bbknn_runs(pca_h5ad_path, temp_h5ad_file):
    """Test that the bbknn command runs successfully."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"BBKNN command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with BBKNN results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.obsp, "BBKNN connectivities not found in obsp"
    assert "distances" in adata.obsp, "BBKNN distances not found in obsp"
    assert "neighbors" in adata.uns, "BBKNN neighbors not found in uns"


def test_bbknn_custom_parameters(pca_h5ad_path, temp_h5ad_file):
    """Test BBKNN with custom parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
        "--neighbors-within-batch",
        "5",
        "--n-pcs",
        "30",
        "--metric",
        "manhattan",
        "--use-annoy",
        "--annoy-n-trees",
        "20",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"BBKNN command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with BBKNN results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.obsp, "BBKNN connectivities not found in obsp"
    assert "distances" in adata.obsp, "BBKNN distances not found in obsp"
    assert "neighbors" in adata.uns, "BBKNN neighbors not found in uns"


def test_bbknn_pynndescent(pca_h5ad_path, temp_h5ad_file):
    """Test BBKNN with PyNNDescent instead of Annoy."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
        "--no-use-annoy",
        "--pynndescent-n-neighbors",
        "40",
        "--pynndescent-random-state",
        "42",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"BBKNN command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with BBKNN results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.obsp, "BBKNN connectivities not found in obsp"
    assert "distances" in adata.obsp, "BBKNN distances not found in obsp"
    assert "neighbors" in adata.uns, "BBKNN neighbors not found in uns"


def test_bbknn_error_handling(pca_h5ad_path, temp_h5ad_file):
    """Test BBKNN error handling with invalid parameters."""
    # Test with invalid batch key
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "nonexistent_batch",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 1, "BBKNN should fail with invalid batch key"
    assert "Error" in result.stderr, "Error message not found in stderr"

    # Test with invalid use_rep
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--use-rep",
        "nonexistent_rep",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 1, "BBKNN should fail with invalid use_rep"
    assert "Error" in result.stderr, "Error message not found in stderr"


def test_bbknn_trim(pca_h5ad_path, temp_h5ad_file):
    """Test BBKNN with trim parameter."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(pca_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--batch-key",
        "batch",
        "--trim",
        "10",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"BBKNN command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with BBKNN results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.obsp, "BBKNN connectivities not found in obsp"
    assert "distances" in adata.obsp, "BBKNN distances not found in obsp"
    assert "neighbors" in adata.uns, "BBKNN neighbors not found in uns"

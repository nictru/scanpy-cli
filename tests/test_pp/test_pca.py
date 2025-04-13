import scanpy as sc
import subprocess
import pickle
import numpy as np


def test_pca_runs(test_h5ad_path, temp_h5ad_file):
    """Test that the pca command runs successfully."""
    cmd = [
        "scanpy-cli",
        "pp",
        "pca",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"PCA command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with PCA results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_pca" in adata.obsm, "PCA results not found in obsm"
    assert "PCs" in adata.varm, "PCA loadings not found in varm"
    assert "pca" in adata.uns, "PCA parameters not found in uns"


def test_pca_pickle_output(test_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the pca command saves the embedding as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "pca_embedding.pkl"

    cmd = [
        "scanpy-cli",
        "pp",
        "pca",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--embedding-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"PCA command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pickle_embedding = pickle.load(f)

    # Check that the pickle file contains the same embedding as in the AnnData object
    assert np.array_equal(
        adata.obsm["X_pca"], pickle_embedding
    ), "Pickle file embedding does not match AnnData embedding"

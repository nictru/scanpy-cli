import subprocess
import scanpy as sc
import numpy as np
import pickle
import pandas as pd


def test_leiden_basic(test_h5ad_path, temp_h5ad_file):
    """Test basic leiden clustering functionality."""
    cmd = [
        "scanpy-cli",
        "tl",
        "leiden",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"

    # Load the output and verify results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert (
        "leiden" in adata.obs.columns
    ), "Leiden clustering results not found in adata.obs"
    assert len(adata.obs["leiden"].unique()) > 1, "All cells are in the same cluster"


def test_leiden_with_resolution(test_h5ad_path, temp_h5ad_file):
    """Test leiden clustering with different resolution parameter."""
    cmd = [
        "scanpy-cli",
        "tl",
        "leiden",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--resolution",
        "0.5",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"

    # Load the output and verify results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "leiden" in adata.obs.columns
    n_clusters_low_res = len(adata.obs["leiden"].unique())

    # Run with higher resolution
    cmd[8] = "2.0"  # Change resolution to 2.0
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"

    # Verify higher resolution gives more clusters
    adata = sc.read_h5ad(temp_h5ad_file)
    n_clusters_high_res = len(adata.obs["leiden"].unique())
    assert (
        n_clusters_high_res > n_clusters_low_res
    ), "Higher resolution should give more clusters"


def test_leiden_with_random_state(test_h5ad_path, temp_h5ad_file):
    """Test leiden clustering with fixed random state for reproducibility."""
    # Run first time with random state
    cmd = [
        "scanpy-cli",
        "tl",
        "leiden",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--random-state",
        "42",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"
    adata1 = sc.read_h5ad(temp_h5ad_file)
    clusters1 = adata1.obs["leiden"].values

    # Run second time with same random state
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"
    adata2 = sc.read_h5ad(temp_h5ad_file)
    clusters2 = adata2.obs["leiden"].values

    # Verify results are identical
    assert np.array_equal(
        np.array(clusters1), np.array(clusters2)
    ), "Results should be identical with same random state"


def test_leiden_with_key_added(test_h5ad_path, temp_h5ad_file):
    """Test leiden clustering with custom key name."""
    custom_key = "custom_clusters"
    cmd = [
        "scanpy-cli",
        "tl",
        "leiden",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--key-added",
        custom_key,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"

    adata = sc.read_h5ad(temp_h5ad_file)
    assert (
        custom_key in adata.obs.columns
    ), f"Custom key {custom_key} not found in adata.obs"


def test_leiden_pickle_output(test_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the leiden command saves the cluster assignments as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "leiden_clusters.pkl"

    cmd = [
        "scanpy-cli",
        "tl",
        "leiden",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--clusters-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Leiden command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pickle_clusters = pickle.load(f)

    # Check that the pickle file contains the same cluster assignments as in the AnnData object
    assert isinstance(
        pickle_clusters, pd.DataFrame
    ), "Pickle file should contain a DataFrame"
    assert (
        "leiden" in pickle_clusters.columns
    ), "Pickle file should contain 'leiden' column"
    assert np.array_equal(
        adata.obs["leiden"].values, pickle_clusters["leiden"].values
    ), "Pickle file cluster assignments do not match AnnData cluster assignments"

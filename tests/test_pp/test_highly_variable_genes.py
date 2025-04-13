import scanpy as sc
import subprocess
import pickle
import pandas as pd
import numpy as np


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


def test_hvg_pickle_output(batch_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the highly-variable-genes command saves the highly variable genes information as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "hvg_info.pkl"

    cmd = [
        "scanpy-cli",
        "pp",
        "highly-variable-genes",
        "--input-file",
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--hvg-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert (
        result.returncode == 0
    ), f"Highly variable genes command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pickle_hvg = pickle.load(f)

    # Check that the pickle file contains the same highly variable genes information as in the AnnData object
    assert isinstance(
        pickle_hvg, pd.DataFrame
    ), "Pickle file should contain a DataFrame"
    assert (
        "highly_variable" in pickle_hvg.columns
    ), "Pickle file should contain 'highly_variable' column"
    assert np.array_equal(
        adata.var["highly_variable"].values, pickle_hvg["highly_variable"].values
    ), "Pickle file highly variable genes information does not match AnnData information"

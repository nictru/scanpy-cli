import scanpy as sc
import subprocess
import pandas as pd


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


def test_scrublet_pickle_output(raw_batch_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the scrublet command saves doublet predictions and scores as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "doublet_predictions.pkl"

    cmd = [
        "scanpy-cli",
        "pp",
        "scrublet",
        "--input-file",
        str(raw_batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--doublet-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Scrublet command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    pickle_df = pd.read_pickle(pickle_path)

    # Check that the pickle file contains the same doublet predictions and scores as in the AnnData object
    assert (
        "scrublet_predicted_doublet" in pickle_df.columns
    ), "Predicted doublets not found in pickle"
    assert (
        "scrublet_doublet_score" in pickle_df.columns
    ), "Doublet scores not found in pickle"
    assert pickle_df.index.equals(
        adata.obs.index
    ), "Pickle index does not match AnnData obs index"
    assert pickle_df["scrublet_predicted_doublet"].equals(
        adata.obs["predicted_doublet"]
    ), "Pickle predicted doublets do not match AnnData predicted doublets"
    assert pickle_df["scrublet_doublet_score"].equals(
        adata.obs["doublet_score"]
    ), "Pickle doublet scores do not match AnnData doublet scores"

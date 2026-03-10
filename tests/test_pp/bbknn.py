import scanpy as sc
import subprocess


def test_bbknn_runs(batch_h5ad_path, temp_h5ad_file):
    """Test that the bbknn command runs successfully."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
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
    assert result.returncode == 0, f"BBKNN command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with BBKNN results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.obsp, "BBKNN connectivities not found in obsp"
    assert "distances" in adata.obsp, "BBKNN distances not found in obsp"
    assert "neighbors" in adata.uns, "BBKNN neighbors not found in uns"


def test_bbknn_custom_parameters(batch_h5ad_path, temp_h5ad_file):
    """Test BBKNN with custom parameters."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(batch_h5ad_path),
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


def test_bbknn_pynndescent(batch_h5ad_path, temp_h5ad_file):
    """Test BBKNN with PyNNDescent instead of Annoy."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(batch_h5ad_path),
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


def test_bbknn_error_handling(batch_h5ad_path, temp_h5ad_file):
    """Test BBKNN error handling with invalid parameters."""
    # Test with invalid batch key
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(batch_h5ad_path),
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
        str(batch_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--use-rep",
        "nonexistent_rep",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 1, "BBKNN should fail with invalid use_rep"
    assert "Error" in result.stderr, "Error message not found in stderr"


def test_bbknn_trim(batch_h5ad_path, temp_h5ad_file):
    """Test BBKNN with trim parameter."""
    cmd = [
        "scanpy-cli",
        "pp",
        "bbknn",
        "--input-file",
        str(batch_h5ad_path),
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

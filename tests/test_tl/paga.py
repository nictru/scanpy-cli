import scipy.sparse
import scanpy as sc
import subprocess
import numpy as np


def test_paga_decimals(test_h5ad_path, temp_h5ad_file):
    """Test that --decimals rounds the PAGA sparse matrices to the specified number of decimal places."""
    cmd = [
        "scanpy-cli",
        "tl",
        "paga",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--neighbors-key",
        "neighbors",
        "--decimals",
        "3",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"PAGA command failed: {result.stderr}"

    adata = sc.read_h5ad(temp_h5ad_file)
    for key in ("connectivities", "connectivities_tree"):
        mat = adata.uns["paga"][key]
        if scipy.sparse.issparse(mat) and mat.data.size > 0:
            assert np.all(mat.data == np.round(mat.data, 3)), (
                f"paga['{key}'] values are not rounded to 3 decimal places"
            )
    """Test that the paga command runs successfully."""
    cmd = [
        "scanpy-cli",
        "tl",
        "paga",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--neighbors-key",
        "neighbors",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"PAGA command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with PAGA results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "connectivities" in adata.uns["paga"], "PAGA connectivities not found"
    assert "connectivities_tree" in adata.uns["paga"], (
        "PAGA connectivities tree not found"
    )

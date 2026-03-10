import scipy.sparse
import scanpy as sc
import subprocess
import numpy as np


def test_neighbors_runs(test_h5ad_path, temp_h5ad_file):
    """Test that the neighbors command runs successfully."""
    cmd = [
        "scanpy-cli",
        "pp",
        "neighbors",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"neighbors command failed: {result.stderr}"
    assert temp_h5ad_file.exists(), "Output file was not created"

    adata = sc.read_h5ad(temp_h5ad_file)
    assert "neighbors" in adata.uns, "neighbors parameters not found in uns"
    assert "distances" in adata.obsp, "distance matrix not found in obsp"
    assert "connectivities" in adata.obsp, "connectivities matrix not found in obsp"


def test_neighbors_decimals(test_h5ad_path, temp_h5ad_file):
    """Test that --decimals rounds the neighbors output to the specified number of decimal places."""
    cmd = [
        "scanpy-cli",
        "pp",
        "neighbors",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"neighbors command failed: {result.stderr}"

    adata = sc.read_h5ad(temp_h5ad_file)
    conn = adata.obsp["connectivities"]
    dist = adata.obsp["distances"]

    assert scipy.sparse.issparse(conn), "connectivities should be sparse"
    assert scipy.sparse.issparse(dist), "distances should be sparse"
    assert np.all(conn.data == np.round(conn.data, 3)), (
        "connectivities values are not rounded to 3 decimal places"
    )
    assert np.all(dist.data == np.round(dist.data, 3)), (
        "distances values are not rounded to 3 decimal places"
    )
    """Test that the neighbors command respects the --n-neighbors parameter."""
    cmd = [
        "scanpy-cli",
        "pp",
        "neighbors",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--n-neighbors",
        "10",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"neighbors command failed: {result.stderr}"

    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata.uns["neighbors"]["params"]["n_neighbors"] == 10

import scanpy as sc
import subprocess


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


def test_neighbors_n_neighbors(test_h5ad_path, temp_h5ad_file):
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

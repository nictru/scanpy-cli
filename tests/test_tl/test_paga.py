import scanpy as sc
import subprocess


def test_paga_runs(test_h5ad_path, temp_h5ad_file):
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
    assert (
        "connectivities_tree" in adata.uns["paga"]
    ), "PAGA connectivities tree not found"

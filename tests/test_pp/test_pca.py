import scanpy as sc
import subprocess


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

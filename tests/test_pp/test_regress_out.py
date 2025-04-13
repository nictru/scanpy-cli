import scanpy as sc
import subprocess


def test_regress_out_single_key(test_h5ad_path, temp_h5ad_file):
    """Test that the regress_out command runs successfully with a single key."""
    # Add a test column to regress out
    adata = sc.read_h5ad(test_h5ad_path)
    adata.obs["test_key"] = [i % 2 for i in range(adata.n_obs)]  # Binary values
    adata.write(test_h5ad_path)

    cmd = [
        "scanpy-cli",
        "pp",
        "regress-out",
        "test_key",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Regress out command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object
    adata_regressed = sc.read_h5ad(temp_h5ad_file)
    assert (
        adata_regressed.X.shape == adata.X.shape
    ), "Data shape changed after regression"


def test_regress_out_multiple_keys(test_h5ad_path, temp_h5ad_file):
    """Test that the regress_out command runs successfully with multiple keys."""
    # Add test columns to regress out
    adata = sc.read_h5ad(test_h5ad_path)
    adata.obs["test_key1"] = [i % 2 for i in range(adata.n_obs)]  # Binary values
    adata.obs["test_key2"] = [i % 3 for i in range(adata.n_obs)]  # Three categories
    adata.write(test_h5ad_path)

    cmd = [
        "scanpy-cli",
        "pp",
        "regress-out",
        "test_key1,test_key2",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Regress out command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object
    adata_regressed = sc.read_h5ad(temp_h5ad_file)
    assert (
        adata_regressed.X.shape == adata.X.shape
    ), "Data shape changed after regression"

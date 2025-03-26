import pytest
import scanpy as sc
import subprocess
import shutil
from pathlib import Path


@pytest.fixture
def processed_h5ad_path(test_h5ad_path, temp_h5ad_file):
    """Create a temporary h5ad file with neighbors computed for UMAP."""
    # Read the test data
    adata = sc.read_h5ad(test_h5ad_path)

    # Compute neighbors (required for UMAP)
    sc.pp.neighbors(adata)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


def test_umap_runs(processed_h5ad_path, temp_h5ad_file):
    """Test that the umap command runs successfully."""
    # Create a copy of the input file since we'll use it in another test
    input_path = Path(str(processed_h5ad_path) + ".umap_input.h5ad")
    shutil.copy(processed_h5ad_path, input_path)

    cmd = [
        "scanpy-cli",
        "tl",
        "umap",
        "--input-file",
        str(input_path),
        "--output-file",
        str(temp_h5ad_file),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"UMAP command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with UMAP results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_umap" in adata.obsm, "UMAP results not found in obsm"

    # Clean up the extra input file
    input_path.unlink()

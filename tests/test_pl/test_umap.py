import pytest
import scanpy as sc
import subprocess
import shutil
from pathlib import Path


@pytest.fixture
def umap_ready_h5ad_path(test_h5ad_path, temp_h5ad_file):
    """Create a temporary h5ad file with UMAP coordinates for plotting."""
    # Read the test data
    adata = sc.read_h5ad(test_h5ad_path)

    # Compute neighbors and UMAP (required for UMAP plotting)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Add a cluster column for coloring
    sc.tl.leiden(adata)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


def test_umap_plot(umap_ready_h5ad_path):
    """Test that the umap plotting command runs successfully."""
    # Create a copy of the input file
    input_path = Path(str(umap_ready_h5ad_path) + ".plot_input.h5ad")
    shutil.copy(umap_ready_h5ad_path, input_path)

    # Create output path for the plot
    output_path = Path(str(umap_ready_h5ad_path) + ".umap_plot.png")

    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(input_path),
        "--output-file",
        str(output_path),
        "--color",
        "leiden",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"

    # Check that the output file exists
    assert output_path.exists(), "Output plot file was not created"

    # Check that the output file is a valid image file
    assert output_path.stat().st_size > 0, "Output plot file is empty"

    # Clean up
    input_path.unlink()
    output_path.unlink()

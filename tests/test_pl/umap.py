import subprocess
import scanpy as sc
import pytest


def test_umap_plot_basic(test_h5ad_path, temp_plot_file):
    """Test basic UMAP plotting functionality."""
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"

    # Check that the output file exists and is valid
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_multiple_colors(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting with multiple color parameters."""
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--color",
        "n_genes",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_with_custom_dimensions(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting with custom dimensions."""
    # First check the available dimensions
    adata = sc.read_h5ad(test_h5ad_path)
    try:
        umap_coords = adata.obsm["X_umap"]
        if umap_coords is None or not hasattr(umap_coords, "shape"):
            pytest.skip("Test data does not have valid UMAP coordinates")

        n_dims = umap_coords.shape[1]
        if n_dims < 2:
            pytest.skip("Test data does not have enough UMAP dimensions")
    except (KeyError, AttributeError):
        pytest.skip("Test data does not have UMAP coordinates")

    # Use valid dimension indices
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--dimensions",
        "0,1",  # Use first two dimensions
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_with_custom_style(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting with custom style parameters."""
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--size",
        "20",
        "--color-map",
        "viridis",
        "--legend-loc",
        "on data",
        "--frameon",
        "--title",
        "Test UMAP Plot",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_with_outline(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting with outline parameters."""
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--add-outline",
        "--outline-color",
        "black,white",
        "--outline-width",
        "0.3,0.05",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_with_custom_figure_size(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting with custom figure size."""
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--figsize",
        "8,6",
        "--dpi",
        "150",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"
    assert temp_plot_file.exists(), "Output plot file was not created"
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"


def test_umap_plot_error_handling(test_h5ad_path, temp_plot_file):
    """Test UMAP plotting error handling."""
    # Test with invalid dimensions
    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_plot_file),
        "--color",
        "louvain",
        "--dimensions",
        "invalid",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Command should fail with invalid dimensions"

    # Test with invalid color
    cmd[8] = "nonexistent_column"
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Command should fail with invalid color column"

    # Test with invalid figure size
    cmd[8] = "louvain"
    cmd.extend(["--figsize", "invalid"])
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Command should fail with invalid figure size"

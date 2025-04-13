import subprocess


def test_umap_plot(test_h5ad_path, temp_plot_file):
    """Test that the umap plotting command runs successfully."""
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

    # Check that the output file exists
    assert temp_plot_file.exists(), "Output plot file was not created"

    # Check that the output file is a valid image file
    assert temp_plot_file.stat().st_size > 0, "Output plot file is empty"

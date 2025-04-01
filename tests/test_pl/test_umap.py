import subprocess
from pathlib import Path


def test_umap_plot(test_h5ad_path):
    """Test that the umap plotting command runs successfully."""
    # Create output path for the plot
    output_path = Path(str(test_h5ad_path) + ".umap_plot.png")

    cmd = [
        "scanpy-cli",
        "pl",
        "umap",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(output_path),
        "--color",
        "louvain",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"UMAP plot command failed: {result.stderr}"

    # Check that the output file exists
    assert output_path.exists(), "Output plot file was not created"

    # Check that the output file is a valid image file
    assert output_path.stat().st_size > 0, "Output plot file is empty"

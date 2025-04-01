import scanpy as sc
import subprocess


def test_rank_genes_groups_runs(test_h5ad_path, temp_h5ad_file):
    """Test that the rank_genes_groups command runs successfully."""
    cmd = [
        "scanpy-cli",
        "tl",
        "rank-genes-groups",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--groupby",
        "louvain",
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Rank genes groups command failed: {result.stderr}"

    # Check that the output file exists
    assert temp_h5ad_file.exists(), "Output file was not created"

    # Check that the output file is a valid AnnData object with rank_genes_groups results
    adata = sc.read_h5ad(temp_h5ad_file)
    assert (
        "rank_genes_groups" in adata.uns
    ), "Rank genes groups results not found in uns"

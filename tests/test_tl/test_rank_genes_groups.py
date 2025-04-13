import scanpy as sc
import subprocess
import pickle
import numpy as np


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


def test_rank_genes_groups_pickle_output(test_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the rank_genes_groups command saves the dictionary as a pickle file when requested."""
    # Create a temporary pickle file path
    pickle_path = tmp_path / "rank_genes_groups.pkl"

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
        "--rank-genes-output",
        str(pickle_path),
    ]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Rank genes groups command failed: {result.stderr}"

    # Check that both output files exist
    assert temp_h5ad_file.exists(), "Output h5ad file was not created"
    assert pickle_path.exists(), "Output pickle file was not created"

    # Load both the AnnData object and the pickle file
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pickle_data = pickle.load(f)

    # Check that the pickle file contains the same rank_genes_groups dictionary as in the AnnData object
    assert isinstance(pickle_data, dict), "Pickle file should contain a dictionary"
    assert "names" in pickle_data, "Gene names not found in pickle"
    assert "scores" in pickle_data, "Scores not found in pickle"
    assert "pvals" in pickle_data, "P-values not found in pickle"
    assert "pvals_adj" in pickle_data, "Adjusted p-values not found in pickle"

    # Compare the dictionaries properly, handling numpy arrays and nested dictionaries
    adata_dict = adata.uns["rank_genes_groups"]
    for key in adata_dict:
        if key == "params":
            # For params dictionary, only compare keys that exist in both
            for param_key in adata_dict[key]:
                if param_key in pickle_data[key]:
                    assert (
                        adata_dict[key][param_key] == pickle_data[key][param_key]
                    ), f"Parameter '{param_key}' does not match"
        elif isinstance(adata_dict[key], np.ndarray):
            assert np.array_equal(
                adata_dict[key], pickle_data[key]
            ), f"Arrays for key '{key}' do not match"
        else:
            assert (
                adata_dict[key] == pickle_data[key]
            ), f"Values for key '{key}' do not match"

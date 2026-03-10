import scanpy as sc
import subprocess
import numpy as np


def test_score_genes_runs(test_h5ad_path, temp_h5ad_file):
    """Use the first few gene names from the dataset."""
    adata = sc.read_h5ad(test_h5ad_path)
    gene_list = ",".join(list(adata.var_names[:5]))

    cmd = [
        "scanpy-cli",
        "tl",
        "score-genes",
        "--gene-list",
        gene_list,
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"score-genes failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata_out = sc.read_h5ad(temp_h5ad_file)
    assert "score" in adata_out.obs.columns


def test_score_genes_custom_score_name(test_h5ad_path, temp_h5ad_file):
    adata = sc.read_h5ad(test_h5ad_path)
    gene_list = ",".join(list(adata.var_names[:5]))

    cmd = [
        "scanpy-cli",
        "tl",
        "score-genes",
        "--gene-list",
        gene_list,
        "--score-name",
        "my_score",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"score-genes failed: {result.stderr}"
    adata_out = sc.read_h5ad(temp_h5ad_file)
    assert "my_score" in adata_out.obs.columns


def test_score_genes_decimals(test_h5ad_path, temp_h5ad_file):
    adata = sc.read_h5ad(test_h5ad_path)
    gene_list = ",".join(list(adata.var_names[:5]))

    cmd = [
        "scanpy-cli",
        "tl",
        "score-genes",
        "--gene-list",
        gene_list,
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"score-genes failed: {result.stderr}"
    adata_out = sc.read_h5ad(temp_h5ad_file)
    scores = adata_out.obs["score"].to_numpy()
    assert np.all(scores == np.round(scores, 3)), "Scores not rounded to 3 decimals"


def test_score_genes_missing_gene(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "tl",
        "score-genes",
        "--gene-list",
        "NOTAREALGENE12345",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail with genes not found in adata"

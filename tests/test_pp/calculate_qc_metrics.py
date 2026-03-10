import scanpy as sc
import subprocess
import numpy as np
import scipy.sparse
import pytest
import anndata


@pytest.fixture
def raw_counts_h5ad(tmp_path):
    """Create a small h5ad with raw integer counts and a mitochondrial gene flag."""
    rng = np.random.default_rng(0)
    n_obs, n_vars = 50, 30
    counts = scipy.sparse.random(
        n_obs, n_vars, density=0.5, random_state=rng, format="csr"
    )
    counts.data = (counts.data * 100 + 1).astype(np.float32)
    import pandas as pd

    var = pd.DataFrame(
        {"mt": [i < 5 for i in range(n_vars)]},
        index=[f"gene_{i}" for i in range(n_vars)],
    )
    adata = anndata.AnnData(counts, var=var)
    path = tmp_path / "raw_counts.h5ad"
    adata.write_h5ad(path)
    return path


def test_calculate_qc_metrics_runs(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "calculate-qc-metrics",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"calculate-qc-metrics failed: {result.stderr}"
    assert temp_h5ad_file.exists()


def test_calculate_qc_metrics_adds_obs_columns(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "calculate-qc-metrics",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"calculate-qc-metrics failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "total_counts" in adata.obs.columns
    assert "n_genes_by_counts" in adata.obs.columns


def test_calculate_qc_metrics_adds_var_columns(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "calculate-qc-metrics",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"calculate-qc-metrics failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "total_counts" in adata.var.columns
    assert "n_cells_by_counts" in adata.var.columns


def test_calculate_qc_metrics_qc_vars(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "calculate-qc-metrics",
        "--qc-vars",
        "mt",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"calculate-qc-metrics failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "pct_counts_mt" in adata.obs.columns


def test_calculate_qc_metrics_decimals(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "calculate-qc-metrics",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"calculate-qc-metrics failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    total = adata.obs["total_counts"].to_numpy()
    assert np.all(total == np.round(total, 3)), "total_counts not rounded to 3 decimals"

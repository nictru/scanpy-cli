import scanpy as sc
import subprocess
import numpy as np
import scipy.sparse
import pytest
import anndata


@pytest.fixture
def raw_counts_h5ad(tmp_path):
    """Create a small h5ad with raw integer counts suitable for normalization."""
    import numpy.random as npr

    rng = npr.default_rng(0)
    n_obs, n_vars = 50, 30
    counts = scipy.sparse.random(
        n_obs, n_vars, density=0.5, random_state=rng, format="csr"
    )
    counts.data = np.round(counts.data * 100).astype(np.float32) + 1
    adata = anndata.AnnData(counts)
    path = tmp_path / "raw_counts.h5ad"
    adata.write_h5ad(path)
    return path


def test_normalize_total_runs(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "normalize-total",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"normalize-total failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata.n_obs > 0


def test_normalize_total_target_sum(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "normalize-total",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--target-sum",
        "1e4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"normalize-total failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    totals = np.asarray(adata.X.sum(axis=1)).flatten()
    assert np.allclose(totals, 1e4, rtol=1e-3), "Cell totals should be ~1e4"


def test_normalize_total_key_added(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "normalize-total",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--key-added",
        "norm_factor",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"normalize-total failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "norm_factor" in adata.obs.columns


def test_normalize_total_decimals(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "normalize-total",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"normalize-total failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    data = (
        adata.X.data
        if scipy.sparse.issparse(adata.X)
        else np.asarray(adata.X).flatten()
    )
    assert np.all(data == np.round(data, 3)), "X values are not rounded to 3 decimals"

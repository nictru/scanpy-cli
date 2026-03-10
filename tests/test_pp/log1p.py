import scanpy as sc
import subprocess
import numpy as np
import scipy.sparse
import pytest
import anndata


@pytest.fixture
def raw_counts_h5ad(tmp_path):
    """Create a small h5ad with raw integer counts suitable for log1p."""
    rng = np.random.default_rng(0)
    n_obs, n_vars = 50, 30
    counts = scipy.sparse.random(
        n_obs, n_vars, density=0.5, random_state=rng, format="csr"
    )
    counts.data = (counts.data * 100 + 1).astype(np.float32)
    adata = anndata.AnnData(counts)
    path = tmp_path / "raw_counts.h5ad"
    adata.write_h5ad(path)
    return path


def test_log1p_runs(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "log1p",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"log1p failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata.n_obs > 0


def test_log1p_values(raw_counts_h5ad, temp_h5ad_file):
    """Values should be approximately log(x+1) of the input."""
    adata_in = sc.read_h5ad(raw_counts_h5ad)
    cmd = [
        "scanpy-cli",
        "pp",
        "log1p",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"log1p failed: {result.stderr}"
    adata_out = sc.read_h5ad(temp_h5ad_file)
    # All non-zero values should be less than the original (log compresses)
    x_in = np.asarray(adata_in.X.todense())
    x_out = np.asarray(adata_out.X.todense())
    mask = x_in > 1
    assert np.all(x_out[mask] < x_in[mask]), "log1p should reduce values > 1"


def test_log1p_decimals(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "log1p",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"log1p failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    data = (
        adata.X.data
        if scipy.sparse.issparse(adata.X)
        else np.asarray(adata.X).flatten()
    )
    assert np.all(data == np.round(data, 3)), "X values are not rounded to 3 decimals"


def test_log1p_base(raw_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "log1p",
        "--input-file",
        str(raw_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--base",
        "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"log1p with --base failed: {result.stderr}"

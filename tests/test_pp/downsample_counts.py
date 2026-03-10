import scanpy as sc
import subprocess
import numpy as np
import scipy.sparse
import pytest
import anndata


@pytest.fixture
def integer_counts_h5ad(tmp_path):
    """Create a small h5ad with integer count data suitable for downsampling."""
    rng = np.random.default_rng(0)
    n_obs, n_vars = 40, 20
    counts = scipy.sparse.random(
        n_obs, n_vars, density=0.6, random_state=rng, format="csr"
    )
    counts.data = np.ceil(counts.data * 200).astype(np.float32)
    adata = anndata.AnnData(counts)
    path = tmp_path / "integer_counts.h5ad"
    adata.write_h5ad(path)
    return path


def test_downsample_counts_per_cell(integer_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "downsample-counts",
        "--counts-per-cell",
        "50",
        "--input-file",
        str(integer_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"downsample-counts failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    totals = np.asarray(adata.X.sum(axis=1)).flatten()
    assert np.all(totals <= 50 + 1e-5), "No cell should exceed counts_per_cell"


def test_downsample_total_counts(integer_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "downsample-counts",
        "--total-counts",
        "500",
        "--input-file",
        str(integer_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"downsample-counts failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    total = np.asarray(adata.X.sum())
    assert total <= 500 + 1e-5, "Total counts should not exceed target"


def test_downsample_requires_one_option(integer_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "downsample-counts",
        "--input-file",
        str(integer_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, (
        "Should fail without --counts-per-cell or --total-counts"
    )


def test_downsample_mutex_options(integer_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "downsample-counts",
        "--counts-per-cell",
        "50",
        "--total-counts",
        "500",
        "--input-file",
        str(integer_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail when both options are provided"


def test_downsample_decimals(integer_counts_h5ad, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "downsample-counts",
        "--counts-per-cell",
        "50",
        "--input-file",
        str(integer_counts_h5ad),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"downsample-counts failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    data = (
        adata.X.data
        if scipy.sparse.issparse(adata.X)
        else np.asarray(adata.X).flatten()
    )
    nonzero = data[data != 0]
    assert np.all(nonzero == np.round(nonzero, 2)), "X values not rounded to 2 decimals"

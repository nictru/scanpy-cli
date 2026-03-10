import scanpy as sc
import subprocess
import numpy as np


def test_scale_runs(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "scale",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"scale failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata.n_obs > 0


def test_scale_unit_variance(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "scale",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"scale failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    X = np.asarray(adata.X)
    std_per_gene = X.std(axis=0)
    # After scaling, std should be ~1 for genes that had variance
    nonzero_std = std_per_gene[std_per_gene > 0]
    assert np.allclose(nonzero_std, 1.0, atol=0.05), "Genes should have unit variance"


def test_scale_max_value(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "scale",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--max-value",
        "10",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"scale failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert np.all(np.asarray(adata.X) <= 10.0 + 1e-5), (
        "Values should be clipped to max_value"
    )


def test_scale_no_zero_center(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "scale",
        "--no-zero-center",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"scale --no-zero-center failed: {result.stderr}"


def test_scale_decimals(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "pp",
        "scale",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"scale failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    X = np.asarray(adata.X).flatten()
    assert np.all(X == np.round(X, 3)), "X values are not rounded to 3 decimals"

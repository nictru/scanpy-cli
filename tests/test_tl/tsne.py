import scanpy as sc
import subprocess
import numpy as np
import pickle


def test_tsne_runs(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "tl",
        "tsne",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"tsne failed: {result.stderr}"
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_tsne" in adata.obsm


def test_tsne_n_components(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "tl",
        "tsne",
        "--n-components",
        "3",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"tsne failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata.obsm["X_tsne"].shape[1] == 3


def test_tsne_key_added(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "tl",
        "tsne",
        "--key-added",
        "X_tsne_custom",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"tsne failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    assert "X_tsne_custom" in adata.obsm


def test_tsne_decimals(test_h5ad_path, temp_h5ad_file):
    cmd = [
        "scanpy-cli",
        "tl",
        "tsne",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--decimals",
        "3",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"tsne failed: {result.stderr}"
    adata = sc.read_h5ad(temp_h5ad_file)
    embedding = np.array(adata.obsm["X_tsne"])
    assert np.all(embedding == np.round(embedding, 3)), (
        "X_tsne not rounded to 3 decimals"
    )


def test_tsne_pickle_output(test_h5ad_path, temp_h5ad_file, tmp_path):
    pickle_path = tmp_path / "tsne_embedding.pkl"
    cmd = [
        "scanpy-cli",
        "tl",
        "tsne",
        "--input-file",
        str(test_h5ad_path),
        "--output-file",
        str(temp_h5ad_file),
        "--embedding-output",
        str(pickle_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"tsne failed: {result.stderr}"
    assert pickle_path.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    with open(pickle_path, "rb") as f:
        pkl_embedding = pickle.load(f)
    assert np.array_equal(adata.obsm["X_tsne"], pkl_embedding)

import anndata
import numpy as np
import scipy.sparse
import pytest
import pandas as pd
from click.testing import CliRunner
from scanpy_cli.io.view import view


@pytest.fixture
def rich_h5ad(tmp_path):
    """Create a small h5ad with all major sections populated."""
    rng = np.random.default_rng(0)
    n_obs, n_vars = 30, 20
    X = scipy.sparse.random(
        n_obs, n_vars, density=0.4, random_state=rng, format="csr", dtype=np.float32
    )
    obs = pd.DataFrame(
        {"n_counts": X.sum(1).A1, "leiden": pd.Categorical(["A", "B"] * 15)},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"highly_variable": rng.choice([True, False], n_vars)},
        index=[f"gene_{i}" for i in range(n_vars)],
    )
    adata = anndata.AnnData(X, obs=obs, var=var)
    adata.obsm["X_pca"] = rng.random((n_obs, 5)).astype(np.float32)
    adata.varm["PCs"] = rng.random((n_vars, 5)).astype(np.float32)
    adata.obsp["connectivities"] = scipy.sparse.eye(
        n_obs, format="csr", dtype=np.float32
    )
    adata.layers["raw"] = X.copy()
    adata.uns["pca"] = {"variance_ratio": rng.random(5).astype(np.float32)}

    path = tmp_path / "test.h5ad"
    adata.write_h5ad(path)
    return path


def test_view_exit_code(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert result.exit_code == 0, result.output


def test_view_shows_dimensions(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "30" in result.output
    assert "20" in result.output


def test_view_shows_obs_columns(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "n_counts" in result.output
    assert "leiden" in result.output


def test_view_shows_var_columns(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "highly_variable" in result.output


def test_view_shows_obsm_keys(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "X_pca" in result.output


def test_view_shows_varm_keys(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "PCs" in result.output


def test_view_shows_obsp_keys(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "connectivities" in result.output


def test_view_shows_layers(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "raw" in result.output


def test_view_shows_uns_keys(rich_h5ad):
    runner = CliRunner()
    result = runner.invoke(view, [str(rich_h5ad)])
    assert "pca" in result.output


def test_view_nonexistent_file():
    runner = CliRunner()
    result = runner.invoke(view, ["nonexistent.h5ad"])
    assert result.exit_code != 0

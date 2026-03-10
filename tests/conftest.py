import os
import anndata
import pytest
import tempfile
import scanpy as sc
import sys
import numpy as np
from pathlib import Path

# Add the src directory to the path so we can import the package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

anndata.settings.allow_write_nullable_strings = True


@pytest.fixture(scope="session")
def test_data_dir():
    """Return the directory with static test data (e.g. committed 10x files)."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_h5ad_path(tmp_path_factory):
    """Return the path to the generated test h5ad file (per-worker temp copy)."""
    path = tmp_path_factory.mktemp("data") / "test.h5ad"
    adata = sc.datasets.pbmc3k_processed()
    adata.write_h5ad(path)
    return path


@pytest.fixture
def raw_h5ad_path(temp_h5ad_file):
    """Create a temporary h5ad file with raw count data for Scrublet."""
    # Create a fresh AnnData object from pbmc3k
    adata = sc.datasets.pbmc3k()

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


@pytest.fixture
def raw_log_h5ad_path(temp_h5ad_file):
    """Create a temporary h5ad file with log-normalized raw data for HVG."""
    # Create a fresh AnnData object from pbmc3k
    adata = sc.datasets.pbmc3k()

    # Log normalize the data
    sc.pp.log1p(adata)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


@pytest.fixture
def batch_h5ad_path(test_h5ad_path, temp_h5ad_file):
    """Create a temporary h5ad file with batch information for batch correction."""
    # Read the test data
    adata = sc.read_h5ad(test_h5ad_path)

    # Set random seed for reproducibility
    np.random.seed(42)

    # Create batch column with two categories
    n_cells = adata.n_obs
    adata.obs["batch"] = np.random.choice(["batch1", "batch2"], size=n_cells)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


@pytest.fixture
def raw_batch_h5ad_path(temp_h5ad_file):
    """Create a temporary h5ad file with raw count data and batch information for Scrublet."""
    # Create a fresh AnnData object from pbmc3k
    adata = sc.datasets.pbmc3k()

    # Set random seed for reproducibility
    np.random.seed(42)

    # Create batch column with two categories
    n_cells = adata.n_obs
    adata.obs["batch"] = np.random.choice(["batch1", "batch2"], size=n_cells)

    # Write to temporary file
    adata.write_h5ad(temp_h5ad_file)

    return temp_h5ad_file


@pytest.fixture
def temp_h5ad_file():
    """Create a temporary h5ad file for testing."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    yield tmp_path

    # Cleanup after test
    if tmp_path.exists():
        tmp_path.unlink()


@pytest.fixture
def temp_output_file(suffix=".png"):
    """Create a temporary output file for testing."""
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp_path = Path(tmp.name)

    yield tmp_path

    # Cleanup after test
    if tmp_path.exists():
        tmp_path.unlink()


@pytest.fixture
def temp_plot_file():
    """Create a temporary plot file for testing."""
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    yield tmp_path

    # Cleanup after test
    if tmp_path.exists():
        tmp_path.unlink()

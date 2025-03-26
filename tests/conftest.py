import os
import pytest
import tempfile
import scanpy as sc
import sys
from pathlib import Path

# Add the src directory to the path so we can import the package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


@pytest.fixture(scope="session")
def test_data_dir():
    """Return the directory with test data."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_h5ad_path(test_data_dir):
    """Return the path to the test h5ad file."""
    path = test_data_dir / "test.h5ad"
    if not path.exists():
        # Create the test data directory if it doesn't exist
        test_data_dir.mkdir(exist_ok=True, parents=True)

        # If the test.h5ad file doesn't exist, create a small test AnnData object
        adata = sc.datasets.pbmc3k_processed()
        # Take a small subset for faster testing
        adata = adata[:100, :100].copy()
        # Save the test data
        adata.write_h5ad(path)

    return path


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

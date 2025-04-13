import pytest
import scanpy as sc
from click.testing import CliRunner
from scanpy_cli.io.read_10x_h5 import read_10x_h5


@pytest.fixture
def test_10x_h5_path(test_data_dir):
    """Return the path to the test 10x h5 file."""
    path = test_data_dir / "filtered_feature_bc_matrix.h5"
    if not path.exists():
        raise FileNotFoundError(
            f"Test 10x h5 file not found at {path}. "
            "Please ensure the test data file exists."
        )
    return path


def test_read_10x_h5_basic(test_10x_h5_path, temp_h5ad_file):
    """Test basic functionality of read_10x_h5 command."""
    runner = CliRunner()

    # Run the command
    result = runner.invoke(read_10x_h5, [str(test_10x_h5_path), str(temp_h5ad_file)])

    # Check command executed successfully
    assert result.exit_code == 0
    assert "successfully saved" in result.output

    # Verify output file exists and is valid
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata is not None
    assert adata.n_obs > 0  # Should have cells
    assert adata.n_vars > 0  # Should have genes


def test_read_10x_h5_gex_only(test_10x_h5_path, temp_h5ad_file):
    """Test read_10x_h5 with gex-only option."""
    runner = CliRunner()

    # Run the command with gex-only option
    result = runner.invoke(
        read_10x_h5, [str(test_10x_h5_path), str(temp_h5ad_file), "--gex-only"]
    )

    # Check command executed successfully
    assert result.exit_code == 0
    assert "successfully saved" in result.output

    # Verify output file exists and is valid
    assert temp_h5ad_file.exists()
    adata = sc.read_h5ad(temp_h5ad_file)
    assert adata is not None


def test_read_10x_h5_nonexistent_input(temp_h5ad_file):
    """Test read_10x_h5 with nonexistent input file."""
    runner = CliRunner()

    result = runner.invoke(read_10x_h5, ["nonexistent.h5", str(temp_h5ad_file)])

    # Check command failed as expected
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_read_10x_h5_invalid_output_path(test_10x_h5_path, test_data_dir):
    """Test read_10x_h5 with invalid output path."""
    runner = CliRunner()

    # Try to write to a non-existent directory
    invalid_output = test_data_dir / "nonexistent_dir" / "output.h5ad"

    result = runner.invoke(read_10x_h5, [str(test_10x_h5_path), str(invalid_output)])

    # Check command executed successfully (should create directory)
    assert result.exit_code == 0
    assert "successfully saved" in result.output
    assert invalid_output.exists()

    # Cleanup
    invalid_output.unlink()
    invalid_output.parent.rmdir()

import scanpy as sc
import subprocess
import numpy as np
import tempfile
from pathlib import Path


def _make_input_with_keys(test_h5ad_path, **obs_cols) -> Path:
    """Copy test data to a temp file and add obs columns for regression."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        input_path = Path(tmp.name)
    adata = sc.read_h5ad(test_h5ad_path)
    for col, values in obs_cols.items():
        adata.obs[col] = values
    adata.write_h5ad(input_path)
    return input_path


def test_regress_out_single_key(test_h5ad_path, temp_h5ad_file):
    """Test that the regress_out command runs successfully with a single key."""
    adata = sc.read_h5ad(test_h5ad_path)
    input_path = _make_input_with_keys(
        test_h5ad_path, test_key=[i % 2 for i in range(adata.n_obs)]
    )

    try:
        cmd = [
            "scanpy-cli",
            "pp",
            "regress-out",
            "test_key",
            "--input-file",
            str(input_path),
            "--output-file",
            str(temp_h5ad_file),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"Regress out command failed: {result.stderr}"
        assert temp_h5ad_file.exists(), "Output file was not created"

        adata_regressed = sc.read_h5ad(temp_h5ad_file)
        assert adata_regressed.X.shape == adata.X.shape, (
            "Data shape changed after regression"
        )
    finally:
        input_path.unlink(missing_ok=True)


def test_regress_out_multiple_keys(test_h5ad_path, temp_h5ad_file):
    """Test that the regress_out command runs successfully with multiple keys."""
    adata = sc.read_h5ad(test_h5ad_path)
    input_path = _make_input_with_keys(
        test_h5ad_path,
        test_key1=[i % 2 for i in range(adata.n_obs)],
        test_key2=[i % 3 for i in range(adata.n_obs)],
    )

    try:
        cmd = [
            "scanpy-cli",
            "pp",
            "regress-out",
            "test_key1,test_key2",
            "--input-file",
            str(input_path),
            "--output-file",
            str(temp_h5ad_file),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"Regress out command failed: {result.stderr}"
        assert temp_h5ad_file.exists(), "Output file was not created"

        adata_regressed = sc.read_h5ad(temp_h5ad_file)
        assert adata_regressed.X.shape == adata.X.shape, (
            "Data shape changed after regression"
        )
    finally:
        input_path.unlink(missing_ok=True)


def test_regress_out_numpy_output(test_h5ad_path, temp_h5ad_file, tmp_path):
    """Test that the regress_out command saves the regressed data as a numpy file when requested."""
    adata = sc.read_h5ad(test_h5ad_path)
    input_path = _make_input_with_keys(
        test_h5ad_path, test_key=[i % 2 for i in range(adata.n_obs)]
    )
    numpy_path = tmp_path / "regressed_data.npy"

    try:
        cmd = [
            "scanpy-cli",
            "pp",
            "regress-out",
            "test_key",
            "--input-file",
            str(input_path),
            "--output-file",
            str(temp_h5ad_file),
            "--regressed-output",
            str(numpy_path),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        assert result.returncode == 0, f"Regress out command failed: {result.stderr}"
        assert temp_h5ad_file.exists(), "Output h5ad file was not created"
        assert numpy_path.exists(), "Output numpy file was not created"

        adata = sc.read_h5ad(temp_h5ad_file)
        numpy_data = np.load(numpy_path)
        assert np.array_equal(adata.X, numpy_data), (
            "Numpy file data does not match AnnData data"
        )
    finally:
        input_path.unlink(missing_ok=True)

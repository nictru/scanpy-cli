import subprocess
import importlib.metadata


def test_cli_launches():
    """Test that the scanpy-cli command launches successfully."""
    cmd = ["scanpy-cli", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"CLI failed to launch: {result.stderr}"

    # Check that the output contains expected text
    assert "Usage:" in result.stdout
    assert "Options" in result.stdout
    assert "Commands" in result.stdout
    assert "pp" in result.stdout
    assert "tl" in result.stdout
    assert "pl" in result.stdout


def test_version():
    """Test that the --version parameter works correctly."""
    cmd = ["scanpy-cli", "--version"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"Version command failed: {result.stderr}"

    # Get the expected version from package metadata
    expected_version = importlib.metadata.version("scanpy-cli")

    # Check that the output contains the correct version
    assert f"scanpy-cli, version {expected_version}" in result.stdout


def test_pp_subcommand():
    """Test that the pp subcommand works."""
    cmd = ["scanpy-cli", "pp", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"pp subcommand failed: {result.stderr}"

    # Check that the output contains expected commands
    output = result.stdout
    assert "Options" in result.stdout
    assert "Commands" in output
    assert "pca" in output
    assert "neighbors" in output
    assert "regress-out" in output


def test_tl_subcommand():
    """Test that the tl subcommand works."""
    cmd = ["scanpy-cli", "tl", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"tl subcommand failed: {result.stderr}"

    # Check that the output contains expected commands
    output = result.stdout
    assert "Options" in result.stdout
    assert "Commands" in output
    assert "umap" in output
    assert "leiden" in output


def test_pl_subcommand():
    """Test that the pl subcommand works."""
    cmd = ["scanpy-cli", "pl", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"pl subcommand failed: {result.stderr}"

    # Check that the output contains expected commands
    output = result.stdout
    assert "Options" in result.stdout
    assert "Commands" in output
    assert "umap" in output

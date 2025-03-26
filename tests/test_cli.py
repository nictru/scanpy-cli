import pytest
import subprocess
import re


def strip_ansi_codes(text):
    """Remove ANSI color codes from text."""
    ansi_escape = re.compile(r'\x1b\[[0-9;]*[a-zA-Z]')
    return ansi_escape.sub('', text)


def test_cli_launches():
    """Test that the scanpy-cli command launches successfully."""
    cmd = ["scanpy-cli", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"CLI failed to launch: {result.stderr}"

    # Strip ANSI color codes from output
    clean_output = strip_ansi_codes(result.stdout)

    # Check that the output contains expected text
    assert "Usage:" in clean_output
    assert "Commands:" in clean_output
    assert "pp" in clean_output
    assert "tl" in clean_output
    assert "pl" in clean_output


def test_pp_subcommand():
    """Test that the pp subcommand works."""
    cmd = ["scanpy-cli", "pp", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"pp subcommand failed: {result.stderr}"

    # Strip ANSI color codes from output
    clean_output = strip_ansi_codes(result.stdout)

    # Check that the output contains expected commands
    assert "Commands:" in clean_output
    assert "pca" in clean_output
    assert "neighbors" in clean_output
    assert "regress-out" in clean_output


def test_tl_subcommand():
    """Test that the tl subcommand works."""
    cmd = ["scanpy-cli", "tl", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"tl subcommand failed: {result.stderr}"

    # Strip ANSI color codes from output
    clean_output = strip_ansi_codes(result.stdout)

    # Check that the output contains expected commands
    assert "Commands:" in clean_output
    assert "leiden" in clean_output
    assert "umap" in clean_output
    assert "tsne" in clean_output


def test_pl_subcommand():
    """Test that the pl subcommand works."""
    cmd = ["scanpy-cli", "pl", "--help"]

    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check that the command was successful
    assert result.returncode == 0, f"pl subcommand failed: {result.stderr}"

    # Strip ANSI color codes from output
    clean_output = strip_ansi_codes(result.stdout)

    # Check that the output contains expected commands
    assert "Commands:" in clean_output
    assert "umap" in clean_output
    assert "tsne" in clean_output
    assert "heatmap" in clean_output

import functools
import logging
import sys

import numpy as np
import scipy.sparse
import rich_click as click
from rich.logging import RichHandler


logger = logging.getLogger("scanpy_cli")


def catch_errors(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(str(e))
            sys.exit(1)

    return wrapper


def setup_logging(level: int = logging.WARNING) -> None:
    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=[RichHandler(show_path=False, show_time=False)],
    )
    logger.setLevel(level)


def round_array(arr, decimals):
    if np.issubdtype(arr.dtype, np.floating):
        return np.round(arr, decimals)
    return arr


def round_sparse(mat, decimals):
    if scipy.sparse.issparse(mat) and np.issubdtype(mat.dtype, np.floating):
        mat = mat.copy()
        mat.data = np.round(mat.data, decimals)
    return mat


def round_uns_dict(d, decimals):
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            if np.issubdtype(v.dtype, np.floating):
                d[k] = np.round(v, decimals)
            elif v.dtype.names:
                for field in v.dtype.names:
                    if np.issubdtype(v.dtype[field], np.floating):
                        v[field] = np.round(v[field], decimals)
        elif isinstance(v, dict):
            round_uns_dict(v, decimals)


decimals_option = click.option(
    "--decimals",
    type=int,
    default=None,
    help="Round newly created floating-point fields to this many decimal places (for reproducible hashing in tests).",
)

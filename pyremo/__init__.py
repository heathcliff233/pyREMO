"""Utilities and metadata for the pyREMO package."""
from .resources import resolve_data_file, DATA_FILENAMES, DATA_ENV_VAR

__all__ = ["resolve_data_file", "DATA_FILENAMES", "DATA_ENV_VAR", "__version__"]
__version__ = "0.1.0"

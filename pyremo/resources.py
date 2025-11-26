"""Helpers for locating packaged REMO parameter files."""
from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path
from typing import Iterator, Optional

try:
    from importlib import resources
except ImportError:  # pragma: no cover
    import importlib_resources as resources  # type: ignore

DATA_ENV_VAR = "PYREMO_DATA_DIR"
DATA_FILENAMES = ("FAMRcomm", "BBdat")
_RESOURCE_CACHE: dict[str, str] = {}


def resolve_data_file(filename: str) -> str:
    """Return an absolute path to the requested REMO data file."""
    for directory in _candidate_directories():
        candidate = directory / filename
        if candidate.exists():
            return str(candidate)

    cached = _extract_packaged_resource(filename)
    if cached:
        return cached

    raise FileNotFoundError(
        f"Could not locate {filename}. Set {DATA_ENV_VAR} or reinstall pyremo."
    )


def _candidate_directories() -> Iterator[Path]:
    env_dir = os.environ.get(DATA_ENV_VAR)
    if env_dir:
        yield Path(env_dir).expanduser()

    # Repo root when running from a checkout (remo.py lives alongside the data)
    repo_root = Path(__file__).resolve().parent.parent
    yield repo_root

    # Packaged copy that ships inside pyremo/data
    pkg_data = Path(__file__).resolve().parent / "data"
    yield pkg_data


def _extract_packaged_resource(filename: str) -> Optional[str]:
    cached = _RESOURCE_CACHE.get(filename)
    if cached and Path(cached).exists():
        return cached

    try:
        data_pkg = resources.files("pyremo.data")
    except (ImportError, ModuleNotFoundError):
        return None

    resource = data_pkg.joinpath(filename)
    if not resource.is_file():
        return None

    cache_dir = Path(tempfile.gettempdir()) / "pyremo-data"
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / filename
    if not target.exists():
        with resources.as_file(resource) as source_path:
            shutil.copyfile(source_path, target)
    _RESOURCE_CACHE[filename] = str(target)
    return _RESOURCE_CACHE[filename]

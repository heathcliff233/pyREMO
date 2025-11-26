"""High level reconstruction helpers shared by the CLI and library users."""
from __future__ import annotations

import io
import os
import sys
import tempfile
from typing import IO, Optional, Tuple, Union

from famr_python.reconstruct import run_reconstruction
from famr_python.topology import generate_hdd, parse_pdb_sequence, read_famr_comm

from .resources import resolve_data_file

TextSource = Union[str, os.PathLike[str], IO[str]]

_COMM_DATA: Optional[dict] = None
_BBDAT_PATH: Optional[str] = None


def load_parameters() -> Tuple[dict, str]:
    """Return cached FAMRcomm contents and the BBdat file path."""
    global _COMM_DATA, _BBDAT_PATH

    if _COMM_DATA is None:
        famr_comm = resolve_data_file("FAMRcomm")
        _COMM_DATA = read_famr_comm(famr_comm)

    if _BBDAT_PATH is None:
        _BBDAT_PATH = resolve_data_file("BBdat")

    return _COMM_DATA, _BBDAT_PATH


def reconstruct_file(
    input_pdb: TextSource,
    output_pdb: Union[str, os.PathLike[str]],
    overwrite: bool = False,
    *,
    quiet: bool = False,
) -> bool:
    """Reconstruct a single CA-trace PDB file into a full-atom model written to `output_pdb`."""
    comm_data, bbdat_file = load_parameters()
    success, message = _process_source(input_pdb, output_pdb, comm_data, bbdat_file, overwrite)
    if not success and message and not quiet:
        print(message)
    return success


def reconstruct(pdb_string: str) -> str:
    """Reconstruct a CA-only PDB string and return the resulting full-atom PDB contents."""
    if not pdb_string.strip():
        raise ValueError("Input PDB string is empty.")

    comm_data, bbdat_file = load_parameters()

    pdb_stream = io.StringIO(pdb_string)
    output_buffer = io.StringIO()

    success, message = _process_source(
        pdb_stream,
        output_buffer,
        comm_data,
        bbdat_file,
        overwrite=True,
    )
    if not success:
        raise RuntimeError(message or "REMO reconstruction failed.")

    return output_buffer.getvalue()


def _process_source(
    pdb_source: TextSource,
    output_target,
    comm_data: dict,
    bbdat_file: str,
    overwrite: bool,
) -> Tuple[bool, Optional[str]]:
    if isinstance(output_target, str):
        target_path = output_target
    elif isinstance(output_target, os.PathLike):
        target_path = os.fspath(output_target)
    else:
        target_path = None

    if isinstance(target_path, str):
        if os.path.exists(target_path) and not overwrite:
            return False, None

        output_dir = os.path.dirname(target_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    parsed_source = pdb_source
    reconstruction_source = pdb_source
    if hasattr(pdb_source, "read"):
        pdb_text = pdb_source.read()
        parsed_source = io.StringIO(pdb_text)
        reconstruction_source = io.StringIO(pdb_text)

    try:
        seq, _ = parse_pdb_sequence(parsed_source)
        if not seq:
            return False, f"No CA atoms found in {_describe_source(pdb_source)}."
    except Exception as exc:  # pragma: no cover
        return False, f"Error parsing {_describe_source(pdb_source)}: {exc}"

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_hdd:
        hdd_path = tmp_hdd.name

    try:
        generate_hdd(seq, comm_data, hdd_path)

        original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        try:
            run_reconstruction(reconstruction_source, hdd_path, bbdat_file, output_target)
        except Exception as exc:  # pragma: no cover
            return False, f"Error reconstructing {_describe_source(pdb_source)}: {exc}"
        finally:
            if sys.stdout is not original_stdout:
                sys.stdout.close()
                sys.stdout = original_stdout
    finally:
        if os.path.exists(hdd_path):
            os.remove(hdd_path)

    return True, None


def _describe_source(source) -> str:
    if isinstance(source, (str, os.PathLike)):
        return os.fspath(source)
    name = getattr(source, "name", None)
    if isinstance(name, str) and name:
        return name
    return "<memory>"

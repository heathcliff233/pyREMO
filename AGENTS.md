# REMO Python Rebuild – Project Context

## Current Snapshot (Nov 2025)
- Legacy Perl/Fortran tooling has been retired. The project now ships a `pyproject.toml` + `pyremo` package, so `pip install .` drops a global `remo` console script (running `python remo.py` from the repo still works for devs) and exposes a Python API (`pyremo.reconstruct` / `pyremo.reconstruct_file`).
- Only native dependencies are `numpy` (required) and `tqdm` (optional via the `progress` extra). A typical install looks like:

  ```bash
  python -m pip install .[progress]
  remo -i <input_ca.pdb> -o <full_atom.pdb> [--overwrite]
  ```

- The canonical parameter/data assets (`FAMRcomm`, `BBdat`) now live solely inside `pyremo/data` (tracked in git and bundled in wheels/sdists). `pyremo.resources.resolve_data_file` checks `PYREMO_DATA_DIR`, then the repo root (for custom drops), then `pyremo/data`/packaged copies.
- Example inputs (`examples/xx1.pdb`, etc.) double as smoke tests for new changes.

## Reconstruction Pipeline
1. **Sequence Extraction** – `parse_pdb_sequence` grabs the Cα sequence/order straight from the input PDB.
2. **Topology Build** – `generate_hdd` recreates the Fortran-era `Hdd` file (atom definitions, bonds, angles, Ramachandran table) using `FAMRcomm`.
3. **Hdd Parsing** – `read_hdd` converts that topology into Python objects (`Protein`, `Residue`, `Atom`).
4. **Coordinate Mapping** – Input Cα coordinates are assigned back to the reconstructed residues.
5. **Backbone Placement** – `backbone.add_backbone` looks up BBdat bins and positions N/C atoms (with a fallback for degenerate geometry and explicit terminal handling).
6. **Side-chain Build** – `sidechain.build_sidechains` places O/OXT, CB, and the remaining heavy atoms using idealized internal coordinates propagated by bond connectivity.
7. **Hydrogen Placement** – `hydrogen.add_hydrogens` adds hydrogens via heuristic templates (tetrahedral for sp³ carbons, planar for peptide nitrogens, etc.).
8. **PDB Output** – `io.write_pdb` writes a fully populated structure; serial numbers wrap at 99,999 to stay within classic PDB limits.

All intermediate Hdd files are generated in a temporary directory and removed after each run.

## Module Breakdown
| Module | Responsibilities |
| --- | --- |
| `remo.py` | CLI parsing, batch orchestration, temp-file management. |
| `pyremo/api.py` | Shared helpers: parameter caching, `reconstruct_file`, and the in-memory `pyremo.reconstruct` API. |
| `pyremo/resources.py` | Resolves paths to `FAMRcomm`/`BBdat`, honoring `PYREMO_DATA_DIR` and packaged data copies. |
| `famr_python/topology.py` | Reads `FAMRcomm`, builds residue-specific atoms/bonds/angles, writes `Hdd`. |
| `famr_python/io.py` | Converts `Hdd` back into Python objects and handles PDB read/write. |
| `famr_python/backbone.py` | BBdat ingestion + backbone coordinate reconstruction (with fallback for edge bins). |
| `famr_python/sidechain.py` | Places side-chain heavy atoms (currently single-rotamer heuristic) and terminal oxygens. |
| `famr_python/hydrogen.py` | Adds hydrogens for backbone and side chains with geometry heuristics. |
| `famr_python/geometry.py` | Distance/angle/dihedral utilities, rotation matrices, atom placement helper. |
| `famr_python/reconstruct.py` | Wires the steps above into a single callable used by the CLI. |

## Differences vs Legacy REMO
- **No Fortran runtime**: avoids macOS `gfortran` incompatibilities and simplifies deployment.
- **Single-mode focus**: Implements Mode 2 (Cα reconstruction). Template-driven and FASTA-driven modes from the Perl script are not ported yet.
- **Deterministic heuristics**: Rotamers and hydrogens use fixed geometries instead of the stochastic optimization present in Fortran; this keeps the rewrite simple but can create steric clashes.
- **Terminal handling**: Python version explicitly sets NH3/NP, CC/OC/OXT coordinates to avoid zero vectors that plagued earlier drafts.

## Known Limitations / Backlog
1. **Side-chain optimization** – No rotamer search or clash resolution; add a Dunbrack lookup or simple Monte Carlo in future.
2. **Hydrogen refinement** – Donor hydrogens are not rotated toward acceptors; optional minimization could improve H-bond geometry.
3. **Secondary-structure inputs** – `seq.dat`/PRHB hints are ignored; only CA geometry drives BBdat lookup.
4. **Advanced modes** – Building from FASTA or mixing in templates remains unimplemented.
5. **Testing** – No automated regression suite yet; coverage currently relies on running `examples/xx1.pdb`.

## Quick Diagnostic Checklist
- After `python -m pip install .[progress]`, running `remo -i examples/xx1.pdb -o /tmp/xx1.full.pdb --overwrite` should complete without warnings and produce a file where no atoms have coordinates `(0,0,0)`. The legacy `python remo.py …` path remains valid from a repo checkout.
- Programmatic smoke: `python - <<'PY'\nimport pyremo, pathlib\nca = pathlib.Path(\"examples/xx1.pdb\").read_text()\nfull = pyremo.reconstruct(ca)\nprint(len(full))\nPY` should print a positive length without errors.
- If you see missing dependencies, reinstall with `python -m pip install .[progress]` (the extra pulls `tqdm`; base install only needs `numpy`).
- Corrupted `FAMRcomm`/`BBdat` manifests as parsing errors—verify the copies under `pyremo/data/` (or whatever directory `PYREMO_DATA_DIR` points to) remain unedited.

This document should be updated whenever we touch the reconstruction logic, data formats, or supported workflows so future contributors can quickly understand the current state of the Python port.

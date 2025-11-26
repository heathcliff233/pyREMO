# REMO Python Rebuild – Project Context

## Current Snapshot (Nov 2025)
- Legacy Perl/Fortran tooling has been retired. `remo.py` is now the only entry point and runs a fully Python reconstruction stack located under `famr_python/`.
- Only native dependencies are `numpy` (required) and `tqdm` (optional). No compilation steps are necessary; the canonical run command is:

  ```bash
  python remo.py -i <input_ca.pdb> -o <full_atom.pdb> [--overwrite]
  ```

- The original parameter/data assets (`FAMRcomm`, `BBdat`) remain untouched and must stay in the repo root because the Python code reads them directly.
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
- `python remo.py -i examples/xx1.pdb -o /tmp/xx1.full.pdb --overwrite` should complete without warnings and produce a file where no atoms have coordinates `(0,0,0)`.
- If you see missing dependencies, install `numpy` (and optionally `tqdm`) in your active environment.
- Corrupted `FAMRcomm`/`BBdat` manifests as parsing errors—verify the files remain unedited.

This document should be updated whenever we touch the reconstruction logic, data formats, or supported workflows so future contributors can quickly understand the current state of the Python port.


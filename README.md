# REMO (Python Edition)

Reconstruct Atomic MOdels (REMO) now runs as a pure Python tool.  
Given a Cα-only PDB file, it rebuilds backbone atoms, side chains, and hydrogens using the
original `FAMRcomm` parameters and `BBdat` backbone library—no Fortran binaries required.

---

## Installation

### Option 1 – use the CLI everywhere

```bash
python -m pip install .[progress]     # add [progress] to pull in tqdm
remo -i examples/xx1.pdb -o /tmp/xx1.full.pdb --overwrite
```

This drops a `remo` console script into your active environment and bundles the `FAMRcomm` / `BBdat`
parameter files with the package. The helper automatically locates them, but you can override the
lookup by pointing `PYREMO_DATA_DIR` to a directory that holds fresh copies.

### Option 2 – run straight from the repo

- Python 3.9+
- `numpy` (required)
- `tqdm` (optional – progress bar)
- Parameter/data files ship inside `pyremo/data` (set `PYREMO_DATA_DIR` to use a custom directory)

```bash
python -m pip install numpy tqdm    # or reuse an existing environment
python remo.py -i examples/xx1.pdb -o /tmp/xx1.full.pdb --overwrite
```

---

## Repository Layout

```
pyproject.toml        # Packaging metadata / entry points
remo.py               # CLI entry point
pyremo/               # Package glue + bundled data copies
    resources.py      # Data resolution helper used by the CLI
    data/             # Packaged FAMRcomm + BBdat (only canonical copy now)
famr_python/          # Python reconstruction package
    backbone.py       # BBdat lookup + backbone placement
    sidechain.py      # Idealized side-chain builder
    hydrogen.py       # Hydrogen placement heuristics
    geometry.py       # Vector math utilities
    io.py             # Hdd/PDB readers & writer
    topology.py       # FAMRcomm parser + Hdd generator
examples/             # Sample Cα inputs + helper files
```

The canonical `FAMRcomm` / `BBdat` files now live in `pyremo/data` and ship with the package.
Set `PYREMO_DATA_DIR=/path/to/dir` if you need to point the CLI at a custom or newer copy.

---

## Quick Start

### Installed CLI (recommended)

```bash
remo \
    --input examples/xx1.pdb \
    --output output_test.pdb \
    --overwrite
```

### Local checkout

```bash
python remo.py \
    --input examples/xx1.pdb \
    --output output_test.pdb \
    --overwrite
```

### CLI Arguments

| Option | Description |
| --- | --- |
| `-i / --input` | Path to a PDB file or directory containing PDB/ENT files. |
| `-o / --output` | Destination file or directory. When `--input` is a directory, a matching tree is created. |
| `--overwrite` | Overwrite existing outputs (otherwise files are skipped). |

The script automatically:
1. Parses the input Cα trace (`parse_pdb_sequence`).
2. Generates a temporary `Hdd` topology from `FAMRcomm`.
3. Reconstructs the structure via `famr_python.reconstruct`.
4. Writes a full-atom PDB to the requested location.

---

## Outputs

Each processed input produces a standard PDB containing:
- All backbone atoms (N/CA/C/O and terminal modifications NH3/OXT).
- Side-chain heavy atoms using a single idealized rotamer per residue.
- Hydrogens placed with fixed bond lengths/angles (sufficient for visualization or downstream refinement).

Temporary `Hdd` files are generated in `/tmp` and deleted once reconstruction completes.

---

## Known Limitations / Future Work

- Side-chain rotamers are not optimized for clashes; they use a single common geometry per residue.
- Hydrogen orientations are heuristic and may need refinement before energy minimization.
- Secondary-structure-specific placement is limited (helix/sheet cues from `BBdat` bins only).
- Only “mode 2” (Cα reconstruction) from the original REMO workflow is covered for now.

Despite those simplifications, the pipeline consistently outputs complete coordinates—no zeroed atoms—and avoids the Fortran runtime issues that motivated this rewrite.

---

## Example

```
remo -i examples/xx1.pdb -o examples/xx1.full.pdb --overwrite
```

Result: `examples/xx1.full.pdb` (full-atom model), using the bundled example Cα file.

---

For questions about the historical REMO method, refer to:

Yunqi Li & Yang Zhang (2009).  
“REMO: a new protocol to generate full atomic protein models from C-alpha traces by optimizing backbone hydrogen-bonding network.” Protein Sci. 18(3): 665–676.

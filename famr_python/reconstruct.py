import os
from os import PathLike
from typing import Any

import numpy as np

from .io import read_pdb_ca, read_hdd, write_pdb
from .backbone import read_bbdat, add_backbone
from .sidechain import build_sidechains
from .hydrogen import add_hydrogens
from .geometry import distance


def run_reconstruction(pdb_source, hdd_file: str, bbdat_file: str, output_target):
    source_desc = _describe_source(pdb_source)
    print(f"Starting reconstruction for {source_desc}")
    
    # 1. Read Hdd (Topology)
    print("Reading topology (Hdd)...")
    protein = read_hdd(hdd_file)
    print(f"Loaded {len(protein.residues)} residues from Hdd.")
    
    # 2. Read PDB (C-alpha trace)
    print("Reading PDB...")
    ca_atoms = read_pdb_ca(pdb_source)
    print(f"Loaded {len(ca_atoms)} CA atoms.")
    
    # 3. Map Coordinates
    print("Mapping coordinates...")
    # Using index-based mapping
    count = 0
    for i, res in enumerate(protein.residues):
        if i < len(ca_atoms):
            name, chain, num, coords = ca_atoms[i]
            # Verify name mismatch?
            # For MVP, assume sequence matches or trust Hdd generator (remo.py)
            # remo.py guarantees Hdd sequence matches PDB input.
            
            ca = next((a for a in res.atoms if a.name == 'CA'), None)
            if ca:
                ca.coords = coords
                res.chain_id = chain if chain.strip() else 'A'
                res.orig_resnum = num
                count += 1
    print(f"Mapped {count} residues.")
    
    # 4. Backbone Reconstruction
    print("Reconstructing backbone...")
    if os.path.exists(bbdat_file):
        bb_data = read_bbdat(bbdat_file)
        add_backbone(protein, bb_data)
    else:
        print("Error: BBdat file not found!")
        return
        
    # 5. Sidechain Reconstruction
    print("Reconstructing sidechains...")
    build_sidechains(protein)
    
    # 6. Hydrogen Addition
    print("Adding hydrogens...")
    add_hydrogens(protein)
    
    # 7. Output
    output_desc = _describe_source(output_target)
    print(f"Writing output to {output_desc}...")
    write_pdb(protein, output_target)
    print("Done.")


def _describe_source(source: Any) -> str:
    if isinstance(source, (str, PathLike)):
        return os.fspath(source)
    name = getattr(source, "name", None)
    if isinstance(name, str) and name:
        return name
    return "<memory>"

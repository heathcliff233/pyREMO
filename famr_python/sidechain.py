import numpy as np
from typing import Dict, List, Tuple
from .data_structures import Protein, Atom, Residue
from .geometry import place_atom, distance

# Idealized torsions (degrees) for MVP
DEFAULT_TORSIONS = {
    'ALL': -60.0, 
}

def build_sidechains(protein: Protein):
    """
    Reconstructs side-chains (and O atom) for all residues.
    Assumes N, CA, C are present.
    """
    
    for i, res in enumerate(protein.residues):
        # 1. Identify Backbone Atoms
        n = next((a for a in res.atoms if a.name in ['N', 'NH3', 'NP']), None)
        ca = next((a for a in res.atoms if a.name == 'CA'), None)
        c = next((a for a in res.atoms if a.name in ['C', 'CC']), None)
        o_main = next((a for a in res.atoms if a.name in ['O', 'OC']), None)
        oxt = next((a for a in res.atoms if a.name == 'OXT'), None)
        
        if not (n and ca and c):
            continue
            
        # 2. Place O
        if o_main and np.linalg.norm(o_main.coords) < 0.1:
            n_next = None
            if i + 1 < len(protein.residues):
                n_next = next((a for a in protein.residues[i+1].atoms if a.name in ['N', 'NH3', 'NP']), None)
            bond_len = 1.23
            if n_next:
                v1 = ca.coords - c.coords
                v2 = n_next.coords - c.coords
                v1 /= np.linalg.norm(v1)
                v2 /= np.linalg.norm(v2)
                v_o = -(v1 + v2)
                if np.linalg.norm(v_o) < 0.1: v_o = np.cross(v1, np.array([0,0,1]))
                v_o /= np.linalg.norm(v_o)
                o_main.coords = c.coords + v_o * bond_len
            else:
                v_ca_c = c.coords - ca.coords
                o_main.coords = c.coords + (v_ca_c / np.linalg.norm(v_ca_c)) * bond_len

        if oxt and np.linalg.norm(oxt.coords) < 0.1 and c and ca:
            if o_main and np.linalg.norm(o_main.coords) > 0.1:
                v = o_main.coords - c.coords
                norm_v = np.linalg.norm(v)
                if norm_v > 0.1:
                    oxt.coords = c.coords - (v / norm_v) * 1.26
            else:
                vec = c.coords - ca.coords
                norm_vec = np.linalg.norm(vec)
                if norm_vec > 0.1:
                    oxt.coords = c.coords + (vec / norm_vec) * 1.26
                
        # 3. Place CB
        cb = next((a for a in res.atoms if a.name == 'CB'), None)
        if cb and np.linalg.norm(cb.coords) < 0.1:
            v_n = n.coords - ca.coords
            v_c = c.coords - ca.coords
            v_n /= np.linalg.norm(v_n)
            v_c /= np.linalg.norm(v_c)
            
            normal = np.cross(v_n, v_c)
            normal /= np.linalg.norm(normal)
            
            # Axis bisecting N-CA-C
            axis = -(v_n + v_c)
            axis /= np.linalg.norm(axis)
            
            # Tetrahedral placement
            # -0.5 * (v_n + v_c) component is along axis * cos(half_angle)?
            # Use combination of normal and axis.
            # 0.866 * normal - 0.5 * axis? (No, axis is bisector)
            
            vec = -0.5 * (v_n + v_c) + 0.866 * normal
            vec /= np.linalg.norm(vec)
            cb.coords = ca.coords + vec * 1.53
            
        # 4. Place Sidechains (BFS)
        # Queue stores (atom, parent, grandparent)
        # Initial: (CB, CA, N)
        
        placed = {n.index, ca.index, c.index}
        if o_main: placed.add(o_main.index)
        if oxt: placed.add(oxt.index)
        if cb: placed.add(cb.index)
        
        queue = []
        if cb:
            queue.append((cb, ca, n))
        
        while queue:
            curr, parent, gp = queue.pop(0)
            
            neighbors = curr.bonds
            for child in neighbors:
                if child.index in placed: continue
                if child.type.startswith('H'): continue
                
                # Place child using gp -> parent -> curr -> child
                if parent and gp:
                    bond_len = 1.54
                    angle_rad = np.radians(109.5)
                    torsion_deg = -60.0
                    if curr.name != 'CB': torsion_deg = 180.0
                    torsion_rad = np.radians(torsion_deg)
                    
                    child.coords = place_atom(gp.coords, parent.coords, curr.coords, bond_len, angle_rad, torsion_rad)
                    
                else:
                    # Not enough ancestors
                    pass
                    
                placed.add(child.index)
                queue.append((child, curr, parent))

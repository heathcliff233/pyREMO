import numpy as np
from typing import List
from .data_structures import Protein, Atom
from .geometry import place_atom
from .naming import normalize_atom_name

N_TERM_H_NAMES = [normalize_atom_name(n) for n in ('1H','2H','3H')]

def add_hydrogens(protein: Protein):
    """
    Places hydrogen atoms.
    """
    for res in protein.residues:
        for atom in res.atoms:
            if not atom.type.startswith('H') and atom.name != 'H':
                # It is a heavy atom. Check if it has H neighbors to place.
                pass
            else:
                continue
                
    # Iterate over heavy atoms and place attached hydrogens
    for res in protein.residues:
        for atom in res.atoms:
            # If atom is H, skip (it will be placed by its parent)
            if atom.type.startswith('H') or atom.name.startswith('H') or atom.name.isdigit():
                continue
                
            # Find unplaced H neighbors
            if atom.name in ['N', 'NH3', 'NP']:
                if atom.name in ['NH3', 'NP']:
                    target_names = N_TERM_H_NAMES
                else:
                    target_names = ['H']
                h_neighbors = [h for h in res.atoms if h.name in target_names and np.linalg.norm(h.coords) < 0.01]
            else:
                h_neighbors = [
                    n for n in atom.bonds
                    if (n.type.startswith('H') or n.name.startswith('H') or n.name[0].isdigit())
                    and np.linalg.norm(n.coords) < 0.01
                ]
            
            if not h_neighbors:
                continue
                
            # Count placed neighbors (Heavy atoms or previously placed H)
            placed_neighbors = [n for n in atom.bonds if np.linalg.norm(n.coords) > 0.01]
            
            # Case 1: sp3 Carbon (4 neighbors total)
            # CT1 (CH), CT2 (CH2), CT3 (CH3).
            if atom.type.startswith('CT') or atom.type == 'CP1':
                if len(placed_neighbors) == 3: # CH
                    # Place H to complete tetrahedron
                    p1, p2, p3 = placed_neighbors[0].coords, placed_neighbors[1].coords, placed_neighbors[2].coords
                    c = atom.coords
                    v1 = p1 - c; v2 = p2 - c; v3 = p3 - c
                    v1 /= np.linalg.norm(v1); v2 /= np.linalg.norm(v2); v3 /= np.linalg.norm(v3)
                    avg_vec = v1 + v2 + v3
                    avg_vec /= np.linalg.norm(avg_vec)
                    
                    h_pos = c - avg_vec * 1.09 # bond len
                    h_neighbors[0].coords = h_pos
                    
                elif len(placed_neighbors) == 2: # CH2
                    # Plane bisector
                    p1, p2 = placed_neighbors[0].coords, placed_neighbors[1].coords
                    c = atom.coords
                    v1 = p1 - c; v2 = p2 - c
                    v1 /= np.linalg.norm(v1); v2 /= np.linalg.norm(v2)
                    
                    axis = v1 + v2
                    axis /= np.linalg.norm(axis)
                    
                    normal = np.cross(v1, v2)
                    normal /= np.linalg.norm(normal)
                    
                    # Two H's: Up and Down relative to plane defined by normal?
                    # No. Tetrahedral geometry relative to C-P1-P2 plane.
                    # Axis bisects P1-C-P2 angle.
                    # H's are in plane perpendicular to (P1-C-P2) containing Axis?
                    # No. H-C-H plane is perpendicular to P1-C-P2 plane.
                    # Bisector of H-C-H is -Axis.
                    
                    # Direction of bisector of H-C-H
                    h_axis = -axis
                    
                    # H positions: Rotate h_axis by +/- theta/2 around normal?
                    # No. H-C-H plane contains h_axis and normal.
                    # Angle H-C-H ~ 109.
                    # Half angle ~ 54.5.
                    
                    # h1 = c + 1.09 * (h_axis * cos(54.5) + normal * sin(54.5))
                    # h2 = c + 1.09 * (h_axis * cos(54.5) - normal * sin(54.5))
                    
                    sin_a = np.sin(np.radians(54.5))
                    cos_a = np.cos(np.radians(54.5))
                    
                    if len(h_neighbors) >= 1:
                        h_neighbors[0].coords = c + 1.09 * (h_axis * cos_a + normal * sin_a)
                    if len(h_neighbors) >= 2:
                        h_neighbors[1].coords = c + 1.09 * (h_axis * cos_a - normal * sin_a)
                        
                elif len(placed_neighbors) == 1: # CH3
                    # Staggered
                    p1 = placed_neighbors[0].coords
                    c = atom.coords
                    v1 = p1 - c
                    v1 /= np.linalg.norm(v1)
                    
                    # Axis is -v1.
                    axis = -v1
                    
                    # Find arbitrary normal
                    if abs(axis[2]) < 0.9: ref = np.array([0,0,1])
                    else: ref = np.array([1,0,0])
                    n1 = np.cross(axis, ref)
                    n1 /= np.linalg.norm(n1)
                    n2 = np.cross(axis, n1)
                    
                    # Cone at 109.5 deg (angle C-C-H).
                    # Complement angle is 180-109.5 = 70.5.
                    # sin(70.5) radius, cos(70.5) along axis.
                    
                    theta = np.radians(180 - 109.5)
                    sin_t = np.sin(theta)
                    cos_t = np.cos(theta)
                    
                    # 3 H's at 0, 120, 240 deg around axis
                    for k, h in enumerate(h_neighbors):
                        phi = np.radians(k * 120.0)
                        h_vec = axis * cos_t + (n1 * np.cos(phi) + n2 * np.sin(phi)) * sin_t
                        h.coords = c + h_vec * 1.09

            # Case 2: Backbone N (Planar or Pyramidal?)
            # Peptide N is usually planar (sp2).
            # N-term N (NH3) is sp3 (tetrahedral).
            elif atom.name in ['N', 'NH3', 'NP']:
                # Check if N-term (NH3) or Peptide (NH)
                if len(h_neighbors) == 1: # Peptide NH
                    # Planar with C(prev), CA(curr).
                    # Bisect C-N-CA outer angle.
                    if len(placed_neighbors) >= 2:
                        p1, p2 = placed_neighbors[0].coords, placed_neighbors[1].coords
                        c = atom.coords
                        v1 = p1 - c; v2 = p2 - c
                        v1 /= np.linalg.norm(v1); v2 /= np.linalg.norm(v2)
                        v_h = -(v1 + v2)
                        v_h /= np.linalg.norm(v_h)
                        h_neighbors[0].coords = c + v_h * 1.01
                elif len(h_neighbors) == 3 and placed_neighbors: # NH3 / NP
                    # Similar to CH3
                    # Use same logic as CH3 but bond length 1.01
                    p1 = placed_neighbors[0].coords
                    c = atom.coords
                    v1 = p1 - c
                    v1 /= np.linalg.norm(v1)
                    axis = -v1
                    if abs(axis[2]) < 0.9:
                        ref = np.array([0, 0, 1])
                    else:
                        ref = np.array([1, 0, 0])
                    n1 = np.cross(axis, ref)
                    n1 /= np.linalg.norm(n1)
                    n2 = np.cross(axis, n1)
                    theta = np.radians(180 - 109.5)
                    sin_t = np.sin(theta)
                    cos_t = np.cos(theta)
                    for k, h in enumerate(h_neighbors):
                        phi = np.radians(k * 120.0)
                        h_vec = axis * cos_t + (n1 * np.cos(phi) + n2 * np.sin(phi)) * sin_t
                        h.coords = c + h_vec * 1.01
                         
            # Other cases (OH, NH2, etc.) - Simplified for MVP
            # Just place them randomly or along bond?
            else:
                # Fallback: place along x-axis to avoid zeros
                for h in h_neighbors:
                    h.coords = atom.coords + np.array([1.0, 0, 0])

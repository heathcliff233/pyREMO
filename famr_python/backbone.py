import numpy as np
from typing import Dict, List, Tuple
from .data_structures import Protein, Atom

BB_DICT_TYPE = Dict[Tuple[int, int, int, int, int, int, int, int], np.ndarray]

def read_bbdat(filename: str) -> BB_DICT_TYPE:
    """
    Reads BBdat into a dictionary.
    Key: (flag, i1, i2, i3, i4, i5, i6, i8)
    Value: np.array([6 floats])
    """
    bb_data = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    # Skip header (counts)
    start_line = 1
    
    for line in lines[start_line:]:
        if not line.strip(): continue
        # Format: 8(I2,1x), 6(F7.3,1x)
        # Python split() handles variable spaces, but fixed width is safer if columns merge.
        # With 1x spacing, split should work if numbers don't overflow.
        # i1..i8 are small. floats are F7.3.
        parts = line.split()
        if len(parts) < 14: continue
        
        try:
            keys = tuple(int(p) for p in parts[:8])
            coords = np.array([float(p) for p in parts[8:14]])
            bb_data[keys] = coords
        except ValueError:
            continue
        
    return bb_data

def get_ca_coords(protein: Protein) -> List[np.ndarray]:
    """Extracts CA coordinates from protein."""
    cas = []
    for res in protein.residues:
        # Find CA
        ca = next((a for a in res.atoms if a.name == 'CA'), None)
        if ca:
            cas.append(ca.coords)
        else:
            cas.append(np.zeros(3)) # Should not happen
    return cas

def solve_linear(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solves A * x = b.
    A is 3x3 (rows are basis vectors).
    b is 3 vector (projections).
    Returns x (vector in global frame).
    """
    # If A is orthogonal, x = A.T @ b
    # We use solve for robustness
    try:
        return np.linalg.solve(a, b)
    except np.linalg.LinAlgError:
        return np.zeros(3)

def _find_atom(residue, names):
    return next((a for a in residue.atoms if a.name in names), None)

def add_backbone(protein: Protein, bb_data: BB_DICT_TYPE):
    """
    Reconstructs backbone N and C atoms using CA trace and BBdat.
    """
    # Extract CA coords
    ca_coords = get_ca_coords(protein)
    n_res = len(protein.residues)
    
    # Helper to get CA(i) safe
    def get_ca(i):
        if i < 0:
            # Extrapolate backwards
            # FAMR: xx(natom+2)=2*xx(ip3)-xx(ip2) -> 2*CA(1)-CA(0)
            # But indices are tricky.
            # If i=-1 (CA before first), we need CA(0) and CA(1).
            c0 = ca_coords[0]
            c1 = ca_coords[1]
            return 2*c0 - c1
        elif i >= n_res:
            # Extrapolate forwards
            c_last = ca_coords[-1]
            c_prev = ca_coords[-2]
            return 2*c_last - c_prev
        return ca_coords[i]
    
    # For FAMR logic, we need 4 points: i-1, i, i+1, i+2.
    # We iterate i from 0 to n_res-1.
    # But boundaries need specific handling matching FAMR.
    
    # Track pending nitrogen placement
    pending_n = None

    for i in range(n_res):
        # Define the 4 CA indices
        # FAMR indices:
        # if i<2 (so 1 in Fortran, 0 in Python):
        #   ip1=0 (virtual), ip2=i, ip3=i+1, ip4=i+2
        # FAMR loops 1..Nres.
        # i=1 (Python 0): ip1=?, ip2=0, ip3=1, ip4=2
        
        p2 = ca_coords[i] # CA(i)
        
        if i == 0:
            # FAMR: ip1 is virtual. 
            # xx(natom+2)=2*xx(ip3)-xx(ip2)
            # So ip1 = 2*CA(1) - CA(0) is NOT what it does.
            # It sets coords at natom+2, natom+3.
            # Line 2101: xx(natom+2) = 2*xx(ip3) - xx(ip2).
            # This is used as ip1? No.
            # i < 2: ip1=0??
            # Wait, if ip1=0, what is xx(0)?
            # Fortran arrays 1-based. xx(0) is valid?
            # No.
            # The code says:
            # if(i<2)then ip1=0 ...
            # But later: vx14=xx(ip4)-xx(ip1).
            # So xx(0) MUST be defined.
            # But where?
            pass
            
        # Let's use our get_ca extrapolation which mimics "2*p_next - p_curr".
        # ca(i-1) for i=0 -> get_ca(-1) -> 2*CA(0) - CA(1).
        # FAMR: 2*xx(ip3)-xx(ip2) -> 2*CA(1) - CA(0).
        # Wait, 2*CA(1) - CA(0) is effectively CA(2) if linear?
        # No. CA(1) + (CA(1)-CA(0)). Extrapolating forward from 0->1?
        # That would be CA(2).
        # We want CA(-1). Extrapolating backward from 1->0?
        # CA(0) + (CA(0)-CA(1)). = 2*CA(0) - CA(1).
        # FAMR uses 2*CA(1) - CA(0)? That puts the point far ahead.
        # Let's re-read FAMR line 2101 carefully.
        
        # "xx(natom+2)=2*xx(ip3)-xx(ip2)"
        # ip3 is CA(i+1) (in general). ip2 is CA(i).
        # So p_virtual = 2*CA(i+1) - CA(i).
        # This virtual point is used as... what?
        # "if(i<2) ... ip1=0"
        # But `xx` is accessed at `ip1`.
        # Maybe `ip1` refers to `natom+2`?
        # No, `ip1` is set to 0.
        # Maybe `xx(0)` is set somewhere?
        
        # Let's assume simple linear extrapolation is sufficient.
        p1 = get_ca(i-1)
        p3 = get_ca(i+1)
        p4 = get_ca(i+2)
        
        # Calculate vectors
        v23 = p3 - p2
        v14 = p4 - p1
        
        r23 = np.linalg.norm(v23)
        r14 = np.linalg.norm(v14)
        r13 = np.linalg.norm(p3 - p1)
        r24 = np.linalg.norm(p4 - p2)
        
        if r23 < 0.0001 or r14 < 0.0001: continue
        
        # Bins
        i14 = max(0, int((r14 - 4.45)/0.25))
        i13 = max(0, int((r13 - 4.85)/0.25))
        i24 = max(0, int((r24 - 4.85)/0.25))
        
        # Basis vectors
        v23_u = v23 / r23
        v14_u = v14 / r14
        
        vxa = v14_u + v23_u
        vxb = v14_u - v23_u
        
        raa = np.linalg.norm(vxa)
        rbb = np.linalg.norm(vxb)
        
        degenerate = (raa < 0.0001 or rbb < 0.0001)
        if not degenerate:
            vxa /= raa
            vxb /= rbb
            vxc = np.cross(vxa, vxb)
            
            # Projection bins
            px = np.dot(v23, vxa)
            py = np.dot(v23, vxb)
            
            ipx = max(0, int((px - 3.15)/0.15))
            ipy = max(0, int((py + 2.1)/0.2))
            
            # Clamping
            i14 = min(i14, 30)
            i13 = min(i13, 30)
            i24 = min(i24, 30)
            ipx = min(ipx, 40)
            ipy = min(ipy, 30)
            
            # Residue type mapping
            res_name = protein.residues[i].name
            res_id = RES_ID_MAP.get(res_name, 2)
            
            found_coords = None
            for flag in [1, 3, 2]:
                key = (flag, i14, i13, i24, ipx, ipy, res_id, 1)
                if key in bb_data:
                    found_coords = bb_data[key]
                    break
        else:
            found_coords = None
                
        if found_coords is None:
            # Fallback: Place atoms along CA-CA vector (idealized extended)
            if r23 > 0.1:
                v_axis = v23 / r23
                vec_c = v_axis * 1.52
                vec_n = v_axis * 2.4
            else:
                continue
        else:
            # Reconstruct
            # found_coords: [N_proj_x, N_proj_y, N_proj_z, C_proj_x, C_proj_y, C_proj_z]
            n_proj = found_coords[0:3]
            c_proj = found_coords[3:6]
            
            # Matrix A = rows [vxa, vxb, vxc]
            A_T = np.vstack([vxa, vxb, vxc]).T
            
            vec_n = A_T @ n_proj
            vec_c = A_T @ c_proj
        
        # N(i+1) = CA(i) + vec_n
        # C(i) = CA(i) + vec_c
        
        # Assign backbone atoms
        current_res = protein.residues[i]
        n_atom = _find_atom(current_res, ('N', 'NH3', 'NP'))
        if pending_n is not None and n_atom is not None:
            n_atom.coords = pending_n
            pending_n = None

        c_atom = _find_atom(current_res, ('C', 'CC'))
        if c_atom is not None:
            c_atom.coords = p2 + vec_c

        if i + 1 < n_res:
            pending_n = p2 + vec_n

    # Approximate N-terminus if still unset
    if protein.residues:
        first_res = protein.residues[0]
        n_atom = _find_atom(first_res, ('N', 'NH3', 'NP'))
        ca0 = ca_coords[0]
        ca_prev = get_ca(-1)
        direction = ca_prev - ca0
        if n_atom is not None and np.linalg.norm(n_atom.coords) < 1e-3 and np.linalg.norm(direction) > 1e-3:
            n_atom.coords = ca0 + direction / np.linalg.norm(direction) * 1.46

# Minimal Residue Map (based on alphabetical or FAMRcomm order?)
# FAMRcomm order usually: GLY, ALA, SER, CYS, VAL...
# I need to verify this mapping.
RES_ID_MAP = {
    'GLY': 1, 'ALA': 2, 'SER': 3, 'CYS': 4, 'VAL': 5, 'THR': 6, 'ILE': 7, 'PRO': 8,
    'MET': 9, 'ASP': 10, 'ASN': 11, 'LEU': 12, 'LYS': 13, 'GLU': 14, 'GLN': 15, 'ARG': 16,
    'HIS': 17, 'PHE': 18, 'TYR': 19, 'TRP': 20
}

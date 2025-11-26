import numpy as np
from typing import List, Dict, Tuple, Optional

class Atom:
    def __init__(self, index: int, name: str, atom_type: str, res_name: str, 
                 charge: float, vdw_e: float, vdw_r: float, hb_da: int):
        self.index = index  # 1-based index from Hdd
        self.name = name.strip()
        self.type = atom_type.strip()
        self.res_name = res_name.strip()
        self.charge = charge
        self.vdw_epsilon = vdw_e
        self.vdw_radius = vdw_r
        self.hbond_type = hb_da
        self.coords = np.zeros(3)  # x, y, z
        self.bonds: List['Atom'] = [] # Connected atoms
        
    def __repr__(self):
        return f"<Atom {self.index} {self.name} {self.res_name}>"

class Bond:
    def __init__(self, atom1_idx: int, neighbors: List[int], params: List[float]):
        self.atom1_idx = atom1_idx
        self.neighbors = neighbors # List of atom indices
        self.params = params
        # params: [k_bond1, r0_1, k_bond2, r0_2, ...]?
        # From remo.py: 8 params. 
        # bond['params'] = [k1, r1, k2, r2, k3, r3, k4, r4] for 4 possible neighbors.

class Angle:
    def __init__(self, a1: int, a2: int, a3: int, k: float, theta0: float, kubo: float, rubo: float):
        self.atoms = (a1, a2, a3)
        self.k = k
        self.theta0 = theta0
        self.kubo = kubo
        self.rubo = rubo

class Protein:
    def __init__(self):
        self.atoms: List[Atom] = []
        self.bonds_data: List[Bond] = [] # Raw bond data from Hdd
        self.angles_data: List[Angle] = [] # Raw angle data from Hdd
        self.hbonds_data: List[str] = [] # Raw H-bond lines
        self.ram_data: List[str] = []    # Ramachandran data
        self.residues: List['Residue'] = []
        
    def get_atom(self, idx: int) -> Optional[Atom]:
        # Hdd indices are 1-based. atoms list is 0-based.
        if 1 <= idx <= len(self.atoms):
            return self.atoms[idx-1]
        return None

class Residue:
    def __init__(self, res_id: int, name: str, start_atom_idx: int, end_atom_idx: int):
        self.id = res_id
        self.name = name
        self.start_atom_idx = start_atom_idx
        self.end_atom_idx = end_atom_idx
        self.atoms: List[Atom] = []
        self.sec_struct = 'C' # Default coil
        self.phi = 0.0
        self.psi = 0.0
        self.chain_id: str = 'A'
        self.orig_resnum: Optional[int] = None

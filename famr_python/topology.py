import math

# Constants and Data
from .naming import normalize_atom_name

AA_HEAVY_COUNT = {
    'GLY': 4, 'ALA': 5, 'SER': 6, 'CYS': 6, 'VAL': 7, 'THR': 7, 'ILE': 8, 'PRO': 7,
    'MET': 8, 'ASP': 8, 'ASN': 8, 'LEU': 8, 'LYS': 9, 'GLU': 9, 'GLN': 9, 'ARG': 11,
    'HIS': 10, 'PHE': 11, 'TYR': 12, 'TRP': 14, 'ASX': 8, 'GLX': 9
}

AA_ALL_COUNT = {
    'GLY': 7, 'ALA': 10, 'SER': 11, 'CYS': 11, 'VAL': 16, 'THR': 14, 'ILE': 19, 'PRO': 14,
    'MET': 17, 'ASP': 12, 'ASN': 14, 'LEU': 19, 'LYS': 22, 'GLU': 15, 'GLN': 17, 'ARG': 24,
    'HIS': 17, 'PHE': 20, 'TYR': 21, 'TRP': 24,
    'NTR': 6, 'CTR': 2 # N-term and C-term additions
}

def parse_pdb_sequence(pdb_source):
    """Parses CA atoms from a PDB path or text stream to get the sequence."""
    seq = []
    res_nums = []

    handle, should_close = _ensure_text_io(pdb_source)
    try:
        for line in handle:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                res_name = line[17:20]
                res_num = line[22:26]
                seq.append(res_name)
                res_nums.append(res_num)
    finally:
        if should_close:
            handle.close()
    return seq, res_nums


def _ensure_text_io(source):
    if hasattr(source, "read"):
        return source, False
    return open(source, "r"), True

def read_famr_comm(comm_file):
    """Reads the FAMRcomm parameter file."""
    with open(comm_file, 'r') as f:
        lines = f.readlines()

    # Parse header
    header = lines[0].split()
    nao = int(header[0])
    nang = int(header[1])
    nhbond = int(header[2])

    atoms = []
    bonds = []
    angles = {}
    hbonds = []
    rams = []

    idx = 1
    # Read Atoms
    entries_read = 0
    while entries_read < nao:
        if idx >= len(lines): break
        line = lines[idx].strip()
        idx += 1
        entries_read += 1
        
        if not line or line.startswith('#'):
            continue
            
        parts = line.split()
        if len(parts) < 8: continue
        
        atom_name = normalize_atom_name(parts[1], parts[3])
        atoms.append({
            'atom': atom_name,
            'type': parts[2],
            'res': parts[3],
            'charge': float(parts[4]),
            'vdwE': float(parts[5]),
            'vdwR': float(parts[6]),
            'hb_da': int(float(parts[7])) 
        })

    nao = len(atoms)

    # Read Bonds
    count = 0
    while count < nao:
        if idx >= len(lines): break
        line = lines[idx].strip()
        idx += 1
        if not line or line.startswith('#'): 
            continue
            
        parts = line.split()
        if len(parts) < 13: continue
        
        bonds.append({
            'connects': [int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])],
            'params': [float(x) for x in parts[5:13]]
        })
        count += 1
    
    # Read Angles  
    for _ in range(nang):
        if idx >= len(lines): break
        line = lines[idx].strip()
        idx += 1
        
        if not line or line.startswith('#'):
            continue
            
        parts = line.split()
        key1 = (parts[0], parts[1], parts[2])
        key2 = (parts[2], parts[1], parts[0])
        val = [float(x) for x in parts[3:7]]
        angles[key1] = val
        angles[key2] = val
        
    # Read HBonds
    count = 0
    while count < nhbond:
        if idx >= len(lines): break
        line = lines[idx].rstrip()
        idx += 1
        hbonds.append(line)
        count += 1
        
    # Skip/Read the comment line before RAM
    if idx < len(lines):
        ram_comment = lines[idx].rstrip()
        idx += 1
    else:
        ram_comment = ""
    
    # Read RAM
    nRAM = 1370
    for _ in range(nRAM):
        if idx < len(lines):
            rams.append(lines[idx].rstrip())
            idx += 1
            
    return {'atoms': atoms, 'bonds': bonds, 'angles': angles, 'hbonds': hbonds, 'rams': rams, 'ram_comment': ram_comment, 'nao': nao}

def generate_hdd(seq, comm_data, output_file='Hdd'):
    """Generates the Hdd topology file."""
    atoms_def = comm_data['atoms']
    bonds_def = comm_data['bonds']
    angles_def = comm_data['angles']
    
    generated_atoms = []
    generated_bonds = []
    
    natom = 0
    
    # Iterate sequence
    for i, res in enumerate(seq):
        res_idx = i + 1
        nheavy = AA_HEAVY_COUNT.get(res, 0)
        
        # Get base definitions
        res_atoms_indices = [idx for idx, a in enumerate(atoms_def) if a['res'] == res]
        
        # Current residue atoms container
        curr_atoms = []
        curr_bonds = []
        
        # 1. Load base atoms
        for comm_idx in res_atoms_indices:
            atom_info = atoms_def[comm_idx].copy()
            bond_info = bonds_def[comm_idx].copy()
            
            # Update bond connections (local relative mapping)
            new_connects = [0]*4
            for k in range(4):
                val = bond_info['connects'][k]
                if val > 0.5:
                    if val < 99990:
                        # Relative mapping
                        offset = val - (comm_idx + 1) # comm_idx is 0-based
                        # curr_atoms contains atoms added so far in this residue
                        current_local_idx = len(curr_atoms) 
                        # natom is global start of residue (0-based count of previous atoms)
                        # current global index = natom + current_local_idx + 1
                        new_connects[k] = (natom + current_local_idx + 1) + offset
                    else:
                        # Inter-residue placeholder
                        new_connects[k] = val
            
            curr_atoms.append(atom_info)
            curr_bonds.append({'connects': new_connects, 'params': bond_info['params']})
            
        # 2. Handle N-term (Residue 1)
        if res_idx == 1:
            if res == 'PRO':
                # PRO N-term
                curr_atoms[0]['atom'] = 'NP'
                curr_atoms[1]['atom'] = 'CP1'
                cd_idx = next((k for k, a in enumerate(curr_atoms) if a['atom'] == 'CD'), -1)
                if cd_idx != -1: curr_atoms[cd_idx]['atom'] = 'CP3'
                
                curr_atoms[0]['charge'] = -0.07
                curr_atoms[1]['charge'] = 0.16
                if cd_idx != -1: curr_atoms[cd_idx]['charge'] = 0.16
                
                # Insert 2 H atoms
                h_insert_idx = nheavy
                h1 = {'atom': 'HC', 'type': 'HC', 'res': res, 'charge': 0.24, 'vdwE': -0.046, 'vdwR': 0.2245, 'hb_da': -1}
                curr_atoms.insert(h_insert_idx, h1)
                curr_atoms.insert(h_insert_idx + 1, h1.copy())
                
                curr_bonds.insert(h_insert_idx, {'connects': [0]*4, 'params': [0.0]*8})
                curr_bonds.insert(h_insert_idx + 1, {'connects': [0]*4, 'params': [0.0]*8})
                
                # Update N bonds
                curr_bonds[0]['connects'][0] = natom + nheavy # CD (approx)
                curr_bonds[0]['connects'][1] = natom + 2 # CA
                curr_bonds[0]['connects'][2] = natom + nheavy + 1 # H1
                curr_bonds[0]['connects'][3] = natom + nheavy + 2 # H2
                curr_bonds[0]['params'] = [320.0, 1.502, 320.0, 1.485, 460.0, 1.006, 460.0, 1.006]
                
                # Update H bonds
                curr_bonds[h_insert_idx]['connects'][0] = natom + 1
                curr_bonds[h_insert_idx]['params'][0] = 460.0; curr_bonds[h_insert_idx]['params'][1] = 1.006
                curr_bonds[h_insert_idx+1]['connects'][0] = natom + 1
                curr_bonds[h_insert_idx+1]['params'][0] = 460.0; curr_bonds[h_insert_idx+1]['params'][1] = 1.006
                
            else:
                # Non-PRO N-term
                curr_atoms[0]['atom'] = 'NH3'
                curr_atoms[0]['charge'] = -0.30
                if res == 'GLY': curr_atoms[1]['charge'] = 0.13
                else: curr_atoms[1]['charge'] = 0.21
                
                h_insert_idx = nheavy
                h_new = {'atom': 'HC', 'type': 'HC', 'res': res, 'charge': 0.33, 'vdwE': -0.046, 'vdwR': 0.2245, 'hb_da': -1}
                
                curr_atoms.insert(h_insert_idx, h_new.copy())
                curr_atoms.insert(h_insert_idx, h_new.copy())
                
                curr_bonds.insert(h_insert_idx, {'connects': [0]*4, 'params': [0.0]*8})
                curr_bonds.insert(h_insert_idx, {'connects': [0]*4, 'params': [0.0]*8})
                
                curr_atoms[h_insert_idx]['atom'] = normalize_atom_name('1H')
                curr_atoms[h_insert_idx+1]['atom'] = normalize_atom_name('2H')
                curr_atoms[h_insert_idx+2]['atom'] = normalize_atom_name('3H')
                
                curr_bonds[0]['connects'][0] = natom + nheavy + 1
                curr_bonds[0]['connects'][1] = natom + 2
                curr_bonds[0]['connects'][2] = natom + nheavy + 2
                curr_bonds[0]['connects'][3] = natom + nheavy + 3
                curr_bonds[0]['params'] = [403.0, 1.04, 200.0, 1.48, 403.0, 1.04, 403.0, 1.04]
                
                for k in range(3):
                    idx_h = h_insert_idx + k
                    curr_bonds[idx_h]['connects'][0] = natom + 1
                    curr_bonds[idx_h]['params'][0] = 403.0; curr_bonds[idx_h]['params'][1] = 1.04

        # 3. Handle C-term (Last Residue)
        if res_idx == len(seq):
            c_idx = 2; o_idx = 3
            curr_atoms[c_idx]['atom'] = 'CC'; curr_atoms[c_idx]['charge'] = 0.34
            curr_atoms[c_idx]['vdwE'] = -0.07; curr_atoms[c_idx]['vdwR'] = 2.00
            curr_atoms[o_idx]['atom'] = 'OC'; curr_atoms[o_idx]['charge'] = -0.67
            
            oxt_idx = nheavy
            oxt = {'atom': 'OXT', 'type': 'OC', 'res': res, 'charge': -0.67, 'vdwE': -0.12, 'vdwR': 1.7, 'hb_da': 1}
            curr_atoms.insert(oxt_idx, oxt)
            curr_bonds.insert(oxt_idx, {'connects': [0]*4, 'params': [0.0]*8})
            
            curr_bonds[c_idx]['connects'][0] = natom + 2 # CA
            curr_bonds[c_idx]['connects'][1] = natom + 4 # O
            curr_bonds[c_idx]['connects'][2] = natom + oxt_idx + 1 # OXT
            
            if res == 'PRO':
                curr_bonds[c_idx]['params'][0] = 250.0; curr_bonds[c_idx]['params'][1] = 1.49
            else:
                curr_bonds[c_idx]['params'][0] = 200.0; curr_bonds[c_idx]['params'][1] = 1.522
            curr_bonds[c_idx]['params'][2] = 525.0; curr_bonds[c_idx]['params'][3] = 1.26
            curr_bonds[c_idx]['params'][4] = 525.0; curr_bonds[c_idx]['params'][5] = 1.26
            
            curr_bonds[oxt_idx]['connects'][0] = natom + 3 # C
            curr_bonds[oxt_idx]['params'][0] = 525.0; curr_bonds[oxt_idx]['params'][1] = 1.26

        # 4. Resolve Connections and Add to Global
        res_start_global = natom
        
        for k, bond in enumerate(curr_bonds):
            final_connects = []
            for c_val in bond['connects']:
                if c_val == 0:
                    final_connects.append(0)
                    continue
                
                if c_val < 99990:
                    # Local connection
                    # Logic fixed in previous turn
                    # c_val is global index but we modified it above.
                    # Wait, logic above: new_connects[k] = (natom + current_local_idx + 1) + offset
                    # That was updating `new_connects` which is `bond['connects']` essentially.
                    # But here we iterate `bond['connects']` again.
                    # The `new_connects` created in Step 1 are ALREADY global indices (roughly).
                    # Wait, let's re-examine Step 1 loop.
                    # new_connects[k] = (natom + current_local_idx + 1) + offset
                    # Yes, this is global index.
                    
                    # However, in Step 4, we check `c_val`.
                    # If we inserted atoms (H atoms), we might need to shift indices?
                    # Step 1 uses `natom` which is start of residue.
                    # If we insert atoms at `h_insert_idx` (Step 2/3), the indices of atoms *after* insertion point shift.
                    # But `new_connects` calculated in Step 1 point to specific atom indices.
                    # If `offset` pointed to an atom that is now shifted, we need to update.
                    
                    # But `offset` is relative to `comm_idx` (original FAMRcomm index).
                    # If `comm_idx` corresponds to an atom that is now shifted, how do we know?
                    # We inserted atoms into `curr_atoms`.
                    # The atom at `local_idx` might have moved to `local_idx + shift`.
                    
                    # Actually, the logic in Step 4 of `remo.py` (original) handled this shift:
                    # if res_idx == 1 (N-term):
                    #   if c_val > natom + shift_start: c_val += shift_amt
                    
                    # So we DO need to apply shift here.
                    
                    # N-term insertions:
                    if res_idx == 1:
                        shift_start = nheavy
                        shift_amt = 2 
                        if c_val > natom + shift_start:
                            c_val += shift_amt
                            
                    # C-term insertions:
                    if res_idx == len(seq):
                        shift_start = nheavy
                        shift_amt = 1
                        if c_val > natom + shift_start:
                            c_val += shift_amt
                            
                    final_connects.append(c_val)
                else:
                    # Inter-residue
                    if k == 0 and curr_atoms[k]['atom'].startswith('N'): # N connecting to Prev C
                        if res_idx > 1:
                            final_connects.append(-1) 
                        else:
                            final_connects.append(0)
                    elif k == 2 and curr_atoms[k]['atom'].startswith('C'): # C connecting to Next N
                         if res_idx < len(seq):
                             final_connects.append(-2)
                         else:
                             final_connects.append(0)
                    else:
                        final_connects.append(0)

            bond['connects'] = final_connects
        
        generated_atoms.extend(curr_atoms)
        generated_bonds.extend(curr_bonds)
        natom += len(curr_atoms)
        
    # Second pass: Fix inter-residue connections
    res_boundaries = [] # (start_idx, end_idx, n_idx, c_idx)
    
    curr_res = None
    start = 0
    n_idx = -1
    c_idx = -1
    
    for i, atom in enumerate(generated_atoms):
        if atom['res'] != curr_res:
            if curr_res is not None:
                res_boundaries.append({'start': start, 'end': i-1, 'n': n_idx, 'c': c_idx})
            curr_res = atom['res']
            start = i
            n_idx = i # N is first
            c_idx = i + 2 # C is usually 3rd (index 2)
            # Verify C
            if generated_atoms[c_idx]['atom'] not in ['C', 'CC']:
                # Search for C
                for k in range(start, min(start+10, len(generated_atoms))):
                    if generated_atoms[k]['atom'] in ['C', 'CC']:
                        c_idx = k
                        break
    # Last res
    res_boundaries.append({'start': start, 'end': len(generated_atoms)-1, 'n': n_idx, 'c': c_idx})
    
    # Fix connections
    for i, bond in enumerate(generated_bonds):
        new_connects = []
        for c in bond['connects']:
            if c == -1: # Prev C
                res_idx = -1
                for r_i, r in enumerate(res_boundaries):
                    if i >= r['start'] and i <= r['end']:
                        res_idx = r_i
                        break
                if res_idx > 0:
                    prev_c = res_boundaries[res_idx-1]['c'] + 1 # 1-based
                    new_connects.append(prev_c)
                else:
                    new_connects.append(0)
            elif c == -2: # Next N
                res_idx = -1
                for r_i, r in enumerate(res_boundaries):
                    if i >= r['start'] and i <= r['end']:
                        res_idx = r_i
                        break
                if res_idx < len(res_boundaries) - 1:
                    next_n = res_boundaries[res_idx+1]['n'] + 1 # 1-based
                    new_connects.append(next_n)
                else:
                    new_connects.append(0)
            else:
                new_connects.append(c)
        bond['connects'] = new_connects

    # 5. Generate Angles
    # Angles generation logic is simple: for each bond pair, check if angle exists in definition
    generated_angles = []
    for i in range(1, natom + 1):
        atom1_idx = i
        atom1_type = generated_atoms[atom1_idx-1]['type']
        bonds1 = generated_bonds[atom1_idx-1]['connects']
        for neighbor in bonds1:
            if neighbor <= 0 or neighbor > natom: continue
            if neighbor == atom1_idx: continue
            atom2_idx = neighbor
            atom2_type = generated_atoms[atom2_idx-1]['type']
            bonds2 = generated_bonds[atom2_idx-1]['connects']
            for neighbor2 in bonds2:
                if neighbor2 <= 0 or neighbor2 > natom: continue
                if neighbor2 == atom2_idx: continue
                if neighbor2 == atom1_idx: continue
                atom3_idx = neighbor2
                atom3_type = generated_atoms[atom3_idx-1]['type']
                if atom3_idx <= atom1_idx: continue
                
                key = (atom1_type, atom2_type, atom3_type)
                if key in angles_def: params = angles_def[key]
                else: params = [0.0, 0.0, 0.0, 0.0]
                
                generated_angles.append({'atoms': [atom1_idx, atom2_idx, atom3_idx], 'params': params})

    # 6. Write Hdd
    with open(output_file, 'w') as f:
        f.write(f"{natom:6d} {len(generated_angles):6d} {len(comm_data['hbonds'])-1:6d}\n")
        for i, atom in enumerate(generated_atoms):
            idx = i + 1
            if idx > 9999: idx = 9999
            atom_name = normalize_atom_name(atom['atom'], atom['res'])
            if len(atom_name) == 1: atom_name = f" {atom_name}  "
            elif len(atom_name) == 2: atom_name = f" {atom_name} "
            elif len(atom_name) == 3:
                if atom_name[0].isdigit(): atom_name = f"{atom_name} "
                else: atom_name = f" {atom_name}"
            f.write(f"{idx:4d} {atom_name:4s} {atom['type']:4s} {atom['res']:4s} "
                    f"{atom['charge']:6.3f} {atom['vdwE']:7.4f} {atom['vdwR']:7.4f} {atom['hb_da']:4d}\n")
        for i, bond in enumerate(generated_bonds):
            connects = bond['connects']
            params = bond['params']
            nbond = sum(1 for c in connects if c > 0)
            idx = i + 1
            if idx > 9999: idx = 9999
            f.write(f"{idx:4d} {connects[0]:6d} {connects[1]:6d} {connects[2]:6d} {connects[3]:6d} "
                    f"{params[0]:8.3f} {params[1]:7.4f} {params[2]:8.3f} {params[3]:7.4f} "
                    f"{params[4]:8.3f} {params[5]:7.4f} {params[6]:8.3f} {params[7]:7.4f} {nbond:2d}\n")
        for i, ang in enumerate(generated_angles):
            idx = i + 1
            if idx > 99999: idx = 99999
            atoms = ang['atoms']
            params = ang['params']
            f.write(f"{idx:5d} {atoms[0]:5d} {atoms[1]:5d} {atoms[2]:5d} "
                    f"{params[0]:7.3f} {params[1]:7.3f} {params[2]:7.3f} {params[3]:7.3f}\n")
        
        for line in comm_data['hbonds'][1:-1]:
            f.write(f"{line}\n")
            
        if 'ram_comment' in comm_data and comm_data['ram_comment']:
            f.write(f"{comm_data['ram_comment']}\n")
            
        for line in comm_data['rams']:
            f.write(f"{line}\n")

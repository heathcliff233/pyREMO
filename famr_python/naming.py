GENERAL_OVERRIDES = {
    '1H': 'H1',
    '2H': 'H2',
    '3H': 'H3',
    '1HA': 'HA1',
    '2HA': 'HA2',
    '3HA': 'HA3',
    '1HB': 'HB1',
    '2HB': 'HB2',
    '3HB': 'HB3',
    '1HD': 'HD1',
    '2HD': 'HD2',
    '3HD': 'HD3',
    '1HG': 'HG1',
    '2HG': 'HG2',
    '3HG': 'HG3',
    '1HE': 'HE1',
    '2HE': 'HE2',
    '3HE': 'HE3',
    '1HZ': 'HZ1',
    '2HZ': 'HZ2',
    '3HZ': 'HZ3',
    '1HH': 'HH1',
    '2HH': 'HH2',
}

RESIDUE_OVERRIDES = {
    'ILE': {
        '1HG1': 'HG12',
        '2HG1': 'HG13',
    }
}

def normalize_atom_name(name: str, residue: str = "") -> str:
    """
    Normalize atom names from FAMRcomm to PDB-friendly versions.
    """
    if not name:
        return name
    trimmed = name.strip()
    res = (residue or "").upper()
    if res in RESIDUE_OVERRIDES and trimmed in RESIDUE_OVERRIDES[res]:
        return RESIDUE_OVERRIDES[res][trimmed]
    if trimmed in GENERAL_OVERRIDES:
        return GENERAL_OVERRIDES[trimmed]
    if trimmed and trimmed[0].isdigit():
        trimmed = trimmed[1:] + trimmed[0]
    return trimmed[:4]

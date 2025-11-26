import numpy as np
import math
from typing import Tuple

def distance(p1: np.ndarray, p2: np.ndarray) -> float:
    """Euclidean distance between two points."""
    return np.linalg.norm(p1 - p2)

def angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Calculate angle (in radians) defined by p1-p2-p3."""
    v1 = p1 - p2
    v2 = p3 - p2
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
        
    cosine = np.dot(v1, v2) / (norm1 * norm2)
    # Clamp for numerical stability
    cosine = max(-1.0, min(1.0, cosine))
    return np.arccos(cosine)

def dihedral(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> float:
    """Calculate dihedral angle (in radians) defined by p1-p2-p3-p4."""
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    
    # Normalize b1 so that it does not influence magnitude of vector rejections
    b1 /= np.linalg.norm(b1)
    
    # Vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    # w = projection of b2 onto plane perpendicular to b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    
    # Angle between v and w
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    
    return np.arctan2(y, x)

def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc - ad), 2 * (bd + ac)],
        [2 * (bc + ad), aa + cc - bb - dd, 2 * (cd - ab)],
        [2 * (bd - ac), 2 * (cd + ab), aa + dd - bb - cc]
    ])

def superimpose_backbone(target_n: np.ndarray, target_ca: np.ndarray, target_c: np.ndarray,
                        mobile_n: np.ndarray, mobile_ca: np.ndarray, mobile_c: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate rotation and translation to superimpose mobile onto target.
    Uses Kabsch algorithm or simpler 3-point alignment since we have N-CA-C.
    Returns (R, t) such that p_aligned = dot(p_mobile, R.T) + t
    """
    # Center both sets of points on CA (or centroid)
    target_points = np.vstack([target_n, target_ca, target_c])
    mobile_points = np.vstack([mobile_n, mobile_ca, mobile_c])
    
    centroid_tgt = np.mean(target_points, axis=0)
    centroid_mob = np.mean(mobile_points, axis=0)
    
    tgt_centered = target_points - centroid_tgt
    mob_centered = mobile_points - centroid_mob
    
    # Covariance matrix
    H = np.dot(mob_centered.T, tgt_centered)
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Rotation matrix
    R = np.dot(Vt.T, U.T)
    
    # Handle reflection
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = np.dot(Vt.T, U.T)
        
    # Translation
    t = centroid_tgt - np.dot(centroid_mob, R.T) # Wait, dot order depends on row/col vectors
    # If p is row vector: p_new = p @ R + t.
    # If p is col vector: p_new = R @ p + t.
    # numpy default is row vectors usually if (N, 3).
    # Here I used (3, 3).
    # mob_centered @ R should approx tgt_centered.
    
    # Let's verify:
    # p' = p @ R
    # H was mob.T @ tgt
    # R should be such that mob @ R ~ tgt.
    
    return R, centroid_tgt - np.dot(centroid_mob, R) 

def place_atom(a1: np.ndarray, a2: np.ndarray, a3: np.ndarray, 
               bond_len: float, angle_rad: float, torsion_rad: float) -> np.ndarray:
    """
    Place atom A4 given A1-A2-A3 and internal coordinates.
    a1: Great-grandparent
    a2: Grandparent
    a3: Parent
    bond_len: dist(A3, A4)
    angle_rad: angle(A2, A3, A4)
    torsion_rad: dihedral(A1, A2, A3, A4)
    """
    # 1. Define local frame at A3
    # v1 = A3 - A2 (normalized)
    v1 = a3 - a2
    v1 /= np.linalg.norm(v1)
    
    # v2 in plane A1-A2-A3, perpendicular to v1
    # v_plane = A2 - A1
    v_plane = a2 - a1
    # Component orthogonal to v1
    v2 = v_plane - np.dot(v_plane, v1) * v1
    if np.linalg.norm(v2) < 1e-6:
        # Collinear A1-A2-A3. Pick arbitrary normal.
        # If v1 is z, pick x.
        if abs(v1[2]) < 0.9:
            v2 = np.array([0, 0, 1.0])
        else:
            v2 = np.array([1.0, 0, 0])
        v2 = v2 - np.dot(v2, v1) * v1
        
    v2 /= np.linalg.norm(v2)
    
    # v3 = v1 x v2
    v3 = np.cross(v1, v2)
    
    # 2. Construct A4 in local frame
    # x = r * sin(angle) * cos(torsion)
    # y = r * sin(angle) * sin(torsion)
    # z = -r * cos(angle)  (Since v1 points A2->A3, and angle is A2-A3-A4.
    # If A4 is "forward", angle is > 90? No.
    # Standard definition: angle is 0 if A4 extends A2-A3 linearly? 
    # No, angle usually A2-A3-A4. 180 is linear extension.
    # If angle is 180, cos=-1. z = r.
    # So z = -r * cos(angle).
    
    # Wait.
    # v1 is A2->A3.
    # Extension is along v1.
    # Angle theta = A2-A3-A4.
    # Deviation from A2-A3 line is (180 - theta).
    # If theta=120 (standard sp3), deviation 60.
    # Component along v1: r * cos(180-theta) = r * (-cos(theta)).
    # Component along v2/v3 plane: r * sin(180-theta) = r * sin(theta).
    
    # So:
    # projection_on_v1 = -bond_len * np.cos(angle_rad)
    # projection_on_plane = bond_len * np.sin(angle_rad)
    
    # In plane components:
    # Torsion is rotation around v1 (A2-A3 bond).
    # Reference vector for torsion?
    # Usually plane A1-A2-A3 is torsion=0 (cis) or 180 (trans).
    # Definition: Dihedral A1-A2-A3-A4.
    # If torsion=0, A4 is cis to A1 (same side).
    # v2 was defined from A2-A1.
    # v2 points roughly A2->A1.
    # If A4 is cis to A1, it projects onto v2 (positive)?
    # Check v2 direction.
    # v_plane = A2 - A1. (Points back to A1).
    # v2 is part of A2-A1 orthogonal to A2-A3.
    # So v2 points towards "inside" of bend A1-A2-A3.
    # If torsion=0 (cis), A4 should point towards A1.
    # So positive v2 component.
    
    x_local = bond_len * np.sin(angle_rad) * np.cos(torsion_rad) # along v2
    y_local = bond_len * np.sin(angle_rad) * np.sin(torsion_rad) # along v3
    z_local = -bond_len * np.cos(angle_rad) # along v1
    
    pos = a3 + x_local * v2 + y_local * v3 + z_local * v1
    return pos 

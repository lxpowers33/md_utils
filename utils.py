import numpy as np
from IPython import embed

####################################################
#####                Basic                     #####
####################################################
def _sliding_mean(data_array, window):
    """Robin's smoothing function"""
    new_list = []
    for i in range(data_array.shape[0]):
        indices = range(max(i - window + 1, 0),
                        min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)

def _sliding_median(data_array, window):
    """Robin's smoothing function"""
    new_list = []
    for i in range(data_array.shape[0]):
        indices = range(max(i - window + 1, 0),
                        min(i + window + 1, len(data_array)))
        avg = np.median(data_array[indices])
        new_list.append(avg)
    return np.array(new_list)

def _vectorize_coords(atomsel):
    """
    Returns N_atoms x 3 np.array
    with columns x, y, z coordinate
    """
    return np.vstack([atomsel.get(direction)
                      for direction in ['x', 'y', 'z']]).T

def _calcdist(atomsel1, atomsel2):
    """
    Calculates the minimum distance between the two atomsels.
    
    Errors if either is empty.
    """
    XYZ1 = _vectorize_coords(atomsel1)
    XYZ2 = _vectorize_coords(atomsel2)
    return min(np.linalg.norm(xyz1-xyz2)
               for xyz1 in XYZ1 for xyz2 in XYZ2)

####################################################
#####                Geometry                  #####
####################################################

def v_projection(a, b):
    """
    The vector projection of a onto b 
    Parameters
    a: (numpy array nx1) the vector to project
    b: (numpy array nx1) the vector to project onto
    Return: scalar
    """
    s = np.sum(a*b)  / np.linalg.norm(b)
    return s


def measure_dihedral_vec(p0, p1, p2, p3):
    """Function from 
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    vectorized by Alex Powers
    each p should be a vector nx3
    returns nx1 array of dihedral angles
    """
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1, axis=1)[:, np.newaxis] 
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - npsumdot(b0, b1)*b1
    w = b2 - npsumdot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = npsumdot(v, w)
    y = npsumdot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def measure_dihedral(p0, p1, p2, p3):
    """Original formula from stackoverflow
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    each p should be a vector 3 by 1"""
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def npsumdot(x, y):
    out = np.sum(x*y, axis=1)
    return out[:, np.newaxis] 

def eindot(x, y):
    return np.einsum('ij,ij->i', x, y)

def test_dihedrals():

    p0 = np.array([24.969, 13.428, 30.692]) # N
    p1 = np.array([24.044, 12.661, 29.808]) # CA
    p2 = np.array([22.785, 13.482, 29.543]) # C
    p3 = np.array([21.951, 13.670, 30.431]) # O
    p4 = np.array([23.672, 11.328, 30.466]) # CB
    p5 = np.array([22.881, 10.326, 29.620]) # CG
    p6 = np.array([23.691,  9.935, 28.389]) # CD1
    p7 = np.array([22.557,  9.096, 30.459]) # CD2

    v_1 = np.vstack((p0, p0, p1, p1))
    v_2 = np.vstack((p1, p1, p4, p4))
    v_3 = np.vstack((p2, p4, p5, p5))
    v_4 = np.vstack((p3, p5, p6, p7))

    out = measure_dihedral_vec(v_1, v_2, v_3, v_4)

    assert(abs(out[0] - (-71.21515)) < 1E-4)
    assert(abs(out[1] - (-171.94319)) < 1E-4)
    assert(abs(out[2] - (60.82226)) < 1E-4)
    assert(abs(out[3] - (-177.63641)) < 1E-4)

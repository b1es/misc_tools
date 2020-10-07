import numpy as np


def get_diam_phys_parametrized(data, value, spacing):
    """
    Get diameter of selected circular plane
    :param data: 3D npy array
    :param value: value the area is marked with
    :param spacing:
    :return: physical diameter
    :return: eccentricity [0-1] where 0 is a circle
    """
    coords = [(np.min(dim), np.max(dim) + 1) for dim in np.where(data == value)]
    if np.sum([1 == (max_ - min_) for min_, max_ in coords]) != 1:
        raise ValueError('Multiple inflow layers')
    inflow = np.squeeze(data[coords[0][0]:coords[0][1], coords[1][0]:coords[1][1], coords[2][0]:coords[2][1]])

    from skimage.measure import regionprops
    props = regionprops(inflow)[0]
    return props.equivalent_diameter * spacing, props.eccentricity


def coords_transform(coordinates, origin, spacing, mode):
    """
    Transform ijk into xyz or vice versa.
    :param coordinates: 1D or 2D array of coordinates
    :param mode: 'ijk_to_xyz' or 'xyz_to_ijk'
    """
    if mode not in ['ijk_to_xyz', 'xyz_to_ijk']:
        raise ValueError('Wrong mode')
        
    if mode == 'ijk_to_xyz':
        data = np.array(origin) + np.array(spacing) * np.array(coordinates).reshape((-1, 3))[:,::-1]
        if np.array(coordinates).ndim == 1:
            data = np.squeeze(data)
        return data
    if mode == 'xyz_to_ijk':
        data = np.rint(
                ((np.array(coordinates).reshape((-1, 3)) - np.array(origin)) / np.array(spacing))[:,::-1]
               ).astype(np.int)
        if np.array(coordinates).ndim == 1:
            data = np.squeeze(data)
        return data
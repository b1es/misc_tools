import numpy as np


def get_diam_phys_parametrized(data, value, spacing):
    """
    Get diameter of selected circular plane
    :param data: 3D npy array
    :param value: value the area is marked with
    :param spacing:
    :return: physical diameter
    """
    coords = [(np.min(dim), np.max(dim) + 1) for dim in np.where(data == value)]
    if np.sum([1 == (max_ - min_) for min_, max_ in coords]) != 1:
        raise ValueError('Multiple inflow layers')
    inflow = np.squeeze(data[coords[0][0]:coords[0][1], coords[1][0]:coords[1][1], coords[2][0]:coords[2][1]])

    from skimage.measure import regionprops
    return regionprops(inflow)[0].equivalent_diameter * spacing

"""Utility functions."""

import numpy as np 

__all__ = ['spherical_to_cartesian', 'search_around_sky']


def spherical_to_cartesian(ra, dec):
    """Convert spherical coordinates into cartesian coordinates on a unit
    sphere.

    Parameters
    ----------
    ra, dec : float or numpy array
        Spherical coordinates.

    Returns
    -------
    x, y, z : float or numpy array
        Cartesian coordinates.
    """
    x = np.cos(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
    y = np.sin(np.deg2rad(ra)) * np.cos(np.deg2rad(dec))
    z = np.sin(np.deg2rad(dec))
    return x, y, z

def search_around_sky(ra, dec, kdtree, rmin, rmax):
    """Cross-match coordinates.

    Parameters
    ----------
    ra, dec : float
        Coordinates around which objects are searched.
    kdtree : scipy.spatial.cKDTree
        KDTree containing the objects to be searched.
    rmin : float, optional
        Minimum radius to search for objects in degrees.
    rmax : float, optional
        Maximum radius to search for objects in degrees.

    Returns
    -------
    idx : numpy array
        Indices of objects in the KDTree that lie in the search area.
    theta : numpy array
        The distance to the search center in degrees.
    """
    # Convert the angular distance into a 3D distance on the unit sphere.
    rmax_3d = np.sqrt(2 - 2 * np.cos(np.deg2rad(rmax)))

    x, y, z = spherical_to_cartesian(ra, dec)
    idx = np.array(kdtree.query_ball_point([x, y, z], rmax_3d))

    # Convert 3D distance back into angular distances.
    if len(idx) > 0:
        dist3d = np.sqrt(np.sum((np.array([x, y, z]) - kdtree.data[idx])**2,
                                axis=1))
    else:
        dist3d = np.zeros(0)
    theta = np.rad2deg(np.arcsin(dist3d * np.sqrt(1 - dist3d**2 / 4)))
    mask = (theta > rmin)

    return idx[mask], theta[mask]

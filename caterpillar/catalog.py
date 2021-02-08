"""Functions to deal with catalogs."""

import os
import copy
import warnings

import numpy as np

from astropy.table import Table, Column, vstack

__all__ = ["remove_is_null", "moments_to_shape"]


def remove_is_null(table, output=None, verbose=True, string='isnull', return_data=True):
    """
    Remove the xxx_isnull columns from the catalog.
    This is an annoying issue with FITS table from HSC database.

    Parameters
    ----------
    table : str or astropy.table object
        Name of the FITS catalog or the data itself
    cat_hdu : int, optional
        The HDU of the catalog data.
        Default: 1
    string : str, optional
        The name of the Null columns.
        Default: 'isnull'
    output : str, optional
        If output is None, will write a new FITS table with '_clean'
        suffix.
        Default: None
    return_data : bool, optional
        Whether return the cleaned data.
        Default : True
    verbose : bool, optional
        Default : True
    """
    if isinstance(table, Table):
        data = table
    else:
        if not os.path.isfile(table):
            raise Exception("# Can not find catalog: %s" % table)
        data = Table.read(table, format='fits')

    if verbose:
        print("Reading the data....")
    col_names = data.colnames
    col_remove = [col for col in col_names if string in col]
    data.remove_columns(col_remove)

    if output is None:
        if verbose:
            print("Saving data to %s ..." % output)
        data.write(table, format='fits', overwrite=True)
    else:
        data.write(output.strip(), format='fits', overwrite=True)

    if return_data:
        return data
    else:
        return None


def moments_to_shape(catalog, shape_type='i_sdss_shape', axis_ratio=False,
                     radian=False, update=True, to_pixel=False):
    """
    Convert the 2nd moments into elliptical shape: radius, ellipticity, position angle.
    Adopted from `unagi.catalog`:
        https://github.com/dr-guangtou/unagi/blob/master/unagi/catalog.py

    Parameters
    ----------
    catalog : Table data
    """
    try:
        xx = catalog["{}_11".format(shape_type)]
        yy = catalog["{}_22".format(shape_type)]
        xy = catalog["{}_12".format(shape_type)]
    except KeyError:
        print("Wrong column name!")
        raise

    e1 = (xx - yy) / (xx + yy)
    e2 = (2.0 * xy / (xx + yy))
    # Get the r50 or determinant radius
    rad = np.sqrt(xx + yy)
    rad = rad / 0.168 if to_pixel else rad
    # Ellipticity or axis ratio
    ell = np.sqrt(e1 ** 2.0 + e2 ** 2.0)
    ell = 1.0 - ell if axis_ratio else ell
    # Position angle in degree or radian
    theta = (-0.5 * np.arctan2(e2, e1))
    theta = (theta * 180. / np.pi) if not radian else theta

    if update:
        rad_col = "{}_r".format(shape_type)
        theta_col = "{}_theta".format(shape_type)
        if axis_ratio:
            ell_col = "{}_ba".format(shape_type)
        else:
            ell_col = "{}_e".format(shape_type)
        if rad_col in catalog.colnames:
            catalog.remove_column(rad_col)
        catalog.add_column(Column(data=rad, name=rad_col))
        if ell_col in catalog.colnames:
            catalog.remove_column(ell_col)
        catalog.add_column(Column(data=ell, name=ell_col))
        if theta_col in catalog.colnames:
            catalog.remove_column(theta_col)
        catalog.add_column(Column(data=theta, name=theta_col))
        return catalog
    return rad, ell, theta

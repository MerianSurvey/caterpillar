"""Functions to deal with catalogs."""

import os
import copy
import distutils.spawn

from shutil import copyfile

import numpy as np
import healpy as hp

from astropy.wcs import WCS
from astropy.table import Table, Column, vstack

__all__ = ["remove_is_null", "moments_to_shape", "filter_through_bright_star_mask",
           "ReferenceCatalog", "PS1_PATTERN"]

PS1_PATTERN = 'ps1-{:05d}.fits'

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
        ell = np.asarray(ell)
        catalog.add_column(Column(data=ell, name=ell_col))
        if theta_col in catalog.colnames:
            catalog.remove_column(theta_col)
        theta = np.asarray(theta)
        catalog.add_column(Column(data=theta, name=theta_col))
        return catalog
    return rad, ell, theta

def filter_through_bright_star_mask(catalog, mask_dir, reg_prefix='new_S18Amask',
                                    filters='grizy', filter_type='outside',
                                    ra='ra', dec='dec', output_suffix='bsm'):
    """Filter the catalog through the .reg files of the bright star masks."""
    # Make the sure venice is installed
    venice = distutils.spawn.find_executable("venice")
    assert venice, "Venice is not installed!"

    # Get the .reg files for the bright star mask
    reg_files = [
        os.path.join(mask_dir, reg_prefix + '_' + band + '.reg') for band in filters]

    # Output catalog
    output_catalogs = [
        catalog.replace('.fits', '_bsm_' + band + '.fits') for band in filters]

    output_final = catalog.replace('.fits', '_%s.fits' % output_suffix)

    # Generate the commands
    for ii, reg_mask in enumerate(reg_files):
        if ii == 0:
            venice_command = (
                venice + ' -m ' + reg_mask + ' -f ' + filter_type + ' -cat ' + catalog +
                ' -xcol ' + ra + ' -ycol ' + dec + ' -o ' + output_catalogs[0]
            )
        else:
            venice_command = (
                venice + ' -m ' + reg_mask + ' -f ' + filter_type + ' -cat ' +
                output_catalogs[ii - 1] + ' -xcol ' + ra + ' -ycol ' + dec +
                ' -o ' + output_catalogs[ii]
            )
        # Execute the command
        _ = os.system(venice_command)

    # Copy the last catalog to the final name
    if not os.path.isfile(output_catalogs[-1]):
        raise Exception("# Something is wrong with the Venice!")
    else:
        _ = copyfile(output_catalogs[-1], output_final)

    # Delete the intermediate catalogs
    for output in output_catalogs:
        try:
            os.remove(output)
        except OSError:
            pass

    return Table.read(output_final)

def add_chunk_id(catalog):
    """Assign chunk ID based on HSC Tract/Patch."""
    chunks = []
    for obj in catalog:
        tract = '{:5d}'.format(obj['tract'])
        patch = '{0:03d}'.format(obj['patch'])
        chunks.append('_'.join([tract, patch[0], patch[2]]))

    catalog.add_column(Column(data=chunks, name='chunk_id'))

    return catalog

class ReferenceCatalog():
    """
    Photometric or astrometric reference catalog.
    """
    def __init__(self, fits_dir, fits_pattern, nside=32, indexing='ring'):
        """
        fits_pattern: string formatter with key "hp", e.g., 'dir/fn-%(hp)05i.fits'
        """
        self.fits_dir = fits_dir
        self.fits_pattern = os.path.join(fits_dir, fits_pattern)
        self.nside = nside
        self.indexing = indexing

    def get_catalogs(self, pixels):
        """
        Get the reference star catalogs.
        """
        ref_cats = []
        for pix in pixels:
            ref_cats.append(Table.read(self.fits_pattern.format(pix)))

        return vstack(ref_cats)

    def get_hp_pixles(self, ra_arr, dec_arr):
        """
        Get the healpix pixels that cover the (RA, Dec) ranges.
        """
        hp_pixels = set()

        # TODO: This is only for ring indexing
        for rr, dd in zip(ra_arr, dec_arr):
            hp_pixels.add(
                hp.pixelfunc.ang2pix(self.nside, rr, dd, lonlat=True, nest=False))

        return hp_pixels

    def get_stars_radec(self, ra_range, dec_range, step=100, margin=0.1):
        """
        Get the reference stars within a (RA, Dec) range.
        """
        ra_min, ra_max = np.min(ra_range) - margin, np.max(ra_range) - margin
        dec_min, dec_max = np.min(dec_range) - margin, np.max(dec_range) - margin

        ra_arr, dec_arr = np.meshgrid(
            np.linspace(ra_min, ra_max, int((ra_max - ra_min) / step)),
            np.linspace(dec_min, dec_max, int((dec_max - dec_min) / step))
        )

        # Get the healpix pixels
        pixels = self.get_hp_pixles(ra_arr, dec_arr)

        # Gather the reference stars
        cat = self.get_catalogs(pixels)

        return cat

    def get_stars_hdu(self, hdu, step=100, margin=10, box_cut=True,
                      ra_col='RA', dec_col='DEC', box_margin=20):
        """
        Get the reference stars based on a FITS HDU.
        """
        # Header of the HDU
        hdu_wcs = WCS(hdu.header)

        # Dimension of the image
        dec_size, ra_size = hdu_wcs.array_shape

        # Build a grid of (X, Y) coordinate
        xx_arr, yy_arr = np.meshgrid(
            np.linspace(1 - margin, ra_size + margin, 2 + int((ra_size + 2 * margin) / step)),
            np.linspace(1 - margin, dec_size + margin, 2 + int((dec_size + 2 * margin) / step)))

        # (RA, Dec) grid
        ra_arr, dec_arr = hdu_wcs.all_pix2world(xx_arr.ravel(), yy_arr.ravel(), 1)

        # Get the healpix pixels
        pixels = self.get_hp_pixles(ra_arr, dec_arr)

        # Gather the reference stars
        cat = self.get_catalogs(pixels)

        # Make a box cut
        if box_cut:
            ra_min, ra_max = ra_arr.min(), ra_arr.max()
            dec_min, dec_max = dec_arr.min(), dec_arr.max()
            ra_sep = (ra_max - ra_min) / box_margin
            dec_sep = (dec_max - dec_min) / box_margin

            cat = cat[
                (cat[ra_col] > ra_min - ra_sep) &
                (cat[ra_col] < ra_max + ra_sep) &
                (cat[dec_col] > dec_min - dec_sep) &
                (cat[dec_col] < dec_max + dec_sep)
            ]

        return cat

"""Generate postamps for Merian galaxies using HSC or DECam images.

Based on Johnny Greco's lsstutils package: https://github.com/johnnygreco/lsstutils
But upgraded to the Gen3 middleware

- 2022-02-17
"""

import os
import shutil
import argparse

from datetime import date

import numpy as np

import astropy.units as u
from astropy.table import Table, QTable

from joblib import Parallel, delayed
from spherical_geometry.polygon import SphericalPolygon

import lsst.log
Log = lsst.log.Log()
Log.setLevel(lsst.log.ERROR)

import lsst.daf.butler as dafButler

import lsst.daf.base
import lsst.geom as geom
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord

import lsst.afw.display
import lsst.afw.display.rgb as afwRgb

MERIAN_REPO = '/projects/MERIAN/repo'
S18A_WIDE_ROOT = '/tigress/HSC/DR/s18a_wide'
PIXEL_SCALE = 0.168 # arcsec / pixel

def _prepare_dataset(root, collection):
    """
    Get the butler for the given dataset and return the skyMap object.
    """
    butler = dafButler.Butler(root, collection=collection)
    return butler, butler.get('skyMap')

def _get_ra_dec_name(id_arr, ra_arr, dec_arr):
    """Get the object name based on ID and (RA, Dec)."""
    return [
        "{:s}_{:s}_{:s}_{:s}".format(
            str(i), "{:8.4f}".format(ra).strip(), "{:8.4f}".format(dec).strip()
            ) for (i, ra, dec) in zip(id_arr, ra_arr, dec_arr)]

def _get_file_prefix(name_arr, band, prefix):
    """Get the prefix of the output files based on the ID."""
    if prefix is None:
        return ["{:s}_{:s}".format(str(name), band) for name in name_arr]
    else:
        return ["{:s}_{:s}_{:s}".format(prefix, str(name), band) for name in name_arr]

def _get_output_dir(output_dir, chunk_arr, name_arr):
    """Get the directory for the output cutout data."""
    # Check the output directory
    if not os.path.isdir(output_dir):
        raise ValueError("Output directory '{:s}' does not exist".format(output_dir))
    
    return [os.path.join(output_dir, str(chunk), str(name)) for (chunk, name) in zip(chunk_arr, name_arr)]
    
def _get_int_chunk(data, n_chunk):
    """Assign integer chunk ID to the data."""
    if n_chunk > len(data):
        raise ValueError("Too many chunks...")
    if n_chunk <= 0:
        raise ValueError("Chunk number has to be larger than 0...")
    
    chunk_arr = np.ones(len(data), dtype=int)
    if n_chunk == 1:
        return chunk_arr
    
    chunk_size = np.ceil(len(data) / n_chunk).astype(int)

    start, end = 0, chunk_size
    for i in np.arange(n_chunk):
        chunk_arr[start: end] = i + 1
        start, end = end, end + chunk_size 
        end = len(data) if end > len(data) else end
    
    return chunk_arr

def _prepare_input_cat(input_cat, half_size, unit, ra_col, dec_col, band, id_col, chunk, 
                       prefix, output_dir, save=True):
    """
    Prepare the input sample for the given dataset.
    
    The cutouts are organized into:
        [output_dir]/[chunk_id]/[galaxy_id]/[file_name].fits
    And the file name prefix is: 
        ([prefix]_[galaxy_id]_[band]
    """
    # Load the input catalog
    if isinstance(input_cat, str):
        input_cat = Table.read(input_cat)

    # Get an array for half size
    if isinstance(half_size, str):
        if half_size.strip() not in input_cat.colnames:
            raise ValueError("Wrong half size column name. [{:s}]".format(half_size))
        half_size_arr = input_cat[half_size]
    else:
        # Using the same size for all objects
        half_size_arr = np.full(len(input_cat), float(half_size))
    
    if np.any(half_size_arr < 0):
        raise ValueError("Negative size value.")
    
    # Add size unit if necessary
    if unit != 'pixel' and half_size_arr.unit is None:
        # Check the half size unit
        if unit.strip() not in ['arcsec', 'arcmin', 'degree', 'pixel']:
            raise ValueError("Wrong size unit. [arcsec, arcmin, degree, pixel]")
        half_size_arr = [s * u.Unit(unit) for s in half_size_arr]

    # Get the RA and DEC arrays
    if ra_col not in input_cat.colnames:
        raise ValueError("Wrong R.A. column name. [{:s}]".format(ra_col))
    if dec_col not in input_cat.colnames:
        raise ValueError("Wrong Dec column name. [{:s}]".format(dec_col))
    ra_arr, dec_arr = input_cat[ra_col], input_cat[dec_col]
    
    # Get the output directory and file name 
    
    # Get the object id or name
    if id_col is None:
        name_arr = _get_ra_dec_name(np.arange(len(ra_arr)) + 1, ra_arr, dec_arr)
    else:
        if id_col not in input_cat.colnames:
            raise ValueError("Wrong ID column name. [{:s}]".format(id_col))
        name_arr = input_cat[id_col]
    
    # Get the output file prefix 
    prefix_arr = _get_file_prefix(name_arr, band, prefix)
    
    # Get the directory of the output file
    if chunk is not None:
        if isinstance(chunk, str):
            if chunk not in input_cat.colnames:
                raise ValueError("Wrong Chunk column name. [{:s}]".format(chunk))
            chunk_arr = input_cat[chunk]
        else:
            chunk_arr = _get_int_chunk(input_cat, int(chunk))
    else:
        chunk_arr = None
        
    # Get the output file directory
    dir_arr = _get_output_dir(output_dir, chunk_arr, name_arr)

    sample = QTable(
        [name_arr, prefix_arr, dir_arr, chunk_arr, list(ra_arr), list(dec_arr), half_size_arr],
        names=('name', 'prefix', 'dir', 'chunk', 'ra', 'dec', 'half_size')
    )

    if save:
        today = date.today()
        prefix = 'postamps' if prefix is None else prefix
        sample.write(
            os.path.join(output_dir, "{:s}-{:4d}-{:02d}-{:02d}.fits".format(
                prefix, today.year, today.month, today.day)), overwrite=True)
    
    return sample


def catalog_cutout(input_cat, root, collection, band, half_size=10, unit='arcsec', psf=True, data_type='deepCoadd_calexp',
                   output_dir='./', prefix=None, ra_col='ra', dec_col='dec', id_col='id', chunk=None):
    """
    Generate cutout/postamps for objects in a catalog using coadd images reduced by `lsstPipe`.

    Parameters
    ----------
    input_cat : `astropy.table.Table` or str
        Input catalog. Should contain at least the (RA, Dec) coordinates of the objects. 
    root : str
        Location of the data repo. 
    collection : str
        The collection sub-directory to the dataset.
    band : str
        Name of the filter.
    half_size : float or str
        Half size of the cutout image. 
        Could be the column name of the half size in the input catalog, or the half size value for the whole catalog.
        Used in combination with the unit.
    unit : str
        Unit of the half size.
    psf : bool
        Whether to generate PSF images.
    data_type : str
        The type of data to be used for the cutout.
    output_dir : str
        The output directory.
    prefix : str
        Prefix of the output file.
    ra_col : str
        Name of the RA column in the input catalog.
    dec_col : str
        Name of the Dec column in the input catalog.
    id_col : str
        Name of the ID column in the input catalog.
    chunk : int
        Number of chunks to split the catalog into.
    """

    # Check the required data type
    if data_type.strip() not in ['deepCoadd', 'deepCoadd_calexp', 'deepCoadd_background', 'deepCoadd_calexp_background']:
        raise ValueError("Wrong coadd type. [deepCoadd, deepCoadd_calexp, deepCoadd_background, deepCoadd_calexp_background]")

    # Get the butler for the given dataset and return the skyMap object.
    butler, skyMap = _prepare_dataset(root, collection)

    # Prepare the input catalog
    sample = _prepare_input_cat(
        input_cat, half_size, unit, ra_col, dec_col, band, id_col, chunk, prefix, output_dir)



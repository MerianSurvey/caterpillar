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
import lsst.daf.persistence as dafPersist
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
    butler = dafPersist.Butler(root, collection=collection)
    return butler, butler.get('skyMap')

def _prepare_input_cat(input_cat, half_size, unit, ra_col, dec_col, id_col, prefix, output_dir):
    """
    Prepare the input sample for the given dataset.
    """
    # Load the input catalog
    if isinstance(input_cat, str):
        input_cat = Table.read(input_cat)

    # Check the half size unit
    if unit.strip() not in ['arcsec', 'arcmin', 'degree', 'pixel']:
        raise ValueError("Wrong size unit. [arcsec, arcmin, degree, pixel]")


    # Check the column names
    if ra_col not in input_cat.colnames:
        raise ValueError("Column '{:s}' not found in the input catalog".format(ra_col))
    if dec_col not in input_cat.colnames:
        raise ValueError("Column '{:s}' not found in the input catalog".format(dec_col))
    if id_col not in input_cat.colnames:
        raise ValueError("Column '{:s}' not found in the input catalog".format(id_col))

    # Check the unit
    if unit is not None:
        if unit not in input_cat.colnames:
            raise ValueError("Column '{:s}' not found in the input catalog".format(unit))
    # Check the half size
    if half_size is not None:
        if half_size < 0:
            raise ValueError("The half size must be positive")
    # Check the output directory
    if not os.path.isdir(output_dir):
        raise ValueError("Output directory '{:s}' does not exist".format(output_dir))
    # Check the prefix
    if prefix is not None:
        if prefix == '':
            raise ValueError("The prefix cannot be empty")
    # Return the prepared catalog
    return input_cat, half_size, unit, ra_col, dec_col, id_col, prefix, output_dir

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
    sample = _prepare_input_cat(input_cat, half_size, unit, ra_col, dec_col, id_col, prefix, output_dir)


    if isinstance(half_size, str):
        half_size = input_cat[half_size]

    if unit == 'arcsec':
        half_size = half_size * u.arcsec
    elif unit == 'pixel':
        half_size = half_size * u.pixel
    else:
        raise ValueError('Unknown unit: {}'.format(unit))

    if prefix is None:
        prefix = '{}_{}'.format(band, half_size)

    if chunk is None:
        chunk = 1
    chunk_size = int(len(input_cat) / chunk)

    for i in range(chunk):
        chunk_cat = input_cat[i*chunk_size:(i+1)*chunk_size]
        chunk_cat = chunk_cat[np.isfinite(chunk_cat[ra_col]) & np.isfinite(chunk_cat[dec_col])]
        chunk_cat = chunk_cat[np.abs(chunk_cat[ra_col]) < 360]
        chunk_cat = chunk_cat[np.abs(chunk_cat[dec_col]) < 180]

        if len(chunk_cat) == 0:
            continue

        # Get the dataIds
        dataIds = []
        for row in chunk_cat:
            dataId = dict(instrument='HSC', detector=band, 
    """
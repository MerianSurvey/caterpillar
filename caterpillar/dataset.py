"""Module for the `caterpillar.dataset` class
"""

import os
import copy
import pathlib

import numpy as np

import fitsio

from astropy.table import Table

__author__ = ['Song Huang']

__all__ = ['Dataset']

# Metadata from the header of the primary HDU
DECAM_PRIM_KEYS = [
    'DATE', 'OBJECT', 'OBSTYPE', 'PROCTYPE', 'PRODTYPE', 'PIXSCAL1', 'PIXSCAL2',
    'EXPTIME', 'OBSID', 'DATE-OBS', 'TIME-OBS', 'MJD-OBS', 'MJD-END', 'PROPID',
    'FILTER', 'CENTRA', 'CENTDEC', 'CORN1RA', 'CORN1DEC', 'CORN2RA', 'CORN2DEC',
    'CORN3RA', 'CORN3DEC', 'CORN4RA', 'CORN4DEC', 'HA', 'ZD', 'AZ', 'DIMMSEE',
    'AIRMASS', 'HUMIDITY', 'DTNSANAM', 'RADESYS', 'RADECSYS', 'EQUINOX', 'PHOTFLAG',
    'SCAMPFLG', 'PHOTREF', 'MAGZERO', 'MAGZPT', 'NPHTMTCH', 'SKYSUB', 'PLVER', 'WCSCAL'
]

DECAM_CCD_KEYS = [
    'NAXIS1', 'NAXIS2', 'EXTNAME', 'DETPOS', 'CCDNUM', 'GAINA', 'GAINB',
    'RDNOISEA', 'RDNOISEB', 'SATURATA', 'SATURATB', 'CRPIX1', 'CRPIX2',
    'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'NSATPIX', 'NBLEED', 'FWHM',
    'ELLIPTIC', 'AVSKY', 'AVSIG', 'ARAWGAIN', 'CENRA1', 'CENDEC1',
    'COR1RA1', 'COR1DEC1', 'COR2RA1', 'COR2DEC1', 'COR3RA1', 'COR3DEC1',
    'COR4RA1', 'COR4DEC1'
]


class Dataset(object):
    '''A base class containing common code for handling dataset from DECam.

    Parameters
    ----------

    data_dir : string
        Path to the directory that contains the dataset.
    suffix : string, optional
        Suffix of the file to search for. Default: 'fits.fz'

    Notes
    -----

    '''
    def __init__(self, data_dir, suffix='fits.fz', verbose=True):
        '''Initialize a dataset classs.'''

        assert os.path.isdir(data_dir), "Wrong path to dataset directory."
        self.path = pathlib.Path(data_dir)
        self.suffix = suffix

        # Recursively find all the data files
        file_list = []
        for path in self.path.rglob('*.{:s}'.format(suffix)):
            file_list.append(str(path.absolute()))
        self.file_list = file_list
        self.n_file = len(file_list)
        if self.n_file > 1000:
            print("# Contains {:d} files! Metatable might take a while".format(
                self.n_file))

        # Get the unique exposures
        # For each exposure there should be three different types of products:
        # `i`: image; `d`: data quality mask; `w`: invariance variance weight map.
        naming_components = [
            os.path.basename(f).replace(self.suffix, '').split('_') for f in self.file_list
        ]

        # Exposure ID: telescope/instrument + YYMMDD + HHMMSS
        exposure_ids = [
            ('_').join(name[:3]) for name in naming_components
        ]
        self.exposures = list(set(exposure_ids))
        self.n_exposure = len(self.exposures)

        # Data products ID: exposure_id + [object_type processing_type product_type]
        # For all our data object_type should be `o`, means "object"
        # Processing type can be: `o`: `InstCal`; `c`: `MasterCal`, `p`: projected,
        #                         `s`: `Stacked`; `k`: `SkySub`, `u`: non of the above
        processing_ids = [
            ('_').join(name[:3] + [name[3][:2]]) for name in naming_components
        ]
        self.processes = list(set(processing_ids))
        self.n_process = len(self.processes)

    def _form_meta_table(self, keys=DECAM_PRIM_KEYS):
        """Form a table to summarize the key information in the header.

        Parameters
        ----------
        keys: list, optional
            List of keywords in the header to recover. Default: DECAM_PRIM_KEYS
        """
        return Table(
            [{k.lower(): fitsio.read_header(f)[k] for k in keys} for f in self.file_list])

    def get_products(self, data_id):
        """Gather the data products belongs to a certain exposure or process.

        Parameters
        ----------

        data_id : string
            Unique exposure or data processing ID.
        """
        prod_list = []
        for path in self.path.rglob('{:s}*.{:s}'.format(data_id, self.suffix)):
            prod_list.append(str(path.absolute()))
        return prod_list

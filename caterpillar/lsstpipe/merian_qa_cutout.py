
import os
import shutil
import argparse
import warnings

from datetime import date

import numpy as np

import astropy.units as u
from astropy.table import Table, QTable

from joblib import Parallel, delayed
from spherical_geometry.polygon import SphericalPolygon

try:
    import lsst.log
    Log = lsst.log.Log()
    Log.setLevel(lsst.log.ERROR)

    import lsst.daf.butler as dafButler
    import lsst.geom as geom
    import lsst.afw.image as afwImage
    import lsst.afw.geom as afwGeom
except ImportError:
    warnings.warn("lsstPipe is not installed. Please install it first.")

MERIAN_REPO = '/projects/MERIAN/repo'
S18A_WIDE_ROOT = '/tigress/HSC/DR/s18a_wide'
PIXEL_SCALE = 0.168 # arcsec / pixel

from cutout import *


root = '/projects/MERIAN/repo'
output = '/projects/MERIAN/poststamps/merian_qa/deep_2202'

n708 = 'DECam/runs/merian/w_2022_02/t9813_deep_N708'
n540 = 'DECam/runs/merian/w_2022_02/t9813_deep_N540'

input_cat = os.path.join(output, 'merian_qa_20220419.fits')

data_type='deepCoadd_calexp'

sample_n708 = cutout_batch(
    input_cat, root, n708, 'N708', njobs=1, psf=False, ready=False, save=True, half_size='half_size',
    unit='arcsec', data_type=data_type, output_dir=output, prefix='merian_qa', chunk='chunk', verbose=True)

sample_n708.write(os.path.join(output, 'merian_qa_deep_2202_n708_done.fits'))

sample_n540 = cutout_batch(
    input_cat, root, n540, 'n540', njobs=1, psf=True, ready=True, save=False, half_size='half_size',
    unit='arcsec', data_type=data_type, output_dir=output, prefix='merian_qa', chunk='chunk', verbose=True)

sample_n540.write(os.path.join(output, 'merian_qa_deep_2202_n540_done.fits'))


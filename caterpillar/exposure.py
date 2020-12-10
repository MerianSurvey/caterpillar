"""Module for the `caterpillar.exposure` class
"""

import os
import copy
import pathlib

import numpy as np

import fitsio

from astropy.table import Table

__author__ = ['Song Huang']

__all__ = ['Exposure']


class Exposure(object):
    '''A base class containing common code for handling CP reduced exposure.

    Here, an "exposure" means a collection of the three data products associated
    with one exposure

    Parameters
    ----------

    Notes
    -----

    '''
    pass

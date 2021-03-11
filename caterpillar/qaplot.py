#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""QA plotting tools."""

import os
import re

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clip

import matplotlib.pyplot as plt
from matplotlib import rc

from . import visual

rc('text', usetex=False)

__all__ = ['header_section_to_index', 'get_ccd_section', 'visual_ccd',
           'visual_exposure_ccds']

def header_section_to_index(section):
    """Convert a header section keyword to four indices."""
    return [int(num) for num in re.split(':|,', section.replace('[', '').replace(']', ''))]


def get_ccd_section(img, hdr, keyword, T=True):
    """Retrieve a section of the CCD image."""
    idx = header_section_to_index(hdr[keyword])
    if T:
        return img[idx[0]:idx[1], idx[2]:idx[3]]
    return img[idx[2]:idx[3], idx[0]:idx[1]]

def visual_ccd(ccd, exp_hdr, fig_dir=None, dpi=100, contrast=0.9,
               cmap='coolwarm', **kwargs):
    """Visualize a single CCD."""

    # Data and header information from CCD
    img = ccd.data.T
    hdr = ccd.header
    ccd_id = hdr['EXTNAME'].strip() + '_' + str(hdr['CCDNUM'])

    # Skip N30 CCD
    if hdr['EXTNAME'].strip() == 'N30':
        return

    # Skip focus CCDs
    if hdr['EXTNAME'].strip()[0] == 'F':
        return

    # Header information from the exposure
    exp_name = exp_hdr['DTNSANAM'].split('.')[0]
    exp_time = int(exp_hdr['EXPTIME'])

    # Get the data and overscan regions of amplifier A & B
    data_a = get_ccd_section(img, hdr, 'DATASECA')
    data_b = get_ccd_section(img, hdr, 'DATASECB')

    bias_a = get_ccd_section(img, hdr, 'BIASSECA')
    bias_b = get_ccd_section(img, hdr, 'BIASSECA')

    gain_a, rdn_a = hdr['GAINA'], hdr['RDNOISEA']
    gain_b, rdn_b = hdr['GAINB'], hdr['RDNOISEB']

    data_a_cor = (data_a * gain_a) - rdn_a
    data_b_cor = (data_b * gain_b) - rdn_b

    bias_a_cor = (bias_a * gain_a) - rdn_a
    bias_b_cor = (bias_b * gain_b) - rdn_b

    # A very naive overscane "correction"
    data_a_cor = data_a_cor - np.nanmedian(bias_a_cor)
    data_b_cor = data_b_cor - np.nanmedian(bias_b_cor)

    mean_a = np.nanmean(
        sigma_clip(data_a_cor.flatten(), sigma=2.5, maxiters=3))
    mean_b = np.nanmean(
        sigma_clip(data_b_cor.flatten(), sigma=2.5, maxiters=3))

    data_b_cor += (mean_a - mean_b)

    if (hdr['EXTNAME'][0] == 'S') or (hdr['EXTNAME'][0] == 'F'):
        img_new = np.vstack([data_b_cor, data_a_cor])
    else:
        img_new = np.vstack([data_a_cor, data_b_cor])

    x_size, y_size = img_new.shape
    x_fig, y_fig = int(x_size / dpi), int(y_size / dpi) + 0.2

    # Make the plot
    fig = plt.figure(figsize=(y_fig, x_fig))
    fig.subplots_adjust(
        left=0.005, bottom=0.005, right=0.995, top=0.96)

    ax = fig.add_subplot(1, 1, 1)
    ax = visual.display_single(
        img_new, contrast=contrast,
        scale_bar=False, cmap=cmap, ax=ax, **kwargs)

    _ = ax.set_title("{:s}  {:s}  EXPTIME: {:d}s".format(
        exp_name, ccd_id, exp_time), fontsize=40)

    fig_name = "{:s}_{:s}.jpg".format(exp_name, ccd_id)

    if fig_dir is not None:
        fig_name = os.path.join(fig_dir, fig_name)

    fig.savefig(fig_name, dpi=dpi)

    plt.close(fig)

def visual_exposure_ccds(exp_file, fig_dir=None, dpi=100, contrast=0.9,
                         cmap='coolwarm', **kwargs):
    """Visual check of an exposure."""
    exp = fits.open(exp_file)
    exp_hdr = exp[0].header

    for ccd_id in np.arange(len(exp))[1:]:
        visual_ccd(exp[ccd_id], exp_hdr, fig_dir=fig_dir,
                   dpi=dpi, contrast=contrast, cmap=cmap, **kwargs)

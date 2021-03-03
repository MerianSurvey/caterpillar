"""Functions about HSC sample selection."""

import numpy as np

from astropy.table import Column

def hsc_broad_cut(catalog, aper=True, size=True, aper_frac=True, aper_color=True,
                  risky=False, verbose=True):
    """Apply broad-cut sample selection to HSC catalog.

    - 2021-02: Based on Song Huang's investigation of COSMOS dwarf sample.
    """

    if aper:
        aper_mask = (
            (catalog['g_undeblended_convolvedflux_2_15_flux'] > 0.) &
            (catalog['r_undeblended_convolvedflux_2_15_flux'] > 0.) &
            (catalog['i_undeblended_convolvedflux_2_15_flux'] > 0.) &
            (catalog['z_undeblended_convolvedflux_2_15_flux'] > 0.) &
            (catalog['y_undeblended_convolvedflux_2_15_flux'] > 0.) &
            (catalog['g_undeblended_convolvedflux_3_20_flux'] > 0.) &
            (catalog['r_undeblended_convolvedflux_3_20_flux'] > 0.) &
            (catalog['i_undeblended_convolvedflux_3_20_flux'] > 0.) &
            (catalog['z_undeblended_convolvedflux_3_20_flux'] > 0.) &
            (catalog['y_undeblended_convolvedflux_3_20_flux'] > 0.)
        )
    else:
        aper_mask = catalog['i_cmodel_flux'] > 0.

    if size:
        size_mask = (
            (np.log10(catalog['i_cmodel_exp_ellipse_r']) >= -0.8) &
            (np.log10(catalog['r_cmodel_exp_ellipse_r']) >= -0.8) &
            (np.log10(catalog['i_cmodel_exp_ellipse_r']) <= 0.75) &
            (np.log10(
                catalog['i_cmodel_exp_ellipse_r']) < -0.25 * catalog['i_cmodel_mag'] + 6.2)
        )
    else:
        size_mask = catalog['i_cmodel_flux'] > 0.

    if aper_frac:
        frac_mask = (
            (catalog['i_undeblended_convolvedflux_2_15_flux'] / catalog['i_cmodel_flux'] <= 1.0) &
            (catalog['i_psf_flux'] / catalog['i_cmodel_flux'] <= 0.90) &
            (catalog['i_undeblended_convolvedflux_2_15_flux'] / catalog['i_cmodel_flux'] >=
             0.13 * catalog['i_cmodel_mag'] - 2.65)
        )
    else:
        frac_mask = catalog['i_cmodel_flux'] > 0.

    rz_aper_cos = (-2.5 * np.log10(
        catalog['r_undeblended_convolvedflux_3_20_flux'] /
        catalog['z_undeblended_convolvedflux_3_20_flux']))
    gi_aper_cos = (-2.5 * np.log10(
        catalog['g_undeblended_convolvedflux_3_20_flux'] /
        catalog['i_undeblended_convolvedflux_3_20_flux']))
    gr_aper_cos = (-2.5 * np.log10(
        catalog['g_undeblended_convolvedflux_3_20_flux'] /
        catalog['r_undeblended_convolvedflux_3_20_flux']))
    iy_aper_cos = (-2.5 * np.log10(
        catalog['i_undeblended_convolvedflux_3_20_flux'] /
        catalog['y_undeblended_convolvedflux_3_20_flux']))

    if aper_color:
        aper_color_mask = (
            (rz_aper_cos <= 0.75) & (rz_aper_cos >= 0.0) &
            (rz_aper_cos >= 0.13 * catalog['i_cmodel_mag'] - 2.8) &
            (gi_aper_cos <= 1.5) & (gi_aper_cos >= 0.0) &
            (gi_aper_cos >= 0.25 * catalog['i_cmodel_mag'] - 5.4) &
            (rz_aper_cos <= 0.63 * gr_aper_cos + 0.25) &
            (iy_aper_cos <= 0.35 * gi_aper_cos + 0.12) &
            (iy_aper_cos >= -0.15)
        )
    else:
        aper_color_mask = catalog['i_cmodel_flux'] > 0.

    mask = size_mask & aper_mask & frac_mask & aper_color_mask

    if risky:
        risky_mask = (
            (catalog['i_cmodel_mag'] <= 23.0) &
            (gi_aper_cos <= 1.20) &
            (rz_aper_cos <= 0.65) &
            (gr_aper_cos <= 0.85)
        )
        mask = mask & risky_mask

    if verbose:
        print("# Select {:d}/{:d} objects".format(mask.sum(), len(catalog)))

    return catalog[mask]

def add_dynamic_poststamp_size(catalog, size_col='i_cmodel_exp_ellipse_r',
                               size_factor=10., min_size=3, max_size=20):
    """Add dynamical poststamp size."""
    s_ang = size_factor * catalog[size_col]
    s_ang = np.where(s_ang <= max_size, s_ang, max_size)
    s_ang = np.where(s_ang >= min_size, s_ang, min_size)

    catalog.add_column(Column(data=s_ang, name='half_size'))

    return catalog

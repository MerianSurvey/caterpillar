#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Galactic extinction related."""

import warnings

import numpy as np

from scipy.interpolate import CubicSpline

from sedpy import observate

__all__ = ['fitzpatrick99', 'get_extinction_coefficient']


def get_extinction_coefficient(band, reference, A_V=1.0, R_V=3.1, filter_dir=None,
                               N_factor=0.8855, **kwargs):
    """
    Estimate the extinction coefficient for a given filter.
    """
    # Filter transmission curve
    filter_obj = observate.Filter(band, directory=filter_dir)
    filter_wave = filter_obj.wavelength

    # Reference spectrum
    ref_wave, ref_flux = reference['wave'], reference['flux']
    if ref_wave[0] > filter_wave[0] or ref_wave[-1] < filter_wave[-1]:
        warnings.warn("# Reference spectrum does not cover the filter!")

    # AB magnitude of the
    ref_mag = filter_obj.ab_mag(ref_wave, ref_flux)

    # E(B-V) value
    e_bv = A_V / R_V

    # Extinction magnitude using Fitzpatrick 1999 extinction curve
    ext_mag = fitzpatrick99(ref_wave, e_bv, R_V, **kwargs)
    #ext_1micron = fitzpatrick99(1e4, e_bv, R_V, **kwargs)

    ref_mag_ext = filter_obj.ab_mag(
        ref_wave, ref_flux * 10.0 ** (-0.4 * (ext_mag * N_factor)))

    return (ref_mag_ext - ref_mag) / e_bv


def fitzpatrick99(wave, ebv, R_V=3.1, flux=None, lmc2_set=False,
                  avglmc_set=False, gamma=0.99, x0=4.596,
                  c1=None, c2=None, c3=3.23, c4=0.41):
    '''
    NAME:
     FM_UNRED
    PURPOSE:
     Deredden a flux vector using the Fitzpatrick (1999) parameterization
    EXPLANATION:
     The R-dependent Galactic extinction curve is that of Fitzpatrick & Massa
     (Fitzpatrick, 1999, PASP, 111, 63; astro-ph/9809387 ).
     Parameterization is valid from the IR to the far-UV (3.5 microns to 0.1
     microns).  UV extinction curve is extrapolated down to 912 Angstroms.
    CALLING SEQUENCE:
     fm_unred( wave, flux, ebv [, 'LMC2', 'AVGLMC', 'ExtCurve', R_V = ,
                                   gamma =, x0=, c1=, c2=, c3=, c4= ])
    INPUT:
      wave - wavelength vector (Angstroms)
      flux - calibrated flux vector, same number of elements as "wave"
      ebv  - color excess E(B-V), scalar.  If a negative "ebv" is supplied,
              then fluxes will be reddened rather than dereddened.
    OUTPUT:
      Unreddened flux vector, same units and number of elements as "flux"
    OPTIONAL INPUT KEYWORDS
      R_V - scalar specifying the ratio of total to selective extinction
               R(V) = A(V) / E(B - V).  If not specified, then R = 3.1
               Extreme values of R(V) range from 2.3 to 5.3
      'AVGLMC' - if set, then the default fit parameters c1,c2,c3,c4,gamma,x0
             are set to the average values determined for reddening in the
             general Large Magellanic Cloud (LMC) field by Misselt et al.
            (1999, ApJ, 515, 128)
      'LMC2' - if set, then the fit parameters are set to the values determined
             for the LMC2 field (including 30 Dor) by Misselt et al.
             Note that neither /AVGLMC or /LMC2 will alter the default value
             of R_V which is poorly known for the LMC.
      The following five input keyword parameters allow the user to customize
      the adopted extinction curve.  For example, see Clayton et al. (2003,
      ApJ, 588, 871) for examples of these parameters in different interstellar
      environments.
      x0 - Centroid of 2200 A bump in microns (default = 4.596)
      gamma - Width of 2200 A bump in microns (default = 0.99)
      c3 - Strength of the 2200 A bump (default = 3.23)
      c4 - FUV curvature (default = 0.41)
      c2 - Slope of the linear UV extinction component
           (default = -0.824 + 4.717 / R)
      c1 - Intercept of the linear UV extinction component
           (default = 2.030 - 3.007 * c2)
    OPTIONAL OUTPUT KEYWORD:
      'ExtCurve' - If this keyword is set, fm_unred will return two arrays.
                  First array is the unreddend flux vector.  Second array is
                  the E(wave-V)/E(B-V) extinction curve, interpolated onto the
                  input wavelength vector.
    EXAMPLE:
       Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
       is altered by a reddening of E(B-V) = 0.1.  Assume an "average"
       reddening for the diffuse interstellar medium (R(V) = 3.1)
       >>> w = 1200 + arange(40)*50       #Create a wavelength vector
       >>> f = w*0 + 1                    #Create a "flat" flux vector
       >>> fnew = fm_unred(w, f, -0.1)    #Redden (negative E(B-V)) flux vector
       >>> plot(w, fnew)
    NOTES:
       (1) The following comparisons between the FM curve and that of Cardelli,
           Clayton, & Mathis (1989), (see ccm_unred.pro):
           (a) - In the UV, the FM and CCM curves are similar for R < 4.0, but
                 diverge for larger R
           (b) - In the optical region, the FM more closely matches the
                 monochromatic extinction, especially near the R band.
       (2)  Many sightlines with peculiar ultraviolet interstellar extinction
               can be represented with the FM curve, if the proper value of
               R(V) is supplied.
    REQUIRED MODULES:
       scipy, np
    REVISION HISTORY:
       Written   W. Landsman        Raytheon  STX   October, 1998
       Based on FMRCurve by E. Fitzpatrick (Villanova)
       Added /LMC2 and /AVGLMC keywords,  W. Landsman   August 2000
       Added ExtCurve keyword, J. Wm. Parker   August 2000
       Assume since V5.4 use COMPLEMENT to WHERE  W. Landsman April 2006
       Ported to Python, C. Theissen August 2012
    '''
    x = 10000. / np.array([wave])  # Convert to inverse microns
    curve = x * 0.

    if lmc2_set:
        x0, gamma = 4.626, 1.05
        c1, c2, c3, c4 = -2.16, 1.31, 1.92, 0.42
    elif avglmc_set:
        x0, gamma = 4.596, 0.91
        c1, c2, c3, c4 = -1.28, 1.11, 2.73, 0.64
    else:
        if c2 is None:
            c2 = -0.824 + 4.717 / R_V
        if c1 is None:
            c1 = 2.030 - 3.007 * c2

    # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and
    # R-dependent coefficients
    xcutuv = 10000.0 / 2700.0
    xspluv = 10000.0 / np.array([2700.0, 2600.0])

    iuv = x >= xcutuv
    iuv_comp = ~iuv

    if len(x[iuv]) > 0:
        xuv = np.concatenate((xspluv, x[iuv]))
    else:
        xuv = xspluv.copy()

    yuv = c1  + c2 * xuv
    yuv = yuv + c3 * xuv**2 / ((xuv ** 2 - x0 ** 2) ** 2 + (xuv * gamma) ** 2)

    filter1 = xuv.copy()
    filter1[xuv <= 5.9] = 5.9

    yuv = yuv + c4 * (0.5392 * (filter1 - 5.9) ** 2 + 0.05644 * (filter1 - 5.9) ** 3)
    yuv = yuv + R_V
    # Save spline points
    yspluv = yuv[0:2].copy()

    # Remove spline points
    if len(x[iuv]) > 0:
        curve[iuv] = yuv[2:len(yuv)]

    # Compute optical portion of A(lambda)/E(B-V) curve
    # using cubic spline anchored in UV, optical, and IR
    xsplopir = np.concatenate(
        ([0], 10000.0 / np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])))
    ysplir = np.array([0.0, 0.26469, 0.82925]) * R_V / 3.1
    ysplop = [
        np.polyval(np.array([2.13572e-04, 1.00270, -4.22809e-01]), R_V),
        np.polyval(np.array([-7.35778e-05, 1.00216, -5.13540e-02]), R_V),
        np.polyval(np.array([-3.32598e-05, 1.00184, 7.00127e-01]), R_V),
        np.polyval(np.array(
            [-4.45636e-05, 7.97809e-04, -5.46959e-03, 1.01707, 1.19456]), R_V)]

    ysplopir = np.concatenate((ysplir, ysplop))

    if len(iuv_comp) > 0:
        cubic = CubicSpline(np.concatenate((xsplopir, xspluv)),
                            np.concatenate((ysplopir, yspluv)), bc_type='natural')
        curve[iuv_comp] = cubic(x[iuv_comp])

    curve = ebv * curve[0]

    # Now apply extinction correction to input flux vector
    if flux is not None:
        flux = flux * 10. ** (0.4 * curve)
        return curve, flux
    return curve

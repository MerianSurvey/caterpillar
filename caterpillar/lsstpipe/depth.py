"""Codes related to estimating the imaging depth."""

import numpy as np 

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.table import Table

import lsst.daf.butler as dafButler
from lsst.daf.base import DateTime

__all__ = ["estimate_photometric_limits"]

def estimate_photometric_limits(butler, visit, ccd, mag_low=25., mag_upp=18., poly_deg=1,
                                collections = 'u/lskelvin/testing',
                                flux_col="base_CircularApertureFlux_6_0_instFlux", 
                                flux_err_col="base_CircularApertureFlux_6_0_instFluxErr", 
                                flux_short='aper6',
                                verbose=False, visual=False, save_fig=False):
    """Estimate the photometric limits of a CCD.
    
    By default, the function will try to estimate the 3-, 5-, and 10-sigma 
    magnitude limits for a certain type of photometry based on:
    1. The relation between S/N and magnitude.
    2. The flux distribution of sky objects.
    
    TODO:
    - This needs to be modified when we have multiband data.
    """
    ccd_summary = {"visit": visit, "ccd": ccd}
    
    dataId = dict(instrument='DECam', visit=int(visit), detector=str(ccd))

    # Get the single exposure source catalog
    try:
        src = butler.get('src', collections=collections, dataId=dataId)
        src_cat = src.asAstropy()
    except:
        print("!!! Failed to load data for visit={:d} ccd={:s}".format(visit, ccd))
        return None 

    # Get the exposure for calibration
    # TODO: is there anyway to get the zeropoint
    calexp = butler.get('calexp', collections=collections, dataId=dataId)
    # Photometric calibration
    calib = calexp.getPhotoCalib()
    
    # Basic info of the exposure
    visit_info = calexp.getInfo().getVisitInfo()
    ccd_summary['time_utc'] = visit_info.getDate().toString(DateTime.Timescale.UTC)
    ccd_summary['alt_degree'] = visit_info.getBoresightAzAlt().getLongitude().asDegrees()
    ccd_summary['az_degree'] = visit_info.getBoresightAzAlt().getLatitude().asDegrees()
    ccd_summary['ha_degree'] = visit_info.getBoresightHourAngle().asDegrees()
    ccd_summary['rot_degree'] = visit_info.getBoresightRotAngle().asDegrees()
    ccd_summary['t_exposure'] = visit_info.getExposureTime()
    ccd_summary['airmass'] = visit_info.getBoresightAirmass()
    ccd_summary['temperature'] = visit_info.getWeather().getAirTemperature()
    ccd_summary['air_pressure'] = visit_info.getWeather().getAirPressure()
    ccd_summary['humidity'] = visit_info.getWeather().getHumidity()

    # Separate the sky objects
    sky_src = src_cat[src_cat['sky_source']]
    obj_src = src_cat[~src_cat['sky_source']]

    obj_mask = ~((obj_src[flux_col] > 0) & obj_src['slot_ApFlux_flag'] & obj_src['slot_PsfFlux_flag'] & 
        obj_src['slot_Centroid_flag'] & obj_src['slot_CalibFlux_flag'] & 
        obj_src['slot_PsfShape_flag'])

    obj_use = obj_src[obj_mask]
    
    if verbose:
        print("\n# Dealing with visit: {:8d} - CCD {:s}".format(visit, ccd))
        print("# There are {:d} detections and {:d} sky objects".format(len(obj_src), len(sky_src)))
        print('# There are {:d} objects with useful photometry'.format(obj_mask.sum()))
    ccd_summary['n_obj'] = len(obj_src)
    ccd_summary['n_obj_use'] = obj_mask.sum()
    ccd_summary['n_sky'] = len(sky_src)
        
    mag, s2n, mag_err = [], [], []
    for obj in obj_use:
        phot = calib.instFluxToMagnitude(obj[flux_col], obj[flux_err_col])
        s2n.append(obj[flux_col] / obj[flux_err_col])
        mag.append(phot.value)
        mag_err.append(phot.error)

    mag = np.asarray(mag)
    s2n = np.asarray(s2n)
    mag_err = np.asarray(mag_err)

    s2n_mask = (s2n > 0) & (mag <= mag_low) & (mag >= mag_upp)

    s2n_poly = np.poly1d(np.polyfit(np.log10(s2n[s2n_mask]), mag[s2n_mask], deg=poly_deg))
    s2n_grid = np.linspace(-1.0, 2.0, 20)

    mag_limits = s2n_poly([np.log10(3.), np.log10(5.), np.log10(10.)])
    ccd_summary['mag_3sig_obj'] = mag_limits[0]
    ccd_summary['mag_5sig_obj'] = mag_limits[1]
    ccd_summary['mag_10sig_obj'] = mag_limits[2]

    if verbose:
        print("# Based on S/N - magnitude relation of real detections:")
        print("   3-sigma limit: {:6.3f}".format(ccd_summary['mag_3sig_obj']))
        print("   5-sigma limit: {:6.3f}".format(ccd_summary['mag_5sig_obj']))
        print("  10-sigma limit: {:6.3f}".format(ccd_summary['mag_10sig_obj']))

    # Estimate the detection limit using sky object
    ccd_summary['mag_3sig_sky'] = calib.instFluxToMagnitude(np.nanstd(sky_src[flux_col]) * 3, 1.0).value
    ccd_summary['mag_5sig_sky'] = calib.instFluxToMagnitude(np.nanstd(sky_src[flux_col]) * 5, 1.0).value
    ccd_summary['mag_10sig_sky'] = calib.instFluxToMagnitude(np.nanstd(sky_src[flux_col]) * 10, 1.0).value

    if verbose:
        print("# Based on flux distribution of sky object:")
        print("   3-sigma limit: {:6.3f}".format(ccd_summary['mag_3sig_sky']))
        print("   5-sigma limit: {:6.3f}".format(ccd_summary['mag_5sig_sky']))
        print("  10-sigma limit: {:6.3f}".format(ccd_summary['mag_10sig_sky']))
    
    if visual:
        fig = plt.figure(figsize=(6.5, 6))
        fig.subplots_adjust(
            left=0.18, bottom=0.12, right=0.99, top=0.93, wspace=0, hspace=0)

        ax1 = fig.add_subplot(1, 1, 1)
        ax1.set_yscale("log", nonpositive='clip')
        ax1.grid(linestyle='--', alpha=0.5)

        ax1.scatter(mag, s2n, marker='o', s=30, facecolor='orangered', alpha=0.6,
                    edgecolor='w')

        ax1.plot(s2n_poly(s2n_grid), 10.0 ** s2n_grid, linewidth=3, c='k', alpha=0.9, linestyle='--')

        ax1.axhline(3, linestyle='--', linewidth=3, color='grey', alpha=0.5)
        ax1.text(16., 3.5, 'S/N=3', fontsize=15)

        ax1.axhline(5.0, linestyle='--', linewidth=3, color='orange', alpha=0.5)
        ax1.text(18., 5.5, 'S/N=5', fontsize=15)

        ax1.axhline(10.0, linestyle='--', linewidth=3, color='forestgreen', alpha=0.5)
        ax1.text(16., 11.0, 'S/N=10', fontsize=15)

        ax1.axvline(ccd_summary['mag_3sig_sky'], linewidth=3, linestyle='-.', color='grey', alpha=0.5)
        ax1.axvline(ccd_summary['mag_5sig_sky'], linewidth=3, linestyle='-.', color='orange', alpha=0.5)
        ax1.axvline(ccd_summary['mag_10sig_sky'], linewidth=3, linestyle='-.', color='forestgreen', alpha=0.5)

        ax1.set_xlabel('Magnitude: {:s}'.format(flux_short), fontsize=22)
        ax1.set_ylabel('log10(S/N)', fontsize=22)

        ax1.set_xlim(ax1.get_xlim()[::-1])

        ax1.tick_params(axis='both', which='major', labelsize=20)

        _ = ax1.set_title('Visit: {:d} - CCD: {:s}'.format(visit, ccd), fontsize=24, pad=8)

        ax1.text(0.05, 0.92, 'Sky', fontsize=18, transform=ax1.transAxes)
        ax1.text(0.04, 0.85, '3sig: {:5.2f}'.format(ccd_summary['mag_3sig_sky']), fontsize=14, transform=ax1.transAxes)
        ax1.text(0.04, 0.80, '5sig: {:5.2f}'.format(ccd_summary['mag_5sig_sky']), fontsize=14, transform=ax1.transAxes)
        ax1.text(0.025, 0.75, '10sig: {:5.2f}'.format(ccd_summary['mag_10sig_sky']), fontsize=14, transform=ax1.transAxes)

        ax1.text(0.75, 0.22, 'Obj', fontsize=18, transform=ax1.transAxes)
        ax1.text(0.74, 0.15, '3sig: {:5.2f}'.format(ccd_summary['mag_3sig_obj']), fontsize=14, transform=ax1.transAxes)
        ax1.text(0.74, 0.10, '5sig: {:5.2f}'.format(ccd_summary['mag_5sig_obj']), fontsize=14, transform=ax1.transAxes)
        ax1.text(0.725, 0.05, '10sig: {:5.2f}'.format(ccd_summary['mag_10sig_obj']), fontsize=14, transform=ax1.transAxes)
        
        if save_fig:
            fig.savefig("s2n_{:d}_{:s}_{:s}.png".format(visit, ccd, flux_short), dpi=100)
            plt.close(fig)
    
    return ccd_summary

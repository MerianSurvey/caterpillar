from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.table import Table, vstack, hstack, Column
from astropy.stats import SigmaClip, sigma_clipped_stats

from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

 

import os, sys
import glob
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.gridspec as gridspec

plt.style.use('mplstyle')

# Time New Roman font doesn't works in tiger, so I use STIXGeneral
plt.rcParams['font.family'] = 'STIXGeneral'


# GAIA query
def query_gaia_refcat(coord, width, height, pm_cut=10):
  
    v = Vizier(columns=["**", "+_r"], catalog='I/355/gaiadr3', row_limit=100000000)   
    result = v.query_region(coord, width=width, height=height, catalog='I/355/gaiadr3')
    
    if len(result) > 0: result = result[0]
    else: return None

    ### GAIA catalog
    gaia_tbl = Table()
    gaia_tbl['RA'] = result['RA_ICRS'].data.data
    gaia_tbl['DEC'] = result['DE_ICRS'].data.data
    gaia_tbl['RA_ERR'] = result['e_RA_ICRS'].data.data
    gaia_tbl['DEC_ERR'] = result['e_DE_ICRS'].data.data
    gaia_tbl['THETA_ERR'] = np.zeros_like(gaia_tbl['RA'].data.data)
    gaia_tbl['PRIMARY'] = (result['APF']==0).data ## set primary to be True
    gaia_tbl['POINT_SOURCE'] = np.ones_like(gaia_tbl['RA'].data.data, dtype='bool')  ## point source is True
    gaia_tbl['pmRA'] = result['pmRA'].data.data
    gaia_tbl['pmDEC'] = result['pmDE'].data.data
    gaia_tbl['pmRA_ERR'] = result['e_pmRA'].data.data
    gaia_tbl['pmDEC_ERR'] = result['e_pmDE'].data.data
    gaia_tbl['Epoch'] = np.ones_like(gaia_tbl['RA'].data.data) * 2016 # According to VIZIER, RA_ICRS are based on J2016
    gaia_tbl['SourceID'] = result['Source'].data.data
    gaia_tbl['RAJ2000'] = result['RAJ2000'].data.data
    gaia_tbl['DECJ2000'] = result['DEJ2000'].data.data
    gaia_tbl['G_MAG'] = result['Gmag'].data.data
    gaia_tbl['G_MAG_ERR'] = result['e_Gmag'].data.data
    gaia_tbl['BP_MAG'] = result['BPmag'].data.data
    gaia_tbl['BP_MAG_ERR'] = result['e_BPmag'].data.data
    gaia_tbl['RP_MAG'] = result['RPmag'].data.data
    gaia_tbl['RP_MAG_ERR'] = result['e_RPmag'].data.data


    # small proper motion cut
    flg = gaia_tbl['pmRA']   < pm_cut  
    flg &=   gaia_tbl['pmDEC']  / np.cos(np.deg2rad(gaia_tbl['DEC'] * u.deg))   < pm_cut 
    gaia_tbl = gaia_tbl[flg]

    return gaia_tbl



def astrometry_plotter(dra_gaia, ddec_gaia, dra_hsc, ddec_hsc,
            dra_median_gaia, ddec_median_gaia, dra_median_hsc, ddec_median_hsc):
    lim = 0.25
    binwidth = 0.01

    fig = plt.figure(figsize=(26,12))
    gs = gridspec.GridSpec(4,9)

    # LEFT PANEL
    ax_main = plt.subplot(gs[1:4, :3])
    ax_xDist = plt.subplot(gs[0, :3],sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:4, 3],sharey=ax_main)
    ax_xDist.tick_params(axis='x', labelbottom=False, labelleft=False)
    ax_yDist.tick_params(axis='y', labelbottom=False, labelleft=False)
    ax_xDist.tick_params(axis='y', labelbottom=False, labelleft=False)
    ax_yDist.tick_params(axis='x', labelbottom=False, labelleft=False)

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_xDist.hist(dra_gaia, bins=bins, color='red', histtype='step',fill =True, alpha=0.2, lw=1)
    ax_yDist.hist(ddec_gaia, bins=bins, color='red', orientation='horizontal',align='mid', histtype='step', fill=True, alpha=0.2, lw=1)
    ax_xDist.axvline(x=0, ymin=0, ymax=1, color='grey', lw=1, ls='--')
    ax_yDist.axhline(y=0, xmin=0, xmax=1, color='grey', lw=1, ls='--')
    ax_xDist.text(0.16, 0.8, r'$\langle\Delta$RA$\rangle$=%.3f'%(dra_median_gaia), horizontalalignment='center', verticalalignment='center', transform=ax_xDist.transAxes, fontsize=25)
    ax_yDist.text(0.5, 0.93, r'$\langle\Delta$Dec$\rangle$=%.3f'%(ddec_median_gaia), horizontalalignment='center', verticalalignment='center', transform=ax_yDist.transAxes, fontsize=25)

    ax_main.scatter(dra_gaia, ddec_gaia, color='r', alpha=0.2)
    ax_main.scatter(dra_median_gaia, ddec_median_gaia, marker='P', s=300, color='k')
    ax_main.axhline(y=0, xmin=0, xmax=1, color='grey', lw=1, ls='--')
    ax_main.axvline(x=0, ymin=0, ymax=1, color='grey', lw=1, ls='--')
    ax_main.set_xlim(-0.2,0.2)
    ax_main.set_ylim(-0.2,0.2)
    ax_main.set_xlabel(r'$\Delta$RA$_\mathrm{Merian-GAIA}$')
    ax_main.set_ylabel(r'$\Delta$Dec$_\mathrm{Merian-GAIA}$')
   

    # RIGHT PANEL
    ax_main2 = plt.subplot(gs[1:4, 5:8])
    ax_xDist2 = plt.subplot(gs[0, 5:8],sharex=ax_main)
    ax_yDist2 = plt.subplot(gs[1:4, 8],sharey=ax_main)
    ax_xDist2.tick_params(axis='x', labelbottom=False, labelleft=False)
    ax_yDist2.tick_params(axis='y', labelbottom=False, labelleft=False)
    ax_xDist2.tick_params(axis='y', labelbottom=False, labelleft=False)
    ax_yDist2.tick_params(axis='x', labelbottom=False, labelleft=False)

    ax_main2.scatter(dra_hsc, ddec_hsc, color='g', alpha=0.2)
    ax_main2.scatter(dra_median_hsc, ddec_median_hsc, marker='P', s=300, color='k')
    ax_main2.axhline(y=0, xmin=0, xmax=1, color='grey', lw=1, ls='--')
    ax_main2.axvline(x=0, ymin=0, ymax=1, color='grey', lw=1, ls='--')
    ax_main2.set_xlim(-0.2,0.2)
    ax_main2.set_ylim(-0.2,0.2)
    ax_main2.set_xlabel(r'$\Delta$RA$_\mathrm{Merian-HSC}$')
    ax_main2.set_ylabel(r'$\Delta$Dec$_\mathrm{Merian-HSC}$')
    
    ax_xDist2.hist(dra_hsc, bins=bins, color='g', histtype='step',fill =True, alpha=0.2, lw=1)
    ax_yDist2.hist(ddec_hsc, bins=bins, color='g', orientation='horizontal',align='mid', histtype='step', fill=True, alpha=0.2, lw=1)
    ax_xDist2.axvline(x=0, ymin=0, ymax=1, color='grey', lw=1, ls='--')
    ax_yDist2.axhline(y=0, xmin=0, xmax=1, color='grey', lw=1, ls='--')
    ax_xDist2.text(0.16, 0.8, r'$\langle\Delta$RA$\rangle$=%.3f'%(dra_median_hsc), horizontalalignment='center', verticalalignment='center', transform=ax_xDist2.transAxes, fontsize=25)
    ax_yDist2.text(0.5, 0.93, r'$\langle\Delta$Dec$\rangle$=%.3f'%(ddec_median_hsc), horizontalalignment='center', verticalalignment='center', transform=ax_yDist2.transAxes, fontsize=25)


    plt.subplots_adjust(wspace=0, hspace=0)
    

    return fig




def astrometry_crossmatch(merian_tab):


    median_ra = np.nanmedian(merian_tab['coord_ra_Merian'])
    median_dec = np.nanmedian(merian_tab['coord_dec_Merian'])

    offset_ra_dcos = merian_tab['coord_ra_Merian'] * np.cos(np.deg2rad(merian_tab['coord_dec_Merian']) )
    offset_ra = offset_ra_dcos.max() - offset_ra_dcos.min()
    offset_dec = merian_tab['coord_dec_Merian'].max() - merian_tab['coord_dec_Merian'].min()
    
    #-----------------------
    # Point-Source selection
    #-----------------------

    # Merian point source selection: point-source in all filters
    ps_flg = (merian_tab['g_extendedness_value_HSCS20A'] == 0) & (merian_tab['i_extendedness_value_HSCS20A']==0) 
    ps_flg &= (merian_tab['r_extendedness_value_HSCS20A'] == 0) & (merian_tab['z_extendedness_value_HSCS20A']==0)
    ps_tab = merian_tab[ps_flg]

    # Merian magnitude cut: g band magnitude < 20 
    gmag = -2.5 * np.log10(ps_tab['g_cModelFlux_Merian']) + 31.4
    ps_tab = ps_tab[(gmag <= 20.)]


    #-----------------------
    # GAIA  matching
    #-----------------------

    # GAIA reference catalog
    coord = SkyCoord(ra=median_ra, dec=median_dec, unit=(u.deg, u.deg), frame='icrs')
    width = u.Quantity(offset_ra, u.deg)
    height = u.Quantity(offset_dec, u.deg)
    gaia_tbl =  query_gaia_refcat(coord, width, height, pm_cut=10)
    gaia_ra = gaia_tbl['RAJ2000']
    gaia_dec = gaia_tbl['DECJ2000']

    # Match GAIA
    coords_gaia = SkyCoord(ra=gaia_ra, dec=gaia_dec, unit=(u.deg, u.deg))
    coords_merian = SkyCoord(ra=ps_tab['coord_ra_Merian'], dec=ps_tab['coord_dec_Merian'], unit=(u.deg, u.deg))

    idx_merian, d2d, d3d = coords_gaia.match_to_catalog_sky(coords_merian)
    matched = d2d < 0.5 * u.arcsec

    matched_gaia = gaia_tbl[matched]
    matched_merian = ps_tab[idx_merian[matched]]

    dra = matched_merian['coord_ra_Merian'] - matched_gaia['RA']
    ddec = matched_merian['coord_dec_Merian'] - matched_gaia['DEC']

    dra_gaia = dra * 3600
    ddec_gaia = ddec * 3600

   
    dra_median_gaia = sigma_clipped_stats(dra_gaia, sigma=3)[1]
    ddec_median_gaia = sigma_clipped_stats(ddec_gaia, sigma=3)[1]

    #-----------------------
    # HSC  matching
    #-----------------------

    # Match HSC
    ps_coord_merian = SkyCoord(ra=ps_tab['coord_ra_Merian'], dec=ps_tab['coord_dec_Merian'], unit=(u.deg, u.deg))
    ps_coord_hsc = SkyCoord(ra=ps_tab['ra_HSCS20A'], dec=ps_tab['dec_HSCS20A'], unit=(u.deg, u.deg))

    idx, d2d, d3d = ps_coord_merian.match_to_catalog_sky(ps_coord_hsc)
    matched = d2d < 0.5 * u.arcsec
    matched &= ps_tab['hsc_match'].astype(bool)
    matched_ps = ps_tab[matched]
    
    # # the 50% brightest are selected
    # g_50, g_45, g_55 = np.nanpercentile(matched_ps['g_gaap2p5Flux_Merian'], [50, 45, 55])
    # matched_ps = matched_ps[ (matched_ps['g_gaap2p5Flux_Merian'] > g_50)]

    # remove invalid values
    flg = (matched_ps['dec_HSCS20A'] !=0 ) & (matched_ps['ra_HSCS20A'] !=0)

    dra_ps = (matched_ps['coord_ra_Merian'] - matched_ps['ra_HSCS20A'])[flg]
    ddec_ps = (matched_ps['coord_dec_Merian'] - matched_ps['dec_HSCS20A'])[flg]

    dra_hsc = dra_ps * 3600
    ddec_hsc = ddec_ps * 3600

    dra_median_hsc = sigma_clipped_stats(dra_hsc, sigma=3)[1]
    ddec_median_hsc = sigma_clipped_stats(ddec_hsc, sigma=3)[1]
    
    return dra_gaia, ddec_gaia, dra_median_gaia, ddec_median_gaia, dra_hsc, ddec_hsc, dra_median_hsc, ddec_median_hsc


def save_astrom_array(outname, 
                      dra_gaia, ddec_gaia, dra_median_gaia, ddec_median_gaia, 
                      dra_hsc, ddec_hsc, dra_median_hsc, ddec_median_hsc):
    
    # save the array
    from astropy.io import fits
    hdul = fits.HDUList([fits.PrimaryHDU()])

    gaia_tab = Table()
    gaia_tab.meta['dra_median'] = dra_median_gaia
    gaia_tab.meta['ddec_median'] = ddec_median_gaia
    # gaia_tab.meta['dra_std_median'] = dra_std / np.sqrt(len(dra_gaia))
    # gaia_tab.meta['ddec_std_median'] = ddec_std / np.sqrt(len(ddec_gaia))
    gaia_tab['dra'] = dra_gaia
    gaia_tab['ddec'] = ddec_gaia
    hdul.append(fits.BinTableHDU(gaia_tab, name='GAIA'))

    hsc_tab = Table()
    hsc_tab.meta['dra_median'] = dra_median_hsc
    hsc_tab.meta['ddec_median'] = ddec_median_hsc
    # hsc_tab.meta['dra_std_median'] = dra_std_ps / np.sqrt(len(dra_ps))
    # hsc_tab.meta['ddec_std_median'] = ddec_std_ps / np.sqrt(len(ddec_ps))
    hsc_tab['dra'] = dra_hsc
    hsc_tab['ddec'] = ddec_hsc
    hdul.append(fits.BinTableHDU(hsc_tab, name='HSC'))

    hdul.writeto(outname, overwrite=True)


if __name__ == '__main__':

  
    cosmos_tract = [9813]
    random_tract = [8280,8525,8764,8769,9010,9078,9088,9099,9102,9131,9135,9227,9313,9319,9344,9378,9456,9459,9470,9556,9564,9574,9589,9798,9800,9812,9828,9939,9946,9953,10040,10048,10061,10183,10293,10427]

    all_tract = cosmos_tract + random_tract

    for tractnum in all_tract:
        tab = Table.read(f'/scratch/gpfs/sd8758/merian/catalog/DR1/{tractnum}/meriandr1_use_{tractnum}_S20A.fits')
       
        if len(tab) == 0: 
            print('Tract %d is empty' % tractnum)
            
            with open('empty_tract.txt', 'a') as f:
                f.write(f'{tractnum}\n')

            continue

        dra_gaia, ddec_gaia, dra_median_gaia, ddec_median_gaia, dra_hsc, ddec_hsc, dra_median_hsc, ddec_median_hsc = astrometry_crossmatch(tab)

        save_astrom_array(f'DR1_astrom_calib_gaia_hsc_{tractnum}.fits', 
                          dra_gaia, ddec_gaia, dra_median_gaia, ddec_median_gaia, 
                          dra_hsc, ddec_hsc, dra_median_hsc, ddec_median_hsc)
        

        fig = astrometry_plotter(dra_gaia, ddec_gaia, dra_hsc, ddec_hsc,
                                    dra_median_gaia, ddec_median_gaia, dra_median_hsc, ddec_median_hsc)
        fig.savefig(f'DR1_astrom_calib_gaia_hsc_{tractnum}.png')

        plt.clf()
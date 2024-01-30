
"""
Utils to plot QA for sky objects for coaads and exposures of tracts.

-----------------------
Author: Xiaojing Lin
Date: 2023-09-27
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import os
from matplotlib.patches import Rectangle

import copy
import pickle
import warnings

import numpy as np
from astropy.table import Table, QTable, Column

import lsst.daf.butler as dafButler
import lsst.afw.image as afwImage

plt.rcParams["font.family"] = "STIXGeneral"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'out'  
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True
plt.rcParams['font.size']=15



plt.rcParams['xtick.bottom'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True



def plot_coadd2d_oneTract(cat,tract = 9813,colname = 'N708_ap09Flux',fname='Tract.png',xlim=None):

    """
    plot statistics in historgram and PATCH layout for coadds of tracts.

    Parameter:
        cat: Table. contained columns: tract, colname
        tract: int. ID of targeted tract
        colname: str. Name of targeted parameter in the catalog.
        fname: str. Name of the saved plot. If fname is not provided, the plot will be shown via plt.show(). Otherwise it would be saved and then closed.
        xlim: None or list. xlim of the histogram.

    Return:
        median: float. median value across the targeted tract.
        std: float. standard deviation across the targeted tract.
        mid_patch: array. median values of each patch in this tract.
        std_patch: array. standard deviation of each patch in this tract.
    """

  
    cat_coadd = cat[cat['tract']==tract]

    flux_raw = cat_coadd[colname]
    mask, flux = sigmaclip(flux_raw,low=3,high=3)
    Neff_aper = np.sum(mask)
    
    
    median = np.nanmedian(flux)
    std = np.sqrt(np.nanmean((flux-median)**2))
    

    vlim_max = median + 3*std
    vlim_min = median - 3*std
    vmax = np.nanmax(flux)
    vmin = np.nanmin(flux)

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(20,8))

    # Left panel: histogram
    ww = (flux<vlim_max) & (flux>vlim_min)
    if xlim is not None:
        ww &= (flux>xlim[0]) &(flux<xlim[1])
    n,bins,_ = ax1.hist(flux[ww],bins='auto',facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5,alpha=0.5)
    ax1.axvline(0,ls='--',color='#FF6666',lw=0.8)
    ax1.axvline(median,ls='--',color='black',lw=0.8)
    ax1.axvline(median-std,ls='--',color='black',lw=0.8)
    ax1.axvline(median+std,ls='--',color='black',lw=0.8)
    ax1.text(0.71,0.95,s='$\mu$=%.2f\n$\sigma$=%.2f\n\nmax=%.2f\nmin=%.2f\n\nNeff/Tot: %d/%d'%(median,std,vmax,vmin,Neff_aper,len(flux_raw)),transform=ax1.transAxes,verticalalignment='top')
    if xlim is None:
        ax1.set_xlim(vlim_min,vlim_max)
    else:
        ax1.set_xlim(xlim[0],xlim[-1])
    ax1.set_xlabel(colname)
    ax1.set_ylabel('Num')
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_title(' '.join((colname, 'Coadd Tract', str(tract))),fontsize=20)

    # Right panel: Patches
    patch_list = np.unique(cat_coadd['patch'])
    nside = int(np.sqrt(np.nanmax(patch_list)+1))
    patch_matrix = np.zeros((nside,nside))
    cmap=plt.cm.seismic
    norm = plt.Normalize(-median-std,+median+std)
    mid_patch = []
    std_patch = []

    for patch_id in patch_list:
        patch_cat = cat_coadd[cat_coadd['patch']==patch_id]
        patch_flux_raw = patch_cat[colname]
        mask, patch_flux = sigmaclip(patch_flux_raw)
        Neff = len(patch_flux)
        
        mid_flux = np.nanmedian(patch_flux)
        std_flux = np.sqrt(np.nanmean((patch_flux-mid_flux)**2))
        
        mid_patch.append(mid_flux)
        std_patch.append(std_flux)
        

        i_row = patch_id//nside
        i_column = patch_id%nside
        patch_matrix[i_row,i_column] += mid_flux

        c = 'black' if abs(norm(mid_flux)-0.5) < 0.1 else 'white'
        t1 = ax2.text(i_column+0.25, i_row+0.25, str(Neff), horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 12.5})
        t2 = ax2.text(i_column-0.25, i_row-0.25, '%.2f'%mid_flux, horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 7.5})
        t3 = ax2.text(i_column+0.25, i_row-0.25, '%.2f'%std_flux, horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 7.5})

    vnorm = np.nanmax(abs(patch_matrix))
    im = ax2.imshow(patch_matrix,cmap='seismic',origin='lower',vmin=-vnorm,vmax=vnorm)
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.set_xlabel('Patch')
    ax2.set_ylabel('Patch')
    ax2.set_title('Coadd Tract %d'%tract)
    plt.minorticks_off()
    fig.colorbar(im, ax=ax2, label='Median')

    if fname is None:
        plt.show()
    else:
        plt.savefig(fname,bbox_inches='tight',dpi=300,facecolor='white', transparent=False)
        plt.close()
        
    return median, std, np.array(mid_patch), np.array(std_patch)
        

def plot_single2d_byCCD(cat_single, exp_id, colname, filter_name = 'N540',width=0.13,height=0.075,fname='Single_byCCD.png',xlim=None):

    """
    plot statistics in historgram and CCD layout for single exposures

    Parameter
        cat_single: Table. contained column:  ccdVisitId,  coord_ra,  coord_dec, colname
        exp_id: int. Visit ID of this exposure.
        colname: str. Name of targeted parameters.
        filter_name: str. Targeted filter.
        width: float. margin width.
        height: float. margin height.
        fname: str. Name of saved file. If fname is not provided, the plot will be shown via plt.show(). Otherwise it would be saved and then closed.
        xlim: None or list. xlim of the histogram.

    Return:
        median: float. median value of this exposure.
        std: float. standard deviation of this exposure
        mid_ccd: array. median values of each CCD in this exposure
        std_ccd: array. standard deviation of each CCD in this exposure
    """
  
    ccd_ids = cat_single['ccdVisitId']-exp_id*100
    ras = cat_single['coord_ra']
    decs = cat_single['coord_dec']
    ccd_list = np.unique(ccd_ids)
    flux_raw = cat_single[colname]
    mask,flux = sigmaclip(flux_raw)
    Neff_aper = np.sum(mask)
    
    median = np.nanmedian(flux)
    std = np.sqrt(np.nanmean((flux-median)**2))
    vlim_max = median + 3*std
    vlim_min = median - 3*std
    vmax = np.nanmax(flux)
    vmin = np.nanmin(flux)



    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(24,8))

    # Left panel: histogram
    n,bins,_ = ax1.hist(flux,bins='auto',facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5,alpha=0.5)
    ax1.axvline(0,ls='--',color='#FF6666',lw=0.8)
    ax1.axvline(median,ls='--',color='black',lw=0.8)
    ax1.axvline(median-std,ls='--',color='black',lw=0.8)
    ax1.axvline(median+std,ls='--',color='black',lw=0.8)
    ax1.text(0.71,0.95,s='$\mu$=%.2f\n$\sigma$=%.2f\n\nmax=%.2f\nmin=%.2f\n\nNeff/Tot: %d/%d'%(median,std,vmax,vmin,Neff_aper,len(flux_raw)),transform=ax1.transAxes,verticalalignment='top')
    if xlim is None:
        ax1.set_xlim(vlim_min,vlim_max)
    else:
        ax1.set_xlim(xlim[0],xlim[-1])
    ax1.set_xlabel(colname)
    ax1.set_ylabel('Num')
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_title(' '.join((filter_name,'Exp',str(exp_id),colname)),fontsize=20)


    # RIGHT panel: CCD layout
    cmap=plt.cm.seismic
    norm = plt.Normalize(-median-std,+median+std)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    ax2.set_xlim(min(ras)-width/2,max(ras)+width/2)
    ax2.set_ylim(min(decs)-height/2,max(decs)+height/2)
    
    mid_ccd = []
    std_ccd = []
    for ccd_id in ccd_list:
        selected_ra = ras[ccd_ids==ccd_id]
        selected_dec = decs[ccd_ids==ccd_id]
        selected_flux_raw = flux[ccd_ids==ccd_id]
        
        mask, selected_flux = sigmaclip(selected_flux_raw)
        Neff = np.sum(mask)
        
        mid_ra = np.nanmean(selected_ra)
        mid_dec = np.nanmean(selected_dec)
        mid_flux = np.nanmedian(selected_flux)
        std_flux = np.sqrt(np.nanmean((selected_flux-mid_flux)**2))
        
        mid_ccd.append(mid_flux)
        std_ccd.append(std_flux)
        

        c = 'black' if abs(norm(mid_flux)-0.5) < 0.1 else 'white'

        trans = ax2.transData.transform((mid_ra,mid_dec))
        trans = ax2.transAxes.inverted().transform(trans)
        x, y = trans[0],trans[1]
        t1 = ax2.text(x-width*0.25, y+0.25*height, str(ccd_id), horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 7.5},transform=ax2.transAxes)

        t2 = ax2.text(x+width*0.2, y+0.2*height, str(Neff), horizontalalignment='center', verticalalignment='center',
                        fontdict={'color':c, 'size': 12.5},transform=ax2.transAxes)
        
        t3 = ax2.text(x-width*0.2, y-0.2*height, str('%.2f'%mid_flux), horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 12.5},transform=ax2.transAxes)
        t4 = ax2.text(x+width*0.25, y-height*0.25, str('%.2f'%std_flux), horizontalalignment='center', verticalalignment='center',
                        fontdict={'color': c, 'size': 7.5},transform=ax2.transAxes)
        
        rect = Rectangle((x - 0.5*width, y - 0.5*height), width=width,height=height, linewidth=1, ec='white',color=cmap(norm(mid_flux)) ,transform=ax2.transAxes)
        ax2.add_patch(rect)

    ax2.set_xlabel('RA [deg]')
    ax2.set_ylabel('DEC [deg]')
    ax2.set_title('Exp%d CCD layout: '%exp_id + colname)
    sm.set_array([])
    fig.colorbar(sm, ax=ax2, label='Median')
    
    if fname is None:
        plt.show()
    else:    
        plt.savefig(fname,bbox_inches='tight',dpi=300,facecolor='white', transparent=False) 
        plt.close()
    
    return median,std, np.array(mid_ccd), np.array(std_ccd)
 



def plot_skyqa_for_Tract(butler, tract_info,
                         skymap='hsc_rings_v1',
                         mode='both',
                         single_savedir = '/home/ml6704/tigress/public_html/singleCCD/catalog/',
                         coadd_savedir =  '/home/ml6704/tigress/public_html/coadd/catalog/',
                         plot_kwargs={'single':{'width':0.13,'height':0.075,'xlim':None},
                                          'coadd':{'xlim':None}}):

    """
    Plot Quality Assessment for sky objects for coaads and exposures of tracts.
    For single exposures,  show statistics in histograms and CCD layout. For coadds, show statistics in histograms and PATCHs layout.

    Parameter:
        butler: butler that initialized to the data repository
        tract_info: Dict. 
                    a dictionary that contains the information of QA targets. Format:
                    {'tract': list,  'visits': list, 'version': str, e.g. '202202',
                     'SingleCollection': str, data repository for single exposures. e.g. 'DECam/runs/merian/w_2022_02',
                     'CoaddCollection': str, data repository for coadds. e.g. 'DECam/runs/merian/w_2022_02/t%d_wideTEST'%tract_id,
                     'colnames': list containing colname recorded in the input catalog, e.g. ['ap09Flux', 'psfFlux'],
                     'filters': list containing filtername recorded in the input catalog, e.g. ['N708','N540']} 
        skymap: str. value of `skymap' being input in butler
        mode: str. 'single', 'coadd', 'both', 'coadd-patch'.
                    
        single_savedir: str. the directory to save QA plots for single exposure.
        coadd_savedir: str. the directory to save QA plots for coadd exposure.
        plot_kwargs: dict. args for plots.  Format (optional):
                     {'single': {'width':float,'height':float,'xlim':float or None},
                       'coadd':{'xlim':float or None}}   

    Return:
        z_visit: dict. Dictionary that contains statistics for single exposures.
        z_coadd: dict. Dictionary that contains statistics for coadds. 
          
    """

    if mode == 'both': run_single = True; run_coadd = True; run_patch = False; 
    if mode == 'single': run_single = True; run_coadd = False; run_patch = False; 
    if mode == 'coadd': run_single = False; run_coadd  = True; run_patch = False; 
    if mode == 'coadd-patch': run_single = False; run_coadd  = True; run_patch = True; 

    z_visit = {}
    z_coadd = {}
    #-------- Single exposure ---------#
    if run_single:
        print('SIngle')
        visits_list = tract_info['visits']
        tract = tract_info['tract']
        savedir = single_savedir + 'v'+ tract_info['version'] + '/' + str(tract) + '/' 
        print('SIngle will be saved in ', savedir)
        if not os.path.isdir(savedir):
            print('create dir:', savedir)
            os.makedirs(savedir)
        
        for c in tract_info['colnames']:
            for ff in tract_info['filters']:
                z_visit[ff+'_'+c] = {'med_single':[],'std_single':[],'med_ccd':[],'std_ccd':[],'num_visit':[]}
                z_coadd[ff+'_'+c] = {'med_coadd':[],'std_coadd':[],'med_patch':[],'std_patch':[]}
        
        for num_visit in visits_list:
            try:
                cat_single = butler.get(
                "sourceTable_visit", instrument="DECam", visit=num_visit, collections=tract_info['SingleCollection'])
                sky_flag = (cat_single['sky_source'] == True)
                skyobj_single = cat_single[sky_flag]

                filter_name = skyobj_single.iloc[0]['band']
            
                for colname in tract_info['colnames']:

                    fname = filter_name + '_' + tract_info['version'] + '-' + str(tract) + '-'+ str(num_visit) + '_' + colname + '.png'
                    print(fname)
                    fname = savedir + fname
                    med_single, std_single, med_ccd, std_ccd = plot_single2d_byCCD(skyobj_single, num_visit, colname, filter_name = filter_name,
                            fname=fname,**plot_kwargs['single'])
                    z_visit[filter_name+'_'+colname]['med_single'].append(med_single)
                    z_visit[filter_name+'_'+colname]['std_single'].append(std_single)
                    z_visit[filter_name+'_'+colname]['med_ccd'].append(med_ccd)
                    z_visit[filter_name+'_'+colname]['std_ccd'].append(std_ccd)
                    z_visit[filter_name+'_'+colname]['num_visit'].append(num_visit)
                    
            
            except Exception as err:
                print(err)

    
    #-------- Coadds ---------#
    if run_coadd:
        print('COadd')
        tract = tract_info['tract']
        savedir = coadd_savedir + 'v'+ tract_info['version'] + '/' + str(tract) + '/' 
        print('COadd will be saved in ', savedir)
        if not os.path.isdir(savedir):
            print('create dir:', savedir)
            os.makedirs(savedir)
            
        cat_coadd = butler.get(
        'objectTable_tract', tract=tract, instrument='DECam',
        skymap=skymap, collections=tract_info['CoaddCollection'])
        skyobj_cat = cat_coadd[cat_coadd.merge_peak_sky==True]
        
    
        for colname in tract_info['colnames']:
            for filter_name in tract_info['filters']:
                try:
                    fname = filter_name + '_' + tract_info['version'] + '-' + str(tract) + '_' + colname + '.png'
                    print(fname)
                    fname = savedir + fname
                    med_coadd,std_coadd, med_patch, std_patch = plot_coadd2d_oneTract(skyobj_cat, tract = tract,colname = filter_name + '_' + colname,
                        fname=fname,**plot_kwargs['coadd'])
                    z_coadd[filter_name+'_'+colname]['med_coadd'].append(med_coadd)
                    z_coadd[filter_name+'_'+colname]['std_coadd'].append(std_coadd)
                    z_coadd[filter_name+'_'+colname]['med_patch'].append(med_patch)
                    z_coadd[filter_name+'_'+colname]['std_patch'].append(std_patch)
                except Exception as err:
                    print(err)

    return z_visit, z_coadd
           

def sigmaclip(a, low=3., high=3.): 
    
    # mask=0 means the value should be discarded
    
    # clip nan
    mask = np.ones(a.shape,dtype=int)
    mask[np.isnan(a)] = 0
    c = a[mask==1]
    
    delta = 1
    while delta:
        c_std = np.sqrt(np.nanmean((c - np.nanmedian(c))**2))
        c_median =  np.nanmedian(c)
        size = c.size
        critlower = c_median - c_std * low
        critupper = c_median + c_std * high
        mask[(a <= critlower) | (a >= critupper)] = 0
        c = c[(c >= critlower) & (c <= critupper)]
        delta = size - c.size

    return mask, a[mask==1]
    



def plot_violin(ax,data,pos,color='red',label=None):
    """
    Plot violin plot for the input data.

    Parameter:
    ax: matplotlib.axes.Axes. Axes to plot on.
    data: array. Data to plot. in shape of (N,).
    pos: float. Position of the violin plot on the x-axis for the input data.
    color: str. Color of the violin plot.
    label: str. Label of the violin plot.

    Return:
    None
    """
    parts = ax.violinplot(data, [pos],vert=True, points=60, widths=0.75, showmeans=False,
                      showextrema=False, showmedians=False, quantiles=[0.16,0.50,0.84],
                       bw_method=0.5)
    for pc in parts['bodies']:
        pc.set_facecolor(color)
        pc.set_edgecolor(color)
        pc.set_alpha(0.1)
    parts['cquantiles'].set_edgecolor(color)
    parts['cquantiles'].set_alpha(0.2)
    parts['cquantiles'].set_linestyle('--')
    
    quartile1, medians, quartile3 = np.percentile(data, [16, 50, 84])
    
    ax.scatter(pos, medians, marker='o', color='white', s=30, zorder=3)
    ax.vlines(pos, quartile1, quartile3, color=color, linestyle='-', lw=5)
    
    if label:
        y = np.percentile(data,99)
        ax.text(pos+0.05,y,s=label,color=color)
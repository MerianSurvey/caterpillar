#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Visulization tools."""

import numpy as np

from astropy.visualization import ZScaleInterval
from astropy.visualization import AsymmetricPercentileInterval

import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('text', usetex=False)

__all__ = ['display_single', 'display_all']


def display_single(img,
                   pixel_scale=0.168,
                   physical_scale=None,
                   xsize=8,
                   ysize=8,
                   ax=None,
                   alpha=1.0,
                   stretch='arcsinh',
                   scale='zscale',
                   zmin=None,
                   zmax=None,
                   contrast=0.25,
                   no_negative=False,
                   lower_percentile=1.0,
                   upper_percentile=99.0,
                   cmap='viridis',
                   scale_bar=True,
                   scale_bar_length=5.0,
                   scale_bar_fontsize=20,
                   scale_bar_y_offset=0.5,
                   scale_bar_color='w',
                   scale_bar_loc='left',
                   color_bar=False,
                   color_bar_loc=1,
                   color_bar_width='75%',
                   color_bar_height='5%',
                   color_bar_fontsize=18,
                   color_bar_color='w',
                   add_text=None,
                   text_fontsize=30,
                   text_color='w'):
    """Display single image.
    Parameters
    ----------
        img: np 2-D array for image
        xsize: int, default = 8
            Width of the image.
        ysize: int, default = 8
            Height of the image.
    """
    if ax is None:
        fig = plt.figure(figsize=(xsize, ysize))
        ax1 = fig.add_subplot(111)
    else:
        ax1 = ax
    ax1.grid(False)

    # Stretch option
    if img.ndim == 3:
        img_scale = img
        vmin, vmax = None, None
    else:
        if stretch.strip() == 'arcsinh':
            img_scale = np.arcsinh(img)
            if zmin is not None:
                zmin = np.arcsinh(zmin)
            if zmax is not None:
                zmax = np.arcsinh(zmax)
        elif stretch.strip() == 'log':
            if no_negative:
                img[img <= 0.0] = 1.0E-10
            img_scale = np.log(img)
            if zmin is not None:
                zmin = np.log(zmin)
            if zmax is not None:
                zmax = np.log(zmax)
        elif stretch.strip() == 'log10':
            if no_negative:
                img[img <= 0.0] = 1.0E-10
            img_scale = np.log10(img)
            if zmin is not None:
                zmin = np.log10(zmin)
            if zmax is not None:
                zmax = np.log10(zmax)
        elif stretch.strip() == 'linear':
            img_scale = img
        else:
            raise Exception("# Wrong stretch option.")

        # Scale option
        if scale.strip() == 'zscale':
            try:
                vmin, vmax = ZScaleInterval(contrast=contrast).get_limits(img_scale)
            except IndexError:
                # TODO: Deal with problematic image
                vmin, vmax = -1.0, 1.0
        elif scale.strip() == 'percentile':
            try:
                vmin, vmax = AsymmetricPercentileInterval(
                    lower_percentile=lower_percentile,
                    upper_percentile=upper_percentile).get_limits(img_scale)
            except IndexError:
                # TODO: Deal with problematic image
                vmin, vmax = -1.0, 1.0
        elif scale.strip() == 'minmax':
            vmin, vmax = np.nanmin(img_scale), np.nanmax(img_scale)
        else:
            vmin, vmax = np.nanmin(img_scale), np.nanmax(img_scale)

        if zmin is not None:
            vmin = zmin
        if zmax is not None:
            vmax = zmax

    show = ax1.imshow(img_scale, origin='lower', cmap=cmap, interpolation='none',
                      vmin=vmin, vmax=vmax, alpha=alpha)

    # Hide ticks and tick labels
    ax1.tick_params(
        labelbottom=False,
        labelleft=False,
        axis=u'both',
        which=u'both',
        length=0)

    # Put scale bar on the image
    if img.ndim == 3:
        img_size_x, img_size_y = img[:, :, 0].shape
    else:
        img_size_x, img_size_y = img.shape

    if physical_scale is not None:
        pixel_scale *= physical_scale

    if scale_bar:
        if scale_bar_loc == 'left':
            scale_bar_x_0 = int(img_size_x * 0.04)
            scale_bar_x_1 = int(img_size_x * 0.04 +
                                (scale_bar_length / pixel_scale))
        else:
            scale_bar_x_0 = int(img_size_x * 0.95 -
                                (scale_bar_length / pixel_scale))
            scale_bar_x_1 = int(img_size_x * 0.95)
        scale_bar_y = int(img_size_y * 0.10)
        scale_bar_text_x = (scale_bar_x_0 + scale_bar_x_1) / 2
        scale_bar_text_y = (scale_bar_y * scale_bar_y_offset)
        if physical_scale is not None:
            scale_bar_text = r'$%d\ \mathrm{kpc}$' % int(scale_bar_length)
        else:
            scale_bar_text = r'$%d^{\prime\prime}$' % int(scale_bar_length)
        scale_bar_text_size = scale_bar_fontsize

        ax1.plot(
            [scale_bar_x_0, scale_bar_x_1], [scale_bar_y, scale_bar_y],
            linewidth=3,
            c=scale_bar_color,
            alpha=1.0)
        ax1.text(
            scale_bar_text_x,
            scale_bar_text_y,
            scale_bar_text,
            fontsize=scale_bar_text_size,
            horizontalalignment='center',
            color=scale_bar_color)
    if add_text is not None:
        text_x_0 = int(img_size_x*0.08)
        text_y_0 = int(img_size_y*0.80)
        ax1.text(
            text_x_0, text_y_0, r'$\mathrm{'+add_text+'}$',
            fontsize=text_fontsize, color=text_color)

    # Put a color bar on the image
    if color_bar:
        ax_cbar = inset_axes(ax1,
                             width=color_bar_width,
                             height=color_bar_height,
                             loc=color_bar_loc)
        if ax is None:
            cbar = plt.colorbar(show, ax=ax1, cax=ax_cbar,
                                orientation='horizontal')
        else:
            cbar = plt.colorbar(show, ax=ax, cax=ax_cbar,
                                orientation='horizontal')

        cbar.ax.xaxis.set_tick_params(color=color_bar_color)
        cbar.ax.yaxis.set_tick_params(color=color_bar_color)
        cbar.outline.set_edgecolor(color_bar_color)
        plt.setp(plt.getp(cbar.ax.axes, 'xticklabels'),
                 color=color_bar_color, fontsize=color_bar_fontsize)
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'),
                 color=color_bar_color, fontsize=color_bar_fontsize)

    if ax is None:
        return fig
    return ax1


def display_all(img_list, n_column=3, img_size=3., hdu_index=None, label_list=None,
                cmap_list=None, label_x=0.1, label_y=0.9, fontsize=20, fontcolor='k',
                hdu_list=False, hdu_start=1, **kwargs):
    """Display a list of images."""
    if not isinstance(img_list, list):
        raise TypeError("Provide a list of image to show or use display_single()")

    # Make a numpy array list if the input is HDUList
    if hdu_list:
        img_list = [img_list[ii].data for ii in np.arange(len(img_list))[hdu_start:]]

    if cmap_list is not None:
        assert len(cmap_list) == len(img_list), "Wrong number of color maps!"

    if label_list is not None:
        assert len(label_list) == len(img_list), "Wrong number of labels!"

    # Number of image to show
    n_img = len(img_list)

    if n_img <= n_column:
        n_col = n_img
        n_row = 1
    else:
        n_col = n_column
        n_row = int(np.ceil(n_img / n_column))

    fig = plt.figure(figsize=(img_size * n_col, img_size * n_row))
    fig.subplots_adjust(left=0., right=1., bottom=0., top=1., wspace=0., hspace=0.)

    gs = gridspec.GridSpec(n_row, n_col)
    gs.update(wspace=0.0, hspace=0.00)

    for ii in range(n_img):
        if hdu_index is None:
            img_show = img_list[ii]
        else:
            img_show = img_list[ii][hdu_index].data

        ax = plt.subplot(gs[ii])
        if cmap_list is not None:
            ax = display_single(img_show, cmap=cmap_list[ii], ax=ax, **kwargs)
        else:
            ax = display_single(img_show, ax=ax, **kwargs)

        if label_list is not None:
            if len(label_list) != n_img:
                print("# Wrong number for labels!")
            else:
                ax.text(label_x, label_y, label_list[ii], fontsize=fontsize,
                        transform=ax.transAxes, color=fontcolor)

    return fig

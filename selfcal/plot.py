import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import mad_std
from astropy.visualization import simple_norm
import matplotlib as mpl
from astropy.io import fits


import os, glob, subprocess, time, logging, math
from typing import Callable, Any
import casatools
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.patches as patches
from matplotlib.ticker import ScalarFormatter
from typing import Union, Tuple, List
from pathlib import Path

import numpy as np
from radio_beam import Beam
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from casatasks import *


def set_rc_params():
    """
    Global configuration for matplotlib.pyplot
    """
    mpl.rcParams.update({
        'font.size': 16,
        "text.usetex": False,  # Disable LaTeX
        "font.family": "sans-serif",
        'mathtext.fontset': 'stix',
        'font.weight': 'medium',
        'xtick.labelsize': 16,
        'figure.figsize': (8,6),
        'ytick.labelsize': 16,
        'axes.labelsize': 16,
        'xtick.major.width': 1,
        'ytick.major.width': 1,
        'axes.linewidth': 1.5,
        'axes.edgecolor': 'orange',
        'lines.linewidth': 2,
        'legend.fontsize': 14,
        'grid.linestyle': '--',
        'axes.grid.which': 'major',
        'axes.grid.axis': 'both',
        'axes.spines.right': True,
        'axes.grid': True,
        'axes.titlesize': 16,
        'legend.framealpha': 1.0
    })



def plot_image_with_beam(imagename):
        
        """
        Plot the FITS image and place the beam at the bottom-left corner.
        """

        set_rc_params()

        figsize=(8,6)
        color='magma_r'
        # Read the FITS file data
        hdu = fits.open(imagename)
        image_data = hdu[0].data[0, 0, :, :] 
        header = hdu[0].header
        w = WCS(header, naxis=2)
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': w})
        

        num_contours = 5
        image_rms = imstat(imagename=imagename)['rms'][0]
        peak_flux = imstat(imagename=imagename)['max'][0]
        extent = None
        vmax = 1.05*peak_flux
        vmin = 1*image_rms

        # color_map = 'magma_r'
        color_map = 'magma_r'

        # # Define the contour color palette
        contour_palette_ = ['#000000', '#444444', '#666666', '#EEEEEE', '#EEEEEE', '#FFFFFF']
        contour_palette = contour_palette_ if '_r' in color_map else contour_palette_[::-1]

        # # Set normalization for plotting
        norm0 = simple_norm(image_data, stretch='linear', max_percent=99.0)
        norm = simple_norm(image_data, stretch='linear', asinh_a=0.02, min_cut=vmin, max_cut=vmax)

        im = ax.imshow(image_data, cmap=color_map, origin='lower', norm=norm)
        ax.coords[0].set_auto_axislabel(True) 
        ax.coords[1].set_auto_axislabel(True) 
        # shape = header['NAXIS1'], header['NAXIS2']
        # bmaj, bmin, pa = get_beam()
        # relative_x = 15
        # relative_y = 15
        # x_pos = (relative_x / 320) * shape[0]  
        # y_pos = (relative_y / 320) * shape[1]  
        # beam_ellipse = patches.Ellipse(
        #     (x_pos,y_pos), width=bmaj, height=bmin, angle=pa, edgecolor='white', facecolor='none', lw=2)
        # ax.add_patch(beam_ellipse)

    

        ax.set_xlabel('RA (J2000)', size=16)
        ax.set_ylabel('Dec (J2000)', size=16)
        ax.tick_params(axis="x", which="both", bottom=True, top=False)
        ax.tick_params(axis="y", which="both", right=False, left=True)
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label('Jy/beam', rotation=90, labelpad=12, size=18)
        cbar.formatter = ScalarFormatter()
        cbar.formatter.set_powerlimits((-3, 3))
        cbar.update_ticks()


        # # Set contour levels
        levels_low = np.asarray([4 * image_rms, 3 * image_rms])
        levels_g = np.geomspace(2.0 * peak_flux, 5 * image_rms, num_contours)
        neg_levels = np.asarray([-3])
        levels_neg = neg_levels * image_rms

        # Draw positive contours
        contour = ax.contour(image_data, levels=levels_g[::-1], colors=contour_palette, linewidths=1.2, extent=extent, alpha=1.0)

        # Draw low-level contours
        contour = ax.contour(image_data, levels=levels_low[::-1], colors='brown', linewidths=1.0, extent=extent, alpha=1.0)

        # Draw negative contours (if any)
        try:
            ax.contour(image_data, levels=levels_neg[::-1], colors='k', linewidths=1.0, extent=extent, alpha=1.0)
        except Exception as e:
            print(f"Error while plotting negative contours: {e}")
        plt.tight_layout()
        # plt.show()
        plt.savefig(imagename.replace('.fits','.pdf'),dpi=300)
        


# fits_filename = '/raid1/scratch/kelvinw/k2_18b/K2-18_split__1_image_dirty.image.tt0.fits'
fits_filename = '/raid1/scratch/kelvinw/k2_18b/23B-307.MJD60286.699115069445.K2-18_sci.C_band.cont.regcal.I.tt0.fits'

plot_image_with_beam(fits_filename)

# hdu = fits.open(fits_filename)
# image_data = hdu[0].data[0, 0, :, :] 
# header = hdu[0].header
# # image_data = np.squeeze(image_data)   # Removes any singleton dimensions

# g=image_data
# # Extract the data from the FITS file (typically in the first HDU)
# g = np.nan_to_num(g)
# # Set up plotting details
# CM = 'magma_r'  # Colormap
# num_contours = 5  # Number of contour levels
# std = mad_std(g)  # Estimate of the standard deviation
# extent = None  # Adjust if you have an extent for the image
# vmax = 0.9 * np.nanmax(g)
# vmin = 3 * std  # Setting minimum value based on MAD

# # Define the contour color palette
# contour_palette_ = ['#000000', '#444444', '#666666', '#EEEEEE', '#EEEEEE', '#FFFFFF']
# contour_palette = contour_palette_ if '_r' in CM else contour_palette_[::-1]

# # Set normalization for plotting
# norm0 = simple_norm(g, stretch='linear', max_percent=99.0)
# norm = simple_norm(g, stretch='sqrt', asinh_a=0.02, min_cut=vmin, max_cut=vmax)

# # Create the plot
# fig, ax = plt.subplots()

# # Plot the image with low-level normalization (transparent)
# im_plot = ax.imshow(g, origin='lower', cmap='gray', norm=norm0, alpha=0.5, extent=extent)

# # Plot the image with higher-level normalization (opaque)
# im_plot = ax.imshow(g, cmap=CM, origin='lower', alpha=1.0, extent=extent, norm=norm)

# # Set contour levels
# levels_low = np.asarray([4 * std, 3 * std])
# levels_g = np.geomspace(2.0 * g.max(), 5 * std, num_contours)
# neg_levels = np.asarray([-3])
# levels_neg = neg_levels * std

# # Draw positive contours
# contour = ax.contour(g, levels=levels_g[::-1], colors=contour_palette, linewidths=1.2, extent=extent, alpha=1.0)

# # Draw low-level contours
# contour = ax.contour(g, levels=levels_low[::-1], colors='brown', linewidths=1.0, extent=extent, alpha=1.0)

# # Draw negative contours (if any)
# try:
#     ax.contour(g, levels=levels_neg[::-1], colors='k', linewidths=1.0, extent=extent, alpha=1.0)
# except Exception as e:
#     print(f"Error while plotting negative contours: {e}")
    
# plt.show()

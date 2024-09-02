import os, subprocess, logging, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import casatools
from casatasks import *
from astropy.constants import c
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.ticker import FormatStrFormatter







def working_dir(filepath):
    """
    Gets the fitsfiles from the path and creates a working directory

    Parameters:
        filepath: path to directory where casa images of wsclean fitsfiles are stored
    """
    
    working_dir = os.path.join(filepath,'stacking_working_dir')

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    os.chdir(working_dir)



def convert_bmaj_to_pix(bmaj_deg, cdelt1):

    """
    Convert the BMAJ from degrees to pixels.
    
    Parameters:
        bmaj (float): Beam major axis in degrees.
        cdelt1 (float): The size of each pixel in degrees.
    
    Returns:
        float: Beam major axis in pixels.
    """

    ### NB: WSClean images with niter=0 have no beam information
    ### Convert the beam major axis to pixels
    
    return bmaj_deg / abs(cdelt1)


def load_fits_data(filepath):
    """
    Load all FITS data and header information.
    
    Parameters:
        fitsfiles (list): List of FITS file paths.
    
    Returns:
        tuple: A tuple containing a numpy array of all data and a list of headers.

    """

    casa_images = [os.path.join(filepath, filename) for filename in glob.glob(os.path.join(filepath, '*.image'))]
    if casa_images:
        print("CASA images exist, will attempt to convert to fits")
        for image in casa_images:
            exportfits(imagename=image,fitsimage=image.replace('.image','.fits'),overwrite=True)
            fitsfiles = glob.glob(os.path.join(filepath,'*.fits'))
    else:
        print("No CASA images found. Using WSClean fitsfiles")
        fitsfiles = glob.glob(os.path.join(filepath,'*-image.fits'))
        print(os.getcwd())


    all_data = []
    headers = []
    for fitsfile in fitsfiles:
        with fits.open(fitsfile) as hdul:
            all_data.append(hdul[0].data)
            headers.append(hdul[0].header)  # Store header info
    # print(all_data)
    return np.array(all_data), headers,fitsfiles

def apply_mask(image_data, crpix1, crpix2, bmaj_pixel):
    """
    Apply a mask to the image data based on beam size.
    
    Parameters:
        image_data (numpy array): 2D array of image data.
        crpix1 (float): The reference pixel along the x-axis.
        crpix2 (float): The reference pixel along the y-axis.
        bmaj_pixel (float): Beam major axis in pixels.
    
    Returns:
        numpy array: Masked image data with NaNs where the mask is applied.
    """
    ny, nx = image_data.shape
    y, x = np.ogrid[:ny, :nx]
    distance = np.sqrt((x - crpix1) ** 2 + (y - crpix2) ** 2)
    mask = distance >= bmaj_pixel
    masked_data = np.where(mask, image_data, np.nan)
    return mask

def calculate_rms_noise(all_data,headers,fitsfiles):
    """
    Calculate RMS noise for each stacked image

    Parameters:
        all_data (numpy array): Array of all FITS data.
        headers (list): List of headers for the FITS files.
        fitsfiles (list): List of FITS file paths.
    Returns:
        list: A list of RMS noise values for each stacked image.
    """
    rms_values = []
    for i in range(1,len(fitsfiles)+1):
        stacked_data = np.median(all_data[:i],axis=0)

        header = headers[0]  # Always use the header from the first original FITS file
        bmaj = header.get('BMAJ')
        cdelt1 = header.get('CDELT1')
        crpix1 = header.get('CRPIX1')
        crpix2 = header.get('CRPIX2')
        """YOU NEED TO PASS BMAJ AS AN ARG OR PUT IT HERE"""
        bmaj = 3.178236333446e-06 
        if bmaj is None or bmaj == 0.0:
            print("BMAJ is either not set or zero. Please enter a manual beam size.")
            # try:
            #     bmaj = float(input("Enter the BMAJ value (beam major axis in degrees): "))
            # except ValueError:
            #     print("Invalid input. Please enter a numerical value for BMAJ.")
            #     continue
        ## Convert to micro janskies
        image_data = stacked_data[0, 0, :, :]*1e6 
        # image_data = stacked_data[0, 0, :, :]
        bmaj_pixel = convert_bmaj_to_pix(bmaj, cdelt1)

        masked_data = apply_mask(image_data, crpix1, crpix2, bmaj_pixel)

        # Calculate the standard deviation (RMS) of the unmasked data
        unmasked_values = image_data[masked_data]
        std_dev = np.std(unmasked_values)
        rms_values.append(std_dev)

        plt.subplot(1, 2, 1)
        plt.title("Original Image Data")
        plt.imshow(image_data, origin='lower', cmap='viridis')
        plt.colorbar()
        
        plt.subplot(1, 2, 2)
        plt.title("Masked Image Data")
        plt.imshow(masked_data, origin='lower', cmap='viridis')
        plt.colorbar()
        
        plt.suptitle(f"Stacked {i} Images")
        plt.savefig(f'masked_image_{i}.pdf', dpi=300)
        plt.show()

    return rms_values

def plot_rms_noise(rms_values, fitsfiles):
    """
    Plot RMS noise versus the number of stacked images.
    
    Parameters:
        rms_values (list): List of RMS noise values.
        fitsfiles (list): List of FITS file paths.
    """
    initial_rms = rms_values[0]
    theoretical_rms = initial_rms / np.sqrt(np.arange(1, len(fitsfiles) + 1))

    plt.figure(figsize=(12, 8))
    plt.plot(range(1, len(fitsfiles) + 1), rms_values, marker='o', label='Observed RMS Noise')
    plt.plot(range(1, len(fitsfiles) + 1), theoretical_rms, linestyle='--', color='red', label='Theoretical RMS Noise ($1/\\sqrt{N}$)')
    plt.xlabel('Number of Stacked Images')
    plt.ylabel('RMS Noise ($\mu$Jy)')
    # plt.title('RMS Noise vs Number of Stacked Images')
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig('rms_noise.pdf',dpi=300)

def plot_fits(fitsname):
    """
    Plot FITS files using astropy.
    
    Parameters:
        fitsname (str): Name of the FITS file.
    """
    fitsfile = fits.open(fitsname)
    image_data = fitsfile[0].data[0, 0, :, :] * 1e6  # convert to micro janskies
    header = fitsfile[0].header
    axis1 = header['NAXIS1']

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    im = ax.imshow(image_data, cmap=plt.get_cmap('viridis'), extent=[-axis1 / 2, axis1 / 2, -axis1 / 2, axis1 / 2])
    cbar = plt.colorbar(im, ax=ax, orientation='vertical')
    cbar.set_label('$\mu$Jy', rotation=90, labelpad=-1)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    visible_ticks = {"top": False, "right": False}
    ax.tick_params(axis="x", which="both", **visible_ticks)
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')

    plt.savefig(fitsname.replace('.fits', '.pdf'))


if __name__ =="__main__":

    filepath = '/raid1/scratch/kelvinw/gv020_working_dir/gv020b_aoflagger_working_dir/pbcor_dir'

    bmaj = 3.178236333446e-06 


    working_dir(filepath)
    all_data, headers,fitsfiles = load_fits_data(filepath)
    rms_values = calculate_rms_noise(all_data,headers,fitsfiles)
    plot_rms_noise(rms_values,fitsfiles)

    # Save the last stacked data to a FITS file and plot it
    stacked_fits_filename = 'stacked_fitsfile.fits'
    fits.writeto(stacked_fits_filename, np.median(all_data, axis=0), overwrite=True)

    plot_fits(stacked_fits_filename)
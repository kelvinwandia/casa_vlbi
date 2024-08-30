import os, subprocess, logging
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import casatools
from casatasks import *
from astropy.constants import c
from astropy.io import fits
from astropy.wcs import WCS

filename ='/raid1/scratch/kelvinw/stacker/example/coordinates.list'
msname = '/raid1/scratch/kelvinw/stacker/example/testdata.ms'
stacking_dir = os.path.join(os.getcwd(),'stacking_dir')

if not os.path.exists(stacking_dir):
    os.mkdir(stacking_dir)

os.chdir(stacking_dir)
print(f"Current working directory: {os.getcwd()}")


def get_imaging_params():

    

    ms = casatools.ms()
    tb = casatools.table()
    ms.open(msname)
    max_uv = ms.getdata('uvdist')['uvdist'].max()
    ms.close()

    tb.open(msname+'/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    highest_freq = chan_freq.max()
    tb.close()

    # 3.6e6 converts the degrees to mas
    # 5 is the sampling

    cell_size = ((c.value/highest_freq)/max_uv)*(180./np.pi)*(3.6e6/5)
    cell_size = np.round(cell_size)
    print("The imaging cell size is:", cell_size)

    return cell_size

def plot_fits(fitsname):
    """
    Plots fitsfiles using astropy
    """
    fitsfile = fits.open(fitsname)
    image_data = fitsfile[0].data[0,0,:,:]
    ny, nx = image_data.shape
    x_center = nx // 2
    y_center = ny // 2
    x_new = np.arange(nx) - x_center
    y_new = np.arange(ny) - y_center

    fig, ax = plt.subplots()

    # image_plot = ax.imshow(image_data, origin='lower', 
    #                    extent=[x_new.min(), x_new.max(), y_new.min(), y_new.max()],cmap='viridis')
    image_plot = ax.imshow(image_data, origin='lower', 
                       extent=[-32, 32, -32, 32],cmap='viridis')
    cbar = plt.colorbar(image_plot,ax=ax,orientation='vertical')
    # ax.set_title(sources_to_image,fontsize=16)
    plt.savefig(fitsname.replace('.fits','.pdf'))

def get_im_stats(imagename):
    
    """
    Gets the statistics for either a 256x256 pix image and writes
    them to a logfile
    """


    rms=imstat(imagename=imagename,box='60,60,580,240')['rms'][0]  # for 640x640 px
    peak=imstat(imagename=imagename,box='300,300,340,340')['max'][0]
    print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))
    
    logfile = 'imstat.txt'
    casa_imstat = imstat(imagename)
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))

        txt_file.write(f"For {imagename}, the maximum pos for imstat is {casa_imstat['maxposf']}\n")


# coordinates = []
# with open(filename, 'r') as infile:
#     for line in infile:
#         parts = line.strip().split(',')
#         cleaned_line = ','.join(parts[:-1])
#         coordinates.append(cleaned_line)


# cell_size = get_imaging_params()
# cell_size = str(cell_size)+'mas'

# for coord in coordinates:
#     ra, dec = coord.split(',')
#     ra = float(ra)
#     dec = float(dec)
#     coord_conv = SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
#     hmsdms = coord_conv.to_string('hmsdms')
#     print(hmsdms)
    
    # phaseshifted_ms = hmsdms.replace(' ','')+'_phaseshifted.ms'
    # subprocess.run(['rm','-r',phaseshifted_ms])
    # os.system(f"rm -r {phaseshifted_ms}*")
    # phasecenter = 'J2000'+ ' '+ hmsdms
    # if not os.path.exists(phaseshifted_ms):
    #     print(f"======>>> Phaseshifting {msname} to {phasecenter}")
    #     phaseshift(
    #         vis = msname, outputvis = phaseshifted_ms, datacolumn='corrected',
    #         phasecenter = phasecenter
    #     )

    # imagename = phasecenter.replace(".ms","").replace(' ','')
    # print(f"======>>> Making {imagename}")
    # if not os.path.exists(imagename):
    #     os.system(f"rm -r {imagename}.*")
    #     logging.info(f"Making image {imagename}")
    #     tclean(
    #         vis = phaseshifted_ms, imagename=imagename, cell=cell_size, niter=0, deconvolver='clark',
    #         imsize=[256,256], phasecenter = phasecenter
    #     )
    #     exportfits(fitsimage=imagename+'.image',imagename=imagename+'.image',overwrite=True)
    #     get_im_stats(imagename+'.image')
    #     plot_fits(imagename+'.fits')


fitsfile = '/raid1/scratch/kelvinw/gv020_working_dir/gv020b_aoflagger_working_dir/selfcal_dir/J2139+1423_selfcal_loop_1-image.fits'

beam_major_axis = []
pixel_size_ra = []

def convert_bmaj_to_pix(bmaj_deg,cdelt1):
    return bmaj_deg/abs(cdelt1)

# def ra_dec_to_pixel(ra,dec,wcs):
#     # ra = np.array([ra])
#     # dec = np.array([dec])
#     # pixel_coords = wcs.world_to_pixel(ra,dec)
#     # return pixel_coords[0][0], pixel_coords[1][0]
#     pixel_values = wcs.all_world2pix(ra,dec,1)
#     # print(pixel_values)
#     return pixel_values

# def extract_ra_dec_wcs(header):
#     # Extracting RA and DEC specific WCS parameters
#     new_header = fits.Header()
    
#     # Copy RA and DEC specific information
#     new_header['CTYPE1'] = header['CTYPE1']
#     new_header['CTYPE2'] = header['CTYPE2']
#     new_header['CRVAL1'] = header['CRVAL1']
#     new_header['CRVAL2'] = header['CRVAL2']
#     new_header['CRPIX1'] = header['CRPIX1']
#     new_header['CRPIX2'] = header['CRPIX2']
#     new_header['CDELT1'] = header['CDELT1']
#     new_header['CDELT2'] = header['CDELT2']
    
#     # Ensure the header is in a correct FITS format
#     new_header['NAXIS'] = 2
#     new_header['NAXIS1'] = header['NAXIS1']
#     new_header['NAXIS2'] = header['NAXIS2']

#      # Create a new WCS object from the modified header
#     new_wcs = WCS(new_header)
    
#     # Verify the new WCS object
#     # print("New WCS Header:")
#     # print(new_header)

#     return new_wcs

with fits.open(fitsfile) as hdul:
    # Access the primary HDU (Header/Data Unit)
    header = hdul[0].header
    bmaj = header.get('BMAJ') #*3.6e6 
    cdelt1 = header.get('CDELT1')
    crpix1 = header.get('CRPIX1')
    crpix2 = header.get('CRPIX2')

    # print(f"BMAJ: {bmaj}, CDELT1: {cdelt1} ")

    # central_ra = header['CRVAL1']
    # central_dec = header['CRVAL2']

    # print(f"CRVAL1: {central_ra}, CRVAL2: {central_dec}")

    data = hdul[0].data
    image_data = data[0,0,:,:]
    ny, nx = image_data.shape
    x_center = nx // 2
    y_center = ny // 2

    bmaj_pixel = convert_bmaj_to_pix(bmaj,cdelt1)   

    # ra_dec_wcs = extract_ra_dec_wcs(header)
    # # print(ra_dec_wcs)
    # central_pixel = ra_dec_to_pixel(central_ra,central_dec,ra_dec_wcs)
    # cx = central_pixel[0]
    # cy = central_pixel[1]
    # # Print or use the new WCS object
    # print(cx,cy)

    ## Create a mask
    y,x = np.ogrid[:ny,:nx]
    distance = np.sqrt((x - crpix1)**2 + (y - crpix2)**2)
    mask = distance >= bmaj_pixel
    # print(mask)
    masked_data = np.where(mask, data[0, 0, :, :], np.nan)  # Mask the data, assuming first channel/layer
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.title("Mask")
    plt.imshow(mask, origin='lower', cmap='viridis')
    plt.colorbar()

    plt.subplot(1, 2, 2)
    plt.title("Masked Data")
    plt.imshow(masked_data, origin='lower', cmap='viridis')  # Displaying masked data
    plt.colorbar()
    plt.savefig('masked_fit.pdf',dpi=300)
    plt.show()
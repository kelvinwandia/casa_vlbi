from casatasks import *
import casatools
import numpy as np
from matplotlib import pyplot as plt
import os, time, subprocess
from astropy.io import fits
# from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
from astropy.wcs import WCS
from astropy.constants import c
import matplotlib
matplotlib.use('Agg') 
import argparse


def time_execution(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        
        if execution_time < 60:  # If execution time is less than a minute
            time_unit = "seconds"
            formatted_time = execution_time
        elif execution_time < 3600:  # If execution time is less than an hour
            time_unit = "minutes"
            formatted_time = execution_time / 60
        else:  # If execution time is an hour or more
            time_unit = "hours"
            formatted_time = execution_time / 3600
            
        print(f"======>>>EXECUTION TIME for {func.__name__}: {formatted_time:.2f} {time_unit}")
        return result
    return wrapper


def get_im_stats(imagename):
    
    """
    Gets the statistics for either a 256x256 pix image and writes
    them to a logfile
    """


    rms=imstat(imagename=imagename,box='51,7,247,76')['rms'][0]  
    peak=imstat(imagename=imagename,box='124,122,133,134')['max'][0]
    print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))

    logfile = 'imstat.txt'
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))
        

def plot_fits(fitsname):
    """
    Plots fitsfiles using astropy
    """
    fitsfile = fits.open(fitsname)
    image_data = fitsfile[0].data[0,0,:,:]
    w = WCS(fitsfile[0].header,naxis=2)
    header = fitsfile[0].header
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=w)
    axis1 = header['NAXIS1']

    im=ax.imshow(image_data, cmap=plt.get_cmap('viridis'),extent=[-axis1/2,axis1/2,-axis1/2,axis1/2]) 
    cbar = plt.colorbar(im,ax=ax,orientation='vertical')
    cbar.set_label('Jy/beam',rotation=90,labelpad=-1)
    cbar.formatter.set_powerlimits((0, 0))

    visible_ticks = {
   "top": False,
   "right": False
        }
    ax.tick_params(axis="x", which="both", **visible_ticks)
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')




def working_dir(filepath):
    """
    Gets the fitsfiles from the path and creates a working directory

    Parameters:
        filepath: path to directory where casa images of wsclean fitsfiles are stored
    """
    
    working_dir = os.path.join(filepath,'imaging_dir')

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    os.chdir(working_dir)


def get_imaging_params(vis):

    

    ms = casatools.ms()
    tb = casatools.table()
    ms.open(vis)
    max_uv = ms.getdata('uvdist')['uvdist'].max()
    ms.close()

    tb.open(vis+'/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    highest_freq = chan_freq.max()
    tb.close()

    # 3.6e6 converts the degrees to mas
    # 5 is the sampling

    cell_size = ((c.value/highest_freq)/max_uv)*(180./np.pi)*(3.6e6/5)
    cell_size = np.round(cell_size)
    print("The imaging cell size is:", cell_size)

    return cell_size


def create_mms(vis,outputvis,field):

    msmd = casatools.msmetadata()
    msmd.open(vis)
    scans = msmd.scansforfield(field=field)
    nspw = msmd.nspw()
    nscans = len(scans)

    print(f"You have {nspw} spectral windows and {nscans} scans")
    print(f"Partitioning measurement set")
    partition(vis=vis, outputvis = outputvis,  createmms=True, separationaxis='scan', numsubms = nscans)





@time_execution
def get_phasecenters(filename):
    """
    Reads the extracted stellar positions form the HR diagram and formats them
    to J2000 HMS DMS
    """

    """
    Image stars from a txt file
    """
    df = pd.read_csv(filename, delimiter='\t')
    phasecenter = []
    for index, row in df.iterrows():
        # ra_str = str(row['RA'])
        # dec_str = str(row['Dec'])
        c = SkyCoord(row['RA']*u.deg,row['Dec']*u.deg,frame='icrs')
        hmsdms = c.to_string('hmsdms')
        phasecenter.append('J2000'+' '+hmsdms)
        # print(phasecenter)
    return phasecenter

@time_execution
def make_map(phasecenter,vis,field):

    ## Partition before imaging
    # outputvis = field+'_partitioned.ms'
    # create_mms(vis,outputvis,field)

    for center in phasecenter[0:1]:
        imagename = center.replace("J2000 ","").replace(" ","")
        os.system(f"rm -r {imagename}.*")
        print(f"Making map {imagename}")

        tclean(
            vis=vis,imagename=imagename, cell='1mas',imsize=[320],niter=0, weighting='natural',deconvolver='clark', phasecenter=center, parallel=True,
            
        )
        fitsname = imagename+'.fits'
        exportfits(imagename=imagename+'.image',fitsimage=fitsname,overwrite=True)
        get_im_stats(imagename=imagename+'.image')
        plot_fits(fitsname)

def main(filename, vis, field, filepath):
    
    # Use the provided arguments
    working_dir(filepath)
    phasecenter = get_phasecenters(filename)
    make_map(phasecenter, vis, field=field)
    get_imaging_params(vis)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some files.")
    parser.add_argument("--filename", required=True, help="Path to the input file with phase center coordinates")
    parser.add_argument("--vis", required=True, help="Path to the visibility dataset (e.g., M15PSRC.ms)")
    parser.add_argument("--field", required=True, help="Field name to use in the map creation")
    parser.add_argument("--filepath", required=True, help="Path to the working directory")

    args = parser.parse_args()
    
    main(args.filename, args.vis, args.field, args.filepath)
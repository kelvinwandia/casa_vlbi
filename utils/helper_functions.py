
import time

import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from config_file import *
import casatools
from casatools import msmetadata, table
import datetime

msmd = casatools.msmetadata()
tb = casatools.table()

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
            
        log_message(f"======>>>EXECUTION TIME for {func.__name__}: {formatted_time:.2f} {time_unit}")
        return result
    return wrapper


def log_message(msg, level="INFO"):
    if level.upper() == "INFO":
        color = "\033[92m"  # Green
    elif level.upper() == "ERROR":
        color = "\033[91m"  # Red
    else:
        color = "\033[0m"   # Default

    casalog.post(msg, priority=level.upper())
    print(f"{color}{msg}\033[0m")

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
    log_message('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))
    
    logfile = 'imstat.txt'
    casa_imstat = imstat(imagename)
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))

        txt_file.write(f"For {imagename}, the maximum pos for imstat is {casa_imstat['maxposf']}\n")


def pybdsf(input_image):

    imagename = input_image
    fitsname = imagename

    img = bdsf.process_image(fitsname,adaptive_rms_box=True, thresh='hard',
                            thresh_isl=True, thresh_pix = detection_threshold, advanced_opts=True,
                            mean_map='map', rms_map =True, group_by_isl=True)
    # adaptive_rms_box=False, spline_rank=4, thresh='hard', thresh_isl=True, thresh_pix = detection_threshold
    # Write out island mask and FITS catalog -- for the large map
    fits_maskfile = imagename.replace('.fits','.maskfile.fits')
    catalog_file = imagename.replace('.fits','.cat')
    img.export_image(outfile=fits_maskfile,img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catalog_file, format='fits', clobber=True, catalog_type ='gaul')
    
    regionfile = imagename.replace('.fits','.casabox')
    ascii_file = imagename.replace('.fits','.ascii')
    rmsfile = imagename.replace('.fits','.rmsfile')

    img.write_catalog(outfile=regionfile,format='casabox',clobber=True,catalog_type='srl')
    img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
    img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

    return regionfile


def run_wsclean(wsclean_sif,command):

    """
    Runs wsclean commands 
    """
    container = wsclean_sif
    if os.path.exists(container):
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname(container)))

    command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + command

    try:
        log_message("Executing: %s", ' '.join(command_to_execute))
        process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        log_message("stdout: %s", stdout)
        log_message("stderr: %s", stderr)

        return_code = process.returncode
        if return_code == 0:
            log_message(f"Strategy executed successfully. Output:\n{stdout}")
        else:
            log_message(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")  

    except Exception as e:
        log_message(f"An error occurred: {e}")


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
    log_message("The imaging cell size is:", cell_size)

    return cell_size




def search_sbd_fringefit_soln(vis,field,refant,minsnr,interval,sbd_search):

    # current_directory = os.getcwd()
    # os.chdir(sbd_search)
    
    msmd = msmetadata()
    msmd.open(vis)
    scanlist = msmd.scansforfield(field)
    log_message(f"Scans for field {field}: {scanlist}")

    # Collect time ranges (in seconds since ref)
    scan_times = []
    for scan in scanlist:
        times = msmd.timesforscan(scan)
        scan_times.append(times)
    msmd.close()

    # Helper: convert MJD seconds to time string HH:MM:SS
    def mjdsec_to_str(t):
        dt = datetime.datetime.utcfromtimestamp(t)
        return dt.strftime('%H:%M:%S')

    # Build list of timeranges
    chunks = []
    for scan in scan_times:
        t0, t1 = scan[0], scan[-1]
        t_edges = np.arange(t0, t1, interval)
        for start, end in zip(t_edges[:-1], t_edges[1:]):
            chunks.append(f"{mjdsec_to_str(start)}~{mjdsec_to_str(end)}")

    log_message(f"Testing {len(chunks)} time ranges...")

    # Run fringefit and record median SNR for each timerange
    tb = table()
    results = []

    for tr in chunks:
        caltable = f"{sbd_search}/fringe_{tr.replace(':','').replace('~','_')}.cal"
        try:
            log_message(f"Running fringefit for timerange {tr}...")
            fringefit(vis=vis,
                    caltable=caltable,
                    field=field,
                    refant=refant,
                    timerange=tr,
                    minsnr=minsnr,
                    )   
            

            # Open caltable and measure SNR
            tb.open(caltable)
            snr = tb.getcol('SNR')
            tb.close()
            median_snr = np.nanmedian(snr)
            results.append((tr, median_snr))
        except Exception as e:
            log_message(f"Failed for {tr}: {e}")
            results.append((tr, 0.0))

    # Select best timerange
    if len(results) > 0:
        best_tr, best_snr = max(results, key=lambda x: x[1])
        log_message(f"\n Best timerange: {best_tr} (median SNR = {best_snr:.2f})")

        plotfile = f"{sbd_search}/fringe_snr_vs_timerange ({field}).png"
        plt.figure(figsize=(8,6))
        plt.plot([r[0] for r in results], [r[1] for r in results], marker='o')
        plt.xticks(rotation=90)
        plt.ylabel('Median SNR')
        plt.savefig(plotfile)
        
        plt.tight_layout()
        
        return best_tr

    else:
        log_message("No valid results found.")
        
    # os.chdir(current_directory)
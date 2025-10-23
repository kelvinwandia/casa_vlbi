
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

def get_number_of_threads():
    try:
        num_threads = os.cpu_count()
        if num_threads is None:
            print("Could not determine the number of threads.")
        else:
            print(f"Number of threads (logical processors) available: {num_threads}")
    except Exception as e:
        print(f"An error occurred while determining the number of threads: {e}")
    
    return int(num_threads)


def run_wsclean(wsclean_sif, command):
    """
    Runs wsclean commands
    """
    container = wsclean_sif
    if os.path.exists(container):
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname(container)))

    command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + command

    try:
        log_message(f"Executing: {' '.join(command_to_execute)}")
        process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        log_message(f"stdout: {stdout}")
        log_message(f"stderr: {stderr}")

        return_code = process.returncode
        if return_code == 0:
            log_message(f"Strategy executed successfully. Output:\n{stdout}")
        else:
            log_message(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")  

    except Exception as e:
        log_message(f"An error occurred: {e}")




def search_sbd_fringefit_soln(vis,field,refant,minsnr,interval,sbd_search):

    # current_directory = os.getcwd()
    # os.chdir(sbd_search)
    
    msmd = msmetadata()
    msmd.open(vis)
    scanlist = msmd.scansforfield(field)
    log_message(f"Scans for field {field}: {scanlist}")

    scan_times = []
    for scan in scanlist:
        times = msmd.timesforscan(scan)
        scan_times.append(times)
    msmd.close()

    def mjdsec_to_str(t):
        dt = datetime.datetime.utcfromtimestamp(t)
        return dt.strftime('%H:%M:%S')

    chunks = []
    for scan in scan_times:
        t0, t1 = scan[0], scan[-1]
        t_edges = np.arange(t0, t1, interval)
        for start, end in zip(t_edges[:-1], t_edges[1:]):
            chunks.append(f"{mjdsec_to_str(start)}~{mjdsec_to_str(end)}")

    log_message(f"Testing {len(chunks)} time ranges...")

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
    
    

def get_msinfo(msname):

    nchan = []
    msmd = casatools.msmetadata()
    msmd.open(msname)
    bandwidth = msmd.bandwidths()
    nspw = len(bandwidth)
    for spw in range(nspw):
        nchan.append(msmd.nchan(spw))
    msmd.close()
    # log_message(f"The measurement set contains {len(nspw)} spectral windows divided into {len(nchan)} channels")
    return nspw,nchan

def get_observing_band(msname: str) -> tuple:

    """
    Identify the frequency band of the data and return relevant frequency information.
    
    Parameters:
        vis (str): Path to the visibility data file.
    
    Returns:
        tuple: Band name, mean frequency, maximum frequency, and minimum frequency (in GHz).
    """

    band_name = None
    freq_ranges = {(1, 2): "L",(2, 4): "S",(4, 8): "C",(8, 12): 
                "X",(12, 18): "U",(18, 26.5): "K", (26.5, 
                    40): "A",(40, 50): "Q",
                            }
    msmd = casatools.msmetadata()
    msmd.open(msname)
    nspw = msmd.nspw()
    
    # Calculate mean frequency for each spectral window
    spws_freq = np.array([np.nanmean(msmd.chanfreqs(spw)) for spw in range(nspw)])
    msmd.done()
    
    # Calculate mean, max, and min frequencies across all spectral windows in GHz
    mean_freq = np.nanmean(spws_freq) * 1e-9
    max_freq = np.nanmax(spws_freq) * 1e-9
    min_freq = np.nanmin(spws_freq) * 1e-9
    
    # Identify band based on mean frequency
    for freq_range, band in freq_ranges.items():
        if freq_range[0] <= mean_freq <= freq_range[1]:
            band_name = band
            break

    log_message(f"Band: {band_name}, Mean Frequency: {mean_freq:.2f} GHz, "
        f"Min Frequency: {min_freq:.2f} GHz, Max Frequency: {max_freq:.2f} GHz")
    
    return band_name, mean_freq, max_freq, min_freq

def get_longest_baseline(msname:str) ->str:
    """
    Calculate the longest baseline in terms of wavelength (lambda).
    
    Parameters:
        vis (str): Path to the visibility data file.
    
    Returns:
        float: Longest baseline in units of wavelength.
    """

    from astropy.constants import c as LIGHT_SPEED

    # Open measurement set and retrieve uvw data
    ms = casatools.ms()
    ms.open(msname)
    # ms.selectinit(datadescid=0)
    # ms.selectinit()
    uvw = ms.getdata('uvw')['uvw']
    ms.close()
    
    # Compute baseline in meters
    uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
    longest_baseline_meters = np.nanmax(uvdist_meters)
    
    # Get frequency data
    band_name, mean_freq, max_freq, min_freq = get_observing_band(msname)
    frequency_hz = max_freq * 1e9
    wavelength_meters = LIGHT_SPEED.value / frequency_hz
    
    # Calculate longest baseline in terms of wavelength
    longest_baseline_lambda = longest_baseline_meters / wavelength_meters
    if longest_baseline_lambda >= 1e6:
        scaled_baseline = longest_baseline_lambda / 1e6
        unit = "Mλ"  # Mega wavelengths
    elif longest_baseline_lambda >= 1e3:
        scaled_baseline = longest_baseline_lambda / 1e3
        unit = "kλ"  # Kilo wavelengths
    else:
        scaled_baseline = longest_baseline_lambda
        unit = "λ"    # Wavelengths

    log_message(f"Longest Baseline: {scaled_baseline:.2f} {unit}")
    return longest_baseline_lambda

def get_imaging_cellsize(msname) -> str:
    """
    Calculate the cell size for imaging based on the longest baseline.

    Parameters:
    ----------
    msname : str
        The name of the measurement set.

    Returns:
    -------
    str
        The size of the imaging cell, either in arcseconds or milliarcseconds, depending on the value.
    """
    # Get the longest baseline in wavelength units, accounting for flags if needed
    longest_baseline_lambda = get_longest_baseline(msname)
    
    # Calculate cell size in arcseconds
    cell_float = (180.0 * 3600 / (np.pi * 5)) * (1.0 / longest_baseline_lambda)
    
    # Convert to mas if the value is very small (e.g., <1 arcsecond)
    if cell_float < 0.01:
        cell_float *= 1000  # convert to mas
        cell = f'{cell_float:.2f} mas'
    else:
        cell = f'{cell_float:.2f} arcsec'
    
    log_message(f"Imaging with a cell size of {cell}")
    return cell
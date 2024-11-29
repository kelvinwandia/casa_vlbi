from casatasks import tclean, imhead
import os, numpy
import casatools
from casatasks import *
import numpy as np
import logging

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

    logging.info(f"Band: {band_name}, Mean Frequency: {mean_freq:.2f} GHz, "
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
    # Open measurement set and retrieve uvw data
    ms = casatools.ms()

    ms.open(msname)
    ms.selectinit(datadescid=0)
    uvw = ms.getdata('uvw')['uvw']
    ms.close()
    
    # Compute baseline in meters
    uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
    longest_baseline_meters = np.nanmax(uvdist_meters)
    # print(longest_baseline_meters)
    
    # Get frequency data
    band_name, mean_freq, max_freq, min_freq = get_observing_band(msname)
    frequency_hz = max_freq * 1e9
    from astropy.constants import c as LIGHT_SPEED

    wavelength_meters =LIGHT_SPEED.value / frequency_hz
    
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

    print(f"Longest Baseline: {scaled_baseline:.2f} {unit}")
    return longest_baseline_lambda


vis = '/raid1/scratch/kelvinw/k2_18b/tapering/K2-18.ms'

def convert_lm_to_uv(theta_lm):
    """
    theta_lm: fwhm in the image domain, in units of arcsec
    """

    theta_uv = ( 4*numpy.log(2)/numpy.pi ) / ( (theta_lm / 3600) * numpy.pi/180.0)
    print("FWHM of %3.3f arcsec maps to a FWHM of %3.3e lambda"%(theta_lm,theta_uv))
    return theta_uv

def convert_uv_to_lm(theta_uv):
    """
    theta_uv : Full width at half maximum, in the uv domain, in units of lambda
    """
    theta_lm  = 3600 * ( 4*numpy.log(2)/numpy.pi ) / (theta_uv * numpy.pi/180.0)
    print("FWHM of %3.3f arcsec maps to a FWHM of %3.3e lambda"%(theta_lm,theta_uv))
    return theta_lm



def dispbeam(beam):
    """
    Print restoring beam info...
    """
    print("Restoring Beam : %3.4e %s  X  %3.4e %s ,  %3.4f %s"%(beam['major']['value'],
                beam['major']['unit'],
                beam['minor']['value'],
                beam['minor']['unit'],
                beam['positionangle']['value'],
                beam['positionangle']['unit'] ))


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
    # print(f"==========> {cell_float}")
    # Convert to mas if the value is very small (e.g., <1 arcsecond)
    if cell_float < 0.01:
        cell_float *= 1000  # convert to arcsec
        cell = f'{cell_float:.2f} mas'
    else:
        cell = f'{cell_float:.2f} arcsec'
    
    print(f"=============>>>Imaging with a cell size of {cell}")
    return cell



def run_im(imnames,uvtapers,weighting='briggs'):
    os.system('rm -rf uvt*')
    
    cell = get_imaging_cellsize(vis)
    for (imname,uvtaper) in zip(imnames,uvtapers):
        print("\n%s : uvtaper = %s"%(imname,uvtaper))
        tclean(vis=vis, spw='', imagename=imname,
                uvtaper=uvtaper,
                weighting=weighting,
                imsize=320, robust=0.5,
                cell=cell,niter=0, 
                #    restoringbeam=['1.4640e-02arcsec','1.4640e-02arcsec','0deg'],
                calcpsf=True, interactive=False, parallel = True,
                restoration=True)
        beam =imhead(imname+'.psf')['restoringbeam']
        dispbeam(beam)

imnames = ['uvt_orig' , 'uvt_taper_im' , 'uvt_taper_uv']




def calc_convolve(theta_orig, theta_taper):
    """
    Calculate the width of a Gaussian resulting from the convolution of two Gaussians.
    This calculation is only for a circular Gaussian.
    Units of inputs : arcsec.
    """
    arcsec_to_radians = (1/3600.0)*numpy.pi/180.0
    sigma_orig = arcsec_to_radians * theta_orig/numpy.sqrt(8*numpy.log(2.0))
    sigma_taper = arcsec_to_radians * theta_taper/numpy.sqrt(8*numpy.log(2.0))

    sigma_new = numpy.sqrt(sigma_orig**2 + sigma_taper**2)
    theta_new = sigma_new * numpy.sqrt(8*numpy.log(2.0)) / arcsec_to_radians

    print("Convolution of FWHMs of %3.4f arcsec and %3.4f arcsec \
          yields %3.4f arcsec"%(theta_orig, theta_taper, theta_new))
    
imtaper = 11.8
uvtapers=['' , '%3.2f arcsec '%(imtaper) , '%3.2elambda'%(convert_lm_to_uv(imtaper)/2.0)]
print("\nSettings for uvtaper in tclean : \n\
[ None,  FWHM in the image domain, HWHM in the uv-domain] ")
print(uvtapers)

# run_im(imnames,uvtapers,weighting='natural')

# convert_uv_to_lm()
# calc_convolve(20,20.896)


vis = '/raid1/scratch/kelvinw/k2_18b/selfcal_d_config/s_band_2/K2_18b_split_target.ms'
remove_ids = ['ea02', 'ea19', 'ea27', 'ea04', 'ea25', 'ea03', 'ea21']

def get_baseline_dist(vis):
     # Get the antenna names and offsets.

     msmd = casatools.msmetadata()

     msmd.open(vis)
     ids = msmd.antennasforscan(msmd.scansforintent("*OBSERVE_TARGET*")[0])
     names = msmd.antennanames(ids)
     # Filter out the unwanted antenna IDs
     ids = [id for id in ids if names[id] not in remove_ids]
     print(f"{len(ids)} antennas : {ids} used for imaging")
     offset = [msmd.antennaoffset(id) for id in ids]
     msmd.close()
     baselines=np.array([])
     for i in range(len(offset)):
        for j in range(i+1,len(offset)):
           baseline = numpy.sqrt((offset[i]["longitude offset"]['value'] -\
             offset[j]["longitude offset"]['value'])**2 + (offset[i]["latitude offset"]\
             ['value'] - offset[j]["latitude offset"]['value'])**2)
           
           baselines=np.append(baselines,np.array([baseline]))
     return baselines


def get_max_uvdist(vis,telescope='VLA'):

    all_baselines=np.array([])
    baselines=get_baseline_dist(vis)
    all_baselines=np.append(all_baselines,baselines)
    # print(all_baselines)
    max_baseline=np.max(all_baselines)
    min_baseline=np.min(all_baselines)
    if 'VLA' in telescope:
        baseline_5=numpy.percentile(all_baselines[all_baselines > 0.05*all_baselines.max()],5.0)
    else: # ALMA
        baseline_5=numpy.percentile(all_baselines,5.0)
    baseline_75=numpy.percentile(all_baselines,75.0)
    # baseline_median=numpy.percentile(all_baselines,50.0)
    # print(max_baseline)
    from astropy.constants import c
    freq = 3.0e9
    wavelength_meters = c.value/freq
    print(wavelength_meters*100)
        
    # Calculate longest baseline in terms of wavelength
    longest_baseline_lambda = max_baseline / wavelength_meters
    print(longest_baseline_lambda)
    if longest_baseline_lambda >= 1e6:
        scaled_baseline = longest_baseline_lambda / 1e6
        unit = "Mλ"  # Mega wavelengths
    elif longest_baseline_lambda >= 1e3:
        scaled_baseline = longest_baseline_lambda / 1e3
        unit = "kλ"  # Kilo wavelengths
    else:
        scaled_baseline = longest_baseline_lambda
        unit = "λ"    # Wavelengths
    print(freq)
    cell_size = (c.value/(freq*longest_baseline_lambda))*(180./np.pi)*(3.6e3/5)
    print(f"Longest Baseline: {scaled_baseline:.2f} {unit}")
    print(f"Cell size: {cell_size:.2f} arcseconds")

get_max_uvdist(vis)


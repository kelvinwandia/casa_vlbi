import casatools
from casaplotms import plotms
from casatasks import *
import numpy as np
import matplotlib.pyplot as plt

import glob, os, logging


msmd = casatools.msmetadata()
tb = casatools.table()
ms = casatools.ms()

cal_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir_d_config'
msname = '/raid1/scratch/kelvinw/k2_18b/working_dir_d_config/23B-307.sb44594812.eb44725045.60239.588568113424.ms'
splitvis = msname.replace('.ms','_calibrated.ms')
working_directory = '/raid1/scratch/kelvinw/k2_18b/pol_cal'


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
    
@staticmethod
def get_longest_baseline(msname:str) ->str:
    """
    Calculate the longest baseline in terms of wavelength (lambda).
    
    Parameters:
        vis (str): Path to the visibility data file.
    
    Returns:
        float: Longest baseline in units of wavelength.
    """
    # Open measurement set and retrieve uvw data
    ms.open(msname)
    ms.selectinit(datadescid=0)
    uvw = ms.getdata('uvw')['uvw']
    ms.close()
    
    # Compute baseline in meters
    uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
    longest_baseline_meters = np.nanmax(uvdist_meters)
    
    # Get frequency data
    band_name, mean_freq, max_freq, min_freq = get_observing_band(msname)
    frequency_hz = max_freq * 1e9
    
    from astropy.constants import c as LIGHT_SPEED
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

    logging.info(f"Longest Baseline: {scaled_baseline:.2f} {unit}")
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
    
    logging.info(f"Imaging with a cell size of {cell}")
    return cell


def set_working_dir(working_directory):

    """
    Creates a working dir if one does not exist
    """

    if not os.path.exists(working_directory):
        os.makedirs(working_directory)

    try:
        os.chdir(working_directory)
        logging.info(f"Changed working directory to {working_directory}")
    except Exception as e:
        logging.error(f"An error occurred while changing directory: {e}")
    
    logging.info(f"Setting logfile in working dir")

    # plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    # if not os.path.exists(plots_dir):
    #     os.makedirs(plots_dir)

def remove_parang_corrections(msname):

    """
    Applies the flagging state before the final applycal stage of the pipeline, then reapplies the calibration 
    with parang disabled 

    Removes the parallactic angle corrections that were applied by the pipeline
    """


    flagmanager(vis=msname, mode='restore',versionname='applycal_5')


    final_cals = glob.glob(os.path.join(cal_directory, '*.final*'))
    prior_cals = glob.glob(os.path.join(cal_directory,'*.hifv_priorcals*'))

    final_cals = [f for f in final_cals if 'finaldelayinitial' not in f and 'finalBPinitialgain' not in f]
    final_cals = final_cals + glob.glob(os.path.join(cal_directory, '*.averagephasegain*'))


    gaintables = prior_cals+final_cals
    gaintables = [gt for gt in gaintables if 'swpow' not in gt]

    interpolation = ['' for _ in gaintables]
    spwmap = [[] for _ in gaintables]
    calwt = [False for _ in gaintables]

    for i, gt in enumerate(gaintables):
        if 'finalBPcal' in gt:
            interpolation[i] = 'linear,linearflag'

    for gt, interp in zip(gaintables, interpolation):
        logging.info(f" Applying gaintable: {gt}, using interpolation: {interp}")


    applycal(
        vis=msname,
        antenna = '*&*',
        gaintable=gaintables,
        interp = interpolation,
        calwt = calwt,
        spwmap = spwmap,
        parang=False,
        flagbackup=False,
        applymode= 'calflagstrict'

    )

def flagdata_split(msname):

    """
    Target flagging, statwt and splitting
    """

    flagdata(vis=msname,
         mode='rflag', correlation='ABS_RR,LL', intent='*CALIBRATE*',
         datacolumn='corrected', ntime='scan', combinescans=False,
         extendflags=False, winsize=3, timedevscale=4.0, freqdevscale=4.0,
         action='apply', flagbackup=True, savepars=True)

    flagdata(vis=msname,
            mode='rflag', correlation='ABS_RR,LL', intent='*TARGET*',
            datacolumn='corrected', ntime='scan', combinescans=False,
            extendflags=False, winsize=3, timedevscale=4.0, freqdevscale=4.0,
            action='apply', flagbackup=True, savepars=True)

    statwt(vis=msname, minsamp=8, datacolumn='corrected')

    split(vis=msname,outputvis=splitvis,datacolumn='corrected',spw='')


def derive_pol_properties(msname,polarisation_calibrator):


    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    
    polarisation_calibrator = '3c286'
    os.system(f"wget -c 'https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/files/modes/{polarisation_calibrator}_2019' -O {polarisation_calibrator}_2019.dat")

    data = np.loadtxt(f"{polarisation_calibrator}_2019.dat")

    _ , mean_freq, _, _ = get_observing_band(msname)
    
    ### Stokes I

    def S(f,S,alpha,beta):
            return S*(f/round(mean_freq,1))**(alpha+beta*np.log10(f/round(mean_freq,1)))

    # Select the frequencies that are strictly greater than mean_freq
    higher_freq_data = data[data[:, 0] > mean_freq]
    sorted_higher_freq_data = higher_freq_data[higher_freq_data[:, 0].argsort()]
    max_freq_to_fit = sorted_higher_freq_data[2:] ## Go two frequencies highter
    max_freq_to_fit_index = np.where(data[:, 0] == max_freq_to_fit[0][0])[0][0]
    popt, pcov = curve_fit(S, data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,1])

    print(f"I@{round(mean_freq,1)}GHz {popt[0]} 'Jy'")
    print('alpha', popt[1])
    print('beta', popt[2])
    print( 'Covariance')
    print(pcov)

    #Clear any plots that may already exist
    plt.close()

    plt.plot(data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,1], 'ro', label='data')
    plt.plot(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), S(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), *popt), 'r-', label='fit')

    plt.title(f'{polarisation_calibrator}')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density (Jy)')
    plt.savefig('FluxvFreq.png')

    ### Polarisation fraction
    def PF(f,a,b,c,d):
        return a+b*((f-round(mean_freq,1))/round(mean_freq,1))+c*((f-round(mean_freq,1))/round(mean_freq,1))**2+d*((f-round(mean_freq,1))/round(mean_freq,1))**3

    # Fit 1 - 5 GHz data points
    popt, pcov = curve_fit(PF,data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,2])
    print("Polfrac Polynomial: ", popt)
    print("Covariance")
    print(pcov)

    #Clear any plots that may already exist
    plt.close()

    plt.plot(data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,2], 'ro', label='data')
    plt.plot(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), PF(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), *popt), 'r-', label='fit')

    plt.title(f'{polarisation_calibrator}')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Lin. Pol. Fraction')
    plt.savefig('LinPolFracvFreq.png')
    plt.show()


set_working_dir(working_directory)
# remove_parang_corrections(msname)
# flagdata_split(msname)

derive_pol_properties(msname = splitvis,polarisation_calibrator='3c286')



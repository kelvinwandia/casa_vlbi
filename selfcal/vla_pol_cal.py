import casatools
from casaplotms import plotms
from casatasks import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import glob, os, logging


msmd = casatools.msmetadata()
tb = casatools.table()
ms = casatools.ms()


logging.basicConfig(level=logging.INFO)

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


def download_polcal_table(polarisation_calibrator='3c286'):

    """ Download NRAO polarisation calibrator tables"""

    polcal_file = f'{polarisation_calibrator}_2019.dat'
    if not os.path.exists(polcal_file):
        logging.info(f"Downloading {polcal_file}")
        os.system(f"wget -c 'https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/files/modes/{polarisation_calibrator}_2019' -O {polcal_file}")
    else:
        logging.info(f"{polcal_file} exists")

    data = np.loadtxt(f"{polarisation_calibrator}_2019.dat")

    return data



def calculate_pol_parameters(msname,polarisation_calibrator,centre_freq=None):

    """ Calculate the spectral index of the pol angle calibrator """

    data = download_polcal_table(polarisation_calibrator)

    if centre_freq is not None: 
        logging.info(f"Using provided centre frequency: {centre_freq}")
        mean_freq = centre_freq
    elif msname:  
        logging.info(f"Getting the centre frequency from the measurement set")
        _, mean_freq, _, _ = get_observing_band(msname)
        mean_freq = mean_freq
    else:  
        raise ValueError("Either msname or centre_freq must be provided.")

    logging.info(f"Using frequency: {mean_freq}")
    
    spectral_index = {}
    pol_fraction = {}
    pol_angle = {}


    # Select the frequencies that are strictly greater than mean_freq
    higher_freq_data = data[data[:, 0] > mean_freq]
    sorted_higher_freq_data = higher_freq_data[higher_freq_data[:, 0].argsort()]
    max_freq_to_fit = sorted_higher_freq_data[2:] ## Go two frequencies higher
    max_freq_to_fit_index = np.where(data[:, 0] == max_freq_to_fit[0][0])[0][0]

    ##### Spectral Index ############
    def S(f,S,alpha,beta):
        return S*(f/round(mean_freq,1))**(alpha+beta*np.log10(f/round(mean_freq,1)))
    popt, pcov = curve_fit(S, data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,1])


    spectral_index['alpha'] = popt[1]
    spectral_index['beta'] = popt[2]

    logging.info(f'The mean flux density @{round(mean_freq, 1)}GHz is {popt[0]} Jy')
    logging.info(f"The source's spectral index is : {spectral_index}")
    # logging.info('Covariance:')
    # logging.info(pcov)
    plt.close()
    plt.plot(data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,1], 'ro', label='data')
    plt.plot(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), S(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), *popt), 'r-', label='fit')

    plt.title(f'{polarisation_calibrator}')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density (Jy)')
    plt.savefig('FluxvFreq.png')

    ######## Polarization fraction ##########

    def PF(f,a,b,c,d):
        return a+b*((f-round(mean_freq,1))/round(mean_freq,1))+c*((f-round(mean_freq,1))/round(mean_freq,1))**2+d*((f-round(mean_freq,1))/round(mean_freq,1))**3


    popt, pcov = curve_fit(PF,data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,2])
    pol_fraction[f'Polfrac polynomial'] = popt
    logging.info(f"The sources polarisation fraction is: {pol_fraction}")
    # logging.info("Covariance")
    # logging.info(pcov)

    plt.close()

    plt.plot(data[0:max_freq_to_fit_index,0], data[0:max_freq_to_fit_index,2], 'ro', label='data')
    plt.plot(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), PF(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), *popt), 'r-', label='fit')

    plt.title(f'{polarisation_calibrator}')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Lin. Pol. Fraction')
    plt.savefig('LinPolFracvFreq.png')

    ########## Polarization angle ###############

    ### Increase max freq to fit
    higher_freq_data = data[data[:, 0] > mean_freq]
    sorted_higher_freq_data = higher_freq_data[higher_freq_data[:, 0].argsort()]
    max_freq_to_fit = sorted_higher_freq_data[4:] ## Go four frequencies higher
    max_freq_to_fit_index = np.where(data[:, 0] == max_freq_to_fit[0][0])[0][0]

    data_greq_two_ghz = data[data[:,0]>2.0]
    sorted_data = data_greq_two_ghz[data_greq_two_ghz[:,0].argsort()[::1]] # ascending
    # print(sorted_data)
    # min_freq_to_fit = sorted_higher_freq_data
    logging.info(f"Using frequency range {data[max_freq_to_fit_index,0]} GHz to fit the pol angle")
    def PA(f,a,b,c,d,e):
        return a+b*((f-round(mean_freq,1))/round(mean_freq,1))+c*((f-round(mean_freq,1))/round(mean_freq,1))**2 + \
            d*((f-round(mean_freq,1))/round(mean_freq,1))**3+ e*((f-round(mean_freq,1))/round(mean_freq,1))**4
    

    # De-rotating the 1.8 GHz point
    data[2,3] = data[2,3]-np.pi
    popt, pcov = curve_fit(PA, data[2:max_freq_to_fit_index,0], data[2:max_freq_to_fit_index,3]) # Fit 2 - 19 GHz data points -- for C/X bands
    pol_angle['Polynomial angle'] = popt
    logging.info(f"The Polangle Polynomial is : {pol_angle}")
    # print("Covariance")
    # print(pcov)
    #Clear any plots that may already exist
    plt.close()

    plt.plot(data[2:max_freq_to_fit_index,0], data[2:max_freq_to_fit_index,3], 'ro', label='data')
    plt.plot(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), PA(np.arange(1,np.round(max_freq_to_fit[0][0]),0.1), *popt), 'r-', label='fit')

    plt.title('3c286')
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Lin. Pol. Angle (rad)')
    plt.savefig('LinPolAnglevFreq.png')


    return spectral_index, pol_fraction, pol_angle






# cal_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir_d_config'
msname = '/home/kelvin/Desktop/vla_calibrated/d_config/selfcal/observation.60292.44055671296/23B-307.sb44616223.eb44905127.60292.44054972222_target.ms'
# splitvis = msname.replace('.ms','_calibrated.ms')
working_directory = '/home/kelvin/Desktop/vla_calibrated/d_config/selfcal/polcal'



set_working_dir(working_directory)
calculate_pol_parameters(msname,polarisation_calibrator='3c286',centre_freq=6.0)

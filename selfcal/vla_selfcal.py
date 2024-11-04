
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
import casatools
import numpy as np
import matplotlib.pyplot as plt
import time 
from radio_beam import Beam
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

import logging
from typing import Callable, Any
msmd = casatools.msmetadata()
tb = casatools.table()
ms = casatools.ms()


msname = '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/23B-307.sb44594812.eb44725045.60239.588568113424/23B-307.sb44594812.eb44725045.60239.588568113424.ms'
# msname = '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/23B-307.sb44594812.eb44691528.60230.613198356485/23B-307.sb44594812.eb44691528.60230.613198356485.ms'


working_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir/23B-307.sb44594812.eb44725045.60239.588568113424' # D
# working_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir/23B-307.sb44594812.eb44691528.60230.613198356485' # A to D

target = 'K2-18'
vis = target+'.ms'
outlierfile = '/home/kelvin/Desktop/Synphly/selfcal/outlier_fields.txt'
basename = os.path.splitext(os.path.basename(vis))[0]

# target


# imaging 

nterms = 2
gridder = 'standard'
deconvolver = 'mtmfs'
weighting='briggs'
robust = 0.5
threshold = ['0.05mJy','0.02mJy','0.001mJy'] # in mJy
wprojplanes = 1
outlier_file = '/home/kelvin/Desktop/Synphly/selfcal/outlier_fields.txt'
pblimit = -0.1 # avoid 1,-1 or 0

imsize = [640,640] 
niter = [1000,10000,100000] # the number of iterations for each loop -- needs to be arbitrarily large

# selfcal
nloops = 3 # number of selfcal loops
calmode = ['p','p','ap']
gaintype= ['G','G','G']
solint = ['60s','30s','180s']
minsnr = [1,1,1]

# # pybdsf
detection_threshold = 5.0

# # final image and peeling
niter_final = 1000000
threshold_final = '0.05mJy'
wsclean_sif= '/raid1/scratch/kelvinw/singularity_containers/wsclean_working.simg'
abs_mem = 128 

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Utils():

    @staticmethod
    def set_working_dir(working_directory: str) -> None:
        """
        Sets the working directory to the specified path, creating it if it does not exist.

        Parameters:
        - working_directory (str): Path to the desired working directory.

        Raises:
        - OSError: If changing the directory fails.
        """

        if not os.path.exists(working_directory):
            logging.info(f"{working_directory} does not exist, making one")
            os.makedirs(working_directory)
        else:
            logging.info(f"Working directory {working_directory} already exists")

        logging.info(f"Changing cwd to {working_directory}")
        os.chdir(working_directory)


    @staticmethod
    def create_plots_directory() -> None:
        """
        Creates a 'plots' directory inside the current working directory.

        Raises:
        ------
        OSError
            If creating the plots directory fails.
        """
        plots_directory = os.path.join(os.getcwd(), 'plots')

        if not os.path.exists(plots_directory):
            try:
                os.makedirs(plots_directory)
                logging.info(f"Created plots directory at {plots_directory}.")
            except OSError as e:
                logging.error(f"Failed to create plots directory: {e}")
                raise
        else:
            logging.info(f"Plots directory {plots_directory} already exists.")

    @staticmethod
    def time_execution(func: Callable) ->Callable:

        """
        Decorator to measure and log the execution time for a function

        Parameters: 
        - func (Callable): The function to time
        Returns:
        - Callable: The wrapped function with timing

        """

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
                
            logging.info(f"EXECUTION TIME for {func.__name__}: {formatted_time:.2f} {time_unit}")
            return result
        return wrapper
    
    @staticmethod
    def run_wsclean(wsclean_sif, command):
        """
        Runs wsclean commands.
        """
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname(wsclean_sif)))
        command_to_execute = ['singularity', 'exec', '-B', singularity_bind, wsclean_sif] + command

        try:
            logging.info(f"Executing: {' '.join(command_to_execute)}")
            process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            logging.info(f"stdout: {stdout}")
            logging.error(f"stderr: {stderr}")

            return_code = process.returncode
            if return_code == 0:
                logging.info(f"Strategy executed successfully. Output:\n{stdout}")
            else:
                logging.error(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")

        except Exception as e:
            logging.exception("An error occurred during wsclean execution: %s", e)

    @staticmethod
    def pybdsf(input_image, detection_threshold):

        # Check if the input image is a FITS file; if not, add .fits
        if not input_image.endswith('.fits'):
            imagename = input_image
            fitsname =input_image + '.fits'
        else:
            # If it is already a FITS file, use it directly
            fitsname = input_image

        # Process the FITS image with pybdsf
        img = bdsf.process_image(fitsname, adaptive_rms_box=True, thresh='hard',
                                thresh_isl=True, thresh_pix=detection_threshold, 
                                advanced_opts=True, mean_map='map', rms_map=True, 
                                group_by_isl=True)

        # Write out island mask and FITS catalog
        img.export_image(outfile=imagename + '.maskfile.fits', img_type='island_mask', img_format='fits', clobber=True)
        img.write_catalog(outfile=imagename + '.cat', format='fits', clobber=True, catalog_type='gaul')

        regionfile = imagename + '.casabox'
        ascii_file = imagename + '.ascii'
        rmsfile = imagename + '.rmsfile'

        img.write_catalog(outfile=regionfile, format='casabox', clobber=True, catalog_type='srl')
        img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
        img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

        return regionfile

class MeasurementSetProcessor:

    """ Class to process easurement set 

    Attributes:
        avgtime (str): Time averaging interval for the split, specified as a string.
                       Defaults to '10s'.
        avgwidth (int): Channel averaging width for the split. Defaults to 4.

    """


    # Default values if split data called without timebin and width
    timebin: str = '10s'
    width: int = 4

    @staticmethod
    @Utils.time_execution
    def split_data(msname: str) -> tuple:

        """
        Splits the measurement set into individual fields

        Parameters:
        - msname (str): Name of the measurement set file

        Returns:
        - tuple: A tuple of the output measurement set names.
        """

        msmd = casatools.msmetadata()
        outputvis_list = []
        try:
            listobs(vis=msname, listfile='listobs.txt',overwrite=True)
            msmd.open(msname)
            field_names = msmd.fieldnames()
            logging.info(f"Field names: {field_names} found in {msname}")
            for field in field_names:
                outputvis = field+'.ms'
                if not os.path.exists(outputvis):
                    logging.info(f"Splitting {msname} to {outputvis}")
                    split(vis=msname,outputvis=outputvis,datacolumn='corrected',timebin=MeasurementSetProcessor.timebin,
                        width=MeasurementSetProcessor.width,field=field)
                    logging.info(f"Finished splitting")
                else:
                    logging.info(f"Split measurement set {outputvis} exists")
                outputvis_list.append(outputvis)
        except Exception as e:
            logging.error(f"Error in split_data: {e}")

        finally:
            msmd.close()

        return tuple(outputvis_list)

class MeasurementSetInfo:

    """ Class to get basic measurement set information """

    from astropy.constants import c as LIGHT_SPEED

    @staticmethod
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
        ms.open(vis)
        ms.selectinit(datadescid=0)
        uvw = ms.getdata('uvw')['uvw']
        ms.close()
        
        # Compute baseline in meters
        uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
        longest_baseline_meters = np.nanmax(uvdist_meters)
        
        # Get frequency data
        band_name, mean_freq, max_freq, min_freq = MeasurementSetInfo.get_observing_band(msname)
        frequency_hz = max_freq * 1e9
        wavelength_meters = MeasurementSetInfo.LIGHT_SPEED.value / frequency_hz
        
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
    
    @staticmethod
    def get_imaging_cellsize(msname) -> str:
        """
        Calculate the cell size for imaging based on the longest baseline.
        
        Returns:
        -------
        str
            The size of the imaging cell in arcseconds.
        """
        longest_baseline_lambda = MeasurementSetInfo.get_longest_baseline(msname)
        cell_float = (180.0 * 3600 / (np.pi * 5)) * (1.0 / longest_baseline_lambda)
        cell = f'{cell_float:.2f} arcsec'
        logging.info(f"Imaging with a cell of size {cell}")
        return cell
    
    
    @staticmethod
    # @Utils.time_execution
    def find_refant(msname_tuple:tuple, field:str):

        """
        
        Find the best reference antenna by calculation the SNR by making a mock calibration table.
        The function will find the tuple with the same name as the field and use that

        Parameters
            msname_tuple (tuple): A tuple of the measuremet set files from split_data

            field (str): The field to self calibrate

        """

        vis = None
        for msname in msname_tuple:
            if field in msname:  # Check if the field name is in the measurement set name
                vis = msname
                break
        
        if vis is None:
            logging.error(f"No measurement set found for the field: {field}")
            return None

        tablename = field + '.refant'
        if not os.path.exists(tablename):
            gaincal(vis=vis, caltable=tablename, field=field, refantmode='flex',
                     solint='inf', minblperant=3, gaintype='G', calmode='p')
            logging.info(f"Gain calibration completed for field: {field}")

        # Read solutions (phases):
        tb = casatools.table()
        
        try:
            tb.open(tablename + '/ANTENNA')
            antenna_names = tb.getcol('NAME')
            tb.close()
            tb.open(tablename)
            antenna_ids = tb.getcol('ANTENNA1')
            flags = tb.getcol('FLAG')
            phases = np.angle(tb.getcol('CPARAM'))
            snrs = tb.getcol('SNR')
            tb.close()
        except Exception as e:
            logging.error(f"Error reading calibration table {tablename}: {e}")
            return None

        # Analyze number of good solutions:
        good_frac = []
        good_snrs = []
        for i, ant_id in enumerate(np.unique(antenna_ids)):
            cond = antenna_ids == ant_id
            f = flags[0, 0, :][cond]
            snr = snrs[0, 0, :][cond]
            frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
            snr_mean = np.nanmean(snr[~f]) if len(snr[~f]) > 0 else np.nan
            good_frac.append(frac)
            good_snrs.append(snr_mean)

        sort_idx = np.argsort(good_frac)[::-1]
        truncation = 3 # truncate the printed antenna to save on logfile size
        logging.info(f'Antennas (showing only {truncation}) sorted by percentage of good solutions:')
        for i in sort_idx[:truncation]:
            logging.info(f"{antenna_names[i]:3}: {good_frac[i]:4.1f}%, <SNR> = {good_snrs[i]:4.1f}")

        if good_frac[sort_idx[0]] < 90:
            logging.warning('Small fraction of good solutions with selected reference antenna!')
            logging.warning('Please inspect antennas to select the optimal reference antenna.')
            logging.warning('Consider using refantmode="flex" in gaincal.')

        pref_ant = antenna_names[sort_idx]
        pref_ant_list = list(pref_ant)
        logging.info(f"The following antennas will be used as the reference antennas: {', '.join(pref_ant_list[:truncation])}")
        return pref_ant_list




class tclean_Imager:

    """
    A class for performing imaging on measurement sets using the tclean algorithm.

    Methods:
    -------
    get_imaging_cellsize() -> str:
        Calculates and returns the cell size for imaging in arcseconds.
    
    imager() -> None:
        Performs the imaging process using the tclean algorithm.
    """

    def __init__(self, msname: str, field:str='', gridder: str = 'standard', deconvolver: str = 'multiscale', 
                 weighting: str = 'natural', robust: float = 0, nterms: int = 1, imsize: int = 640,
                 niter: int = 0, threshold: str = None, wprojplanes: int = 1, mask: str = '',
                 imagename: str = f"{msname.replace('.ms','_image')}",
                 usemask: str = 'user', pblimit: float = 0.1,phasecenter: str = '', overwrite:bool = False, 
                 use_pybdsf:bool = True, pybdsf_threshold: int = 5):
        """
        Initializes the tclean_Imager instance with specified parameters.

        Parameters:
        ----------
        msname : str
            The name of the measurement set file.
        field: str,optional
            The name of the field to image. Default is an empty string
        gridder : str, optional
            The type of gridder to use (default is 'standard').
        deconvolver : str, optional
            The type of deconvolver to use (default is 'multiscale').
        weighting : str, optional
            The weighting scheme to use (default is 'natural').
        imagename: str, optional
            Name of the images to be created
        robust : float, optional
            The robust parameter for weighting (default is 0.5).
        nterms : int, optional
            The number of Taylor terms to use (default is 1).
        imsize : int, optional
            The size of the image (default is 640).
        niter : int, optional
            The number of iterations for the tclean process (default is 0).
        threshold : str, optional
            The threshold for stopping the deconvolution (default is None).
        wprojplanes : int, optional
            The number of w-projection planes (default is 1).
        mask : str, optional
            The mask to apply during the imaging process (default is empty string).
        usemask : str, optional
            The mask usage strategy (default is 'user').
        pblimit : float, optional
            The limit for the primary beam (default is 0.1).
        phasecenter: str, optional
            The field phase center (default is an empty string).
        overwrite: bool, optional
            If True, deletes existing images and creates new ones (default is False).
        use_pybdsf: bool, optional
            Calls pybdsf and uses it for masking
        pybdsf_threshold: int, optional
            Masking threshold

        """
        self.msname = msname  
        self.field = field
        self.gridder = gridder
        self.imagename = imagename
        self.deconvolver = deconvolver
        self.weighting = weighting
        self.robust = robust
        self.nterms = nterms
        self.imsize = imsize  
        self.niter = niter
        self.threshold = threshold
        self.wprojplanes = wprojplanes
        self.mask = mask
        self.usemask = usemask
        self.pblimit = pblimit
        self.phasecenter = phasecenter
        self.overwrite = overwrite
        self.use_pybdsf = use_pybdsf
        self.pybdsf_threshold = pybdsf_threshold


    @Utils.time_execution
    def imager(self) -> None:
        """
        Perform the imaging process using the tclean algorithm.

        This method calculates the cell size for imaging, determines the output image name based on the number of iterations,
        and calls the tclean function to generate the dirty image.

        Returns:
        -------
        None
        """
        cell = MeasurementSetInfo.get_imaging_cellsize(msname)
        self.imagename = self.imagename+'_dirty' if self.niter == 0 else self.imagename

    
        matching_files = glob.glob(f"{self.imagename}.*")

        # Check if any matching files exist
        if matching_files:
            logging.info(f"Image exists: {matching_files}")

            # If the image exists and overwrite is requested, proceed with deletion
            if self.overwrite:
                try:
                    for file in matching_files:
                        command = f"rm -rf {file}" 
                        subprocess.run(command, shell=True, check=True)  
                    # logging.info(f"Overwriting requested. Deleted existing image files: {matching_files}")
                    logging.info(f"Overwriting requested. Deleted existing image files")

                except OSError as e:
                    logging.error(f"Error deleting existing images: {e}")
                except Exception as e:
                    logging.error(f"Unexpected error: {e}")
        
            
        
        logging.info(f"Imaging {self.msname}")

        tclean(
            vis=self.msname,
            imagename=self.imagename,
            imsize=[self.imsize, self.imsize], 
            cell=cell,
            gridder=self.gridder,
            deconvolver=self.deconvolver,
            weighting=self.weighting,
            robust=self.robust,
            niter=self.niter,
            nterms=self.nterms,
            threshold=self.threshold,
            wprojplanes=self.wprojplanes,
            mask=self.mask,
            usemask=self.usemask,
            pblimit=self.pblimit,
            field = self.field,
            interactive=False,
        )

        logging.info(f"Finished imaging {self.msname}, created image: {self.imagename}")

        # if self.deconvolver == 'mtmfs':
        #     image_ext = '.image.tt0'
        # else:
        #     image_ext = '.image'

        # exportfits(imagename=self.imagename+image_ext,fitsimage=self.imagename+ image_ext+'.fits',overwrite=True)
        # # try:
        # #     logging.info(f"Running pybdsf on {self.imagename}...")
        # #     if self.use_pybdsf:
        # #         Utils.pybdsf(self.imagename+image_ext,self.pybdsf_threshold)
        # #         logging.info(f"Successfully ran pybdsf on {imagename}.")
        # #     else:
        # #         logging.info(f"Masking using PYBDSF not requested")
        # # except Exception as e:
        # #     logging.error(f"Failed to run pybdsf on {self.imagename}: {e}")


class wsclean_Imager:
    """
    A class for performing imaging on measurement sets using the WSClean algorithm.

    Methods:
    -------
    get_imaging_cellsize() -> str:
        Calculates and returns the cell size for imaging in arcseconds.
    
    imager() -> None:
        Performs the imaging process using the WSClean algorithm.

    """

    def __init__(self, msname: str, imsize: int = 640, niter: int = 0, threshold: str = None, deconvolution:str = None,
                 overwrite: bool = False, use_pybdsf: bool = True, pybdsf_threshold: int = 5, mgain: float = 0.8):
        """
        Initializes the wsclean_Imager instance with specified parameters.

        Parameters:
        ----------
        msname : str
            The name of the measurement set file.
        imsize : int, optional
            The size of the image (default is 640).
        niter : int, optional
            The number of iterations for the WSClean process (default is 0).
        threshold : str, optional
            The threshold for stopping the deconvolution (default is None).
        overwrite: bool, optional
            If True, deletes existing images and creates new ones (default is False).
        use_pybdsf: bool, optional
            Calls pybdsf and uses it for masking.
        pybdsf_threshold: int, optional
            Masking threshold (default is 5).
        mgain : float, optional
            The gain for minor cycles in WSClean (default is 0.8).
        """
        self.msname = msname  
        self.imsize = imsize  
        self.niter = niter
        self.threshold = threshold
        self.deconvolution = deconvolution
        self.overwrite = overwrite
        self.use_pybdsf = use_pybdsf
        self.pybdsf_threshold = pybdsf_threshold
        self.mgain = mgain

    
    
    @Utils.time_execution
    def imager(self) -> None:
        """
        Perform the imaging process using the WSClean algorithm.

        Returns:
        -------
        None
        """
        cell = MeasurementSetInfo.get_imaging_cellsize(msname)
        imagename = self.msname.replace('.ms', '_image')

        
        command = [
            "wsclean",
            "-log-time",
            "-size", str(self.imsize), str(self.imsize),
            "-reorder",      
            "-name", imagename,
            "-scale", cell,
            "-mgain", str(self.mgain),
            "-niter", str(self.niter),
        ]

        if self.deconvolution == "multiscale":
            command.append("-multiscale")

        command.append(self.msname)
        
        # Run the WSClean command
        try:
            Utils.run_wsclean(wsclean_sif,command)
            logging.info(f"Finished imaging {self.msname}, created image: {imagename}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error during WSClean imaging: {e}")

        # If using PYBDSF for masking, call it here
        # if self.use_pybdsf:
        #     try:
        #         logging.info(f"Running pybdsf on {imagename}...")
        #         Utils.pybdsf(imagename + '-image.fits', self.pybdsf_threshold)
        #         logging.info(f"Successfully ran pybdsf on {imagename}.")
        #     except Exception as e:
        #         logging.error(f"Failed to run pybdsf on {imagename}: {e}")
        # else:
        #     logging.info("Masking using PYBDSF not requested.")




class SelfCalibration(tclean_Imager):

    def __init__(self, msname, nloops, thresholds, **kwargs):
        """
        Initializes SelfCalibration, inheriting from tclean_Imager.

        Parameters:
        - msname: Name of the measurement set.
        - nloops: Number of self-calibration loops.
        - thresholds: List of threshold values for each loop.
        - kwargs: Additional parameters for tclean_Imager.
        """
        super().__init__(msname=msname, **kwargs)  # Call the parent constructor
        self.nloops = nloops
        self.thresholds = thresholds  

    def selfcal(self) -> None:
        """
        Perform self-calibration on the measurement set.
        """
        Utils.create_plots_directory()
        regionfile = self.msname + '_first_masking.casabox'
        # cell = MeasurementSetInfo.get_imaging_cellsize(self.msname)

        logging.info("Deleting model column before selfcal")
        delmod(vis=self.msname, otf=True)

        for selfcal_loop in range(self.nloops):
            caltable = f'{self.msname.replace(".ms", "_caltable_loop_")}_{selfcal_loop}.gcal'
            prev_caltables = sorted(glob.glob('*.gcal'))
            if len(prev_caltables) > 0 and calmode[selfcal_loop] != '':
                applycal(vis=self.msname, gaintable=prev_caltables, parang=False)

            # Set niter to 0 for the first loop to create a dirty map and run PYBSDF
            if selfcal_loop == 0:
                """
                Set niter to 0 and run PYBDSF
                """
                self.niter = 0  # Set niter for dirty map
                imagename_dirty =  f'{self.msname.replace(".ms", "_dirty")}'

                if self.deconvolver == 'mtmfs':
                    image_ext = '.image.tt0'
                else:
                    image_ext = '.image'
                exportfits(imagename=imagename_dirty+image_ext,fitsimage=imagename_dirty+ image_ext+'.fits',overwrite=True)

                try:
                    logging.info(f"Running pybdsf on {imagename_dirty}...")
                    if self.use_pybdsf:
                        Utils.pybdsf(imagename_dirty+image_ext,self.pybdsf_threshold)
                        logging.info(f"Successfully ran pybdsf on {imagename_dirty}.")
                    else:
                        logging.info(f"Masking using PYBDSF not requested")
                except Exception as e:
                    logging.error(f"Failed to run pybdsf on {imagename_dirty}: {e}")

            else:
                """
                Set niter to what is set in the imager instance
                also set PYBDSF to false
                """
                self.niter = self.niter 
                self.use_pybdsf = False
                imagename = f'{self.msname.replace(".ms", "_selfcal_loop")}_{selfcal_loop}'


            # Initialize the imager_instance with required parameters
            imager_instance = tclean_Imager(
                msname=self.msname,
                imagename = imagename,
                nterms = self.nterms,
                imsize=self.imsize,
                niter=self.niter,
                deconvolver=self.deconvolver,
                use_pybdsf=self.use_pybdsf,
                pybdsf_threshold=self.pybdsf_threshold,
                threshold=self.thresholds[selfcal_loop],
                mask=regionfile ,
                overwrite = self.overwrite
            )
            
            logging.info("Adding modelcolumn to data")
            # model images from the MTMFS images,
            if self.deconvolver == 'mtmfs':
                ft(vis = msname, model=[imagename+'.model.tt0',imagename+'.model.tt1'], nterms=2,usescratch=True)
            else:
                ft(vis = msname, model=imagename+'.model', usescratch=True)

            logging.info("Plotting the model column")
            plotms(
                vis=msname, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
            )

        #         gaincal( vis =vis, caltable = caltable, refant = refant, solint = solint[selfcal_loop],
        #                 gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
        #                 calmode = calmode[selfcal_loop], append=False, parang=False
        #                 )
        #         coloraxis = ['corr','spw']
        #         for color in coloraxis:
        #             if calmode[selfcal_loop] =='p':
        #                 plotms(
        #                     vis = caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
        #                     iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
        #                     plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750,
        #                 )
        #             else:
        #                 plotms(
        #                         vis = caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
        #                         iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
        #                         plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750
        #                     )

        #         if selfcal_loop == nloops-1:
        #             prev_caltables = sorted(glob.glob('*.gcal'))
        #             print("Applying the caltable derived from last gaincal iteration")
        #             applycal(vis=vis, gaintable = prev_caltables, parang=False )
            
            # ### Get the last imagename from the loop and generate a final mask
            
#         imagename = basename +f'_{nloops-1}'+'.final'
#         ##  tclean here to make the final image
#         print("Make final image with all selfcal corrections applied")
#         tclean(
#             vis = vis, imagename = imagename, imsize=imsize, cell=cell, gridder=gridder,
#             wprojplanes = wprojplanes, deconvolver = deconvolver, weighting = weighting,
#             robust = robust, niter=niter_final, threshold = threshold_final, nterms=nterms,
#             pblimit=pblimit, interactive=False, usemask = 'user', mask=regionfile,
#         )
        
        

Utils.set_working_dir(working_directory)
MeasurementSetProcessor.timebin='6s'
MeasurementSetProcessor.width = 4
msname_tuple = MeasurementSetProcessor.split_data(msname)
band_info = MeasurementSetInfo.get_observing_band(msname)
longest_baseline = MeasurementSetInfo.get_longest_baseline(msname)
refant = MeasurementSetInfo.find_refant(msname_tuple,field=target)

def main():
    # Define your parameters here
    msname = msname_tuple[0]  #
    print(msname)
    """First loop is a dirty map for masking """
    nloops = 4
    thresholds = ['', '0.05mJy', '0.01mJy', '0.005mJy']  # Example thresholds for each loop

    # Create an instance of SelfCalibration for the first loop (dirty map)
    self_calibration_instance = SelfCalibration(
        msname=msname,
        nloops=nloops,
        thresholds=thresholds,
        imsize=640,
        niter=1000,  
        nterms = 2,
        deconvolver='mtmfs',
        weighting='briggs',
        robust='0.5',
        use_pybdsf=True,
        pybdsf_threshold=5,
        overwrite = True
    )

    self_calibration_instance.selfcal()

if __name__ == "__main__":
    main()


#     # def applycal_target():

#     #     """
#     #     Applies cal to the target field
#     #     """
#     #     prev_caltables = sorted(glob.glob('*.gcal'))
#     #     applycal(
#     #         vis = vis_tocal, gaintable = prev_caltables, parang=False
#     #     )



#     # def peeling():

#     #     """
#     #     Subtract the bright sources in the field so you are left with the emission from the star

#     #     Use pybsdf casa region file -- you can change the threshold to only select the bright sources you want

#     #     Run pybsdf on the final image
        
#     #     """

#     #     # Using the final imagename 

#     #     imagename = basename +f'_{nloops-1}'+'.final'
#     #     # region_to_peel = pybdsf(input_image=imagename+'.image.tt0')

#     #     ## NB: You had masked the bright source from the dirty map -- so you can just subtract
#     #     uvsub(vis=vis)

#     #     ## Make a final map without the sources
#     #     cell = get_imaging_cellsize()
#     #     peeled_map = basename+'_peeled_map'
#     #     if not os.path.exists(peeled_map):
#     #         print(f"Making {peeled_map}")
#     #         tclean(
#     #             vis = vis, imagename=peeled_map, imsize=imsize, cell=cell,
#     #             gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
#     #             weighting = weighting, robust = robust, niter=0, # threshold = '0.5mJy',
#     #             nterms = nterms, pblimit = pblimit
#     #         )






# for msname in msname_tuple:
#     # logging.info(f"Processing measurement set {msname}")
#     try:
#         imager = tclean_Imager(
#             msname=msname_tuple[0],
#             gridder='standard',  
#             deconvolver='mtmfs', 
#             weighting='briggs',  
#             robust=0.5,  
#             nterms=2,  
#             imsize=640,  
#             niter=0,  
#             threshold='0.1mJy',  
#             wprojplanes=1,  
#             pblimit=0.1,  
#             phasecenter='',  
#             overwrite=True, 
#             use_pybdsf = False,
#             pybdsf_threshold=5 
#         )
#         imager.imager()  # Call the imaging process

    

#     except Exception as e:
        
#         print(f"An error occurred while processing {msname}: {e}")
#         logging.error(f"An error occurred while processing {msname}: {e}")

# logging.info("Imaging process completed for all measurement sets.")









#     #  def get_im_stats(imagename):
        
#     #     """
#     #     Gets the statistics for either a 256x256 pix image and writes
#     #     them to a logfile
#     #     """


#     #     rms=imstat(imagename=imagename,box='60,60,580,240')['rms'][0]  # for 640x640 px
#     #     peak=imstat(imagename=imagename,box='300,300,340,340')['max'][0]
#     #     print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
#     #                 (imagename, peak*1e3, rms*1e3, peak/rms))
        
#     #     logfile = 'imstat.txt'
#     #     casa_imstat = imstat(imagename)
#     #     with open(logfile,"a") as txt_file:
#     #         txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
#     #                     (imagename, peak*1e3, rms*1e3, peak/rms))

#     #         txt_file.write(f"For {imagename}, the maximum pos for imstat is {casa_imstat['maxposf']}\n")





# class PlottingRoutines():

#     def plot_fits(imagename,fig_width =8,fig_height=6):
#         """
#         Plots fitsfiles using astropy
#         """

#         fitsname = imagename+'.fits'
#         exportfits(imagename=imagename+'.image.tt0',fitsimage=fitsname,overwrite=True)

#         plt.figure(figsize=(fig_width, fig_height))  
#         hdul = fits.open(fitsname)
#         w = WCS(hdul[0].header, naxis=2)
#         w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
#         ax = plt.subplot(projection=w)

#         # Disable automatic labelling
#         ax.coords[0].set_auto_axislabel(True) 
#         ax.coords[1].set_auto_axislabel(True) 

#         # Extract image data
        
#         image_data = hdul[0].data[0,0,:,:]

#         # Display the image with a color map
#         im = ax.imshow(image_data, cmap=plt.get_cmap('cividis'))

#         # Automatically set pixel scale based on WCS header information
#         pixscale = abs(w.wcs.cdelt[0]) * u.deg.to(u.arcsec) * u.arcsec  # Convert from degrees to arcseconds

#         # Define and add the beam ellipse
#         array_beam = imhead(imagename+'.psf.tt0')['restoringbeam']
#         major_axis = array_beam['major']['value']
#         minor_axis = array_beam['minor']['value']
#         pos_angle = array_beam['positionangle']['value']

#         major_axis = major_axis*u.arcsec  # Convert to arcseconds if needed 
#         minor_axis = minor_axis*u.arcsec  # Convert to arcseconds if needed
#         pos_angle = pos_angle*u.deg

#         my_beam = Beam(major_axis, minor_axis, pos_angle)
#         ycen_pix, xcen_pix = 15, 15
#         ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
#         _ = ax.add_artist(ellipse_artist)

#         ax.set_xlabel('RA (J2000)',size=14)
#         ax.set_ylabel('Dec (J2000)',size=14)  

#         ax.tick_params(axis = "x", which = "both", bottom = True, top = False)
#         ax.tick_params(axis = "y", which = "both", right = False, left = True)

#         # ra = ax.coords[0]
#         # dec = ax.coords[1]

#         # ra.set_ticklabel(size=12)
#         # dec.set_ticklabel(size=12)

#         cbar = plt.colorbar(im,extend='both')
#         cbar.ax.tick_params(labelsize=16)
#         cbar.set_label('Jy/beam',rotation=90, labelpad=12,size=18)

#         # ax.contour(image_data,levels=[-3*0.136e-3,3*0.136e-3,5*0.136e-3,10*0.136e-3,15*0.136e-3], colors='white',
#         #         linewidths=0.5)

#         plt.savefig(fitsname.replace('.fits','_1.pdf'),dpi=300)

   



# class SelfCalibration():

#     def selfcal():

#         if os.path.exists(outlier_file) and open(outlier_file).read() == '':
#             outlierfile = ''

#         regionfile = basename + '_first_masking.casabox'
#         cell = get_imaging_cellsize()
#         print("Deleting model column before selfcal")
#         delmod(vis=vis,otf=True)

#         for selfcal_loop in range(nloops):
#             caltable = f'caltable_{selfcal_loop}.gcal'
#             prev_caltables = sorted(glob.glob('*.gcal'))
#             if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
#                 applycal(vis=vis, gaintable = prev_caltables, parang=False )
        
#             imagename = f'target_selfcal_{selfcal_loop}'
#             if os.path.exists(imagename):
#                 print("Continuing to the next image")
            
#             else:
#                 # imagename = f'target_selfcal_{selfcal_loop}'
#                 print(f"Making image {imagename}")
#                 tclean(
#                     vis = vis, imagename=imagename, imsize=imsize, cell=cell,
#                     parallel=False,
#                     gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
#                     weighting = weighting, robust = robust, niter=niter[selfcal_loop], threshold = threshold[selfcal_loop],
#                     nterms = nterms, pblimit = -1,interactive=False, usemask='user', mask=regionfile
#                 )

#                 ## NB: The problem was niter -- there was a space in the list []

#                 print("Adding modelcolumn to data")
#                 # model images from the MTMFS images,
#                 ft(vis = vis, model=[imagename+'.model.tt0',imagename+'.model.tt1'], nterms=2,usescratch=True)

#                 # plot the model column
#                 plotms(
#                     vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
#                     showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
#                 )

#                 gaincal( vis =vis, caltable = caltable, refant = refant, solint = solint[selfcal_loop],
#                         gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
#                         calmode = calmode[selfcal_loop], append=False, parang=False
#                         )
#                 coloraxis = ['corr','spw']
#                 for color in coloraxis:
#                     if calmode[selfcal_loop] =='p':
#                         plotms(
#                             vis = caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
#                             iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
#                             plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750,
#                         )
#                     else:
#                         plotms(
#                                 vis = caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
#                                 iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
#                                 plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750
#                             )

#                 if selfcal_loop == nloops-1:
#                     prev_caltables = sorted(glob.glob('*.gcal'))
#                     print("Applying the caltable derived from last gaincal iteration")
#                     applycal(vis=vis, gaintable = prev_caltables, parang=False )
            
#             # ### Get the last imagename from the loop and generate a final mask
            
#         imagename = basename +f'_{nloops-1}'+'.final'
#         ##  tclean here to make the final image
#         print("Make final image with all selfcal corrections applied")
#         tclean(
#             vis = vis, imagename = imagename, imsize=imsize, cell=cell, gridder=gridder,
#             wprojplanes = wprojplanes, deconvolver = deconvolver, weighting = weighting,
#             robust = robust, niter=niter_final, threshold = threshold_final, nterms=nterms,
#             pblimit=pblimit, interactive=False, usemask = 'user', mask=regionfile,
#         )
        
        


#     # def applycal_target():

#     #     """
#     #     Applies cal to the target field
#     #     """
#     #     prev_caltables = sorted(glob.glob('*.gcal'))
#     #     applycal(
#     #         vis = vis_tocal, gaintable = prev_caltables, parang=False
#     #     )



#     # def peeling():

#     #     """
#     #     Subtract the bright sources in the field so you are left with the emission from the star

#     #     Use pybsdf casa region file -- you can change the threshold to only select the bright sources you want

#     #     Run pybsdf on the final image
        
#     #     """

#     #     # Using the final imagename 

#     #     imagename = basename +f'_{nloops-1}'+'.final'
#     #     # region_to_peel = pybdsf(input_image=imagename+'.image.tt0')

#     #     ## NB: You had masked the bright source from the dirty map -- so you can just subtract
#     #     uvsub(vis=vis)

#     #     ## Make a final map without the sources
#     #     cell = get_imaging_cellsize()
#     #     peeled_map = basename+'_peeled_map'
#     #     if not os.path.exists(peeled_map):
#     #         print(f"Making {peeled_map}")
#     #         tclean(
#     #             vis = vis, imagename=peeled_map, imsize=imsize, cell=cell,
#     #             gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
#     #             weighting = weighting, robust = robust, niter=0, # threshold = '0.5mJy',
#     #             nterms = nterms, pblimit = pblimit
#     #         )








#     # delete model column to be d

#     def peeling():

#         """
#         Subtract all the sources in the field such that you are left with a blank

#         Use wsclean and pybdsf casa region from the final iterations 

#         Then perform uvsub in CASA

#         """

#         # Make a region file of the final self calibrated image and use it to peel the sources
#         imagename = basename +f'_{nloops-1}'+'.final.image.tt0'
#         regionfile_to_peel = pybdsf(input_image=imagename)
#         fitsmask = imagename.replace('.image.tt0','')+'.maskfile.fits'
#         model_fits = imagename.replace('.final.image.tt0','.final.image.tt0-model.fits')
#         os.rename(fitsmask,model_fits )

#         cell = get_imaging_cellsize()


#         threshold_cmd = ['wsclean', '-auto-threshold','3', '-size', f'{imsize[0]}', f'{imsize[1]}','-scale', f'{cell}',\
#                         '-mgain', '0.8', '-niter', '0',f'{vis}']
        
#         predict_cmd = ['wsclean', '-log-time', '-predict', '-field', '', '-reorder' ,'-name', f'{imagename}', '-abs-mem',f'{abs_mem}', vis]


#         run_wsclean(predict_cmd)

#         ## NB: wsclean needs to find an image named my-image-model.fits or reg 
#         ## works by replacing model column with model for the problem sources using

        
        
#         ## Subtract the models put in the model column from the data and make an image
            
#         print("Running uvsub")
#         uvsub(vis=vis)

#         ## Run wsclean to check if the subtraction has been successful -- make dirty map

#         run_wsclean(threshold_cmd)


# def main():

#     set_working_dir()
#     split_data()
#     check_band()
#     global refant
#     refant = find_refant(vis)
#     make_dirty_map()
#     selfcal()
#     # applycal_target()
#     peeling()

#     # imagename = 'K2-18_dirty'
#     # plot_fits(imagename)



# main()

        


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
from radio_beam import Beam

from casatasks import *
from casaplotms import *
from casatools import componentlist

msmd = casatools.msmetadata()
tb = casatools.table()
ms = casatools.ms()


wsclean_sif = '/raid1/scratch/kelvinw/singularity_containers/wsclean_working.simg'



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
    def run_wsclean(command):
            """
            Runs wsclean commands and logs the output and errors in real-time.
            """
            # Get the path to bind the singularity container
            singularity_bind = os.path.join(os.path.dirname(os.path.dirname(wsclean_sif)))
            # Form the command to execute with Singularity
            command_to_execute = ['singularity', 'exec', '-B', singularity_bind, wsclean_sif] + command
            print(f"Executing: {' '.join(command_to_execute)}")
            try:
                # logging.info(f"Executing: {' '.join(command_to_execute)}")
                process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

                # Read stdout and stderr as they are produced
                for stdout_line in iter(process.stdout.readline, ""):
                    # Log stdout in real time
                    logging.info(stdout_line.strip())
                for stderr_line in iter(process.stderr.readline, ""):
                    # Log stderr in real time
                    logging.error(stderr_line.strip())

                # Wait for the process to finish and get the return code
                return_code = process.wait()

                if return_code == 0:
                    logging.info("Strategy executed successfully.")
                else:
                    logging.error(f"Error executing strategy. Return code: {return_code}")

            except Exception as e:
                logging.exception("An error occurred during wsclean execution: %s", e)


    @staticmethod
    def get_im_stats(imagename):
        
        #### Not properly integrated into the code yet -- this is useful for EVN images

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




    # @staticmethod
    # def pybdsf(input_image, detection_threshold):

    #     # Check if the input image is a FITS file; if not, add .fits
    #     if not input_image.endswith('.fits'):
    #         input_image = input_image
    #         fitsname =input_image + '.fits'
    #     else:
    #         # If it is already a FITS file, use it directly
    #         fitsname = input_image

    #     # Process the FITS image with pybdsf
    #     img = bdsf.process_image(fitsname, adaptive_rms_box=True, thresh='hard',
    #                             thresh_isl=True, thresh_pix=detection_threshold, 
    #                             advanced_opts=True, mean_map='map', rms_map=True, 
    #                             group_by_isl=True)

    #     # Write out island mask and FITS catalog
    #     img.export_image(outfile=input_image + '.maskfile.fits', img_type='island_mask', img_format='fits', clobber=True)
    #     img.write_catalog(outfile=input_image + '.cat', format='fits', clobber=True, catalog_type='gaul')

    #     regionfile = input_image + '.casabox'
    #     ascii_file = input_image + '.ascii'
    #     rmsfile = input_image + '.rmsfile'

    #     img.write_catalog(outfile=regionfile, format='casabox', clobber=True, catalog_type='srl')
    #     img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
    #     img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

    #     return regionfile
    
    @staticmethod
    def make_mask(fits_file, rms, threshold):
        """
        This function reads a FITS file, masks all values above the given RMS threshold,
        and writes the mask as a new FITS file, retaining the WCS and header information.

        Parameters:
        -----------
        fits_file : str
            The path to the input FITS file.
        rms : float
            The RMS value for scaling the threshold.
        threshold : float
            The multiplicative factor for the RMS threshold to create the mask.

        Returns:
        --------
        output_file : str
            The path to the output FITS mask file.
        """
        with fits.open(fits_file) as hdul:
            data = hdul[0].data
            header = hdul[0].header

        if data is None:
            raise ValueError("No image data found in the FITS file.")

        mask = data > (threshold * rms)
        mask = mask.astype(np.int16)
        hdu = fits.PrimaryHDU(data=mask, header=header)

        output_file = fits_file.replace('.fits', '_masking.fits')
        hdu.writeto(output_file, overwrite=True)
        logging.info(f"Mask saved to {output_file}")

        return output_file

        


class MeasurementSetProcessor:

    """ Class to process easurement set 

    Attributes:
        avgtime (str): Time averaging interval for the split, specified as a string.
                       Defaults to '10s'.
        avgwidth (int): Channel averaging width for the split. Defaults to 4.

    """


    # Default values if no specific `avgtime` or `width` is provided
    default_timebin: str = '10s'
    default_width: int = 4

    @staticmethod
    @Utils.time_execution
    def prepare_measurement_set(msname: str, fieldname: str = None, split_required: bool = True, avgtime: str = None, width: int = None) -> Union[Tuple[str], Tuple[List[str]]]:
        """
        Prepares the measurement set by splitting it into individual fields if necessary.

        Parameters:
        - msname (str): Path to the measurement set file. Can be a full path or a relative path.
        - split_required (bool): Whether to split the measurement set into individual fields.
                                 Defaults to True.
        - avgtime (str, optional): Custom time averaging interval for splitting.
                                   Uses `default_timebin` if not provided.
        - width (int, optional): Custom channel averaging width for splitting.
                                 Uses `default_width` if not provided.
        - fieldname (str, optional): Specific fieldname to check and split (if desired).

        Returns:
        - tuple: If splitting is not required, returns msname directly.
                If splitting is required and a fieldname is provided, splits by that field.
                Otherwise, returns a tuple of the output measurement set names.
        """
        msmd = casatools.msmetadata()
        outputvis_list = []

        try:
            # Ensure msname is a full path for consistent processing
            msname = os.path.abspath(msname)
            listobs(vis=msname, listfile='listobs.txt', overwrite=True)
            msmd.open(msname)
            field_names = msmd.fieldnames()
            num_fields = len(field_names)
            logging.info(f"Field names: {field_names} found in {msname}")

            # If splitting is not required, return msname directly
            if not split_required:
                logging.info("Splitting not required.")
                return msname  # Return as a single-item tuple

            # Otherwise, proceed with splitting even if only one field exists
            timebin = avgtime if avgtime is not None else MeasurementSetProcessor.default_timebin
            split_width = width if width is not None else MeasurementSetProcessor.default_width

            current_working_dir = os.getcwd()
            if isinstance(fieldname, str):
                fieldname = [fieldname]  
            if fieldname:
                # Split the specified fields (target, phasecalibrator, etc.)
                for field in fieldname:
                   
                    if field in field_names:
                        outputvis = os.path.join(current_working_dir, f"{field}_split_{timebin}_{split_width}.ms")
                        if not os.path.exists(outputvis):
                            logging.info(f"Splitting {msname} to {outputvis} with timebin={timebin} and width={split_width}")
                            split(vis=msname, outputvis=outputvis, datacolumn='corrected',
                                timebin=timebin, width=split_width, field=field)
                            listobs(vis=outputvis, listfile=outputvis.replace('.ms', '_listobs.txt'), overwrite=True)
                            logging.info("Finished splitting")
                        else:
                            logging.info(f"Split measurement set {outputvis} exists")
                        outputvis_list.append(outputvis)
                    else:
                        logging.warning(f"Field {field} not found in measurement set.")
                
            else:
                # If no specific field is given, split all fields
                for field in field_names:
                    outputvis = os.path.join(current_working_dir, f"{field}_split_{timebin}_{split_width}.ms")
                    if not os.path.exists(outputvis):
                        logging.info(f"Splitting {msname} to {outputvis} with timebin={timebin} and width={split_width}")
                        split(vis=msname, outputvis=outputvis, datacolumn='corrected',
                            timebin=timebin, width=split_width, field=field)
                        listobs(vis=outputvis, listfile=outputvis.replace('.ms', '_listobs.txt'), overwrite=True)
                        logging.info("Finished splitting")
                    else:
                        logging.info(f"Split measurement set {outputvis} exists")
                    outputvis_list.append(outputvis)

        except Exception as e:
            logging.error(f"Error in prepare_measurement_set: {e}")

        finally:
            msmd.close()

        return outputvis


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
        ms.open(msname)
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
        longest_baseline_lambda = MeasurementSetInfo.get_longest_baseline(msname)
        
        # Calculate cell size in arcseconds
        cell_float = (180.0 * 3600 / (np.pi * 5)) * (1.0 / longest_baseline_lambda)
        
        # Convert to mas if the value is very small (e.g., <1 arcsecond)
        if cell_float < 1.0:
            cell_float *= 1000  # convert to mas
            cell = f'{cell_float:.2f} mas'
        else:
            cell = f'{cell_float:.2f} arcsec'
        
        logging.info(f"Imaging with a cell size of {cell}")
        return cell

        
    
    @staticmethod
    # @Utils.time_execution
    def find_refant(msname):

        """
        
        Find the best reference antenna by calculation the SNR by making a mock calibration table.
        The function will find the tuple with the same name as the field and use that

        Parameters
            msname: Measurement set file
            field (str): The field to self calibrate

        """
        msmd = casatools.msmetadata()
        msmd.open(msname)
        field_names = msmd.fieldnames()
        # field = msname.replace('.ms', '')
        field = field_names[0]
        if field in field_names:
            logging.info(f"Phase solutions for {field} will be used to select the reference antenna.")
        else:
            logging.warning(f"Field '{field}' not found in field names. Available fields: {field_names}")
        msmd.close()

        tablename = field+'.refant'
        if not os.path.exists(tablename):
            gaincal(vis= msname, caltable=tablename, field=field, refantmode='flex',
                     solint='inf', minblperant=3, gaintype='G', calmode='p')
            logging.info(f"Phase solutions to select refant generated")

        # Read solutions (phases):
        tb = casatools.table()
        try:
            tb.open(tablename + '/ANTENNA')
            antenna_names = tb.getcol('NAME')
            tb.close()
            tb.open(tablename)
            antenna_ids = tb.getcol('ANTENNA1')
            num_antennas = len(antenna_ids)
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

    @staticmethod
    def get_msinfo(msname):

        nchan = []
        msmd = casatools.msmetadata()
        msmd.open(msname)
        bandwidth = msmd.bandwidths()
        nspw = len(bandwidth)
        for spw in range(nspw):
            nchan.append(msmd.nchan(spw))
        msmd.close()

        return nspw,nchan



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
                 niter: int = 0, threshold: str = None, wprojplanes: int = 1, mask: str = '', imagename: str = None,
                 usemask: str = 'user', pblimit: float = 0.1,phasecenter: str = '', overwrite:bool = False, 
                 use_pybdsf:bool = False , pybdsf_threshold: int = 5,cell:list = None,parallel:bool = False,
                 outlierfile:str = None, pbcorrect:bool = False):
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
        parallel: bool (default is False)
            Use MPI for imaging
        outlierfile: str (default is empty)
            "Outlier fields to image"
        pbcorrect: bool 
            "Correct the primary beam using widebandpbcor

        """
        self.msname = msname  
        self.field = field
        self.gridder = gridder
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
        self.parallel = parallel
        self.pbcorrect = pbcorrect
        
        self.outlierfile = outlierfile if outlierfile else ""  
        self.imagename = imagename if imagename else f"{self.msname.replace('.ms', '_image')}"
        
        logging.info(f"Using outlierfile: {self.outlierfile if self.outlierfile else 'None'}")

        if cell is not None:
            self.cell = cell
            logging.info(f"Manually supplied imaging cell size: {self.cell}")
        else:
            self.cell = MeasurementSetInfo.get_imaging_cellsize(self.msname)
            logging.info(f"Calculated imaging cell size based on MS: {self.cell}")


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
        cell = MeasurementSetInfo.get_imaging_cellsize(self.msname)
        self.imagename = self.imagename

    
        matching_files = glob.glob(f"{self.imagename}.*")

        # Check if any matching files exist
        if matching_files:
            logging.info(f"Image exists: {matching_files}")

                # Skip the function if overwrite is False
            if not self.overwrite:
                logging.info(f"Overwrite is set to False. Skipping imaging for {self.imagename}")
                return  # Exit the function without proceeding

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
    
        logging.info(f"Imaging {self.msname} to make: {self.imagename}")


        ### NB: The threshold here is an rms based threshold (nsigma) as opposed to an absolute one (threshold)

        if isinstance(self.threshold, str) and self.threshold == '':
            self.threshold = 0.0
        elif isinstance(self.threshold, int):
           self.threshold = float(self.threshold)
        tclean(
            vis=self.msname,
            imagename=self.imagename,
            imsize=[self.imsize, self.imsize], 
            cell=self.cell,
            gridder=self.gridder,
            deconvolver=self.deconvolver,
            weighting=self.weighting,
            robust=self.robust,
            niter=self.niter,
            nterms=self.nterms,
            nsigma=self.threshold,
            wprojplanes=self.wprojplanes,
            mask=self.mask,
            usemask=self.usemask,
            pblimit=self.pblimit,
            field = self.field,
            parallel = self.parallel,
            outlierfile = self.outlierfile,
            interactive=False,
            
        )

        logging.info(f"Finished imaging {self.msname}, created image: {self.imagename}")



    def pbcorr(self) -> None:
        """
        Perform primary beam correction using widebandpbcor.

        This method applies primary beam correction to the imaged data using the widebandpbcor function.
        
        Returns:
        -------
        None
        """
        widebandpbcor(
            vis=self.msname,
            imagename=self.imagename,
            nterms=2,
            action='pbcor',
        )
        logging.info(f"Applied primary beam correction to {self.imagename}.")

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



class WSClean_Imager:
    """
    A class for performing imaging on measurement sets using the WSClean algorithm.

    Methods:
    -------
    get_imaging_cellsize() -> str:
        Calculates and returns the cell size for imaging in arcseconds.
    
    imager() -> None:
        Performs the imaging process using the WSClean algorithm.

    """

    def __init__(self, msname: str, imsize: int = 640, niter: int = 0, threshold: str = 0.0, deconvolution:str = None,
                 overwrite: bool = False, use_pybdsf: bool = False , pybdsf_threshold: int = 5, mgain: float = 0.8,
                imagename: str = None, wsclean_sif: str = None, maskfile: str = '', cell:list = None):
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
        self.maskfile = maskfile

        self.imagename = imagename if imagename else f"{self.msname.replace('.ms', '_image')}"

         # Use the manually supplied cell or calculate it if not provided
        self.cell = cell if cell is not None else MeasurementSetInfo.get_imaging_cellsize(self.msname)

        logging.info(f"Using an imaging cell size: {cell}")

    @Utils.time_execution
    def imager(self) -> None:
        """
        Perform the imaging process using the WSClean algorithm.

        Returns:
        -------
        None
        """
                
        command = [
            "wsclean",
            "-log-time",
            "-size", str(self.imsize), str(self.imsize),
            "-reorder",      
            "-name", self.imagename,
            "-scale", self.cell,
            "-mgain", str(self.mgain),
            "-niter", str(self.niter),
            "-auto-threshold", str(self.threshold),
            "-fits-mask", str(self.maskfile),
        ]
        if self.maskfile == '':
            logging.info(f"Masking not requested.")
        else:
            logging.info(f"Using maskfile: {self.maskfile}")
        if self.deconvolution == "multiscale":
            command.append("-multiscale")

        command.append(self.msname)
        
        try:
            Utils.run_wsclean(command)
            logging.info(f"Finished imaging {self.msname}, created image: {self.imagename}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error during WSClean imaging: {e}")

    @Utils.time_execution
    def predict(self) -> None:
        """
        Perform the prediction process using the WSClean algorithm.

        Parameters:
        ----------
        phasecal_ms : str
            The path to the measurement set for the prediction process.
        """
        predict_cmd = [
            'wsclean',
            '-log-time',
            '-predict',
            '-reorder',
            '-name', self.imagename,
            self.msname
        ]

        # Run the WSClean prediction command
        try:
            Utils.run_wsclean(predict_cmd)
            logging.info(f"Prediction completed successfully with imagename: {self.imagename}")
        except Exception as e:
            logging.error(f"Error during prediction: {e}")


class PlottingRoutines:

    """
    Functions to plot fits files
    """


    def __init__(self, imagename,color='magma',figsize=(10,8)):
        """
        Initialize the ImageProcessor class.

        :param imagename: Name of the FITS file for the image (excluding file extension).

        """

        self.imagename = imagename
        self.color = color
        self.figsize = figsize

    @property
    def imagename(self):
        """Getter for the imagename property."""
        return self._imagename
    
    @imagename.setter
    def imagename(self, new_imagename):
        """Setter for the imagename property."""

        if not isinstance(new_imagename, str):
            raise ValueError("Imagename must be a string")
        self._imagename = new_imagename
        print(f"Imagename updated to: {self._imagename}")

    def get_beam(self):
        """
        Get the beam in arcsec.
        """
        imagename_header = fits.getheader(self.imagename)
        imaging_beam = Beam.from_fits_header(imagename_header)

        ## beam in arcsec
        bmaj = imaging_beam.major.to(u.arcsec).value
        bmin = imaging_beam.minor.to(u.arcsec).value
        pa = imaging_beam.pa.to(u.deg).value  

        return (bmaj, bmin, pa)

        
    def plot_image_with_beam(self):
            """
            Plot the FITS image and place the beam at the bottom-left corner.
            """
            # Read the FITS file data
            hdu = fits.open(self.imagename)
            image_data = hdu[0].data[0, 0, :, :] 
            header = hdu[0].header
            w = WCS(header, naxis=2)
            fig, ax = plt.subplots(figsize=self.figsize, subplot_kw={'projection': w})
            im = ax.imshow(image_data, cmap=self.color, origin='lower', interpolation='none')
            ax.coords[0].set_auto_axislabel(True) 
            ax.coords[1].set_auto_axislabel(True) 
            shape = header['NAXIS1'], header['NAXIS2']
            bmaj, bmin, pa = self.get_beam()
            relative_x = 15
            relative_y = 15
            x_pos = (relative_x / 320) * shape[0]  
            y_pos = (relative_y / 320) * shape[1]  
            beam_ellipse = patches.Ellipse(
                (x_pos,y_pos), width=bmaj, height=bmin, angle=pa, edgecolor='white', facecolor='none', lw=2)
            ax.add_patch(beam_ellipse)
            ax.set_xlabel('RA (J2000)', size=14)
            ax.set_ylabel('Dec (J2000)', size=14)
            ax.tick_params(axis="x", which="both", bottom=True, top=False)
            ax.tick_params(axis="y", which="both", right=False, left=True)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=16)
            cbar.set_label('Jy/beam', rotation=90, labelpad=12, size=18)
            cbar.formatter = ScalarFormatter()
            cbar.formatter.set_powerlimits((-3, 3))
            cbar.update_ticks()
            # plt.show()
            plt.savefig(self.imagename.replace('.fits','.pdf'),dpi=300)


class SelfCalibrationWSClean(WSClean_Imager):

    def __init__(self, msname, nloops, thresholds, calmode, gaintype, solint, minsnr, refant=None, final_image: bool = False, pbcorrect:bool = False, **kwargs):
        super().__init__(msname=msname, **kwargs)
        self.nloops = nloops
        self.thresholds = thresholds
        self.calmode = calmode
        self.gaintype = gaintype
        self.solint = solint
        self.minsnr = minsnr
        self.final_image = final_image
        self.pbcorrect = pbcorrect


        self.refant = refant if refant is not None else str(MeasurementSetInfo.find_refant(self.msname))

        # print(self.refant)

        

    @Utils.time_execution
    def selfcal(self) -> None:
        Utils.create_plots_directory()

        logging.info("Deleting model column before selfcal")
        delmod(vis=self.msname, otf=True)

        for selfcal_loop in range(self.nloops):
            caltable = f'{self.msname.replace(".ms", "_caltable_loop")}_{selfcal_loop}.gcal'
            prev_caltables = sorted(glob.glob('*.gcal'))

            # Apply previous calibration tables if they exist and if it's not the dirty map loop
            if selfcal_loop > 0 and len(prev_caltables) > 0 and self.calmode[selfcal_loop] != '':
                applycal(vis=self.msname, gaintable=prev_caltables, parang=False)
            
            
            imagename = f'{self.msname.replace(".ms", "_selfcal_loop")}_{selfcal_loop}'

            # Set niter to 0 for the first loop to create a dirty map and run PYBDSF
            if selfcal_loop == 0:
                # self.niter = 0
                imagename_dirty = self.imagename+'_dirty'

                imager_instance = WSClean_Imager(
                    msname = self.msname, 
                    imagename = imagename_dirty,
                    imsize = self.imsize,  
                    threshold  =  self.threshold, ## threshold here can be 0.0 which is the default, 
                    overwrite = self.overwrite,
                    use_pybdsf = self.use_pybdsf,
                    pybdsf_threshold = self.pybdsf_threshold,
                    mgain = self.mgain,
                    cell = self.cell,
                    niter = 1, # niter will ensure you dont hit a thresh of 0.0, also note niter=0 will fail in pybdsf
                    )
                imager_instance.imager()

                plotter = PlottingRoutines(imagename = imagename_dirty + '-image.fits')
                plotter.plot_image_with_beam()

                # If using PYBDSF for masking, call it here
                if self.use_pybdsf:
                    try:
                        logging.info(f"Running pybdsf on {imagename_dirty}...")
                        Utils.pybdsf(imagename_dirty + '-image.fits', self.pybdsf_threshold)
                        self.maskfile = imagename_dirty + '-image.fits.maskfile.fits'
                        self.maskfile.replace('.fits','') + '.fits' ## removing the repeated .fits name
                        logging.info(f"PYBDSF output file: {self.maskfile} will be used for masking")
                        logging.info(f"Successfully ran pybdsf on {imagename_dirty}.")
                    except Exception as e:
                        logging.error(f"Failed to run pybdsf on {imagename_dirty}: {e}")
                else:
                    logging.info("Masking using PYBDSF not requested.")
                    self.maskfile = ''

            else:
                # Set niter for subsequent loops and disable pybdsf
                self.use_pybdsf = False
                
                imager_instance = WSClean_Imager(
                    msname = self.msname, 
                    imagename = imagename,
                    imsize = self.imsize,  
                    niter = self.niter,
                    threshold = self.thresholds[selfcal_loop],
                    overwrite = self.overwrite,
                    use_pybdsf = self.use_pybdsf,
                    pybdsf_threshold = self.pybdsf_threshold,
                    mgain = self.mgain,
                    cell = self.cell,
                    maskfile = self.maskfile
                    )
                imager_instance.imager()

                plotter = PlottingRoutines(imagename = imagename + '-image.fits')
                plotter.plot_image_with_beam()

                ###TODO logging.info not working properly, using print
                ### Put image model in measurement set  MODEL COLUMN -- similar to ft in CASA
                # logging.info(f"Making imagename: {imagename} for selfcal loop: {selfcal_loop} using niter: {self.niter}")
                print(f"Making imagename: {imagename} for selfcal loop: {selfcal_loop} using niter: {self.niter}")
                ## WSClean must find an image named --model.fits (IN CWD) in order to predict !
                model_fits = imagename.replace('-image.fits','-model.fits')
                logging.info(f"Adding modelcolumn to data. Using {model_fits} to predict")
                print(f"======>>>Adding modelcolumn to data. Using {model_fits} to predict")


                imager_instance.predict()

                ### Make plot to verify that predict is working
                nspw,nchan = MeasurementSetInfo.get_msinfo(self.msname)
                logging.info("Plotting the model column")
                try:
                    plotms(
                        vis=self.msname, xaxis='UVwave', yaxis='amp', ydatacolumn='model', avgchannel=str(nspw), avgtime='300',
                        showgui=False, plotfile=imagename + '_modelcolumn.png', overwrite=True, width=1500, height=750,
                    )
                except Exception as e:
                    logging.error(f"Error plotting model column: {e}")


                ## Run gaincal
                
                if self.calmode[selfcal_loop] == 'p':
                    minblperant = 3
                else:
                    minblperant = 4

                # try:
                gaincal(vis=self.msname,
                        caltable=caltable,
                        refant=self.refant,
                        solint=self.solint[selfcal_loop],
                        gaintype=self.gaintype[selfcal_loop],
                        gaintable=prev_caltables,
                        minsnr=self.minsnr[selfcal_loop],
                        calmode=self.calmode[selfcal_loop],
                        minblperant = minblperant,
                        # combine ='spw',
                        append=False,
                        parang=False)
                # except Exception as e:
                #     logging.error(f"Error during gain calibration: {e}")

                gridcols = 7
                gridrows = 4  
                # Loop over the coloraxis values (corr and spw)
                coloraxis = ['corr', 'spw']
                for color in coloraxis:
                    if self.calmode[selfcal_loop] == 'p':
                        plotms(
                            vis=caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
                            iteraxis='antenna', coloraxis=color, showgui=False, overwrite=True,
                            plotfile=caltable.replace('.gcal', f'_{color}.png'), dpi=300, width=1500, height=750,
                        )
                    else:
                        plotms(
                            vis=caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
                            iteraxis='antenna', coloraxis=color, showgui=False, overwrite=True,
                            plotfile=caltable.replace('.gcal', f'_{color}.png'), dpi=300, width=1500, height=750,
                        )
                # Apply calibration tables after the last self-calibration loop
                if selfcal_loop == self.nloops - 1:
                    prev_caltables = sorted(glob.glob('*.gcal'))
                    logging.info("Applying the caltable derived from last gaincal iteration")
                    applycal(vis=self.msname, gaintable=prev_caltables, parang=False)

        if self.final_image:
            ## Generate a final mask of sources to peel -- optional
            imagename_final = self.imagename+'_final_clean'
            logging.info("Making final image with all selfcal corrections applied")
            
            # self.niter = 1000000  # Can be modified to set a new value if needed
            # self.threshold = '0.001mJy'
            self.use_pybdsf = False

            imager_instance = WSClean_Imager(
                msname = self.msname, 
                imagename = imagename_final,
                imsize = self.imsize,  
                niter = self.niter,
                threshold = self.threshold, # attempt to go to 0.0
                overwrite = self.overwrite,
                use_pybdsf = self.use_pybdsf,
                pybdsf_threshold = self.pybdsf_threshold,
                mgain = self.mgain,
                maskfile = self.maskfile,
                cell = self.cell,
                )
            
            imager_instance.imager()
            plotter = PlottingRoutines(imagename = imagename_final + '-image.fits')
            plotter.plot_image_with_beam()


        

    # def applycal_target():

    #     """
    #     Applies cal to the target field
    #     """
    #     prev_caltables = sorted(glob.glob('*.gcal'))
    #     applycal(
    #         vis = vis_tocal, gaintable = prev_caltables, parang=False
    #     )




class SelfCalibrationTclean(tclean_Imager):

    def __init__(self, msname, nloops, thresholds, calmode, gaintype, solint, minsnr, make_final_image:bool = True, refant=None, masking_threshold:int = 5,**kwargs):
        super().__init__(msname=msname, **kwargs)
        self.nloops = nloops
        self.thresholds = thresholds
        self.calmode = calmode
        self.gaintype = gaintype
        self.solint = solint
        self.minsnr = minsnr
        self.make_final_image = make_final_image
        self.masking_threshold = masking_threshold

        self.refant = refant if refant is not None else str(MeasurementSetInfo.find_refant(self.msname))
        # print(self.refant)

    @Utils.time_execution
    def selfcal(self) -> None:
        Utils.create_plots_directory()
         
        
        logging.info("Deleting model column before selfcal")
        delmod(vis=self.msname, otf=True)

        for selfcal_loop in range(self.nloops):
            caltable = f'{self.msname.replace(".ms", "_caltable_loop_")}_{selfcal_loop}.gcal'
            prev_caltables = sorted(glob.glob('*.gcal'))

            # Apply previous calibration tables if they exist and if it's not the dirty map loop
            if selfcal_loop > 0 and len(prev_caltables) > 0 and self.calmode[selfcal_loop] != '':
                applycal(vis=self.msname, gaintable=prev_caltables, parang=False)

            imagename = f'{self.msname.replace(".ms", "_selfcal_loop")}_{selfcal_loop}'

            # Set niter to 0 for the first loop to create a dirty map and run PYBDSF
            if selfcal_loop == 0:
                # self.niter = 0
                imagename_dirty = self.imagename+'_dirty'

                mask = ''

                # Run tclean to generate the dirty map
                imager_instance = tclean_Imager(
                    msname=self.msname,
                    imagename=imagename_dirty, # imagename for dirty defined in imagr class
                    nterms=self.nterms,
                    imsize=self.imsize,
                    cell=self.cell,
                    niter=0,
                    deconvolver=self.deconvolver,
                    threshold=self.thresholds[selfcal_loop],
                    mask=mask,
                    weighting = self.weighting,
                    robust = self.robust,
                    overwrite=self.overwrite,
                    outlierfile = self.outlierfile,
                )
                logging.info(f"Imaging {self.msname} to make: {self.imagename}")
                imager_instance.imager()

                # Check if primary beam correction is requested and call the method
                if self.pbcorrect:
                    logging.info(f"Primary beam correction requested, applying pbcorr.")
                    imager_instance.pbcorr()  # Call the pbcorr method

                # Export dirty map to FITS after it has been created
                image_ext = '.image.tt0' if self.deconvolver == 'mtmfs' else '.image'
                fitsimage = imagename_dirty+ image_ext + '.fits'
                exportfits(imagename=imagename_dirty + image_ext, fitsimage=fitsimage, overwrite=True)
                
                """ Make a masking file """
                image_rms = imstat(fitsimage)['rms'][0]
                global masking_file
                masking_file = Utils.make_mask(fits_file=fitsimage, rms=image_rms, threshold=self.masking_threshold)
                # (Regrid the mask to match the image if necessary)
                masking_image = masking_file.replace('.fits', '.im')
                importfits(fitsimage=masking_file,imagename=masking_image,overwrite=True)
                makemask(
                    mode='copy',
                    inpimage=masking_image,
                    inpmask=masking_image,
                    output=imagename_dirty + '.mask', 
                )
                masking_file = imagename_dirty + '.mask'
                                

                # # Run PYBDSF if requested
                # try:
                #     logging.info(f"Running pybdsf on {imagename_dirty}...")
                #     if self.use_pybdsf:
                #         Utils.pybdsf(imagename_dirty+ image_ext, self.pybdsf_threshold)
                #         mask = imagename_dirty+ image_ext + '.casabox'
                #         logging.info(f"Successfully ran pybdsf on {imagename_dirty}.")
                #     else:
                #         logging.info(f"Masking using PYBDSF not requested")
                #         mask = ''
                # except Exception as e:
                #     logging.error(f"Failed to run pybdsf on {imagename_dirty}: {e}")

            else:
                # Set niter for subsequent loops and disable pybdsf
                self.niter = self.niter  # Can be modified to set a new value if needed
                # self.niter = 1 # the issue is that this was not getting initialised properly
                self.use_pybdsf = False
                logging.info(f"Making imagename: {imagename} for selfcal loop: {selfcal_loop} using niter: {self.niter}")
                # Initialize imager_instance with required parameters for self-calibration
                imager_instance = tclean_Imager(
                    msname=self.msname,
                    imagename=imagename,
                    nterms=self.nterms,
                    imsize=self.imsize,
                    niter=self.niter,
                    cell=self.cell,
                    deconvolver=self.deconvolver,
                    threshold=self.thresholds[selfcal_loop],
                    weighting = self.weighting,
                    robust = self.robust,
                    mask=masking_file,
                    overwrite=self.overwrite,
                    outlierfile=self.outlierfile,
                )
                imager_instance.imager()
                

                # Perform gain calibration only if this is not the dirty map loop
                logging.info("Adding model column to data")
                try:
                    if self.deconvolver == 'mtmfs':
                        ft(vis=self.msname, model=[imagename + '.model.tt0', imagename + '.model.tt1'], nterms=2, usescratch=True)
                    else:
                        ft(vis=self.msname, model=imagename + '.model', usescratch=True)
                except Exception as e:
                    logging.error(f"Error adding model column: {e}")
                nspw,nchan = MeasurementSetInfo.get_msinfo(self.msname)
                logging.info("Plotting the model column")
                try:
                    plotms(
                        vis=self.msname, xaxis='UVwave', yaxis='amp', ydatacolumn='model', avgchannel=str(nspw), avgtime='300',
                        showgui=False, plotfile=imagename + '_modelcolumn.png', overwrite=True, width=1500, height=750,
                    )
                except Exception as e:
                    logging.error(f"Error plotting model column: {e}")

                

                # Perform gain calibration
                logging.info(f"Running gain calibration. Writing caltable: {caltable}")
                
                if self.calmode[selfcal_loop] == 'p':
                    minblperant = 3
                    solnorm = False
                else:
                    minblperant = 4
                    solnorm = True
                # try:
                gaincal(vis=self.msname,
                        caltable=caltable,
                        refant=self.refant,
                        solint=self.solint[selfcal_loop],
                        gaintype=self.gaintype[selfcal_loop],
                        gaintable=prev_caltables,
                        minsnr=self.minsnr[selfcal_loop],
                        calmode=self.calmode[selfcal_loop],
                        minblperant = minblperant,
                        solnorm = solnorm, 
                        append=False,
                        parang=False)
                # except Exception as e:
                #     logging.error(f"Error during gain calibration: {e}")

                gridcols = 7
                gridrows = 4  
                # Loop over the coloraxis values (corr and spw)
                coloraxis = ['corr', 'spw']
                for color in coloraxis:
                    if self.calmode[selfcal_loop] == 'p':
                            plotms(
                                vis=caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
                                iteraxis='antenna', coloraxis=color, showgui=False, overwrite=True,
                                plotfile=caltable.replace('.gcal', f'_phase_{color}_.png'), dpi=300, width=3000, height=1500,
                            )

                    else:
                          plotms(
                                    vis=caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
                                    iteraxis='antenna', coloraxis=color, showgui=False, overwrite=True,
                                    plotfile=caltable.replace('.gcal', f'_amp_{color}_.png'), dpi=300, width=3000, height=1500,
                                )

                # Apply calibration tables after the last self-calibration loop
                if selfcal_loop == self.nloops - 1:
                    prev_caltables = sorted(glob.glob('*.gcal'))
                    logging.info("Applying the caltable derived from last gaincal iteration")
                    applycal(vis=self.msname, gaintable=prev_caltables, parang=False)

        if self.make_final_image:
            ## Generate a final mask of sources to peel -- optional
            imagename_final = self.imagename+'_final_clean'
            logging.info("Making final image with all selfcal corrections applied")
            
            self.niter = 1000000  # Can be modified to set a new value if needed
            # self.threshold = '0.001mJy'
            self.use_pybdsf = False

            # Initialize imager_instance with required parameters for self-calibration
            imager_instance = tclean_Imager(
                msname=self.msname,
                imagename=imagename_final,
                nterms=self.nterms,
                imsize=self.imsize,
                niter=self.niter,
                deconvolver=self.deconvolver,
                threshold=self.thresholds[selfcal_loop],
                mask=masking_file,
                cell = self.cell,
                weighting = self.weighting,
                robust = self.robust,
                overwrite=self.overwrite
            )
            imager_instance.imager()

            # Export dirty map to FITS after it has been created
            image_ext = '.image.tt0' if self.deconvolver == 'mtmfs' else '.image'
            exportfits(imagename=imagename_final + image_ext, fitsimage=imagename_final+ image_ext + '.fits', overwrite=True)

            # # Run PYBDSF if requested
            # try:
            #     logging.info(f"Running pybdsf on {imagename_final}...")
            #     if self.use_pybdsf:
            #         Utils.pybdsf(imagename_final+ image_ext, self.pybdsf_threshold)
            #         mask = imagename_final+ image_ext+ '.casabox'
            #         logging.info(f"Successfully ran pybdsf on {imagename_final}.")
            #     else:
            #         logging.info(f"Masking using PYBDSF not requested")
            #         mask = ''
            # except Exception as e:
            #     logging.error(f"Failed to run pybdsf on {imagename_final}: {e}")


    # def applycal_target():

    #     """
    #     Applies cal to the target field
    #     """
    #     prev_caltables = sorted(glob.glob('*.gcal'))
    #     applycal(
    #         vis = vis_tocal, gaintable = prev_caltables, parang=False
    #     )


class ImageProcessor:

    def __init__(self, imagename, msname, box_size=10, threshold_factor=5, max_iterations=None):
        """
        Initialize the ImageProcessor class.

        :param imagename: Name of the FITS file for the image (excluding file extension).
        :param msname: Name of the measurement set (visibility dataset) file.
        :param box_size: Half-width of the box around the peak for fitting.
        :param threshold_factor: Factor of the noise level to set as the source threshold.
        :param max_iterations: Maximum number of iterations to peel sources. 
                                If None, the peeling will stop when the flux falls below the threshold.
        """
        self.imagename = imagename
        self.msname = msname  # Set the measurement set name
        self.box_size = box_size
        self.threshold_factor = threshold_factor
        self.max_iterations = max_iterations  # Max iterations or None
        self.cl = casatools.componentlist()
    
    

    def peel_sources(self):
        """
        Main function to peel sources iteratively.

        ## NB: To run more than one iteration, you need to create an image for each iteration
        """
        # Get initial noise level from the image statistics
        stats = imstat(self.imagename)
        noise_level = stats['rms'][0]
        threshold = 5 * noise_level  # Threshold set to 5 times the RMS

        print(f"Using threshold: {threshold} Jy")  # Debug print

        ## call get_beam from plotting routines
        bmaj, bmin, pa = PlottingRoutines.get_beam(self)
        # print(bmaj,bmin,pa)

        ## Set region to ignore
        ignore_radius = 2 * bmaj

        # Get max_iterations value
        max_iterations = self.max_iterations
        print(f"Peeling will run {max_iterations} iterations")  

        # Iterate using a for loop
        for iteration in range(max_iterations):
            print(f"Starting iteration {iteration + 1}")  # Debug print

            # Get image statistics for max flux
            image_data = imstat(self.imagename)
            max_flux = image_data['max'][0]
            max_x, max_y = image_data['maxpos'][0], image_data['maxpos'][1]

            # If flux is below the threshold, exit the loop
            if max_flux < threshold:
                print(f"Flux is below threshold {threshold} Jy. Stopping after {iteration + 1} iterations.")
                break

            # Ignore the center of the image (a few beam factors)
            header = imhead(self.imagename)
            shape = header['shape']
            center_x, center_y = shape[0] / 2, shape[1] / 2

            if ((max_x - center_x) ** 2 + (max_y - center_y) ** 2) ** 0.5 < ignore_radius:
                print("Skipping central region of the image.")
                break

            # Define the box around the peak position
            xmin = max(0, max_x - self.box_size)
            xmax = max(0, max_x + self.box_size)
            ymin = max(0, max_y - self.box_size)
            ymax = max(0, max_y + self.box_size)

            imfit_box = f"{xmin},{ymin},{xmax},{ymax}"
        
            print(f"Automatically determined imfit box: {imfit_box}")

            # Fit a Gaussian to the source and extract position and flux
            fit_result = imfit(self.imagename, box=imfit_box)
            if 'component0' not in fit_result['results']:
                print("No fit found for this source, skipping.")
                break

            peak_flux = fit_result['results']['component0']['peak']['value']
            ra = fit_result['deconvolved']['component0']['shape']['direction']['m0']['value']
            dec = fit_result['deconvolved']['component0']['shape']['direction']['m1']['value']

            # Convert RA and Dec to standard strings
            sky_coord = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame='icrs')
            ra_hms = sky_coord.ra.to_string(unit=u.hour, sep=':')
            dec_dms = sky_coord.dec.to_string(unit=u.deg, sep='.', pad=True)

            print(f"Right Ascension (hms): {ra_hms}")
            print(f"Declination (dms): {dec_dms}")
            print(f"Peak flux: {peak_flux} Jy")

            # Define the component list filename using RA and Dec
            clname = f"{ra_hms}_{dec_dms}.cl"
            if os.path.exists(clname):
                os.system(f"rm -r {clname}")

            # Add component to the list with extracted RA, Dec, and flux
            self.cl.addcomponent(flux=peak_flux, fluxunit='Jy', shape='point', dir=f"J2000 {ra_hms} {dec_dms}")
            self.cl.rename(clname)
            self.cl.done()
            print(f"Component with flux {peak_flux} Jy at J2000 {ra_hms} {dec_dms} saved to {clname}")

            # Delete the original model column
            print("Deleting model column.")
            delmod(vis=self.msname, otf=True)
            plotms(vis=self.msname, xaxis='frequency', yaxis='amp', ydatacolumn='model', plotfile='empty_model.png', showgui=False, overwrite=True)

            # Add model to the MODEL column and subtract
            print(f"Adding {clname} to model column using ft")
            ft(vis=self.msname, complist=clname, incremental=False, usescratch=True)
            plotms(vis=self.msname, xaxis='frequency', yaxis='amp', ydatacolumn='model', plotfile='added_model.png', showgui=False, overwrite=True)

            # Perform UV subtraction
            print("Performing UVSUB to remove the source.")
            uvsub(vis=self.msname, reverse=False)

            print(f"Finished peeling source at {ra_hms}, {dec_dms} with flux {peak_flux} Jy.\n")

            # print(f"Making new image to use for peeling")
            # # Re-image to update `self.imagename`
            # new_imagename = f"{self.imagename}_iteration_{iteration + 1}"
            # WSClean_Imager(
            #     msname=self.msname, 
            #     imsize=320,
            #     niter=1, # requires at least one to get a clean beam
            #     use_pybdsf=False,
            #     pybdsf_threshold=5,
            #     overwrite=False,
            #     imagename=new_imagename,
            #     # cell='4.6arcsec'
            # ).imager()

            # Update `self.imagename` to the new image name
            # self.imagename = new_imagename + "-image.fits"
            # print(f"Updated imagename to {self.imagename}")







def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def configure_parameters():

    target = 'K2-18'
    phase_calibrator =''

    fieldnames = [name for name in [target, phase_calibrator] if name]

    """Configure the parameters for self-calibration."""
    return {
        'working_directory': Path('/raid1/scratch/kelvinw/k2_18b/selfcal'),
        'nloops': 6,
        'thresholds': ['', 4, 4,4,4,4 ],
        'calmode': ['', 'p','p','p','ap','ap'],
        'gaintype': ['', 'G' ,'G','G','G','G'],
        'solint': ['', '96s','48s','12s','192s','96s'],
        'minsnr': ['', 2, 2, 2, 2, 2],
        'avgtime': '',
        'width': 1,
        'fieldname':fieldnames,
        'outlierfile': '/raid1/scratch/kelvinw/casa_vlbi/selfcal/outlier.txt',
        # 'msname':'/raid1/scratch/kelvinw/k2_18b/selfcal/K2-18_split__1_phaseshifted.ms',
        # 'msname': '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/s_band_d_config/23B-307.sb44594812.eb44725045.60239.588568113424/K2-18.ms'  
        # 'msname': '/raid1/scratch/kelvinw/gv020_working_dir/gv020b_working_dir/gv020b_3.ms'
        # 'msname': '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/s_band_d_config/23B-307/pipeline.60623.88275462948/23B-307.sb44616223.eb44871184.60286.71989133102.ms'
        'msname': '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/x_band_d_config/23B-307/pipeline.60625.53603009274/23B-307.sb44672076.eb44857900.60279.378077800924.ms'
    }




def prepare_data(working_directory, msname, avgtime, width,fieldname):
    """Set the working directory and prepare the measurement set for multiple fields."""
    Utils.set_working_dir(working_directory)
    msnames = []  # To store the processed measurement sets for each field

    for field in fieldname:
        msname_split = MeasurementSetProcessor.prepare_measurement_set(
            msname=msname, 
            fieldname=field, 
            split_required= True, 
            avgtime=avgtime, 
            width=width
        )
        msnames.append(msname_split)

    return msnames

def perform_selfcalibration(vis, parameters):

    """Perform self-calibration with WSClean."""
    self_calibration_wsclean = SelfCalibrationWSClean(
        msname=vis, 
        nloops=parameters['nloops'],
        thresholds=parameters['thresholds'],
        calmode=parameters['calmode'],
        gaintype=parameters['gaintype'],
        solint=parameters['solint'],
        minsnr=parameters['minsnr'],
        imsize=640,
        niter=10000000,
        use_pybdsf=True,
        pybdsf_threshold=10,
        overwrite=False,
        final_image=True,
        # refant = 'EF' ,
        # cell = '3arcsec' 
    )
    # self_calibration_wsclean.selfcal()

    self_calibration_tclean = SelfCalibrationTclean(
        msname=vis,
        nloops=parameters['nloops'],
        thresholds=parameters['thresholds'],
        calmode=parameters['calmode'],
        gaintype=parameters['gaintype'],
        solint=parameters['solint'],
        minsnr=parameters['minsnr'],
        imsize=640,
        niter=1000,  
        nterms=2,
        deconvolver='mtmfs',
        weighting='briggs',
        robust=0.5,
        use_pybdsf=False, # leave as False -- mask from pybdsf is weird
        masking_threshold=8,
        pybdsf_threshold=5,
        overwrite=False,
        parallel = True,
        # outlierfile = parameters['outlierfile'],
        make_final_image=False,
        pbcorrect = False
    )
    self_calibration_tclean.selfcal()



def check_image_exists(imagename):
    """Check if the image already exists at the given path."""
    if os.path.exists(imagename):
        logging.info(f"Image found at {imagename}")
        return True
    else:
        logging.warning(f"Image {imagename} not found.")
        return False


def process_and_plot_image(vis, parameters):

    """Generate, process, and plot the self-calibration image."""
    # Generate imagename with full path
    imagename = f"{vis.replace('.ms', '')}_selfcal_loop_3"
    imagename = f"{parameters['working_directory']}/{imagename}-image.fits"

    # Check if image exists before proceeding
    if not check_image_exists(imagename):
        return

    # Initialize ImageProcessor
    processor = ImageProcessor(imagename=imagename, msname=vis, box_size=10, max_iterations=1)
    processor.peel_sources()
    peeled_image = imagename.replace('-image.fits', '_peeled')

    # Image the peeled measurement set
    peeled_source_image = WSClean_Imager(
        msname=vis, 
        imsize=320,
        niter=1,
        use_pybdsf=False,
        pybdsf_threshold=5,
        overwrite=False,
        imagename=peeled_image
    )
    peeled_source_image.imager()

    # Plotting and further processing
    # processor.imagename = peeled_image + '-image.fits'
    # plotter = PlottingRoutines(imagename=peeled_image + '-image.fits')
    # plotter.plot_image_with_beam()


def main():
    """Main entry point for the script."""
    setup_logging()

    # Configure parameters
    parameters = configure_parameters()

    # fieldname = parameters.get('fieldnames', ['target', 'phasecalibrator'])
    fieldname = parameters.get('fieldname') 
    msname = prepare_data(parameters['working_directory'], parameters['msname'], 
                          parameters['avgtime'], parameters['width'],fieldname=fieldname)

    vis = msname[0]

    # # Perform self-calibration
    perform_selfcalibration(vis, parameters)

    # # Process the image, peel sources, and plot
    # process_and_plot_image(vis, parameters)


if __name__ == "__main__":
    main()








        
        




import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import zipfile
import logging
import shutil
import glob
import sys
from astropy.io import fits
from natsort import natsorted

import os, glob, subprocess, time, bdsf, logging, math
from typing import Callable, Any
import casatools, casalogger
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


import os
import json
import time
import logging
import subprocess
from functools import wraps
from typing import Any, Callable

class Utils:
    state_file = 'calibration_state.json'

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
    def time_execution(func: Callable) -> Callable:
        """
        Decorator to measure and log the execution time for a function

        Parameters:
        - func (Callable): The function to time
        Returns:
        - Callable: The wrapped function with timing
        """
        @wraps(func)
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
    def run_wsclean(command: list) -> None:
        """
        Runs wsclean commands and logs the output and errors in real-time.
        """
        # Get the path to bind the singularity container
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname('wsclean_sif')))
        # Form the command to execute with Singularity
        command_to_execute = ['singularity', 'exec', '-B', singularity_bind, 'wsclean_sif'] + command
        logging.info(f"Executing: {' '.join(command_to_execute)}")
        
        try:
            process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

            # Read stdout and stderr as they are produced
            for stdout_line in iter(process.stdout.readline, ""):
                logging.info(stdout_line.strip())
            for stderr_line in iter(process.stderr.readline, ""):
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
    def load_state() -> dict:
        """Load the calibration state from the JSON file."""
        if not os.path.exists(Utils.state_file):
            return {"tables": {}}
        with open(Utils.state_file, 'r') as f:
            return json.load(f)

    @staticmethod
    def save_state(step: str, result: Any) -> None:
        """Save the result of a calibration step into the state file."""
        state = Utils.load_state()
        state["tables"][step] = result
        with open(Utils.state_file, 'w') as f:
            json.dump(state, f, indent=4)

    @staticmethod
    def state_tracking(func: Callable) -> Callable:
        """Decorator to track if a function's step has already been completed."""
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Determine the step name from the function name (or other logic)
            step_name = func.__name__
            
            # Load the state to check if this step has been completed
            state = Utils.load_state()
            
            # If this step has been completed, skip it
            if step_name in state["tables"]:
                logging.info(f"Skipping {step_name} since it has already been completed: {state['tables'][step_name]}")
                return state["tables"][step_name]
            
            # Otherwise, execute the function and save the state
            try:
                result = func(*args, **kwargs)
                Utils.save_state(step_name, result)
                return result
            except Exception as e:
                logging.error(f"Error during {step_name}: {e}")
                return None
        
        return wrapper





class AttachMetadata:
    helper_scripts: str = 'casa-vlbi-master.zip'
    helper_scripts_dir: str = 'casa-vlbi'
    repo_url: str = 'https://github.com/jive-vlbi/casa-vlbi/archive/refs/heads/master.zip'

    @staticmethod
    def download_helper_scripts() -> None:
        """Download and extract JIVE helper scripts if not already available."""
        if os.path.exists(AttachMetadata.helper_scripts_dir):
            logging.info(f"Helper scripts already available: {AttachMetadata.helper_scripts_dir}")
            return

        if not os.path.exists(AttachMetadata.helper_scripts):
            subprocess.run(['wget', '-c', AttachMetadata.repo_url, '-O', AttachMetadata.helper_scripts], check=True)
        
        with zipfile.ZipFile(AttachMetadata.helper_scripts, 'r') as zip_ref:
            zip_ref.extractall()
            logging.info("Downloaded and extracted helper scripts.")

        extracted_dir = AttachMetadata.helper_scripts.replace('.zip', '')
        if os.path.exists(extracted_dir):
            shutil.move(extracted_dir, AttachMetadata.helper_scripts_dir)
            logging.info("Helper scripts directory renamed to 'casa-vlbi'.")

    @staticmethod
    def import_casa_vlbi_tools() -> bool:
        """Import CASA VLBI tools from the downloaded helper scripts."""
        sys.path.append(AttachMetadata.helper_scripts_dir)
        try:
            global append_tsys, append_gc, convert_flags, convert_gaincurve
            from casavlbitools.fitsidi import append_tsys, append_gc, convert_flags
            from casavlbitools.casa import convert_gaincurve
            logging.info("Successfully imported CASA VLBI tools.")
            return True
        except ImportError as e:
            logging.error(f"Import error for CASA VLBI tools: {e}")
            return False

    @staticmethod
    def convert_flags(uvflg_file: str, fits_files: list, experiment: str) -> None:
        """Convert flags if a UVFLG file is provided."""
        if not os.path.exists(uvflg_file):
            logging.info("No UVFLG file found.")
            return

        try:
            convert_flags(infile=uvflg_file, idifiles=fits_files, outfile=f'{experiment}_apriori.flag')
            logging.info("Flag conversion completed.")
        except Exception as e:
            logging.error(f"Error during flag conversion: {e}")

    @staticmethod
    def remove_extensions(fits_files: list, extensions=('GAIN_CURVE', 'SYSTEM_TEMPERATURE')) -> None:
        """Remove specified extensions from the FITS files if present."""
        for filename in fits_files:
            with fits.open(filename, mode='update') as hdul:
                indices = [i for i, ext in enumerate(hdul) if ext.header.get('EXTNAME') in extensions]
                if indices:
                    for i in reversed(indices):
                        del hdul[i]
                    hdul.flush()
                    logging.info(f"Removed extensions {extensions} from {filename}.")
                else:
                    logging.info(f"No extensions {extensions} found in {filename}.")

    @staticmethod
    def attach_tsys(fits_files: list, antab_file: str) -> None:
        """Attach SYSTEM_TEMPERATURE (TSYS) table to FITS files if missing."""
        for filename in fits_files:
            with fits.open(filename) as hdul:
                if any(ext.header.get('EXTNAME') == 'SYSTEM_TEMPERATURE' for ext in hdul):
                    logging.info(f"SYSTEM_TEMPERATURE already exists in {filename}.")
                    continue

            try:
                append_tsys(antab_file, idifiles=filename)
                logging.info(f"Attached SYSTEM_TEMPERATURE to {filename}.")
            except Exception as e:
                logging.error(f"Error attaching TSYS to {filename}: {e}")

    @staticmethod
    def attach_gain_curve(fits_file: str, antab_file: str) -> None:
        """Attach GAIN_CURVE table to the first FITS file in the series."""
        try:
            append_gc(antab_file, fits_file)
            logging.info(f"Attached GAIN_CURVE to {fits_file}.")
        except Exception as e:
            logging.error(f"Error attaching GAIN_CURVE to {fits_file}: {e}")

    @staticmethod
    def convert_gaincurve(antab_file: str, experiment: str) -> None:
        """Create a gain curve table for the experiment."""
        gc_table = f'{experiment}.gc'
        if os.path.exists(gc_table):
            os.remove(gc_table)

        try:
            convert_gaincurve(antab_file, gc_table, min_elevation=0.0, max_elevation=90.0)
            logging.info(f"Generated gain curve table: {gc_table}")
        except Exception as e:
            logging.error(f"Error generating gain curve table: {e}")

    @staticmethod
    def check_missing_tsys(fits_files: list) -> None:
        """Log any FITS files missing the SYSTEM_TEMPERATURE extension."""
        missing = [
            filename for filename in fits_files
            if not any(fits.getheader(filename, ext).get('EXTNAME') == 'SYSTEM_TEMPERATURE' for ext in range(len(fits.open(filename))))
        ]
        if missing:
            logging.warning(f"Files missing SYSTEM_TEMPERATURE: {missing}")
        else:
            logging.info("All files contain SYSTEM_TEMPERATURE.")

    @staticmethod
    @Utils.time_execution
    def attach_tsys_gc(experiment: str, antab_file: str, idifitsfiles: str, uvflg_file: str) -> None:
        """Main method to orchestrate TSYS and GAIN_CURVE table attachment."""
        AttachMetadata.download_helper_scripts()
        if not AttachMetadata.import_casa_vlbi_tools():
            return

        fits_files = natsorted(glob.glob(os.path.join(idifitsfiles, f'{experiment}_1.IDI*')))

        AttachMetadata.convert_flags(uvflg_file, fits_files, experiment)
        AttachMetadata.remove_extensions(fits_files)
        AttachMetadata.attach_tsys(fits_files, antab_file)
        AttachMetadata.attach_gain_curve(fits_files[0], antab_file)
        AttachMetadata.convert_gaincurve(antab_file, experiment)
        AttachMetadata.check_missing_tsys(fits_files)


if __name__ == "__main__":
    experiment = 'gv020b_3'
    idifitsfiles = '/raid1/scratch/kelvinw/gv020_fitsfiles/gv020b_fitsfiles'
    antab_file = '/raid1/scratch/kelvinw/casa_vlbi/data/gv020b_1.antab'
    uvflg_file = '/raid1/scratch/kelvinw/casa_vlbi/data/gv020b_1.uvflg'  # Leave as an empty string if not used

    AttachMetadata.attach_tsys_gc(experiment, antab_file, idifitsfiles, uvflg_file)

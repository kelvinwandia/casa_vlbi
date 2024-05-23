import os, glob, subprocess, json, re,sys, logging, traceback
from natsort import natsorted
from datetime import datetime
import numpy as np
from AIPS import AIPS
from AIPSTV import AIPSTV
from AIPSData import AIPSUVData
from AIPSTask import AIPSTask
from Wizardry.AIPSData import AIPSImage


import configparser
config = configparser.ConfigParser()
config.read('aips_config.ini')

exec(open("./aips_functions.py").read())
exec(open("../utils/helper_functions.py").read())
exec(open("./casa_functions.py").read())

print("Creating log dir")
log_dir = os.path.join(os.getcwd(), 'logs')

if not os.path.exists(log_dir):
    print(f"Log directory '{log_dir}' does not exist, creating it...")
    os.makedirs(log_dir, exist_ok=True)
else:
    print(f"Log directory '{log_dir}' exists")

try:
    console
except:
    logfile_name = datetime.now().strftime('parseltongue_%H:%M:%S_%d:%m:%Y.log')
    filename = os.path.join(log_dir, logfile_name)
    
    logging.basicConfig(filename=filename, level=logging.DEBUG)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)



# import logging
# from datetime import datetime

# class ColorizedArgsFormatter(logging.Formatter):
#     arg_colors = {
#         logging.DEBUG: "\x1b[38;21m",    # Grey
#         logging.INFO: "\x1b[1;32m",       # Green
#         logging.WARNING: "\x1b[33;21m",    # Yellow
#         logging.ERROR: "\x1b[31;21m",      # Red
#         logging.CRITICAL: "\x1b[31;1m",    # Bold Red
#     }
#     reset = "\x1b[0m"  # Reset color

#     def format(self, record):
#         level_color = self.arg_colors.get(record.levelno, self.reset)
#         msg = super().format(record)
#         return f"{level_color}{msg}{self.reset}"

# # Configure logging
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)

# # Create console handler and set formatter
# ch = logging.StreamHandler()
# ch.setFormatter(ColorizedArgsFormatter())

# # Create file handler and set formatter
# logfile_name = datetime.now().strftime('parseltongue_%H_%M_%S_%d_%m_%Y.log')
# fh = logging.FileHandler(filename=logfile_name)
# fh.setFormatter(ColorizedArgsFormatter())

# # Add handlers to logger
# logger.addHandler(ch)
# logger.addHandler(fh)


# globals
userno = config.getint('globals','set_user')
experiment = config.get('globals','experiment_name')
working_dir = config.get('globals','working_directory')
start_tv = config.getboolean('globals','start_tv')
indisk = config.getint('globals','indisk')
inseq = config.getint('globals','inseq')
inclass = config.get('globals','inclass')

# load data
load_data = config.getboolean('load_data','load_data')
zap_data = config.getboolean('load_data','zap_data')
fitsfiles_dir = config.get('load_data','fitsfiles_directory')
file_extension_for_flagging = config.get('load_data','file_extension_for_flagging')
file_extension_for_cal = config.get('load_data','file_extension_for_cal')
tasav_file = config.get('load_data','tasav_file')
pointing = config.get('load_data','pointing')
integration_time=config.getfloat('load_data','integration_time')

# sources
phase_calibrator = config.get('sources','phase_calibrator')
fringe_finder = config.get('sources','fringe_finder')
target = config.get('sources','target')

# flagging

do_flagging = config.getboolean('flagging','do_flagging')
aoflagger_sif = config.get('flagging','aoflagger_sif')
flagging_strategy = config.get('flagging','flagging_strategy')

# calibrate data
refant = config.get('calibrate','refant')
searchants=config.get('calibrate','searchants').split(',')
fring_solint = config.getfloat('calibrate','fring_solint')
fring_timerange = config.get('calibrate','fring_timerange')
do_amp_parang_correction = config.getboolean('calibrate','do_amp_parang_correction')
do_tec_correction = config.getboolean('calibrate','do_tec_correction')
do_singleband_fring = config.getboolean('calibrate','do_singleband_fring')
apply_singleband_corrections = config.getboolean('calibrate','apply_singleband_corrections')

do_global_fring = config.getboolean('calibrate','do_global_fring')
fring_snr = config.getfloat('calibrate','fring_snr')
apply_global_fring_corrections = config.getboolean('calibrate','apply_global_fring_corrections')

do_bandpass = config.getboolean('calibrate','do_bandpass')
apply_bpass_corrections = config.getboolean('calibrate','apply_bpass_corrections')

apply_all_calibrations = config.getboolean('calibrate','apply_all_calibrations')


# split
do_splat = config.getboolean('split','do_splat')
write_fits = config.getboolean('split','write_fits')



logging.info(f"AIPS user no.:{userno}")
AIPS.userno = userno
if AIPS.userno == 0:
    raise ValueError("Please set AIPS userno")


# logging.info("Setting indata")
# set_indata()
# logging.info("Indata set")

if not os.path.exists(working_dir):
    logging.info(f"{working_dir} does not exist, making one")
    try:
        set_working_dir(working_dir)
    except Exception as e:
        logging.error(f"An error occured while creating the working directory: {e}")
else:
    logging.info(f"Working directory {working_dir} already exists")
    os.chdir(working_dir)

if start_tv==True:
    try:
        TV(start_tv)
    except Exception as e:
        logging.error(f"An error occurred while controlling TV: {e}")

if zap_data==True:
    try:
        cleanup()
    except Exception as e:
        logging.error(f"{e}")

if load_data==True:
    try:
        logging.info("Loading fitsfiles")
        load_fitsfiles(file_extension_for_flagging)
        load_tasav()
    except Exception as e:
        logging.error(f"An error occurred while loading fitsfiles: {e}")

if do_flagging == True:

    """
    Check the data sorting here, export to UVFITS, convert to ms, flag and then load the data again to AIPS for cal
    """
    unflagged_fitsfile = experiment+'.'+file_extension_for_flagging
    output_unflagged_file_ext = 'FITS'
    flagged_fitfile = f"{experiment}_{pointing}_1.{file_extension_for_cal}"
    vis= experiment+'.ms'
    try:
        if not os.path.exists(unflagged_fitsfile):
            runfittp(output_unflagged_file_ext) 
        else:
            logging.info(f"UVFITS file {unflagged_fitsfile} exists. New one will not be written")
    except Exception as e:
        logging.info(f"An error occurred: {e}")
    
    try:    
        logging.info(f"Making measurement {vis} for flagging")
        makems(vis,unflagged_fitsfile)
    except Exception as e:
        logging.info(f"An error {e} occured")
    try:
        logging.info("Flagging data")
        # plot_check_baddata(save_as="_before_flagging")
        # execute_aoflagger_strategy()
        # plot_check_baddata(save_as="_after_flagging")
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

    try:
        logging.info(f"Exporting flagged {vis} to {flagged_fitfile} ")
        if not os.path.exists(flagged_fitfile):
            makeuvfits(vis,flagged_fitfile)
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

    try:
        logging.info(f"Loading the flagged data")
        load_fitsfiles(file_extension_for_cal)
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

if do_amp_parang_correction == True:
    try:  
        logging.info("Setting indata")
        set_indata()
        logging.info("Indata set")
        logging.info("Copying CL2 from TASAV file to UVDATA")
        runtacop()
    except Exception as e:
        logging.error(f"An error occurred while copying from TASAV: {e}")

if do_tec_correction == True:
    try: 
        logging.info("Performing ionospheric corrections")
        runtecor()
    except Exception as e:
        logging.error(f"An error occurred: {e}")

if do_singleband_fring == True:
    try:
        logging.info("Running single band fring using")
        fring_instr(phase_calibrator)
    except Exception as e:
        logging.error(f"An error occurred: {e}")

if apply_singleband_corrections == True:
    try:
        logging.info("Running CLCAL")
        apply_solutions(phase_calibrator)
    except Exception as e:
        logging.info(f"An error occurred: {e}")

if do_global_fring == True:
    try:
        logging.info("Running global fring")
        global_fring(phase_calibrator)
    except Exception as e:
        logging.info(f"An error occurred: {e}")

if apply_global_fring_corrections == True:
    try:
        logging.info("Running CLCAL")
        apply_solutions(phase_calibrator)
    except Exception as e:
        logging.info(f"An error occurred: {e}")

if do_bandpass == True:
    try:
        logging.info("Running bandpass")
        runbpass(phase_calibrator)
    except Exception as e:
        logging.info(f"An error occurred: {e}")    


if do_splat == True:
    try:
        sources = target+phase_calibrator+fringe_finder
        logging.info(f"Running SPLAT for {sources}")
        runsplat(target,phase_calibrator,fringe_finder)
    except Exception as e:
        logging.info(f"An error occurred: {e}")    

    
if write_fits == True:
    try:
        runfittp()
    except Exception as e:
        logging.info(f"An error occurred: {e}")
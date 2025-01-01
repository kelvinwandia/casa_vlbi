import os, glob, re, logging
from datetime import datetime
import casatasks, casatools
import casaplotms
import numpy as np
import subprocess
import matplotlib
# matplotlib.use('Agg')  
import time
from natsort import natsorted
import zipfile
import shutil


"""
RUNNING THE SINGULARITY CONTAINER
singularity exec ../singularity/casa6_wsclean2_aoflagger3.sif env XDG_RUNTIME_DIR=$(dirname "$0") xvfb-run --auto-servernum python main.py

"""

msmd = casatools.msmetadata()
tb = casatools.table()


# import configparser
# # Load the configuration file
# current_dir = os.path.dirname(os.path.abspath(__file__))
# print(current_dir)
# config_path = os.path.join(current_dir, 'config.ini')
# if not os.path.exists(config_path):
#     raise FileNotFoundError(f"Configuration file not found: {config_path}")

# config = configparser.ConfigParser()
# config.read(config_path)

import sys

current_dir = os.path.dirname(os.path.abspath(__file__))

# Add the parent directory to sys.path
sys.path.append(current_dir)
sys.path.append(os.path.join(current_dir, 'utils'))
sys.path.append(os.path.join(current_dir, 'calibration'))
sys.path.append(os.path.join(current_dir, 'selfcal'))
sys.path.append(os.path.join(current_dir, 'data'))








### This is important for running the sif container
### it allows the logfiles to be placed in the same dir as main.py
### do not touch !!
# script_dir = os.path.dirname(os.path.abspath(__file__))
# os.environ['XDG_RUNTIME_DIR'] = script_dir

logging.info(f"Current working directory: {os.getcwd()}")

print("Creating log dir")
log_dir = os.path.join(os.getcwd(), 'logs')

if not os.path.exists(log_dir):
    print(f"Log directory '{log_dir}' does not exist, creating it...")
    os.makedirs(log_dir, exist_ok=True)
else:
    print(f"Log directory '{log_dir}' exists")

logging.info(f"Current working directory: {os.getcwd()}")


logfile_name = datetime.now().strftime('casa_vlbi_%H_%M_%S_%d_%m_%Y.log')  # Replace colons with underscores
filename = os.path.join(log_dir, logfile_name)
logging.basicConfig(filename=filename, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')




# # globals
# load_data = config.getboolean('globals','load_data')
# experiment = config.get('globals','experiment_name')
# working_directory = config.get('globals', 'working_directory')
# uvfits_file = config.get('globals','uvfits_file')
# aoflagger_path = config.get('globals','aoflagger_path')
# use_singularity = config.getboolean('globals','use_singularity')
# # singularity_bind = config.get('globals','singularity_bind')
# wsclean_sif = config.get('globals','wsclean_sif')
# singularity_container = config.get('globals','singularity_container')

# # apriori
# use_casa = config.getboolean('apriori','use_casa')
# attach_metadata = config.getboolean('apriori','attach_metadata')
# idifitsfiles = config.get('apriori','idifitsfiles')
# antab_file = config.get('apriori','antab_file')
# uvflg_file = config.get('apriori','uvflg_file')
# do_apriori_cal = config.getboolean('apriori','do_apriori_cal')
# apply_apriori_cal = config.getboolean('apriori','apply_apriori_cal')

# target = config.get('basic','target')
# phase_calibrator = config.get('basic','phase_calibrator')
# fringe_finder = config.get('basic','fringe_finder')
# integration_time = config.getfloat('basic','integration_time')
# do_split = config.getboolean('basic','do_split')
# timebin = config.get('basic','timebin')
# width = config.getint('basic','width')
# verbosity = config.getboolean('basic','verbosity')

# # flag
# do_flagging = config.getboolean('flagging','do_flagging')
# # aoflagger_strategies = config.get('flagging','aoflagger_strategies')
# # faint_source_strategy = config.get('flagging','faint_source_strategy')
# # bright_source_strategy = config.get('flagging','bright_source_strategy')
# manual_file = config.get('flagging','manual_file')
# edge_channel_fraction = config.getfloat('flagging','edge_channel_fraction')
# use_aoflagger = config.getboolean('flagging','use_aoflagger')
# flagging_strategy = config.get('flagging','flagging_strategy')

# flag_antenna = config.getboolean('flagging','flag_antenna')
# antenna_to_flag = config.get('flagging','antenna_to_flag')
# export_uvfits = config.getboolean('flagging','export_uvfits')


# # calibrate
# do_tec_corrections = config.getboolean('calibrate','do_tec_corrections')
# do_sbd_fringe = config.getboolean('calibrate','do_sbd_fringe')
# apply_sbd = config.getboolean('calibrate','apply_sbd')
# refant = config.get('calibrate','refant')
# timerange = config.get('calibrate','timerange')
# timerange_ar = config.get('calibrate','timerange_ar')
# snr_sbd = config.getfloat('calibrate','snr_sbd')
# do_mbd_fringe = config.getboolean('calibrate','do_mbd_fringe')
# apply_mbd = config.getboolean('calibrate','apply_mbd')
# snr_mbd = config.getfloat('calibrate','snr_mbd')
# solint = config.getfloat('calibrate','solint')
# do_bpass = config.getboolean('calibrate','do_bpass')
# apply_bpass = config.getboolean('calibrate','apply_bpass')

# make_dirty_map = config.getboolean('calibrate','make_dirty_map')
# split_calibrated = config.getboolean('calibrate','split_calibrated')
# imsize= [int(part) for part in config.get('calibrate', 'imsize').split(',')]
# detection_threshold = config.getfloat('selfcal','detection_threshold')

# # selfcal
# do_selfcal = config.getboolean('selfcal','do_selfcal')
# use_tclean = config.getboolean('selfcal','use_tclean')
# use_wsclean = config.getboolean('selfcal','use_wsclean')

# pybdsf_threshold = config.get('selfcal','pybdsf_threshold')
# pybdsf_niter = config.getint('selfcal','pybdsf_niter')

# weighting = config.get('selfcal','weighting')
# robust = config.getfloat('selfcal','robust')

# nloops = config.getint('selfcal','nloops')
# calmode = config.get('selfcal','calmode').split(',')
# gaintype = config.get('selfcal','gaintype').split(',')
# cell =  config.get('selfcal', 'cell')
# threshold = config.get('selfcal','threshold').split(',')
# minsnr = [float(part) for part in config.get('selfcal', 'minsnr').split(',')]
# # imsize= [int(part) for part in config.get('selfcal', 'imsize').split(',')]
# # niter = [int(part) for part in config.get('selfcal', 'niter').split(',')]
# niter = int(config.get('selfcal','niter'))
# niter_final = int(config.getint('selfcal','niter_final'))
# threshold_final = config.getint('selfcal','threshold_final')
# robust = config.getfloat('selfcal','robust')

# tclean_threshold = config.get('selfcal','tclean_threshold')
# solint_selfcal = config.get('selfcal','solint_selfcal').split(',')
# apply_to_target = config.getboolean('selfcal','apply_to_target')
# detect_sources = config.getboolean('selfcal','detect_sources')

# # pb corrections
# do_pbcor = config.getboolean('pbcor','do_pbcor')
# pb_file = config.get('pbcor','pb_file')



# Import the modules
from utils.helper_functions import *
from selfcal import *
from calibration.calibrate import *
from calibration.pbcor import *

from config_file import *

# from config_file import configuration_file
# # Print all sections and options in the config
# for section in configuration_file.sections():
#     print(f"Section: {section}")
#     for option in configuration_file.options(section):
#         print(f"  {option} = {configuration_file.get(section, option)}")

if not os.path.exists(working_directory):
    print(working_directory)
    logging.info(f"{working_directory} does not exist, making one")
    try:
        set_working_dir()
    except Exception as e:
        logging.error(f"An error occured while creating the working directory: {e}")
else:
    logging.info(f"Working directory {working_directory} already exists")
    os.chdir(working_directory)

if load_data == True:
    if use_casa == True and attach_metadata == True:
        logging.info(f"Attaching tsys and gc to fitsfiles")
        attach_tsys_gc()
    try:
        vis = experiment + '.ms'
        splitvis = None
        logging.info("Making measurement set")
        makems(vis)

        if do_split:
            if timebin == '':
                splitvis = vis.replace('.ms', f'_split_{width}_chan.ms')
            else:
                splitvis = vis.replace('.ms', f'_split_{timebin}_{width}_chan.ms')
            makems(vis,splitvis)
            vis = splitvis

        logging.info("Getting fields")
        getfields()
        logging.info("getfields completed successfully")
        print(f"The current vis is {vis}")
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

if do_flagging == True:
    try:
        logging.info("Flagging data")
        # plot_check_baddata(save_as="_before_flagging")
        if use_aoflagger == True:
            execute_aoflagger_strategy()
        flagging()
        # flag_edge_channels()
        if flag_antenna == True:
            antenna_flag(antenna_to_flag)
        # plot_check_baddata(save_as="_after_flagging")
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

if use_casa == True:
    logging.info(f"Using CASA to do apriori cal")
    if do_apriori_cal == True:
        try:
            logging.info("Doing amplitude calibration using TSYS and GC tables")
            gencal_tsys_gc()
        except Exception as e:
            logging.warning(f"Encountered error {e}")

    if apply_apriori_cal == True:
        try:
            logging.info("Doing amplitude calibration using TSYS and GC tables")
            applycal_tsys_gc()
        except Exception as e:
            logging.warning(f"Encountered error {e}")

if export_uvfits == True:
    try:
        export_to_uvfits(vis)
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_tec_corrections == True:
    try:
        logging.info(f"Calculating ionospheric corrections")
        tec_corrections()
    except Exception as e:
        logging.warning(f"Encountered error {e}")


if do_sbd_fringe == True:
    try:
        logging.info(f"Calculating instrumental delay corrections")
        sbd_fringefit()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_sbd == True:
    try:
        logging.info("Applying instrumental delay corrections")
        applycal_sbd_fringe()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_mbd_fringe == True:
    try:
        logging.info("Running multiband fringefit")
        mbd_fringefit()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_mbd == True:
    try:
        logging.info("Applying multiband corrections")
        applycal_mbd_fringe()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_bpass == True:
    print("Executing this section")
    try:
        logging.info("Calculating bandpass solutions")
        bpass()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_bpass == True:
    try:
        logging.info("Applying bandpass corrections")
        applycal_bpass()
        after_cal_plots()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

# if make_dirty_map == True:
#     try:
#         logging.info("Making dirty map")
#         dirty_map(target)
#         dirty_map(phase_calibrator)
#     except Exception as e:
#         logging.warning(f"Encountered error {e}")

if split_calibrated == True:
    try:
        logging.info("Making dirty map")
        split_calibrated_ms(phase_calibrator,target)
    except Exception as e:
        logging.warning(f"Encountered error {e}")


if do_selfcal == True:
    try:
        logging.info("Self calibrating the data")
        
        selfcal_dir = os.path.join(working_directory,'selfcal_dir')
        logging.info(f"Making and switching to {selfcal_dir}")
        if not os.path.exists(selfcal_dir):
            os.makedirs(selfcal_dir)
        os.chdir(selfcal_dir)
        split_selfcal()
        selfcal_part1()
        selfcal_part2()
        os.chdir(working_directory)
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_to_target == True:
    selfcal_dir = os.path.join(working_directory,'selfcal_dir')
    try:
        logging.info(f"Switching to {selfcal_dir}")
        os.chdir(selfcal_dir)
        logging.info("Applying calibrations to science target")
        applycal_target()
        os.chdir(working_directory)
    except Exception as e:
        logging.warning(f"Encountered error {e}")

# if detect_sources == True:

#     selfcal_dir = os.path.join(working_directory,'selfcal_dir')
#     target_ms = os.path.join(selfcal_dir,target+'.ms')
#     detected_sources = os.path.join(working_directory,'detected_sources')
#     logging.info(f"Making and switching to {detected_sources}")
#     if not os.path.exists(detected_sources):
#         os.makedirs(detected_sources)
 
#     try:
#         logging.info(f"Switching to {detected_sources}")
#         os.chdir(detected_sources)
#         logging.info("Detecting sources")
#         m15_sources()
#         os.chdir(working_directory)

#     except Exception as e:
#         logging.warning(f"Encountered error {e}")

if do_pbcor == True:
    try:
        logging.info(f"You are about to perform pb corrections")
        pbcor_dir = os.path.join(working_directory,'pbcor_dir')
        logging.info(f"Making and switching to {pbcor_dir}")
        selfcal_dir = os.path.join(working_directory,'selfcal_dir')
        target_ms = os.path.join(selfcal_dir,target+'.ms')
        if not os.path.exists(pbcor_dir):
            os.makedirs(pbcor_dir)
        os.chdir(pbcor_dir)
        gencal_pb_table()
        os.chdir(working_directory)
    except Exception as e:
        logging.warning(f"Encountered error {e}")
import os, glob, subprocess, json, re,sys, logging, traceback
from natsort import natsorted
from datetime import datetime
import numpy as np
from  AIPS  import AIPS
from AIPSTV import AIPSTV
from AIPSData import AIPSUVData
from AIPSTask import AIPSTask
from Wizardry.AIPSData import AIPSImage


import configparser
config = configparser.ConfigParser()
config.read('aips_config.ini')

exec(open("./aips_functions.py").read())
# exec(open("../utils/helper_functions.py").read())
exec(open("./casa_functions.py").read())
exec(open("./selfcal.py").read())

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



# globals
userno = config.getint('globals','set_user')
experiment = config.get('globals','experiment_name')
working_dir = config.get('globals','working_directory')
start_tv = config.getboolean('globals','start_tv')
indisk = config.getint('globals','indisk')
inseq = config.getint('globals','inseq')
inclass = config.get('globals','inclass')

# load data
load_data = config.getboolean('process_data','load_data')
zap_data = config.getboolean('process_data','zap_data')
fitsfiles_dir = config.get('process_data','fitsfiles_directory')
file_extension =  config.get('process_data','file_extension')
make_antab = config.getboolean('process_data','make_antab')
export_data = config.getboolean('process_data','export_data')

# file_extension_for_flagging = config.get('load_data','file_extension_for_flagging')
# file_extension_for_cal = config.get('load_data','file_extension_for_cal')
# tasav_file = config.get('load_data','tasav_file')
# pointing = config.get('load_data','pointing')
# integration_time=config.getfloat('load_data','integration_time')

# # sources
# phase_calibrator = config.get('sources','phase_calibrator').split(',')
# fringe_finder = config.get('sources','fringe_finder').split(',')
# target = config.get('sources','target').split(',')

# # flagging

# do_flagging = config.getboolean('flagging','do_flagging')
# use_aoflagger = config.getboolean('flagging','use_aoflagger')
# aoflagger_sif = config.get('flagging','aoflagger_sif')
# flagging_strategy = config.get('flagging','flagging_strategy')
# manual_file = config.get('flagging','manual_file')

# # calibrate data
# refant = config.get('calibrate','refant')
# searchants=config.get('calibrate','searchants').split(',')
# fring_solint = config.getfloat('calibrate','fring_solint')
# fring_timerange = config.get('calibrate','fring_timerange')
# fring_timerange_ar_down = config.get('calibrate','fring_timerange_ar_down')
# fring_timerange_ar_up = config.get('calibrate','fring_timerange_ar_up')
# do_amp_parang_correction = config.getboolean('calibrate','do_amp_parang_correction')
# do_tec_correction = config.getboolean('calibrate','do_tec_correction')
# do_singleband_fring = config.getboolean('calibrate','do_singleband_fring')
# apply_singleband_corrections = config.getboolean('calibrate','apply_singleband_corrections')

# do_global_fring = config.getboolean('calibrate','do_global_fring')
# fring_snr = config.getfloat('calibrate','fring_snr')
# apply_global_fring_corrections = config.getboolean('calibrate','apply_global_fring_corrections')

# do_bandpass = config.getboolean('calibrate','do_bandpass')
# apply_bpass_corrections = config.getboolean('calibrate','apply_bpass_corrections')



# # split
# do_splat = config.getboolean('split','do_splat')
# write_fits = config.getboolean('split','write_fits')

# # selfcal
# split_selfcal = config.getboolean('selfcal','split_selfcal')
# make_dirty_map = config.getboolean('selfcal','make_dirty_map')
# do_selfcal = config.getboolean('selfcal','do_selfcal')
# use_tclean = config.getboolean('selfcal','use_tclean')
# use_wsclean = config.getboolean('selfcal','use_wsclean')
# wsclean_sif = config.get('selfcal','wsclean_sif')

# pybdsf_threshold = config.get('selfcal','pybdsf_threshold')
# pybdsf_niter = config.getint('selfcal','pybdsf_niter')
# imsize= [int(part) for part in config.get('selfcal', 'imsize').split(',')]
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
# niter = config.get('selfcal','niter')
# niter_final = config.getint('selfcal','niter_final')
# threshold_final = config.getint('selfcal','threshold_final')
# robust = config.getfloat('selfcal','robust')
# detection_threshold = config.getfloat('selfcal','detection_threshold')
# tclean_threshold = config.get('selfcal','tclean_threshold')
# solint_selfcal = config.get('selfcal','solint_selfcal').split(',')
# apply_to_target = config.getboolean('selfcal','apply_to_target')
# detect_sources = config.getboolean('selfcal','detect_sources')


logging.info(f"AIPS user no.:{userno}")
AIPS.userno = userno
if AIPS.userno == 0:
    raise ValueError("Please set AIPS userno")


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

if load_data==True:
    try:
        logging.info("Loading fitsfiles")
        load_fitsfiles(fitsfiles_dir,file_extension)
        # load_tasav()
    except Exception as e:
        logging.error(f"An error occurred while loading fitsfiles: {e}")

if make_antab==True:
    try:
        logging.info("Making antenna tables")
        set_indata()
        create_antab_file(working_dir)
    except Exception as e:
        logging.error(f"An error occurred while making antenna tables {e}")

if export_data == True:
    try:
        logging.info("Exporting data to .UVFITS")
        runfittp()
    except Exception as e:
        logging.error(f"An error occurred while exporting data {e}")


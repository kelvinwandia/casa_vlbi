import os, glob, re, logging
from datetime import datetime
import casatasks, casatools
import casaplotms
import numpy as np
import subprocess
import matplotlib
# matplotlib.use('Agg')  


import configparser
config = configparser.ConfigParser()
config.read('config_casa.ini')

logfile_name = datetime.now().strftime('casa_cal_%H_%M_%S_%d_%m_%Y.log')
logging.basicConfig(filename=logfile_name,level=logging.DEBUG)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)


load_data = config.getboolean('globals','load_data')
experiment = config.get('globals','experiment_name')
working_directory = config.get('globals', 'working_directory')
uvfits_file = config.get('globals','uvfits_file')
singularity_path = config.get('globals','singularity_path')
use_singularity = config.getboolean('globals','use_singularity')


target = config.get('basic','target')
phase_calibrator = config.get('basic','phase_calibrator')
fringe_finder = config.get('basic','fringe_finder')
integration_time = config.getfloat('basic','integration_time')
do_split = config.getboolean('basic','do_split')
timebin = config.get('basic','timebin')
width = config.getint('basic','width')

# flag
do_flagging = config.getboolean('flagging','do_flagging')
aoflagger_sif = config.get('flagging','aoflagger_sif')
# aoflagger_strategies = config.get('flagging','aoflagger_strategies')
faint_source_strategy = config.get('flagging','faint_source_strategy')
bright_source_strategy = config.get('flagging','bright_source_strategy')
manual_file = config.get('flagging','manual_file')

# calibrate
do_sbd_fringe = config.getboolean('calibrate','do_sbd_fringe')
apply_sbd = config.getboolean('calibrate','apply_sbd')
refant = config.get('calibrate','refant')
timerange = config.get('calibrate','timerange')
snr_sbd = config.getfloat('calibrate','snr_sbd')
do_mbd_fringe = config.getboolean('calibrate','do_mbd_fringe')
apply_mbd = config.getboolean('calibrate','apply_mbd')
snr_mbd = config.getfloat('calibrate','snr_mbd')
solint = config.getfloat('calibrate','solint')
do_bpass = config.getboolean('calibrate','do_bpass')
apply_bpass = config.getboolean('calibrate','apply_bpass')

make_dirty_map = config.getboolean('calibrate','make_dirty_map')
imsize = config.get('calibrate','imsize')
detection_threshold = config.getfloat('selfcal','detection_threshold')

# selfcal
do_selfcal = config.getboolean('selfcal','do_selfcal')


weighting = config.get('selfcal','weighting')
robust = config.getfloat('selfcal','robust')

nloops = config.getint('selfcal','nloops')
calmode = config.get('selfcal','calmode').split(',')
gaintype = config.get('selfcal','gaintype').split(',')
cell =  config.get('selfcal', 'cell')
threshold = config.get('selfcal','threshold').split(',')
minsnr = [float(part) for part in config.get('selfcal', 'minsnr').split(',')]
imsize= [int(part) for part in config.get('selfcal', 'imsize').split(',')]
niter = [int(part) for part in config.get('selfcal', 'niter').split(',')]

niter_final = config.getint('selfcal','niter_final')
threshold_final = config.get('selfcal','threshold_final')
solint_selfcal = config.get('selfcal','solint_selfcal').split(',')
apply_to_target = config.getboolean('selfcal','apply_to_target')

# pb corrections
do_pbcor = config.getboolean('pbcor','do_pbcor')

exec(open(".calibration/calibrate.py").read())
exec(open(".selfcal/selfcal.py").read())
exec(open(".calibration/pbcor.py"))


if not os.path.exists(working_directory):
    logging.info(f"{working_directory} does not exist, making one")
    try:
        set_working_dir(working_directory)
    except Exception as e:
        logging.error(f"An error occured while creating the working directory: {e}")
else:
    logging.info(f"Working directory {working_directory} already exists")
    os.chdir(working_directory)

if load_data == True:
    try:
        vis = experiment + '.ms'
        splitvis = None
        logging.info("Running CASA task importuvfits")
        makems(vis)

        if do_split:
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
        flagging()
        run_aoflagger_container()
    except Exception as e:
        logging.critical(f"Exception {e} occurred")

if do_sbd_fringe == True:
    try:
        logging.info(f"Deriving corrections for instrumental")
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
        logging.info("Apply multiband corrections")
        applycal_mbd_fringe()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_bpass == True:
    try:
        logging.info("Calculating bandpass solutions")
        bpass()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_bpass == True:
    try:
        logging.info("Applying bandpass corrections")
        applycal_bpass()
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if make_dirty_map == True:
    try:
        logging.info("Imaging")
        mytclean(phase_calibrator,niter=0)
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_selfcal == True:
    try:
        logging.info("Self calibrating the data")
        # selfcal_part1(phase_calibrator)
        selfcal_part2(phase_calibrator)
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if apply_to_target == True:
    try:
        logging.info("Applying calibrations to science targe")
        pass
    except Exception as e:
        logging.warning(f"Encountered error {e}")

if do_pbcor == True:
    pass

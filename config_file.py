# config.py
import configparser
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(current_dir, 'config.ini')

if not os.path.exists(config_path):
    raise FileNotFoundError(f"Configuration file not found: {config_path}")

configuration_file = configparser.ConfigParser()
configuration_file.read(config_path)


# globals
load_data = configuration_file.getboolean('globals','load_data')
experiment = configuration_file.get('globals','experiment_name')
working_directory = configuration_file.get('globals', 'working_directory')
uvfits_file = configuration_file.get('globals','uvfits_file')
aoflagger_path = configuration_file.get('globals','aoflagger_path')
use_singularity = configuration_file.getboolean('globals','use_singularity')
telescope = configuration_file.get('globals','telescope')
# singularity_bind = configuration_file.get('globals','singularity_bind')
wsclean_sif = configuration_file.get('globals','wsclean_sif')
singularity_container = configuration_file.get('globals','singularity_container')

# apriori
use_casa = configuration_file.getboolean('apriori','use_casa')
attach_metadata = configuration_file.getboolean('apriori','attach_metadata')
fitsfiles_dir = configuration_file.get('apriori','fitsfiles_dir')
antab_file = configuration_file.get('apriori','antab_file')
uvflg_file = configuration_file.get('apriori','uvflg_file')
do_apriori_cal = configuration_file.getboolean('apriori','do_apriori_cal')
apply_apriori_cal = configuration_file.getboolean('apriori','apply_apriori_cal')

target = configuration_file.get('basic','target')
phase_calibrator = configuration_file.get('basic','phase_calibrator')
fringe_finder = configuration_file.get('basic','fringe_finder')
integration_time = configuration_file.getfloat('basic','integration_time')
do_split = configuration_file.getboolean('basic','do_split')
timebin = configuration_file.get('basic','timebin')
width = configuration_file.getint('basic','width')
verbosity = configuration_file.getboolean('basic','verbosity')

# flag
do_flagging = configuration_file.getboolean('flagging','do_flagging')
# aoflagger_strategies = configuration_file.get('flagging','aoflagger_strategies')
# faint_source_strategy = configuration_file.get('flagging','faint_source_strategy')
# bright_source_strategy = configuration_file.get('flagging','bright_source_strategy')
manual_file = configuration_file.get('flagging','manual_file')
edge_channel_fraction = configuration_file.getfloat('flagging','edge_channel_fraction')
use_aoflagger = configuration_file.getboolean('flagging','use_aoflagger')
flagging_strategy = configuration_file.get('flagging','flagging_strategy')

flag_antenna = configuration_file.getboolean('flagging','flag_antenna')
antenna_to_flag = configuration_file.get('flagging','antenna_to_flag')
export_uvfits = configuration_file.getboolean('flagging','export_uvfits')


# calibrate
do_tec_corrections = configuration_file.getboolean('calibrate','do_tec_corrections')
do_sbd_fringe = configuration_file.getboolean('calibrate','do_sbd_fringe')
apply_sbd = configuration_file.getboolean('calibrate','apply_sbd')
refant = configuration_file.get('calibrate','refant')
timerange = configuration_file.get('calibrate','timerange')
snr_sbd = configuration_file.getfloat('calibrate','snr_sbd')
do_mbd_fringe = configuration_file.getboolean('calibrate','do_mbd_fringe')
apply_mbd = configuration_file.getboolean('calibrate','apply_mbd')
snr_mbd = configuration_file.getfloat('calibrate','snr_mbd')
solint = configuration_file.getfloat('calibrate','solint')
do_bpass = configuration_file.getboolean('calibrate','do_bpass')
apply_bpass = configuration_file.getboolean('calibrate','apply_bpass')
make_dirty_map = configuration_file.getboolean('calibrate','make_dirty_map')


# split_calibrated = configuration_file.getboolean('calibrate','split_calibrated')
# imsize= [int(part) for part in configuration_file.get('calibrate', 'imsize').split(',')]
# detection_threshold = configuration_file.getfloat('selfcal','detection_threshold')

# # selfcal
# do_selfcal = configuration_file.getboolean('selfcal','do_selfcal')
# use_tclean = configuration_file.getboolean('selfcal','use_tclean')
# use_wsclean = configuration_file.getboolean('selfcal','use_wsclean')

# pybdsf_threshold = configuration_file.get('selfcal','pybdsf_threshold')
# pybdsf_niter = configuration_file.getint('selfcal','pybdsf_niter')

# weighting = configuration_file.get('selfcal','weighting')
# robust = configuration_file.getfloat('selfcal','robust')

# nloops = configuration_file.getint('selfcal','nloops')
# calmode = configuration_file.get('selfcal','calmode').split(',')
# gaintype = configuration_file.get('selfcal','gaintype').split(',')
# cell =  configuration_file.get('selfcal', 'cell')
# threshold = configuration_file.get('selfcal','threshold').split(',')
# minsnr = [float(part) for part in configuration_file.get('selfcal', 'minsnr').split(',')]
# # imsize= [int(part) for part in configuration_file.get('selfcal', 'imsize').split(',')]
# # niter = [int(part) for part in configuration_file.get('selfcal', 'niter').split(',')]
# niter = int(configuration_file.get('selfcal','niter'))
# niter_final = int(configuration_file.getint('selfcal','niter_final'))
# threshold_final = configuration_file.getint('selfcal','threshold_final')
# robust = configuration_file.getfloat('selfcal','robust')

# tclean_threshold = configuration_file.get('selfcal','tclean_threshold')
# solint_selfcal = configuration_file.get('selfcal','solint_selfcal').split(',')
# apply_to_target = configuration_file.getboolean('selfcal','apply_to_target')
# detect_sources = configuration_file.getboolean('selfcal','detect_sources')

# # pb corrections
# do_pbcor = configuration_file.getboolean('pbcor','do_pbcor')
# pb_file = configuration_file.get('pbcor','pb_file')


vis = experiment+'.ms'
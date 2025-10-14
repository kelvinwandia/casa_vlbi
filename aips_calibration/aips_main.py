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


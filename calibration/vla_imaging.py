import os, glob, re, subprocess
import numpy as np
from astropy.constants import c
import casatools
from casatasks import *


working_directory = '/raid1/scratch/kelvinw/gcs/m15/measurement_sets/working_dir'
path_to_vis = '/raid1/scratch/kelvinw/gcs/m15/measurement_sets'
vis = os.path.join(path_to_vis,'*.ms')
target = 'M 15 X-2' 

def set_working_dir():

    if not os.path.exists(working_directory):
        # logging.info(f"{working_directory} does not exist, making one")
        os.makedirs(working_directory)
    else:
        # logging.info(f"Working directory {working_directory} already exists")
        pass

    os.chdir(working_directory)


def split():

    for msname in vis:
        outputvis = msname.replace('.ms',f'_{target}_split.ms')
        if not os.path.exists(outputvis):
            split(vis=outputvis,field=target,datacolumn='corrected')
        
        listobs(vis=outputvis,listfile=outputvis.replace('.ms','.txt'),overwrite=True)


def getimaging_params():
    outputvis = msname.replace('.ms',f'_{target}_split.ms')
    for msname in vis:
        ms = casatools.ms()
        tb = casatools.table()
        ms.open(outputvis)
        max_uv = ms.getdata('uvdist')['uvdist'].max()
        ms.close()

        tb.open(vis + '/SPECTRAL_WINDOW')
        chan_freq = tb.getcol('CHAN_FREQ')
        highest_freq = chan_freq.max()
        tb.close()

        # 3.6e6 converts the degrees to mas
        # 5 is the sampling
        cell_size = ((c.value / highest_freq) / max_uv) * (180. / np.pi) * (3.6e6 / 5)
        cell_size = np.round(cell_size)
        logging.info("You are using a cell size of:", cell_size)

 

    return cell_size




import os, glob, re, subprocess
import numpy as np
from astropy.constants import c
import casatools
from casaplotms import *
from casatasks import *


working_directory = '/raid1/scratch/kelvinw/gcs/m15/measurement_sets/working_dir'
path_to_vis = '/raid1/scratch/kelvinw/gcs/m15/measurement_sets'
vis = glob.glob(os.path.join(path_to_vis, '*.ms'))

# vis = '/raid1/scratch/kelvinw/gcs/m15/working_directory_observation.55707.477620949074/11A-269_sb4158083_1.55707.47758597222_hs.ms'
target = 'M 15 X-2' 

def set_working_dir():
    if not os.path.exists(working_directory):
        print(f"{working_directory} does not exist, creating it.")
        os.makedirs(working_directory)
    else:
        print(f"Working directory {working_directory} already exists.")

    os.chdir(working_directory)
    print(f"Current working directory: {os.getcwd()}")


def split_target():
    for msname in vis:
        outputvis = msname.replace(' ','_').replace('.ms',f'_{target}_split.ms')
        print(outputvis)
        if not os.path.exists(outputvis):
            split(vis=msname,outputvis=outputvis,field=target,datacolumn='corrected')
        
        listobs(vis=outputvis,listfile=outputvis.replace('.ms','.txt'),overwrite=True)


def plotuv_cov():
    for msname in vis:
        plotfile = msname.replace('.ms','_uvcov.png')
        if not os.path.exists(plotfile):
            plotms(vis=msname, xaxis='Uwave',yaxis='Vwave',plotfile=plotfile,field=target,
                    height=750,width=1500,showgui=False,overwrite=True)



def getimaging_params():

    
    imaging_params = {}
    for msname in vis:
        outputvis = msname.replace(' ','_').replace('.ms',f'_{target}_split.ms')
        ms = casatools.ms()
        tb = casatools.table()
        ms.open(outputvis)
        max_uv = ms.getdata('uvdist')['uvdist'].max()
        ms.close()

        tb.open(outputvis + '/SPECTRAL_WINDOW')
        chan_freq = tb.getcol('CHAN_FREQ')
        highest_freq = chan_freq.max()
        tb.close()

        # 3.6e6 converts the degrees to mas
        # 5 is the sampling
        cell_size = ((c.value / highest_freq) / max_uv) * (180. / np.pi) * (3.6e6 / 5)
        cell_size = np.round(cell_size)
        imaging_params[msname] = cell_size
        print(imaging_params)


set_working_dir()
split_target()
plotuv_cov()
getimaging_params()



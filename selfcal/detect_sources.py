
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from utils.helper_functions import *





# phasecenter = ['21h30m01.203493s +12d10m38.1592s']
# phasecenter = ['J2000 21h39m01.309269s +14d23m35.99221s']

### correlation centres 
# phasecenter=['21h29m58.3500s +12d10m01.500s','21h29m58.3120s +12d10m02.679s','21h30m01.2034s +12d10m38.160s',
#         '21h29m51.9025s +12d10m17.132s','21h29m56.3050s +12d11m01.500s','21h29m56.3050s +12d09m11.500s',
#         '21h30m02.4410s +12d09m11.500s']

### detected sources in order (from Kirsten 2015) -- M15A, M15C, AC211,S1,S2

phasecenter = ['21h29m58.246512s +12d10m01.2339s','21h30m01.203493s +12d10m38.1592s','21h29m58.312403s +12d10m02.6740s',
               '21h29m51.9034555d +12d10m17.13240s','21h30m02.085700s +12d09m04.2203s']



@time_execution
def m15_sources():

    """
    Image sources detected in Kirsten et.al 2015
    """
    target_ms = target+'.ms'

    msmd = casatools.msmetadata()
    msmd.open(target_ms)
    scans = msmd.scansforfield(field=target)
    nscans = len(scans)


    for coord in range(len(phasecenter)):
        phaseshifted_ms = phasecenter[coord].replace(" ","").replace("+","_").replace("J2000","")+'_phaseshifted.ms'

        if not os.path.exists(phaseshifted_ms):
            phaseshift(vis=target_ms,outputvis=phaseshifted_ms,datacolumn='corrected',phasecenter="J2000"+" "+phasecenter[coord])
        else:
            logging.info(f"{phaseshifted_ms} exists. A new one will not be created")

        transformed_ms = phasecenter[coord].replace(" ","").replace("+","_").replace("J2000","")+'_transformed.ms'
        if not os.path.exists(transformed_ms):
            mstransform(vis=phaseshifted_ms,outputvis=transformed_ms,datacolumn='data',createmms=False,
                        separationaxis='scan',numsubms=msmd.nscans(),timeaverage=True,chanaverage=True,
                        timebin='20s',chanbin=16)
        else:
            logging.info(f"{transformed_ms} exists. A new one will not be created")
        
        os.system(f"rm -r {phaseshifted_ms}")

        imagename =  transformed_ms.replace(".ms","")
        if not os.path.exists(imagename):
            logging.info(f"Making image {imagename}")

            wsclean_cmd = ['wsclean', '-log-time', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                    '-scale', f'{cell}','-mgain', '0.8', '-niter', '0', f'{transformed_ms}']

            run_wsclean(wsclean_sif,wsclean_cmd)

            wsclean_fitsfile = imagename+'-image.fits'
            get_im_stats(wsclean_fitsfile)
            plot_fits(wsclean_fitsfile)





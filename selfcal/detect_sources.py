from casatasks import *
import casatools
import subprocess, os

from utils.helper_functions import *


# phasecenter=['21h29m58.3500s +12d10m01.500s','21h29m58.3120s +12d10m02.679s','21h30m01.2034s +12d10m38.160s',
#         '21h29m51.9025s +12d10m17.132s','21h29m56.3050s +12d11m01.500s','21h29m56.3050s +12d09m11.500s',
#         '21h30m02.4410s +12d09m11.500s']


phasecenter = ['J2000 21h39m01.309269s +14d23m35.99221s']



@time_execution
def m15_sources():

    """
    Image sources detected in Kirsten et.al 2015
    """
    target_ms = phase_calibrator+'.ms'

    msmd = casatools.msmetadata()
    msmd.open(target_ms)
    scans = msmd.scansforfield(field='J2139+1423')
    nscans = len(scans)


    for coord in range(len(phasecenter)):
        phaseshifted_ms = phasecenter[coord].replace(" ","").replace("+","_").replace("J2000","")+'_phaseshifted.ms'

        if not os.path.exists(phaseshifted_ms):
            phaseshift(vis=target_ms,outputvis=phaseshifted_ms,datacolumn='corrected',phasecenter=phasecenter[coord])
        else:
            logging.info(f"{phaseshifted_ms} exists. A new one will not be created")

        transformed_ms = phasecenter[coord].replace(" ","").replace("+","_").replace("J2000","")+'_transformed.ms'
        if not os.path.exists(transformed_ms):
            mstransform(vis=phaseshifted_ms,outputvis=transformed_ms,datacolumn='data',createmms=False,
                        separationaxis='scan',numsubms=msmd.nscans(),timeaverage=True,chanaverage=True,
                        timebin='20s',chanbin=16)
        else:
            logging.info(f"{transformed_ms} exists. A new one will not be created")
            
        imagename =  transformed_ms.replace('.ms','')
        print(f"Making image {imagename}")

        wsclean_cmd = ['wsclean', '-log-time', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                '-scale', f'{cell}','-mgain', '0.8', '-niter', f'{niter}', f'{transformed_ms}']

        run_wsclean(wsclean_cmd)

        wsclean_fitsfile = imagename+'-image.fits'
        get_im_stats(wsclean_fitsfile)
        plot_fits(wsclean_fitsfile)






import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from utils.helper_functions import *

def make_mms():

    # Use mms to speed up casa imaging -- this needs a flag
    # such that it can be easily disabled
    ## split for the field -- e

    # msmd = casatasks.msmetadata()
    # msmd.open(vis)
    # scans = msmd.scansforfield(field=field)
    # nscans = len(scans)

    # sources = [phase_calibrator,target]
    # for source in sources:
    #     mstransform(
    #         vis=vis, outputvis=source+'.ms', datacolumn='corrected',field=source, 
    #         createmms=True, separationaxis = 'scan', numsubms = msmd.nscans()
    #     )
    # vis = '/home/kelvin/Desktop/gv020_working_dir/gv020b/J2139+1423.ms'
    # phase_calibrator = 'J2139+1423'
    # imagename = f'imagename_{phase_calibrator}'
    # # logging.info(f"Making {imagename}")
    
    # os.system(f"rm -r {imagename}.*")
    
    # imsize = [640,640]
    # cell='1mas'
    
    # logging.info("Running tclean")
    # casatasks.tclean(vis=vis,imsize=imsize,imagename=imagename,cell=cell,
    #     niter=0, deconvolver='clark',interactive=False, gridder='standard',field=phase_calibrator,
    #     parallel=True
    #     )  

    pass

@time_execution
def split_selfcal():
    
    sources = [phase_calibrator,target]

    for source in sources:
        outputvis = source+'.ms'
        if not os.path.exists(outputvis):
            print(f"Splitting {vis} to {outputvis}")
            #TODO : CHECK DATA COLUMN CAREFULLY - USING DATA IF FULLY CALIBRATED IN AIPS 
            split(
                vis = vis, outputvis = outputvis, datacolumn='data') 
        else:
            print(f"{outputvis} exists. Will not make a new one")

    global phasecal_ms, target_ms
    phasecal_ms = phase_calibrator+'.ms'
    target_ms = target+'.ms'


def run_pybdsf(input_image):

    # The input image is a casa .image that then gets exported to a FITS
    imagename = input_image

    img = bdsf.process_image(imagename,adaptive_rms_box=True, thresh='hard',
                            thresh_isl=True, thresh_pix = detection_threshold, advanced_opts=True,
                            mean_map='map', rms_map =True, group_by_isl=True)
    # adaptive_rms_box=False, spline_rank=4, thresh='hard', thresh_isl=True, thresh_pix = detection_threshold
    # Write out island mask and FITS catalog -- for the large map
    img.export_image(outfile=imagename.replace('.fits','')+'.maskfile.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=imagename+'.cat', format='fits', clobber=True, catalog_type ='gaul')
    
    regionfile = imagename+'.casabox'
    ascii_file = imagename+'.ascii'
    rmsfile = imagename+'.rmsfile'

    img.write_catalog(outfile=regionfile,format='casabox',clobber=True,catalog_type='srl')
    img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
    img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

    return regionfile

@time_execution
def selfcal_part1():

    """
    Creates an (a large) an image that is used to create a casa region file using pybdsf 
    for masking
    """
    msmd.open(vis)
    source = phase_calibrator
    field_id = msmd.fieldsforname(source)[0]
    msmd.close()
    # global first_part_imagename
    pybdsf_imagename = source.replace('.ms','')+'_pybdsf'
    if not os.path.exists(pybdsf_imagename+'-image.fits'):
        print(f"Making {pybdsf_imagename}")
        # wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{pybdsf_threshold}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{pybdsf_imagename}','-scale', f'{cell}',\
        #                     '-mgain', '0.8', '-niter', f'{pybdsf_niter}', f'{source}']
        wsclean_cmd = ['wsclean', '-log-time','-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{pybdsf_imagename}','-scale', f'{cell}',\
                            '-mgain', '0.8', '-niter', f'{pybdsf_niter}' , '-field',f'{field_id}',f'{vis}']
        
        run_wsclean(wsclean_cmd)

    regionfile = run_pybdsf(input_image=pybdsf_imagename+'-image.fits')



@time_execution
def selfcal_part2():

    msmd.open(vis)
    source = phase_calibrator
    field_id = msmd.fieldsforname(source)[0]
    msmd.close()

    pybdsf_imagename = source.replace('.ms','')+'_pybdsf'
    maskfile = pybdsf_imagename+ '-image.maskfile.fits'

    print("Deleting model column before selfcal")
    delmod(vis=vis,otf=True,field=str(field_id))

    for selfcal_loop in range(nloops):
        caltable = f'caltable_{selfcal_loop}.tb'
        prev_caltables = sorted(glob.glob('*.tb'))
        if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
            applycal(vis=vis, gaintable = prev_caltables, field=str(field_id), parang=False )
    
        imagename = source.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
        if os.path.exists(imagename):
            print("Continuing to the next image")
        
        else:
            imagename =  source.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
            print(f"Making image {imagename}")

            wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{threshold[selfcal_loop]}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                        '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
                        '-mgain', '0.8', '-niter', f'{niter}', '-field',f'{field_id}', f'{vis}']

            run_wsclean(wsclean_cmd)

            wsclean_fitsfile = imagename+'-image.fits'
            get_im_stats(wsclean_fitsfile)
            plot_fits(wsclean_fitsfile)
          
            model_fits = imagename.replace('-image.fits','-model.fits')

            print(f"Adding modelcolumn to data. Using {model_fits} to predict")
            predict_cmd = ['wsclean', '-log-time', '-predict', '-reorder' ,'-field',f'{field_id}','-name', f'{imagename}', vis]

            # Predicting
            run_wsclean(predict_cmd)

            # Plot the model column
            plotms(
                vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
                field = str(field_id)
            )
            if calmode[selfcal_loop] == 'p':
                minblperant = 3
            else:
                minblperant = 4

            gaincal(vis = vis, caltable = caltable, refant = refant, solint = solint_selfcal[selfcal_loop],
                    gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
                    calmode = calmode[selfcal_loop], append=False, parang=False, minblperant=minblperant, field=str(field_id)
                    )
            coloraxis = ['corr','spw']
            for color in coloraxis:
                if calmode[selfcal_loop] =='p':
                    plotms(
                        vis = caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
                        iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
                        plotfile=caltable.replace('.tb',f'_{color}.png'), dpi=300, width=1500, height=750,
                    )
                else:
                    plotms(
                            vis = caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
                            iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
                            plotfile=caltable.replace('.tb',f'_{color}.png'), dpi=300, width=1500, height=750
                        )

            if selfcal_loop == nloops-1:
                prev_caltables = sorted(glob.glob('*.tb'))
                print("Applying the caltable derived from last gaincal iteration")
                applycal(vis=vis, gaintable = prev_caltables,field=str(field_id), parang=False )
        
        # ### Get the last imagename from the loop and generate a final mask
        
    imagename_final = source.replace('.ms','')+f'final_map_loop_{nloops-1}'
    ##  tclean here to make the final image
    print("Make final image with all selfcal corrections applied")
    wsclean_cmd_final = ['wsclean', '-log-time', '-auto-threshold',f'{threshold_final}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename_final}', \
            '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
            '-mgain', '0.8', '-niter', f'{niter_final}','-field',f'{field_id}', f'{vis}']

    run_wsclean(wsclean_cmd_final)

    wsclean_fitsfile = imagename_final+'-image.fits'
    get_im_stats(wsclean_fitsfile)
    plot_fits(wsclean_fitsfile)


@time_execution
def applycal_target():

    """
    Applies cal to the target field
    """
    prev_caltables = sorted(glob.glob('*.tb'))
    print(f"Applying calibration tables {prev_caltables} to {target_ms}")
    applycal(
        vis = target_ms, gaintable = prev_caltables, parang=False
    )




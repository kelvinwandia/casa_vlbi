
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from utils.helper_functions import *



@time_execution
def split_selfcal():

    sources = [phase_calibrator,target]
    vis_to_split = os.path.join(working_directory,vis)
    
    width = 4; timebin='2s'

    for source in sources:
        outputvis = source+'.ms'
        if not os.path.exists(outputvis):
            logging.info(f"======>>>Splitting {vis} to {outputvis}")
            logging.info(f"Averaging to {width} channels and {timebin} ")
            #TODO : CHECK DATA COLUMN CAREFULLY - USING DATA IF FULLY CALIBRATED IN AIPS 
            # split(vis = vis_to_split, outputvis = outputvis, datacolumn='data',field=source) 
            split(vis = vis_to_split, outputvis = outputvis, datacolumn='corrected',field=source,
            width=width,timebin=timebin) 
        else:
            logging.info(f"======>>>{outputvis} exists. Will not make a new one")

       

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
    # msmd.open(vis)
    source = phase_calibrator
    # field_id = msmd.fieldsforname(source)[0]
    # msmd.close()
    # global first_part_imagename
    pybdsf_imagename = source.replace('.ms','')+'_pybdsf'
    if not os.path.exists(pybdsf_imagename+'-image.fits'):
        logging.info(f"======>>>Making {pybdsf_imagename}")

        
        if use_tclean == True:

            tclean(vis= phasecal_ms, imagename=pybdsf_imagename,imsize=imsize, cell=cell,
                gridder='standard',weighting='briggs',robust=robust,niter=0)
            fitsname = pybdsf_imagename+'.fits'
            exportfits(imagename=pybdsf_imagename+'.image',fitsimage=fitsname,overwrite=True)
            get_im_stats(fitsname)
            plot_fits(fitsname)

            regionfile = run_pybdsf(input_image=fitsname)

        elif use_wsclean == True:

            # wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{pybdsf_threshold}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{pybdsf_imagename}','-scale', f'{cell}',\
            #                     '-mgain', '0.8', '-niter', f'{pybdsf_niter}', f'{source}']
            wsclean_cmd = ['wsclean', '-log-time','-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{pybdsf_imagename}','-scale', f'{cell}',\
                                '-mgain', '0.8', '-niter', f'{pybdsf_niter}' ,f'{phasecal_ms}']
            
            run_wsclean(wsclean_sif,wsclean_cmd)

            regionfile = run_pybdsf(input_image=pybdsf_imagename+'-image.fits')
        
        elif use_tclean and use_wsclean == True:
            logging.critical("You cannot use both imagers at once, check the config file")

        else:
            logging.info(f"Imager not selected")



@time_execution
def selfcal_part2():


    # msmd.open(vis)
    source = phase_calibrator
    # field_id = msmd.fieldsforname(source)[0]
    # msmd.close()

    pybdsf_imagename = source.replace('.ms','')+'_pybdsf'
    regionfile = pybdsf_imagename+'.casabox'
    maskfile = pybdsf_imagename+ '-image.maskfile.fits'


    logging.info("======>>>Deleting model column before selfcal")
    # delmod(vis=phasecal_ms,otf=True,field=str(field_id))
    delmod(vis=phasecal_ms,otf=True)

    for selfcal_loop in range(nloops):
        caltable = f'caltable_{selfcal_loop}.tb'
        prev_caltables = sorted(glob.glob('*.tb'))
        if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
            logging.info(f"======>>>Applying {prev_caltables}")
            # applycal(vis=phasecal_ms, gaintable = prev_caltables, field=str(field_id), parang=False )
            applycal(vis=phasecal_ms, gaintable = prev_caltables, parang=False )
    
        imagename = source.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
        if os.path.exists(imagename):
            logging.info("Continuing to the next image")
        
        else:
            imagename =  source.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
            logging.info(f"======>>>Making image {imagename}")

            if use_tclean == True:
                tclean(
                    vis= phasecal_ms, imagename=imagename,imsize=imsize, cell=cell,
                    gridder='standard',weighting='briggs',robust=robust,niter=niter,
                    threshold = threshold[selfcal_loop],usemask='user',mask=regionfile
                    )
                
                fitsname = imagename+'.fits'
                exportfits(imagename=imagename+'.image',fitsname=fitsname,overwrite=True)
                get_im_stats(fitsname)
                plot_fits(fitsname)

                logging.info("======>>>Adding modelcolumn to data")
                ft(vis = vis, model=imagename+'.image',usescratch=True)

            elif use_wsclean == True:

                # wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{threshold[selfcal_loop]}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                #             '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
                #             '-mgain', '0.8', '-niter', f'{niter}', '-field',f'{field_id}', f'{phasecal_ms}']
                
                wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{threshold[selfcal_loop]}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                            '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
                            '-mgain', '0.8', '-niter', f'{niter}', f'{phasecal_ms}']

                run_wsclean(wsclean_sif,wsclean_cmd)

                wsclean_fitsfile = imagename+'-image.fits'
                get_im_stats(wsclean_fitsfile)
                plot_fits(wsclean_fitsfile)
            
                model_fits = imagename.replace('-image.fits','-model.fits')

                logging.info(f"======>>>Adding modelcolumn to data. Using {model_fits} to predict")
                # predict_cmd = ['wsclean', '-log-time', '-predict', '-reorder' ,'-field',f'{field_id}','-name', f'{imagename}', phasecal_ms]
                predict_cmd = ['wsclean', '-log-time', '-predict', '-reorder' ,'-name', f'{imagename}', phasecal_ms]
                # Predicting
                run_wsclean(wsclean_sif,predict_cmd)

            elif use_tclean and use_wsclean == True:
                logging.critical("You cannot use both imagers at once, check the config file")

            else:
                logging.info(f"Imager not selected")




            # Plot the model column
            plotfile = f"{imagename}_modelcolumn.png"
            plotms(
                vis=phasecal_ms, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=plotfile, overwrite=True, width=1500, height=750,
                # field = str(field_id)
            )
            if calmode[selfcal_loop] == 'p':
                minblperant = 3
            else:
                minblperant = 4

            gaincal(vis = phasecal_ms, caltable = caltable, refant = refant, solint = solint_selfcal[selfcal_loop],gainfield=source,
                    gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
                    calmode = calmode[selfcal_loop], append=False, parang=False, minblperant=minblperant,
                    #   field=str(field_id)
                    )
            coloraxis = ['corr','spw']
            for color in coloraxis:
                plotfile = f"{caltable.replace('.tb', f',f_{color}.png')}"
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
                logging.info("======>>>Applying the caltable derived from last gaincal iteration")
                # applycal(vis=phasecal_ms, gaintable = prev_caltables,field=str(field_id), parang=False )
                applycal(vis=phasecal_ms, gaintable = prev_caltables, parang=False )

        
        # ### Get the last imagename from the loop and generate a final mask
        
    imagename_final = source.replace('.ms','')+f'_final_map_loop_{nloops-1}'
    ##  tclean here to make the final image
    logging.info("======>>>Make final image with all selfcal corrections applied")

    if use_tclean == True:
        tclean(
                vis= phasecal_ms, imagename=imagename_final,imsize=imsize, cell=cell,
                gridder='standard',weighting='briggs',robust=robust,niter=niter,
                threshold = '1mJy',usemask='user',mask=regionfile
                )
        fitsname = imagename_final+'.fits'
        exportfits(imagename=imagename_final+'.image',fitsname=fitsname,overwrite=True)
        get_im_stats(fitsname)
        plot_fits(fitsname)

        logging.info("======>>>Adding modelcolumn to data")
        ft(vis = vis, model=imagename_final+'.image',usescratch=True)


    if use_wsclean == True:
        wsclean_cmd_final = ['wsclean', '-log-time', '-auto-threshold',f'{threshold_final}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename_final}', \
                '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
                '-mgain', '0.8', '-niter', f'{niter_final}', f'{phasecal_ms}']

        run_wsclean(wsclean_sif,wsclean_cmd_final)

        wsclean_fitsfile = imagename_final+'-image.fits'
        get_im_stats(wsclean_fitsfile)
        plot_fits(wsclean_fitsfile)


@time_execution
def applycal_target():

    """
    Applies cal to the target field
    """
    target_ms = target + '.ms'
    prev_caltables = sorted(glob.glob('*.tb'))

    logging.info(f"======>>>Applying calibration tables {prev_caltables} to {target_ms}")
    applycal(
        vis = target_ms, gaintable = prev_caltables, parang=False)

    imagename = target_ms.replace('.ms','_target_dirty')
    logging.info(f"Imaging target")

    if use_wsclean == True:
        wsclean_target_dirty= ['wsclean', '-log-time', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                '-scale', f'{cell}','-mgain', '0.8', '-niter', '0', f'{target_ms}']
        
        run_wsclean(wsclean_sif,wsclean_target_dirty)

        wsclean_fitsfile = imagename+'-image.fits'
        get_im_stats(wsclean_fitsfile)
        plot_fits(wsclean_fitsfile)






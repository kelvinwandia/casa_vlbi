
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np

def plot_fits(fitsname):
    """
    Plots fitsfiles using astropy
    """
    fitsfile = fits.open(fitsname)
    image_data = fitsfile[0].data[0,0,:,:]
    ny, nx = image_data.shape
    x_center = nx // 2
    y_center = ny // 2
    x_new = np.arange(nx) - x_center
    y_new = np.arange(ny) - y_center

    fig, ax = plt.subplots()

    # image_plot = ax.imshow(image_data, origin='lower', 
    #                    extent=[x_new.min(), x_new.max(), y_new.min(), y_new.max()],cmap='viridis')
    image_plot = ax.imshow(image_data, origin='lower', 
                       extent=[-32, 32, -32, 32],cmap='viridis')
    cbar = plt.colorbar(image_plot,ax=ax,orientation='vertical')
    # ax.set_title(sources_to_image,fontsize=16)
    plt.savefig(fitsname.replace('.fits','.pdf'))

def get_im_stats(imagename):
    
    """
    Gets the statistics for either a 256x256 pix image and writes
    them to a logfile
    """


    rms=casatasks.imstat(imagename=imagename,box='60,60,580,240')['rms'][0]  # for 640x640 px
    peak=casatasks.imstat(imagename=imagename,box='300,300,340,340')['max'][0]
    print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))
    
    logfile = 'imstat.txt'
    casa_imstat = imstat(imagename)
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))

        txt_file.write(f"For {imagename}, the maximum pos for imstat is {casa_imstat['maxposf']}\n")

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

def pybdsf(input_image):

    # The input image is a casa .image that then gets exported to a FITS
    imagename = input_image.replace('.image','')
    fitsname = imagename+'.fits'
    exportfits(imagename = input_image, fitsimage=fitsname, overwrite=True)

    img = bdsf.process_image(fitsname,adaptive_rms_box=True, thresh='hard',
                            thresh_isl=True, thresh_pix = detection_threshold, advanced_opts=True,
                            mean_map='map', rms_map =True, group_by_isl=True)
    # adaptive_rms_box=False, spline_rank=4, thresh='hard', thresh_isl=True, thresh_pix = detection_threshold
    # Write out island mask and FITS catalog -- for the large map
    img.export_image(outfile=imagename+'.maskfile.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=imagename+'.cat', format='fits', clobber=True, catalog_type ='gaul')
    
    regionfile = imagename+'.casabox'
    ascii_file = imagename+'.ascii'
    rmsfile = imagename+'.rmsfile'

    img.write_catalog(outfile=regionfile,format='casabox',clobber=True,catalog_type='srl')
    img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
    img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

    return regionfile



def selfcal_part1(field):

    """
    Creates an (a large) an image that is used to create a casa region file using pybdsf 
    for masking
    """
    
    # global first_part_imagename
    first_part_imagename = field+'_pybdsf_masking'

    if not os.path.exists(first_part_imagename):
        print(f"Making {first_part_imagename}")
        tclean(
            vis = vis, imagename=first_part_imagename, imsize=imsize, cell=cell, field=field,
            gridder = 'standard', weighting = weighting, robust = robust, niter=pybdsf_niter, threshold = pybdsf_threshold,   
        )

    regionfile = pybdsf(input_image=first_part_imagename+'.image')


def selfcal_part2(field):


    regionfile = field+'_pybdsf_masking.casabox' 

    print("Deleting model column before selfcal")
    delmod(vis=vis,otf=True)

    for selfcal_loop in range(nloops):
        caltable = f'caltable_{selfcal_loop}.tb'
        prev_caltables = sorted(glob.glob('*.tb'))
        if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
            applycal(vis=vis, gaintable = prev_caltables, parang=False )
    
        imagename = f'{field}_{selfcal_loop}'
        if os.path.exists(imagename):
            print("Continuing to the next image")
        
        else:
            imagename = f'target_selfcal_{selfcal_loop}'
            print(f"Making image {imagename}")
            tclean(
                vis = vis, imagename=imagename, imsize=imsize, cell=cell,
                parallel=False,
                gridder = 'standard', weighting = weighting, robust = robust, niter=niter[selfcal_loop], threshold = threshold[selfcal_loop],
                interactive=False, usemask='user', mask=regionfile, field=field
            )
            exportfits(imagename=imagename+'.image',fitsimage=imagename+'.fits',overwrite=True)
            get_im_stats(imagename+'.image')
            plot_fits(imagename+'.fits')
          

            print("Adding modelcolumn to data")
            # model images from the MTMFS images,
            ft(vis = vis, model= imagename+'.model',usescratch=True)

            # plot the model column
            plotms(
                vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
            )

            gaincal(vis =vis, caltable = caltable, refant = refant, solint = solint_selfcal[selfcal_loop],
                    gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
                    calmode = calmode[selfcal_loop], append=False, parang=False
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
                applycal(vis=vis, gaintable = prev_caltables, parang=False )
        
        # ### Get the last imagename from the loop and generate a final mask
        
    imagename = field +f'_{nloops-1}'+'.final'
    ##  tclean here to make the final image
    print("Make final image with all selfcal corrections applied")
    tclean(
        vis = vis, imagename = imagename, imsize=imsize, cell=cell, weighting = weighting,
        robust = robust, niter=niter_final, threshold = threshold_final, field=field,
        interactive=False, usemask = 'user', mask=regionfile,
    )


### Use the output here to peel -- wsclean predict should work
## implement using wsclean -- also no need to create a large image


# def applycal_target():

#     """
#     Applies cal to the target field
#     """
#     prev_caltables = sorted(glob.glob('*.tb'))
#     applycal(
#         vis = vis_tocal, gaintable = prev_caltables, parang=False
#     )





import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


# vis = '/home/kelvin/Desktop/gv020_working_dir/gv020b/J2139+1423.ms'

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


def split_selfcal():
    
    sources = [phase_calibrator,target]

    for source in sources:
        outputvis = source+'.ms'
        if not os.path.exists(outputvis):
            print(f"Splitting {vis} to {outputvis}")
            split(
                vis = vis, outputvis = outputvis, datacolumn='corrected')
        else:
            print(f"{outputvis} exists. Will not make a new one")

    global phasecal_ms, target_ms
    phasecal_ms = phase_calibrator+'.ms'
    target_ms = target+'.ms'

def pybdsf(input_image):

    imagename = input_image
    fitsname = imagename

    img = bdsf.process_image(fitsname,adaptive_rms_box=True, thresh='hard',
                            thresh_isl=True, thresh_pix = detection_threshold, advanced_opts=True,
                            mean_map='map', rms_map =True, group_by_isl=True)
    # adaptive_rms_box=False, spline_rank=4, thresh='hard', thresh_isl=True, thresh_pix = detection_threshold
    # Write out island mask and FITS catalog -- for the large map
    fits_maskfile = imagename.replace('.fits','.maskfile.fits')
    catalog_file = imagename.replace('.fits','.cat')
    img.export_image(outfile=fits_maskfile,img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catalog_file, format='fits', clobber=True, catalog_type ='gaul')
    
    regionfile = imagename.replace('.fits','.casabox')
    ascii_file = imagename.replace('.fits','.ascii')
    rmsfile = imagename.replace('.fits','.rmsfile')

    img.write_catalog(outfile=regionfile,format='casabox',clobber=True,catalog_type='srl')
    img.write_catalog(outfile=ascii_file, format='ascii', clobber=True, catalog_type='gaul')
    img.export_image(outfile=rmsfile, img_type='rms', img_format='fits', clobber=True)

    return regionfile

def run_wsclean(command):

    """
    Runs wsclean commands 
    """

    container = wsclean_sif
    if os.path.exists(container):
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname(wsclean_sif)))

    command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + command
    try:
        print("Executing: %s", ' '.join(command_to_execute))
        process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        print("stdout: %s", stdout)
        print("stderr: %s", stderr)

        return_code = process.returncode
        if return_code == 0:
            print(f"Strategy executed successfully. Output:\n{stdout}")
        else:
            print(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")  

    except Exception as e:
        print(f"An error occurred: {e}")

def selfcal_part1():

    """
    Creates an (a large) an image that is used to create a casa region file using pybdsf 
    for masking
    """
    
    # global first_part_imagename
    pybdsf_imagename = phasecal_ms.replace('.ms','')+'_pybdsf'
    if not os.path.exists(pybdsf_imagename+'-image.fits'):
        print(f"Making {pybdsf_imagename}")
        wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{pybdsf_threshold}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{pybdsf_imagename}','-scale', f'{cell}',\
                            '-mgain', '0.8', '-niter', f'{pybdsf_niter}', f'{phasecal_ms}']
        
        run_wsclean(wsclean_cmd)

    regionfile = pybdsf(input_image=pybdsf_imagename+'-image.fits')



def selfcal_part2():

    pybdsf_imagename = phasecal_ms.replace('.ms','')+'_pybdsf'
    maskfile = pybdsf_imagename+ '-image.maskfile.fits'

    print("Deleting model column before selfcal")
    delmod(vis=phasecal_ms,otf=True)

    for selfcal_loop in range(nloops):
        caltable = f'caltable_{selfcal_loop}.tb'
        prev_caltables = sorted(glob.glob('*.tb'))
        if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
            applycal(vis=phasecal_ms, gaintable = prev_caltables, parang=False )
    
        imagename = phasecal_ms.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
        if os.path.exists(imagename):
            print("Continuing to the next image")
        
        else:
            imagename =  phasecal_ms.replace('.ms','')+f'_selfcal_loop_{selfcal_loop}'
            print(f"Making image {imagename}")

            wsclean_cmd = ['wsclean', '-log-time', '-auto-threshold',f'{threshold[selfcal_loop]}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
                        '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
                        '-mgain', '0.8', '-niter', f'{niter}', f'{phasecal_ms}']

            run_wsclean(wsclean_cmd)

            wsclean_fitsfile = imagename+'-image.fits'
            get_im_stats(wsclean_fitsfile)
            plot_fits(wsclean_fitsfile)
          
            model_fits = imagename.replace('-image.fits','-model.fits')

            print(f"Adding modelcolumn to data. Using {model_fits} to predict")
            predict_cmd = ['wsclean', '-log-time', '-predict', '-reorder' ,'-name', f'{imagename}', phasecal_ms]

            # Predicting
            run_wsclean(predict_cmd)

            # Plot the model column
            plotms(
                vis=phasecal_ms, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
            )
            if calmode[selfcal_loop] == 'p':
                minblperant = 3
            else:
                minblperant = 4

            gaincal(vis =phasecal_ms, caltable = caltable, refant = refant, solint = solint_selfcal[selfcal_loop],
                    gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
                    calmode = calmode[selfcal_loop], append=False, parang=False, minblperant=minblperant
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
                applycal(vis=phasecal_ms, gaintable = prev_caltables, parang=False )
        
        # ### Get the last imagename from the loop and generate a final mask
        
    imagename_final = phasecal_ms.replace('.ms','')+f'final_map_loop_{nloops-1}'
    ##  tclean here to make the final image
    print("Make final image with all selfcal corrections applied")
    wsclean_cmd_final = ['wsclean', '-log-time', '-auto-threshold',f'{threshold_final[0]}', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename_final}', \
            '-scale', f'{cell}', '-fits-mask', f'{maskfile}',\
            '-mgain', '0.8', '-niter', f'{niter_final}', f'{phasecal_ms}']

    run_wsclean(wsclean_cmd_final)

    wsclean_fitsfile = imagename_final+'-image.fits'
    get_im_stats(wsclean_fitsfile)
    plot_fits(wsclean_fitsfile)



# def applycal_target():

#     """
#     Applies cal to the target field
#     """
#     prev_caltables = sorted(glob.glob('*.tb'))
#     applycal(
#         vis = vis_tocal, gaintable = prev_caltables, parang=False
#     )




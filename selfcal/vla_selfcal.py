
import os, glob, subprocess
from casatasks import *
from casaplotms import *
import bdsf
import casalogger
import casatools
import numpy as np
import matplotlib.pyplot as plt
import time 
from radio_beam import Beam
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

msmd = casatools.msmetadata()
tb = casatools.table()
ms = casatools.ms()


msname = '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/23B-307.sb44594812.eb44725045.60239.588568113424/23B-307.sb44594812.eb44725045.60239.588568113424.ms'
# msname = '/raid1/scratch/kelvinw/k2_18b/official_pipe_cal/23B-307.sb44594812.eb44691528.60230.613198356485/23B-307.sb44594812.eb44691528.60230.613198356485.ms'


working_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir/23B-307.sb44594812.eb44725045.60239.588568113424' # D
# working_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir/23B-307.sb44594812.eb44691528.60230.613198356485' # A to D

target = 'K2-18'
vis = target+'.ms'
outlierfile = '/home/kelvin/Desktop/Synphly/selfcal/outlier_fields.txt'
basename = os.path.splitext(os.path.basename(vis))[0]

# target


# imaging 

nterms = 2
gridder = 'standard'
deconvolver = 'mtmfs'
weighting='briggs'
robust = 0.5
threshold = ['0.05mJy','0.02mJy','0.001mJy'] # in mJy
wprojplanes = 1
outlier_file = '/home/kelvin/Desktop/Synphly/selfcal/outlier_fields.txt'
pblimit = -0.1 # avoid 1,-1 or 0

imsize = [640,640] 
niter = [1000,10000,100000] # the number of iterations for each loop -- needs to be arbitrarily large

# selfcal
nloops = 3 # number of selfcal loops
# loop = 0 # large image for selfcal part 1
calmode = ['p','p','ap']
gaintype= ['G','G','G']
solint = ['60s','30s','180s']
minsnr = [1,1,1]

# # pybdsf
detection_threshold = 5.0

# # final image and peeling
niter_final = 1000000
threshold_final = '0.05mJy'
wsclean_sif= '/raid1/scratch/kelvinw/singularity_containers/wsclean_working.simg'
# singularity_bind = '/share/nas/kelvinw/'

# spw = 17 # wsclean chan out
abs_mem = 128 # mem to use in GB

def set_working_dir():

    if not os.path.exists(working_directory):
        print(f"{working_directory} does not exist, making one")
        os.makedirs(working_directory)
    else:
        print(f"Working directory {working_directory} already exists")

    print(f"Changing cwd to {working_directory}")
    os.chdir(working_directory)

def split_data():
    
    listobs(vis=msname, listfile='listobs.txt',overwrite=True)
    msmd.open(msname)
    field_names = msmd.fieldnames()
    print(f"Field names: {field_names} found in {msname}")
    for field in field_names:
        outputvis = field+'.ms'
        if not os.path.exists(outputvis):
            print(f"Splitting {msname} to {outputvis}")
            split(vis=msname,outputvis=outputvis,datacolumn='corrected',timebin='10s',width=4,field=field)
            print(f"Finished splitting")
    msmd.close()


def time_execution(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        
        if execution_time < 60:  # If execution time is less than a minute
            time_unit = "seconds"
            formatted_time = execution_time
        elif execution_time < 3600:  # If execution time is less than an hour
            time_unit = "minutes"
            formatted_time = execution_time / 60
        else:  # If execution time is an hour or more
            time_unit = "hours"
            formatted_time = execution_time / 3600
            
        print(f"======>>>EXECUTION TIME for {func.__name__}: {formatted_time:.2f} {time_unit}")
        return result
    return wrapper



def check_longest_baseline():
    """
    Calculate the longest baseline in terms of wavelength (lambda).
    
    Parameters:
        vis (str): Path to the visibility data file.
    
    Returns:
        float: Longest baseline in units of wavelength.
    """
    # Open measurement set and retrieve uvw data
    ms.open(vis)
    ms.selectinit(datadescid=0)
    uvw = ms.getdata('uvw')['uvw']
    ms.close()
    
    # Compute baseline in meters
    uvdist_meters = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
    longest_baseline_meters = np.nanmax(uvdist_meters)
    
    # Get frequency data
    band_name, mean_freq, max_freq, min_freq = check_band()
    frequency_hz = max_freq * 1e9
    speed_of_light = 3e8  # Speed of light in m/s
    wavelength_meters = speed_of_light / frequency_hz
    
    # Calculate longest baseline in terms of wavelength
    longest_baseline_lambda = longest_baseline_meters / wavelength_meters
    return longest_baseline_lambda


def check_band():
    """
    Identify the frequency band of the data and return relevant frequency information.
    
    Parameters:
        vis (str): Path to the visibility data file.
    
    Returns:
        tuple: Band name, mean frequency, maximum frequency, and minimum frequency (in GHz).
    """
    band_name = None
    freq_ranges = {(1, 2): "L",(2, 4): "S",(4, 8): "C",(8, 12): 
                   "X",(12, 18): "U",(18, 26.5): "K", (26.5, 
                    40): "A",(40, 50): "Q",
                            }
    
    # Open measurement set metadata
    vis=target+'.ms'
    msmd.open(vis)
    nspw = msmd.nspw()
    
    # Calculate mean frequency for each spectral window
    spws_freq = np.array([np.nanmean(msmd.chanfreqs(spw)) for spw in range(nspw)])
    msmd.done()
    
    # Calculate mean, max, and min frequencies across all spectral windows in GHz
    mean_freq = np.nanmean(spws_freq) * 1e-9
    max_freq = np.nanmax(spws_freq) * 1e-9
    min_freq = np.nanmin(spws_freq) * 1e-9
    
    # Identify band based on mean frequency
    for freq_range, band in freq_ranges.items():
        if freq_range[0] <= mean_freq <= freq_range[1]:
            band_name = band
            break

    print(f"Band: {band_name}, Mean Frequency: {mean_freq:.2f} GHz, "
          f"Max Frequency: {max_freq:.2f} GHz, Min Frequency: {min_freq:.2f} GHz")
    
    return band_name, mean_freq, max_freq, min_freq


def get_imaging_cellsize():

    longest_baseline_lambda = check_longest_baseline()
    cell_float = (180.0 * 3600 / (np.pi * 5)) * (1.0 /longest_baseline_lambda)
    cell = f'{cell_float:.2f}arcsec'
    print(f"Imaging with a cell of size {cell}")
    return cell


def plot_fits(imagename,fig_width =8,fig_height=6):
    """
    Plots fitsfiles using astropy
    """

    fitsname = imagename+'.fits'
    exportfits(imagename=imagename+'.image.tt0',fitsimage=fitsname,overwrite=True)

    plt.figure(figsize=(fig_width, fig_height))  
    hdul = fits.open(fitsname)
    w = WCS(hdul[0].header, naxis=2)
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    ax = plt.subplot(projection=w)

    # Disable automatic labelling
    ax.coords[0].set_auto_axislabel(True) 
    ax.coords[1].set_auto_axislabel(True) 

    # Extract image data
    
    image_data = hdul[0].data[0,0,:,:]

    # Display the image with a color map
    im = ax.imshow(image_data, cmap=plt.get_cmap('cividis'))

    # Automatically set pixel scale based on WCS header information
    pixscale = abs(w.wcs.cdelt[0]) * u.deg.to(u.arcsec) * u.arcsec  # Convert from degrees to arcseconds

    # Define and add the beam ellipse
    array_beam = imhead(imagename+'.psf.tt0')['restoringbeam']
    major_axis = array_beam['major']['value']
    minor_axis = array_beam['minor']['value']
    pos_angle = array_beam['positionangle']['value']

    major_axis = major_axis*u.arcsec  # Convert to arcseconds if needed 
    minor_axis = minor_axis*u.arcsec  # Convert to arcseconds if needed
    pos_angle = pos_angle*u.deg

    my_beam = Beam(major_axis, minor_axis, pos_angle)
    ycen_pix, xcen_pix = 15, 15
    ellipse_artist = my_beam.ellipse_to_plot(xcen_pix, ycen_pix, pixscale)
    _ = ax.add_artist(ellipse_artist)

    ax.set_xlabel('RA (J2000)',size=14)
    ax.set_ylabel('Dec (J2000)',size=14)  

    ax.tick_params(axis = "x", which = "both", bottom = True, top = False)
    ax.tick_params(axis = "y", which = "both", right = False, left = True)

    # ra = ax.coords[0]
    # dec = ax.coords[1]

    # ra.set_ticklabel(size=12)
    # dec.set_ticklabel(size=12)

    cbar = plt.colorbar(im,extend='both')
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Jy/beam',rotation=90, labelpad=12,size=18)

    # ax.contour(image_data,levels=[-3*0.136e-3,3*0.136e-3,5*0.136e-3,10*0.136e-3,15*0.136e-3], colors='white',
    #         linewidths=0.5)

    plt.savefig(fitsname.replace('.fits','_1.pdf'),dpi=300)

# def get_im_stats(imagename):
    
#     """
#     Gets the statistics for either a 256x256 pix image and writes
#     them to a logfile
#     """


#     rms=imstat(imagename=imagename,box='60,60,580,240')['rms'][0]  # for 640x640 px
#     peak=imstat(imagename=imagename,box='300,300,340,340')['max'][0]
#     print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
#                 (imagename, peak*1e3, rms*1e3, peak/rms))
    
#     logfile = 'imstat.txt'
#     casa_imstat = imstat(imagename)
#     with open(logfile,"a") as txt_file:
#         txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
#                     (imagename, peak*1e3, rms*1e3, peak/rms))

#         txt_file.write(f"For {imagename}, the maximum pos for imstat is {casa_imstat['maxposf']}\n")


# def make_dirty_map():
#     cell = get_imaging_cellsize()
#     msmd.open(vis)
#     field_names = msmd.fieldnames()
#     for field in field_names:
#         imagename = field+'_dirty'
#         vis_to_image = field+'.ms'
#         print(f"Imaging {vis_to_image}")
#         tclean(
#             vis = vis_to_image, imagename = imagename, imsize=[5120,5120], cell=cell, gridder=gridder,
#             deconvolver = deconvolver, weighting = weighting,
#             robust = robust, niter=0, nterms=nterms, interactive=False, 
#         )
#         print(f"Finished imaging {vis_to_image}")


def find_refant(vis, field=target, tablename=target+'.refant'):
    """
    This function comes from the e-MERLIN CASA Pipeline.
    https://github.com/e-merlin/eMERLIN_CASA_pipeline/blob/master/functions/eMCP_functions.py#L1501

    It finds the best reference antenna for calibration.

    Parameters
    ----------
    msfile : str
        The measurement set file
    field : str
        The field to calibrate
    tablename : str
        The name of the calibration table
    """
    # Find phase solutions per scan:
    if not os.path.exists(tablename):
        gaincal(vis=vis,caltable=tablename,field=field,refantmode='flex', solint='inf',minblperant=3,gaintype='G', calmode='p')
    # find_casa_problems()
    # Read solutions (phases):

    tb.open(tablename + '/ANTENNA')
    antenna_names = tb.getcol('NAME')
    tb.close()
    tb.open(tablename)
    antenna_ids = tb.getcol('ANTENNA1')
    # times  = tb.getcol('TIME')
    flags = tb.getcol('FLAG')
    phases = np.angle(tb.getcol('CPARAM'))
    snrs = tb.getcol('SNR')
    tb.close()
    # Analyse number of good solutions:
    good_frac = []
    good_snrs = []
    for i, ant_id in enumerate(np.unique(antenna_ids)):
        cond = antenna_ids == ant_id
        # t = times[cond]
        f = flags[0, 0, :][cond]
        p = phases[0, 0, :][cond]
        snr = snrs[0, 0, :][cond]
        frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
        snr_mean = np.nanmean(snr[~f])
        good_frac.append(frac)
        good_snrs.append(snr_mean)
    sort_idx = np.argsort(good_frac)[::-1]
    print(' ++==>> Antennas sorted by % of good solutions:')
    for i in sort_idx:
        print(' ++==>> {0:3}: {1:4.1f}, <SNR> = {2:4.1f}'.format(antenna_names[i],
                                                         good_frac[i],
                                                         good_snrs[i]))
    if good_frac[sort_idx[0]] < 90:
        print(' ++==>> Small fraction of good solutions with selected refant!')
        print(' ++==>> Please inspect antennas to select optimal refant')
        print(' ++==>> You may want to use refantmode="flex" in gaincal')
    pref_ant = antenna_names[sort_idx]
    pref_ant_list = ','.join(list(pref_ant))
    print(f"The following will be used as the reference antennas: {pref_ant_list}")
    return pref_ant_list

def pybdsf(input_image):

    # The input image is a casa .image that then gets exported to a FITS
    imagename = input_image.replace('.image.tt0','')
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

def make_dirty_map():

    """
    Creates an (a large) an image that is used to create a casa region file using pybdsf 
    for masking
    """
    cell = get_imaging_cellsize()
    dirty_map = basename + '_first_masking'

    if not os.path.exists(dirty_map):
        print(f"Making {dirty_map}")
        tclean(
            vis = vis, imagename=dirty_map, imsize=imsize, cell=cell,
            gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
            weighting = weighting, robust = robust, niter=0, # threshold = '0.5mJy',
            nterms = nterms, pblimit = pblimit
        )

    regionfile = pybdsf(input_image=dirty_map+'.image.tt0')


def selfcal():

    if os.path.exists(outlier_file) and open(outlier_file).read() == '':
        outlierfile = ''

    regionfile = basename + '_first_masking.casabox'
    cell = get_imaging_cellsize()
    print("Deleting model column before selfcal")
    delmod(vis=vis,otf=True)

    for selfcal_loop in range(nloops):
        caltable = f'caltable_{selfcal_loop}.gcal'
        prev_caltables = sorted(glob.glob('*.gcal'))
        if len(prev_caltables) >0 and calmode[selfcal_loop] !='':
            applycal(vis=vis, gaintable = prev_caltables, parang=False )
    
        imagename = f'target_selfcal_{selfcal_loop}'
        if os.path.exists(imagename):
            print("Continuing to the next image")
        
        else:
            # imagename = f'target_selfcal_{selfcal_loop}'
            print(f"Making image {imagename}")
            tclean(
                vis = vis, imagename=imagename, imsize=imsize, cell=cell,
                parallel=False,
                gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
                weighting = weighting, robust = robust, niter=niter[selfcal_loop], threshold = threshold[selfcal_loop],
                nterms = nterms, pblimit = -1,interactive=False, usemask='user', mask=regionfile
            )

            ## NB: The problem was niter -- there was a space in the list []

            print("Adding modelcolumn to data")
            # model images from the MTMFS images,
            ft(vis = vis, model=[imagename+'.model.tt0',imagename+'.model.tt1'], nterms=2,usescratch=True)

            # plot the model column
            plotms(
                vis=vis, xaxis='UVwave', yaxis='amp', ydatacolumn='model',avgchannel='64',avgtime='300',
                showgui=False, plotfile=imagename+'_modelcolumn.png', overwrite=True, width=1500, height=750,
            )

            gaincal( vis =vis, caltable = caltable, refant = refant, solint = solint[selfcal_loop],
                    gaintype = gaintype[selfcal_loop], gaintable=prev_caltables,  minsnr = minsnr[selfcal_loop],
                    calmode = calmode[selfcal_loop], append=False, parang=False
                    )
            coloraxis = ['corr','spw']
            for color in coloraxis:
                if calmode[selfcal_loop] =='p':
                    plotms(
                        vis = caltable, xaxis='time', yaxis='phase', gridcols=3, gridrows=3,
                        iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
                        plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750,
                    )
                else:
                    plotms(
                            vis = caltable, xaxis='time', yaxis='amp', gridcols=3, gridrows=3,
                            iteraxis='antenna', coloraxis = color, showgui=False, overwrite=True,
                            plotfile=caltable.replace('.gcal',f'_{color}.png'), dpi=300, width=1500, height=750
                        )

            if selfcal_loop == nloops-1:
                prev_caltables = sorted(glob.glob('*.gcal'))
                print("Applying the caltable derived from last gaincal iteration")
                applycal(vis=vis, gaintable = prev_caltables, parang=False )
        
        # ### Get the last imagename from the loop and generate a final mask
        
    imagename = basename +f'_{nloops-1}'+'.final'
    ##  tclean here to make the final image
    print("Make final image with all selfcal corrections applied")
    tclean(
        vis = vis, imagename = imagename, imsize=imsize, cell=cell, gridder=gridder,
        wprojplanes = wprojplanes, deconvolver = deconvolver, weighting = weighting,
        robust = robust, niter=niter_final, threshold = threshold_final, nterms=nterms,
        pblimit=pblimit, interactive=False, usemask = 'user', mask=regionfile,
    )
    
    


# def applycal_target():

#     """
#     Applies cal to the target field
#     """
#     prev_caltables = sorted(glob.glob('*.gcal'))
#     applycal(
#         vis = vis_tocal, gaintable = prev_caltables, parang=False
#     )

def run_wsclean(command):

    """
    Runs wsclean commands 
    """

    container = wsclean_sif
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



# def peeling():

#     """
#     Subtract the bright sources in the field so you are left with the emission from the star

#     Use pybsdf casa region file -- you can change the threshold to only select the bright sources you want

#     Run pybsdf on the final image
    
#     """

#     # Using the final imagename 

#     imagename = basename +f'_{nloops-1}'+'.final'
#     # region_to_peel = pybdsf(input_image=imagename+'.image.tt0')

#     ## NB: You had masked the bright source from the dirty map -- so you can just subtract
#     uvsub(vis=vis)

#     ## Make a final map without the sources
#     cell = get_imaging_cellsize()
#     peeled_map = basename+'_peeled_map'
#     if not os.path.exists(peeled_map):
#         print(f"Making {peeled_map}")
#         tclean(
#             vis = vis, imagename=peeled_map, imsize=imsize, cell=cell,
#             gridder = gridder, wprojplanes = wprojplanes, deconvolver = deconvolver,
#             weighting = weighting, robust = robust, niter=0, # threshold = '0.5mJy',
#             nterms = nterms, pblimit = pblimit
#         )








    # delete model column to be d

def peeling():

    """
    Subtract all the sources in the field such that you are left with a blank

    Use wsclean and pybdsf casa region from the final iterations 

    Then perform uvsub in CASA

    """

    # Make a region file of the final self calibrated image and use it to peel the sources
    imagename = basename +f'_{nloops-1}'+'.final.image.tt0'
    regionfile_to_peel = pybdsf(input_image=imagename)
    fitsmask = imagename.replace('.image.tt0','')+'.maskfile.fits'
    model_fits = imagename.replace('.final.image.tt0','.final.image.tt0-model.fits')
    os.rename(fitsmask,model_fits )

    cell = get_imaging_cellsize()


    threshold_cmd = ['wsclean', '-auto-threshold','3', '-size', f'{imsize[0]}', f'{imsize[1]}','-scale', f'{cell}',\
                    '-mgain', '0.8', '-niter', '0',f'{vis}']
    
    predict_cmd = ['wsclean', '-log-time', '-predict', '-field', '', '-reorder' ,'-name', f'{imagename}', '-abs-mem',f'{abs_mem}', vis]


    run_wsclean(predict_cmd)

    ## NB: wsclean needs to find an image named my-image-model.fits or reg 
    ## works by replacing model column with model for the problem sources using

    
    
    ## Subtract the models put in the model column from the data and make an image
        
    print("Running uvsub")
    uvsub(vis=vis)

    ## Run wsclean to check if the subtraction has been successful -- make dirty map

    run_wsclean(threshold_cmd)


def main():

    set_working_dir()
    split_data()
    check_band()
    global refant
    refant = find_refant(vis)
    make_dirty_map()
    selfcal()
    # applycal_target()
    peeling()

    # imagename = 'K2-18_dirty'
    # plot_fits(imagename)



main()

        

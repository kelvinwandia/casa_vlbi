import os,subprocess,time
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from astropy.io import fits
import casatasks
import casatools
from astropy.coordinates import SkyCoord
import astropy.units as u
import logging

from casatools import componentlist, image, quanta, table
cl = componentlist()
ia = image()
qa = quanta()
tb = table()

working_dir = os.path.join('/raid1/scratch/kelvinw/gv020_working_dir/gv020a_working_dir/simulations')
os.chdir(working_dir)

msname = 'j2139+1423_split.ms'

pointing_centre = ['21h39m01.309269s +14d23m35.99221s']
phasecenter =['21h39m05.309269s +14d20m34.99221s']
# phasecenter = pointing_centre

# phasecenter = ['21h30m25.2034s +12d13m18.160s'] #  7.2 arcmin from the pointing centre
starting_freq='1602.5056MHz' # central channel
freq_increment = '31.25kHz'
num_chan=8*512
imsize=[256,256]
cellsize='1mas'



def run_wsclean(wsclean_sif,command):

    """
    Runs wsclean commands 
    """
    container = wsclean_sif
    if os.path.exists(container):
        singularity_bind = os.path.join(os.path.dirname(os.path.dirname(container)))

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


def get_number_of_threads():
    try:
        num_threads = os.cpu_count()
        if num_threads is None:
            print("Could not determine the number of threads.")
        else:
            print(f"Number of threads (logical processors) available: {num_threads}")
    except Exception as e:
        print(f"An error occurred while determining the number of threads: {e}")
    
    return num_threads

def angsep(offset_source):

    """
    Calculates the angular separation between provided coordinates and the pointing centre and
    converts the coordinates from equatorial to rad using astropy
    # follow this link for the reference eqn
    # https://link.springer.com/referenceworkentry/10.1007/978-0-387-35973-1_761

    Args:
        offset_source: the equatorial coordinates (J2000) of the source

    Returns:
        separation_angle: the angular separation
    """

    separation_angle =  []

    ra_pointing_centre, dec_pointing_centre = pointing_centre[0].split()
    
    pointing_centre_skycoord_obj = SkyCoord(ra_pointing_centre,dec_pointing_centre, unit=(u.hourangle,u.deg),frame='icrs')
    
    for source_coords in offset_source:
        
        source_coords_skycoord_obj = SkyCoord(source_coords.split()[0],source_coords.split()[1],unit=(u.hourangle,u.deg),frame='icrs')
        
        ## This formular needs fixing -- has been fixed
        vector_dot_product =  np.cos(pointing_centre_skycoord_obj.dec.radian)*np.cos(source_coords_skycoord_obj.dec.radian)*\
            np.cos( pointing_centre_skycoord_obj.ra.radian- source_coords_skycoord_obj.ra.radian) + \
            np.sin(pointing_centre_skycoord_obj.dec.radian)*np.sin(source_coords_skycoord_obj.dec.radian)
       
        angle = np.arccos(vector_dot_product)
        # print(f"The offset angle is {angle*180*60/np.pi} arcmin")

        offset_angles_dict = {
            'coordinates':source_coords_skycoord_obj.to_string('hmsdms'),
            'angular_separation_arcmin': angle*180*60/np.pi, # arcmin
            'angular_separation_radians': angle
        }

        separation_angle.append(offset_angles_dict)
    print(separation_angle)
    return separation_angle



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
    image_plot = ax.imshow(image_data, origin='lower', 
                       extent=[-32, 32, -32, 32],cmap='viridis')
    cbar = plt.colorbar(image_plot,ax=ax,orientation='vertical')
    plt.savefig(fitsname.replace('.fits','.pdf'))

def get_im_stats(imagename):
    
    """
    Gets the statistics for either a 256x256 pix image and writes
    them to a logfile
    """


    rms=casatasks.imstat(imagename=imagename,box='51,7,247,76')['rms'][0]  # for 256x256 px
    peak=casatasks.imstat(imagename=imagename,box='124,122,133,134')['max'][0]
    print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))

    logfile = 'imstat.txt'
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))

    

def set_working_dir():

    current_dir = os.getcwd()
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    else:
        print(f"Working dir {working_dir} exists")

    os.chdir(working_dir)
    # return current_dir



def direction_string(ra, dec, frame):

    """helper function for often needed string"""

    return ' '.join([frame, ra, dec])


def addCompandPredictVis():


    ra_dir,dec_dir = phasecenter[0].split(' ')
    
    clname = f"{ra_dir}_{dec_dir.replace('+',' ')}_model.cl"
    os.system(f"rm -r {clname}")
    if not os.path.exists(clname):
        cl.addcomponent(
            dir=direction_string(ra_dir,dec_dir,frame='J2000'),flux=1.0,fluxunit='Jy',
            freq=starting_freq,shape='point'
        )
        cl.rename(clname)
        cl.done()

        print("Adding model to data using ft")
        casatasks.ft(
            vis=msname, complist=clname,incremental=False,usescratch=True
            )

        print("Copying model to DATA column")
        tb.open(msname,nomodify=False)
        moddata = tb.getcol(columnname='MODEL_DATA')
        tb.putcol(columnname='DATA',value=moddata)
        moddata.fill(0.0)
        tb.putcol(columnname='MODEL_DATA',value=moddata)
        tb.close()


def transform_and_image():

    ra_dir,dec_dir = phasecenter[0].split(' ')

    outputms = msname.replace('.ms','_added_source')+'.ms'
    subprocess.run(['rm','-r',outputms])
    if not os.path.exists(outputms):
        print(f"Phase shifting the data to {outputms}")
        casatasks.phaseshift(
            vis=msname,outputvis=outputms,datacolumn='data',
            phasecenter=direction_string(ra_dir,dec_dir,frame='J2000')
        )

    transformed_ms = 'transformed_ms'+'.ms'
    subprocess.run(['rm','-r',transformed_ms])
    if not os.path.exists(transformed_ms):
        print("Transforming the measurement set")
        casatasks.mstransform(
            vis=outputms,outputvis=transformed_ms,
            datacolumn='data',timeaverage=True, timebin='20s',
            createmms = True)
        # casatasks.split(vis=outputms,outputvis=transformed_ms,datacolumn='data',
        #       timebin='20s')
    else:
        print(f"{transformed_ms} exists. No action required")

    ### Imaging
    imagename=f"{ra_dir}_{dec_dir.replace('+','')}"
    print(f"Making image {imagename}")
    if not os.path.exists(imagename):
        os.system(f"rm -r {imagename}.*")

    # print(f"Imaging {transformed_ms}")
    # casatasks.tclean(
    #     vis=msname, imagename=imagename,cell=cellsize, niter=0,
    #     imsize=imsize,parallel=False,
    #     # phasecenter=direction_string(ra_dir,dec_dir,frame='J2000')
    #     )

    num_threads = get_number_of_threads()
    wsclean_cmd = ['wsclean', '-log-time','-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}','-scale', f'{cellsize}',\
                                '-mgain', '0.8', '-j', f'{num_threads}','-log-time', '-niter', '1' ,'-v',f'{transformed_ms}']
    
    wsclean_sif= '/raid1/scratch/kelvinw/singularity_containers/wsclean_working.simg'

    run_wsclean(wsclean_sif,wsclean_cmd)
    # casatasks.exportfits(imagename=imagename+'.image',fitsimage=imagename+'.fits',overwrite=True)
    # get_im_stats(imagename+'.image')
    get_im_stats(imagename+'-image.fits')
    # plot_fits(imagename+'.fits') # dont plot fits - imshow requires X environemt

    return imagename

def summarize(imagename):

    # imagename = imagename+'.image'
    imagename=imagename+'-image.fits'
    x = casatasks.imstat(imagename)
    imfit_box = '124,122,133,134'
    y = casatasks.imfit(imagename,box=imfit_box)

    poserr = y['deconvolved']['component0']['shape']['direction']['error']
    print(f"The simulated source position is {phasecenter}")
    print(f"The coordinates of max position from imstat are {x['maxposf']}")
    cl.fromrecord(y['deconvolved'])
    rd = cl.getrefdir(0)
    cl.done()
    prec=5
    ra_err = qa.div(
    qa.div(qa.quantity(poserr['longitude']), 15),
    qa.cos(qa.quantity(rd['m1'])))
    ra_err['unit'] = 's'
    dec_err = qa.quantity(poserr['latitude'])
    print(
        "fitted position from imfit",
        qa.time(qa.totime(qa.quantity(rd['m0'])), prec=6+prec, form='hms')[0], '\u00b1',
        qa.tos(ra_err, prec=prec),
        qa.angle(qa.totime(qa.quantity(rd['m1'])), prec=6+prec)[0], '\u00b1',
        qa.tos(dec_err, prec=prec),
    )

    summary_file = 'summary_pos.txt'
    with open(summary_file, 'a') as txtfile:
        txtfile.write(
            "fitted position from imfit " +
            qa.time(qa.totime(qa.quantity(rd['m0'])), prec=6+prec, form='hms')[0] + '\u00b1' +
            qa.tos(ra_err, prec=prec) + ' ' +
            qa.angle(qa.totime(qa.quantity(rd['m1'])), prec=6+prec)[0] + '\u00b1' +
            qa.tos(dec_err, prec=prec) + '\n'
        )


start = time.time()
separation_angle = angsep(phasecenter)
addCompandPredictVis()
imagename=transform_and_image()
end = time.time()
print(f"It takes {((end-start)/3600):.2f} hours to simulate,predict and image source")
summarize(imagename)

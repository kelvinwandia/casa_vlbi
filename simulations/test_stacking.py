import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import os
from casatools import (simulator, image, table, coordsys, measures,
                componentlist, quanta, ctsys)
from casatasks import (tclean, ft, imhead, listobs, exportfits, flagdata,
                bandpass,imstat, applycal)
from casatasks.private import simutil
import os, re, sys, json
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import casatasks
import casatools


def generate_random_coordinates(ra_center, dec_center, min_separation_arcmin, num_points, seed=None):

    if seed is not None:
        np.random.seed(seed)

    # Convert minimum separation to degrees
    min_separation_deg = min_separation_arcmin / 60.0

    # Convert the radius to degrees
    max_radius_deg = 1 / 60.0

    # Convert the center coordinates to a SkyCoord object
    center_coord = SkyCoord(ra=ra_center, dec=dec_center, unit=(u.deg, u.deg), frame='icrs')

    # List to hold the generated coordinates
    coordinates = []

    # Generate the coordinates ensuring sufficient separation and within the radius
    while len(coordinates) < num_points:
        # Generate random points within a circle of radius 1 arcmin
        r = max_radius_deg * np.sqrt(np.random.uniform(0, 1))
        theta = np.random.uniform(0, 2 * np.pi)
        delta_ra_deg = r * np.cos(theta)
        delta_dec_deg = r * np.sin(theta)

        # Create the new coordinate by applying the offsets
        new_coord = SkyCoord(ra=center_coord.ra.deg + delta_ra_deg,
                             dec=center_coord.dec.deg + delta_dec_deg,
                             unit=u.deg,
                             frame='icrs')

        # Check the separation with all existing points
        sufficient_separation = True
        for coord in coordinates:
            sep = new_coord.separation(coord).deg
            if sep < min_separation_deg:
                sufficient_separation = False
                break

        # If the point has sufficient separation, add it to the list
        if sufficient_separation:
            coordinates.append(new_coord)

    return coordinates

def plot_coordinates(coordinates, ra_center_deg, dec_center_deg):
    # Extract RA and Dec from SkyCoord objects
    ra_degrees = [coord.ra.deg for coord in coordinates]
    dec_degrees = [coord.dec.deg for coord in coordinates]

    # Create the plot
    plt.figure(figsize=(8, 8))
        # Plot the circle with 1 arcmin radius
    circle = plt.Circle((ra_center_deg, dec_center_deg), 1 / 60.0, color='green', fill=False, label='1 arcmin radius')
    plt.gca().add_patch(circle)
    plt.scatter(ra_degrees, dec_degrees, color='blue', marker='o', label='Random Points')
    plt.scatter(ra_center_deg, dec_center_deg, color='red', marker='x', label='Center')



    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    plt.title('Random Coordinates with Minimum Separation')
    plt.legend()
    plt.grid(True)
    plt.gca().invert_xaxis()  
    plt.tight_layout()
    plt.savefig('coords.pdf')
    # plt.show()

# Example usage
ra_center = 322.49304  # RA in degrees
dec_center = 12.16700  # Dec in degrees
min_separation_arcmin = 0.1  # Minimum separation in arcminutes
num_points = 6
seed = 42

# coordinates = generate_random_coordinates(ra_center, dec_center, min_separation_arcmin, num_points, seed)

# Print the generated coordinates
# for coord in coordinates:
#     print(coord.to_string('hmsdms'))

# Plot the coordinates with the center and the circle
# plot_coordinates(coordinates, ra_center, dec_center)



sm = simulator()
ia = image()
tb = table()
cs = coordsys()
me = measures()
qa = quanta()
cl = componentlist()
mysu = simutil.simutil()


observatory = 'EVN'
integration_time = '0.5s'
starting_freq = '1602.5056MHz'
# channel_width ='31.250kHz'
channel_width='32.5kHz'
freq_resolution = channel_width
freq_increment = channel_width
num_channels=8*512
antenna_cfg_path = '/raid1/scratch/kelvinw/gv020/simulations/config_files/'
# antenna_cfg_path = './config_files/'
cfg_file = observatory.lower()+'.cfg'
antennalist = antenna_cfg_path+cfg_file

msname = observatory.lower()+'_'+integration_time+'_'+channel_width+'.ms'

pointing_centre = '21h29m58.3500s +12d10m01.500s'
phasecenter=['21h29m58.3500s +12d10m01.500s']
ra_dir, dec_dir = pointing_centre.split(' ')

starttime='0h'
stoptime='1h'
position_observatory = me.observatory(observatory.lower())
imsize = [256,256]
cellsize = '1mas'

msname = 'trial.ms'

def direction_string(ra, dec, frame):

    """helper function for often needed string"""

    return ' '.join([frame, ra, dec])


def makeEmptyMS():


    os.system('rm -rf ' + msname)
    os.system('rm -rf *.ms*')

    sm.open(ms=msname)
    (x,y,z,diam,antname,padnames,telescope,posobs) = mysu.readantenna(antennalist)

    antname =['JB','WB','EF','ON','MC','TR','NT','AR','GB']
    padnames = antname
    sm.setconfig(
        telescopename=telescope, x=x, y=y, z=z, dishdiameter=diam.tolist(),
        mount='alt-az',antname=antname,padname=padnames,coordsystem='global',
        referencelocation=position_observatory,
        )

    sm.setfeed(mode='perfect R L',pol=[''])

    sm.setspwindow(
        spwname='Spw 0', freq=starting_freq, deltafreq=freq_increment,
        freqresolution=freq_resolution,refcode='GEO',nchannels=num_channels,stokes='RR LL',
        )
    sm.setfield(
        sourcename=pointing_centre, sourcedirection=me.direction(rf='J2000',v0=ra_dir,v1=dec_dir))
    sm.setlimits(shadowlimit=0.001,elevationlimit='10.0deg')
    sm.setauto(autocorrwt=0.0)
    sm.settimes(
            integrationtime=integration_time,usehourangle=True,
            referencetime=me.epoch('UTC','2010/11/09/18:00:00'),
        )
    sm.observe(sourcename=pointing_centre,spwname='Spw 0',starttime=starttime,stoptime=stoptime)
    sm.close()
    listobs_file = observatory.lower()+'.listobs.txt'

    listobs(vis=msname,listfile=listobs_file,overwrite=True)

    file = open(listobs_file)
    for line in file.readlines():
        print(line.replace('\n',''))
    file.close()

def addCompandPredictVis():

    clname = 'model_component.cl'

    os.system('rm -r *.cl')

    ra_dir,dec_dir = phasecenter[0].split(' ')
    print(f"The coordinates are {ra_dir} and {dec_dir}")

    print(f"Making component {clname}")
    cl.addcomponent(
        dir=direction_string(ra_dir,dec_dir,frame='J2000'),flux=5.0,fluxunit='Jy',
        freq=starting_freq,shape='point'
    )
    cl.rename(clname)
    cl.done()


    casatasks.ft(
        vis=msname,complist=clname,incremental=False,usescratch=True
    )

    print("Copying model to DATA column")
    tb.open(msname,nomodify=False)
    moddata = tb.getcol(columnname='MODEL_DATA')
    tb.putcol(columnname='DATA',value=moddata)
    moddata.fill(0.0)
    tb.putcol(columnname='MODEL_DATA',value=moddata)
    tb.close()



    #### Offset source
    # clname2 = 'offset_source.cl'
    # cl.addcomponent(
    #     dir=direction_string('21h29m45s','12d10m01.30s',frame='J2000'),flux=1.0,fluxunit='Jy',
    #     freq=starting_freq,shape='point'
    # )
    # cl.rename(clname2)
    # cl.done()

    # casatasks.ft(
    #     vis=msname,complist=clname2,incremental=False, usescratch=True
    # )
    # print("Copying model to DATA column")
    # tb.open(msname,nomodify=False)
    # moddata = tb.getcol(columnname='MODEL_DATA')
    # tb.putcol(columnname='DATA',value=moddata)
    # moddata.fill(0.0)
    # tb.putcol(columnname='MODEL_DATA',value=moddata)
    # tb.close()

def get_im_stats(imagename):
    
    """
    Gets the statistics for either a 256x256 pix image and writes
    them to a logfile
    """


    rms=imstat(imagename=imagename,box='51,7,247,76')['rms'][0]  
    peak=imstat(imagename=imagename,box='124,122,133,134')['max'][0]
    print('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                (imagename, peak*1e3, rms*1e3, peak/rms))

    logfile = 'imstat.txt'
    with open(logfile,"a") as txt_file:
        txt_file.write('For %s, the peak %.3f mJy/beam, rms %.3f mJy/beam, S/N %6.0f\n\n' %
                    (imagename, peak*1e3, rms*1e3, peak/rms))
        

def plot_fits(fitsname):
    """
    Plots fitsfiles using astropy
    """
    fitsfile = fits.open(fitsname)
    image_data = fitsfile[0].data[0,0,:,:]
    w = WCS(fitsfile[0].header,naxis=2)
    header = fitsfile[0].header
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=w)
    axis1 = header['NAXIS1']

    im=ax.imshow(image_data, cmap=plt.get_cmap('viridis'),extent=[-axis1/2,axis1/2,-axis1/2,axis1/2]) 
    cbar = plt.colorbar(im,ax=ax,orientation='vertical')
    cbar.set_label('Jy/beam',rotation=90,labelpad=-1)
    cbar.formatter.set_powerlimits((0, 0))

    visible_ticks = {
   "top": False,
   "right": False
        }
    ax.tick_params(axis="x", which="both", **visible_ticks)
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')

    plt.savefig(fitsname.replace('.fits','.pdf'))

## Add Gaussian random noise
def addNoiseSim():
    sm.openfromms(msname)
    sm.setseed(1234)
    sm.setnoise(mode='simplenoise',simplenoise='1Jy')
    sm.corrupt()
    sm.close()

def makeImage():
    imagename = 'image_trial'
    fitsname = imagename+'.fits'
    os.system(f"rm -r {imagename}.*")
    tclean(vis=msname,cell=cellsize,niter=0,imagename=imagename,deconvolver='clark',
           weighting='natural',imsize=[256,256],datacolumn='data')
    exportfits(imagename=imagename+'.image',fitsimage=fitsname,overwrite=True)
    get_im_stats(imagename=imagename+'.image')
    plot_fits(fitsname)
    
makeEmptyMS()
addCompandPredictVis()
addNoiseSim()
makeImage()
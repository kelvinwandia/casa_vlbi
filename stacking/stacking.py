import os, subprocess, logging, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import casatools
from casatasks import *
from astropy.constants import c
from astropy.io import fits
from astropy.wcs import WCS


filepath = '/raid1/scratch/kelvinw/gv020_working_dir/gv020b_aoflagger_working_dir/pbcor_dir'
fitsfiles = glob.glob(os.path.join(filepath,'*-image.fits'))
working_dir = os.path.join(filepath,'stacking_working_dir')

if not os.path.exists(working_dir):
    os.mkdir(working_dir)

os.chdir(working_dir)

casa_images = [os.path.join(filepath, filename) for filename in glob.glob(os.path.join(filepath, '*.image'))]
if casa_images:
    print("CASA images exist, will attempt to convert to fits")
    for image in casa_images:
        exportfits(imagename=image,fitsimage=image.replace('.image','.fits'),overwrite=True)
else:
    print("No CASA images found")


### NB: WSClean images with niter=0 have no beam information
### Convert the beam major axis to pixels
def convert_bmaj_to_pix(bmaj_deg, cdelt1):
    return bmaj_deg / abs(cdelt1)





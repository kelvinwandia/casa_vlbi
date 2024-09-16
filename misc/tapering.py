from casatasks import tclean, imhead
import os, numpy

def convert_lm_to_uv(theta_lm):
    """
    theta_lm: fwhm in the image domain, in units of arcsec
    """

    theta_uv = ( 4*numpy.log(2)/numpy.pi ) / ( (theta_lm / 3600) * numpy.pi/180.0)
    print("FWHM of %3.3f mas maps to a FWHM of %3.3e lambda"%(theta_lm*1000,theta_uv))
    return theta_uv

def convert_uv_to_lm(theta_uv):
    """
    theta_uv : Full width at half maximum, in the uv domain, in units of lambda
    """
    theta_lm  = 3600 * ( 4*numpy.log(2)/numpy.pi ) / (theta_uv * numpy.pi/180.0)
    print("FWHM of %3.3f mas maps to a FWHM of %3.3e lambda"%(theta_lm*1000,theta_uv))
    return theta_lm



def dispbeam(beam):
    """
    Print restoring beam info...
    """
    print("Restoring Beam : %3.4e %s  X  %3.4e %s ,  %3.4f %s"%(beam['major']['value'],
                beam['major']['unit'],
                beam['minor']['value'],
                beam['minor']['unit'],
                beam['positionangle']['value'],
                beam['positionangle']['unit'] ))



def run_im(imnames,uvtapers,weighting='briggs'):
    os.system('rm -rf uvt*')
    vis = '/raid1/scratch/kelvinw/n14c3_working_dir/selfcal_dir/J1849+3024.ms'

    for (imname,uvtaper) in zip(imnames,uvtapers):
        print("\n%s : uvtaper = %s"%(imname,uvtaper))
        tclean(vis=vis, spw='', imagename=imname,
                uvtaper=uvtaper,
                weighting=weighting,
                imsize=200, robust=0.5,
                cell='0.24mas',niter=0, 
                #    restoringbeam=['1.4640e-02arcsec','1.4640e-02arcsec','0deg'],
                calcpsf=True, interactive=False, 
                restoration=True)
        beam =imhead(imname+'.psf')['restoringbeam']
        dispbeam(beam)

imnames = ['uvt_orig' , 'uvt_taper_im' , 'uvt_taper_uv']




def calc_convolve(theta_orig, theta_taper):
    """
    Calculate the width of a Gaussian resulting from the convolution of two Gaussians.
    This calculation is only for a circular Gaussian.
    Units of inputs : arcsec.
    """
    arcsec_to_radians = (1/3600.0)*numpy.pi/180.0
    sigma_orig = arcsec_to_radians * theta_orig/numpy.sqrt(8*numpy.log(2.0))
    sigma_taper = arcsec_to_radians * theta_taper/numpy.sqrt(8*numpy.log(2.0))

    sigma_new = numpy.sqrt(sigma_orig**2 + sigma_taper**2)
    theta_new = sigma_new * numpy.sqrt(8*numpy.log(2.0)) / arcsec_to_radians

    print("Convolution of FWHMs of %3.4f mas and %3.4f mas \
          yields %3.4f mas"%(theta_orig*1000, theta_taper*1000, theta_new*1000))
    
imtaper = 5e-3
uvtapers=['' , '%3.2f mas '%(imtaper*1000) , '%3.2elambda'%(convert_lm_to_uv(imtaper)/2.0)]
print("\nSettings for uvtaper in tclean : \n\
[ None,  FWHM in the image domain, HWHM in the uv-domain] ")
print(uvtapers)

run_im(imnames,uvtapers,weighting='natural')

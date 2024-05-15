from casatasks import *
import casatools
import subprocess, os

# from utils.helper_functions import *


phasecenter=['21h29m58.3500s +12d10m01.500s','21h29m58.3120s +12d10m02.679s','21h30m01.2034s +12d10m38.160s',
        '21h29m51.9025s +12d10m17.132s','21h29m56.3050s +12d11m01.500s','21h29m56.3050s +12d09m11.500s',
        '21h30m02.4410s +12d09m11.500s']


target_ms = '/home/kelvin/Desktop/gv020_working_dir/gv020b/J2139+1423.ms'

# for coord in range(len(phasecenter)):

#     ra_dir, dec_dir = phasecenter[coord].split(' ')
#     transformed_ms = 'transformed_ms_'+phasecenter[coord].replace(" ","")+'.ms'
#     subprocess.run(f"rm -r {transformed_ms}")


msmd = casatools.msmetadata()
msmd.open(target_ms)
scans = msmd.scansforfield(field='J2139+1423')
nscans = len(scans)

outputvis = 'transformed.ms'

if not os.path.exists(outputvis):
    mstransform(vis=target_ms, outputvis=outputvis, datacolumn='corrected', 
        createmms=True, separationaxis = 'scan', numsubms = msmd.nscans(),
        timeaverage=True, timebin='20s', chanaverage=True, chanbin=8,
        phasecenter='J2000 21h39m01.309269s +14d23m35.99221s', field='J2139+1423',
    )

def run_wsclean(command):

    """
    Runs wsclean commands 
    """

    container = '/home/kelvin/Desktop/singularity/wsclean-v3.3-no-cuda.sif'
    wsclean_sif = container
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



niter = 0
imsize=[640,640]
cell='1mas'

imagename =  target_ms.replace('.ms','')
print(f"Making image {imagename}")

wsclean_cmd = ['wsclean', '-log-time', '-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}', \
        '-scale', f'{cell}',\
        '-mgain', '0.8', '-niter', f'{niter}', f'{outputvis}']

run_wsclean(wsclean_cmd)

# wsclean_fitsfile = imagename+'-image.fits'
# get_im_stats(wsclean_fitsfile)
# plot_fits(wsclean_fitsfile)
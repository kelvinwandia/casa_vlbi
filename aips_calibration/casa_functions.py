import os, re, glob, subprocess
from casatasks import *
import casatools
import casaplotms


"""

The functions below are not useful in the CASA pipe -- they are already included in the scripts.
They have been added to aid with the flagging in the aips pipeline -- this ensure that the main
aips script for calibration remains pristine ie just parseltongue
    
"""


# @time_execution
def makems(vis,fitsfile):

    if not os.path.exists(vis):
        logging.info(f"======>>>Making {vis}")
        importuvfits(vis=vis, fitsfile=fitsfile)
        listfile = vis.replace(".ms","_listobs.list")
        
        listobs(vis = vis, listfile = listfile, overwrite=True)

def getfields():
        
        """get list of field names in the ms """

        msmd = casatools.msmetadata()

        msmd.open(vis)  
        fieldnames = msmd.fieldnames()
        msmd.done()
        fields = {}

        for index, item in enumerate(fieldnames):
            if any(char.isdigit() for char in item):
                fields[index] = item
       
        logging.info(f"{fields} found in measurement set")
        
        return fields

def report_flag(summary, axis):
    # logging.info("REPORTING FLAGGING STATS")
    try:
        for id, stats in summary[axis].items():
            logging.info('%s %s: %5.1f percent flagged' % (axis, id, 100. * stats['flagged'] / stats['total']))
    except Exception as e:
        logging.info(f"Exception {e} while reporting flags")


def execute_aoflagger_strategy():

    """
    Flags using aoflagger
    """

    try:
        # aoflagger_sif = os.path.join(aoflagger_path.rstrip('/'), aoflagger_sif)
        print(aoflagger_sif)
        print(f"Checking for aoflagger_sif at: {aoflagger_sif}")
        if os.path.exists(aoflagger_sif):
            logging.info(f"Found {aoflagger_sif}")
            singularity_bind = os.path.join(os.path.dirname(os.path.dirname(aoflagger_sif)))
            print(singularity_bind)
            logging.info(f"You are binding singularity to {singularity_bind}")
        else:
            print(f"{aoflagger_sif} not found")
    except FileNotFoundError:
        logging.critical(f"Singularity aoflagger_sif not found")

    # fields  = getfields()

    aoflagger_cmds = ['aoflagger', '-v', '-indirect-read', '-strategy', flagging_strategy, vis]

    logging.info(f"Using strategy {flagging_strategy}")
    command_to_execute = ['singularity', 'exec', '-B', singularity_bind, aoflagger_sif] + aoflagger_cmds

    try:
        logging.info("Executing: %s", ' '.join(command_to_execute))
        process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        logging.info("stdout: %s", stdout)
        logging.info("stderr: %s", stderr)

        return_code = process.returncode
        if return_code == 0:
            logging.info(f"Strategy executed successfully. Output:\n{stdout}")
        else:
            logging.critical(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")

    except Exception as e:
        logging.critical(f"An error occurred: {e}")

    flagmanager(vis=vis, mode='save', versionname="after_automatic_flagging")

    aoflagger_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after automatic flagging")
    report_flag(aoflagger_flagging_summary, 'field')



def get_msinfo():

    nchan = []
    msmd = casatools.msmetadata()
    msmd.open(vis)
    bandwidth = msmd.bandwidths()
    nspw = len(bandwidth)
    for spw in range(nspw):
        nchan.append(msmd.nchan(spw))
    msmd.close()

    return nspw,nchan

def plot_check_baddata(save_as=None):
    """
    Plots the vis over each spectral window to check the effect before and after flagging

    Parameters:
        save_as (str): Name to save the plot file as. If None, default naming will be used.
    """
    nspw, _ = get_msinfo()

    plots_dir = os.path.join(working_dir).rstrip('/') + '/' + 'plots'
    flags_dir = os.path.join(plots_dir,'vis_before_after_flagging')

    if not os.path.exists(flags_dir):
        os.makedirs(flags_dir)

    logging.info("======>>> Plot visibilities to check bad data")

    for spw in range(0,nspw):
        plotfile = f"{flags_dir}/spw_{spw}.png" if save_as is None else f"{flags_dir}/{save_as}_spw_{spw}.png"
        plotms(vis=vis, xaxis='channel', yaxis='amp', field=phase_calibrator, iteraxis='antenna', gridcols=3, 
            spw=str(spw),gridrows=3, plotfile=plotfile, width=1500, height=750, dpi=300, showgui=False,
            overwrite=True)

    logging.info("======>>> Finished plotting the visibilities")



def makeuvfits(vis,fitsfile):
 
    logging.info(f"======>>>Making {fitsfile}")
    exportuvfits(vis=vis, fitsfile=fitsfile)

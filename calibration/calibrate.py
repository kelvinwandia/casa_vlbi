
import logging, zipfile, os, sys
from utils.helper_functions import *
import os, glob, re, logging
from datetime import datetime
import casatasks, casatools
import casaplotms
import numpy as np
import subprocess
import matplotlib
# matplotlib.use('Agg')  
import time
from natsort import natsorted
import zipfile
import shutil
from casatools import msmetadata as msmdtool
from casatools import table as tbtool

from config_file import *

global cal_tables_dict
cal_tables_dict = {}

msmd = casatools.msmetadata()
tb = casatools.table()


try:
   from casampi.MPIEnvironment import MPIEnvironment   
   parallel=MPIEnvironment.is_mpi_enabled
except:
   parallel=False
   



def set_working_dir():

    """
    Creates a working dir if one does not exist
    """

    os.makedirs(working_directory)

    try:
        os.chdir(working_directory)
        log_message(f"Changed working directory to {working_directory}")
    except Exception as e:
        logging.error(f"An error occurred while changing directory: {e}")
    
    log_message(f"Setting logfile in working dir")

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

@time_execution
def attach_tsys_gc():


    """
    Downloads JIVE helper scripts and 
    checks if system temperatures have been appended to fitsfiles
    select only the first or a random fits file
    check only tsys - problematic if already added - breaks the code

    check GAIN_CURVE
    requires only appending to one of the FITS-IDI files since they're the same throughout the
    observation
    """

    """ Use the casavlbi tools inside the package
        Saves a lot of time with copying and re-copying the edited fitsidi.py file
    """
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    helper_scripts_dir = helper_scripts_dir = os.path.join(base_dir, 'casa-vlbi')
    log_message(f"Using JIVE helper scripts: {helper_scripts_dir}")

    helper_scripts = 'casa-vlbi-master.zip'
    helper_scripts_dir = 'casa-vlbi'
    repo_url = 'https://github.com/jive-vlbi/casa-vlbi/archive/refs/heads/master.zip'
    os.system(f"rm -r {helper_scripts}")
    # Download and extract the helper scripts if not already done
    if not os.path.exists(helper_scripts_dir):
        if not os.path.exists(helper_scripts):
            subprocess.run(['wget', '-c', '--no-check-certificate', repo_url, '-O', helper_scripts], check=True)

        with zipfile.ZipFile(helper_scripts, 'r') as zip_ref:
            zip_ref.extractall()
            log_message("Zipped file extracted")

        # Move the extracted folder to the desired directory
        extracted_dir = helper_scripts.strip('.zip')
        if os.path.exists(extracted_dir):
            shutil.move(extracted_dir, helper_scripts_dir)
            log_message("renamed to casa-vlbi")
        else:
            logging.error("Extraction failed or the directory was not found.")
    else:
        log_message("JIVE helper scripts already downloaded: %s", helper_scripts_dir)


    # Append the helper scripts directory to sys.path
    sys.path.append(helper_scripts_dir)
    log_message(f"Appended {helper_scripts_dir} to sys.path")


    # Import necessary functions from the casa-vlbi package
    try:
        from casavlbitools.fitsidi import append_tsys, append_gc, convert_flags
        from casavlbitools.casa import convert_gaincurve
        log_message("imported casa-vlbi tools.")
    except ImportError as e:
        logging.error(f"importing casa-vlbi tools: {e}")
        sys.exit(1)

    # Check for UVFLG file and log_message
    if os.path.exists(uvflg_file):
        log_message(f"UVFLG File: {uvflg_file}")
    else:
        logging.warning("No uvflg file found")


    valid_extensions = ('IDI', 'idifits')
    fitsfiles = [os.path.join(fitsfiles_dir, f) 
                for f in os.listdir(fitsfiles_dir) 
                    if any(f.split('.')[-1].lower().startswith(ext.lower()) for ext in valid_extensions)]

    fitsfiles = natsorted(fitsfiles)

    log_message(f"Found FITS files: {fitsfiles}")

    if not fitsfiles:
        logging.error(f"No FITS files with extensions {valid_extensions} found in {fitsfiles_dir}")
        raise FileNotFoundError("No FITS files to load!")


    # Convert flags
    try:
        convert_flags(infile=uvflg_file, idifiles=fitsfiles, outfp=sys.stdout, outfile='{}_apriori.flag'.format(experiment))
        log_message("Flag conversion completed.")
    except Exception as e:
        logging.warning(f"Error during flag conversion: {e}")


    """ To remove all GC and TSYS"""
    extension_names = ['GAIN_CURVE', 'SYSTEM_TEMPERATURE']
    for filename in fitsfiles:
        log_message(f"Processing file: {filename}")
        with fits.open(filename, mode='update') as hdul:
            # Find extensions to remove
            extensions_to_remove = [i for i, ext in enumerate(hdul) if ext.header.get('EXTNAME') in extension_names]
            
            if extensions_to_remove:
                log_message(f"Extensions {', '.join(extension_names)} exist in the FITS file. Removing them.")
                # Remove the extensions
                for i in reversed(extensions_to_remove):
                    del hdul[i]
                hdul.flush()
                log_message("Extensions removed.")
            else:
                log_message(f"Extensions {', '.join(extension_names)} do not exist in the FITS file.")

        extension_name = 'SYSTEM_TEMPERATURE'

        if any(extension_name == ext.header.get('EXTNAME') for ext in hdul):
            log_message(f"'{extension_name}' exists in the FITS file.")
            hdul.close()
        else:
            log_message(f"Extension '{extension_name}' does not exist in the FITS file.")
        
            hdul.close()

    log_message("Attaching TSYS table")
    for i in fitsfiles:
        append_tsys(antab_file,idifiles=i)
    log_message("Finished attaching TSYS table")

    log_message("Attaching GAIN_CURVE table.")
    append_gc(antab_file, fitsfiles[0])  # Append the GAIN_CURVE table
    log_message("Finished attaching GAIN_CURVE table.")

    # Convert gain curves and flag files
    log_message("Converting gain curves")
    gc_table = f'{experiment}.gc'
    os.system(f"rm -r {gc_table}")
    if not os.path.exists(gc_table):
        log_message(f"Writing gaincurve input table {gc_table}")
        convert_gaincurve(antab_file, gc_table, min_elevation=0.0, max_elevation=90.0)

    # Find missing system temperature extension
    extension_name = 'SYSTEM_TEMPERATURE'
    missing_extensions = []
    for filename in fitsfiles:
        hdul = fits.open(filename)
        # checks the extension through the ext.header.get ... loop ext through entire hdul
        if any(extension_name == ext.header.get('EXTNAME') for ext in hdul):
            log_message(f"======>>>{extension_name}' exists in the FITS file.")
        else:
            log_message(f"======>>>Extension '{extension_name}' does not exist in the {filename} file.")
            missing_extensions.append(filename)

        # Close the FITS file
        hdul.close()
        
    # # log_message the filenames that do not contain the 'SYSTEM_TEMPERATURE' extension
    # log_message("Filenames with missing 'SYSTEM_TEMPERATURE' extension:", missing_extensions)


def check_pols(vis):

    tb = casatools.table()
    tb.open(f"{vis}/FEED", nomodify=False)
    feeds = tb.getcol("POLARIZATION_TYPE")

    # log_message(f"Original POLARIZATION_TYPE: {feeds}")
    feeds = np.where(feeds=="l","L",feeds)
    feeds = np.where(feeds=="?","R",feeds)
    tb.putcol("POLARIZATION_TYPE",feeds)
    tb.close()
    
@time_execution
def makems(vis,splitvis=None):

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)


    valid_extensions = ('IDI', 'idifits')

    fitsfiles = [os.path.join(fitsfiles_dir, f) 
    for f in os.listdir(fitsfiles_dir) 
        if any(f.split('.')[-1].lower().startswith(ext.lower()) for ext in valid_extensions)]


    fitsfiles = natsorted(fitsfiles)

    log_message(f"Found FITS files: {fitsfiles}")

    if not fitsfiles:
        logging.error(f"No FITS files with extensions {valid_extensions} found in {fitsfiles_dir}")
        raise FileNotFoundError("No FITS files to load!")

    if use_casa == True:
        log_message("use CASA has been requested")
        log_message("Assuming TSYS and GC already attached to fitsfiles")
        if not os.path.exists(vis):
            log_message(f"Making {vis}")
            importfitsidi(
                vis= vis, fitsidifile=fitsfiles,scanreindexgap_s=15.0,constobsid=True)
            listfile = vis.replace(".ms","_listobs.list")
            listobs(vis = vis, listfile = listfile, overwrite=True)

            check_pols(vis)
        else:
            check_pols(vis)

    else:
        log_message("Using UVFITS from AIPS")
        log_message(f"======>>>Using {uvfits_file}")
        if not os.path.exists(vis):
            log_message(f"======>>>Making {vis}")
            importuvfits(vis=vis, fitsfile=uvfits_file)
            listfile = vis.replace(".ms","_listobs.list")
            listobs(vis = vis, listfile = listfile, overwrite=True)

    if splitvis and not os.path.exists(splitvis):
        if not os.path.exists(splitvis):
            log_message(f"Averaging to {timebin} and {width} channels")
            split(
                vis = vis, outputvis = splitvis, timebin=timebin, width=width,
                datacolumn='data'
            )

            listfile = splitvis.replace(".ms","_split_listobs.list")
            listobs(vis = splitvis, listfile = listfile, overwrite=True)


def export_to_uvfits(vis):
    """
    Exports the flagged measurement set to UVFITS for calibration in AIPS
    """
    uvfitsfile = vis.replace('.ms','.FITS')
    log_message(f"Exporting {vis} to AIPS compatible {uvfitsfile}")
    if not os.path.exists(uvfitsfile):
        exportuvfits(
            vis = vis, fitsfile=uvfitsfile, writesyscal=False, overwrite=False
        )
    else:
        log_message(f"{uvfitsfile} exists.")
                
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
       
        log_message(f"{fields} found in measurement set")
        
        return fields


def report_flag(summary, axis):
    # log_message("REPORTING FLAGGING STATS")
    try:
        for id, stats in summary[axis].items():
            log_message('%s %s: %5.1f percent flagged' % (axis, id, 100. * stats['flagged'] / stats['total']))
    except Exception as e:
        log_message(f"Exception {e} while reporting flags",level="ERROR")
    
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


@time_execution
def plot_check_baddata(save_as=None):
    """
    Plots the vis over each spectral window to check the effect before and after flagging

    Parameters:
        save_as (str): Name to save the plot file as. If None, default naming will be used.
        avgtime (int): averaging time
    """
    nspw, _ = get_msinfo()

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    flags_dir = os.path.join(plots_dir,'vis_before_after_flagging')

    if not os.path.exists(flags_dir):
        os.makedirs(flags_dir)


    log_message("Plot visibilities to check bad data")

    sources = [phase_calibrator,target]

    for spw in range(0,nspw):
        for source in sources:
            plotfile = f"{flags_dir}/spw_{spw}.png" if save_as is None else f"{flags_dir}/{save_as}_{source}_spw_{spw}.png"
            plotms(vis=vis, xaxis='channel', yaxis='amp', field=source, iteraxis='antenna', gridcols=3, 
                spw=str(spw),gridrows=3, plotfile=plotfile, width=1500, height=750, dpi=300, showgui=False, 
                overwrite=True)

    log_message("Finished plotting the visibilities")

def get_number_of_threads():
    try:
        num_threads = os.cpu_count()
        if num_threads is None:
            log_message("Could not determine the number of threads.")
        else:
            log_message(f"Number of threads (logical processors) available: {num_threads}")
    except Exception as e:
        log_message(f"An error occurred while determining the number of threads: {e}")
    
    return num_threads/8


def flag_autocorr():
    log_message("Flagging the auto-correlations")
    flagdata(
            vis = vis, autocorr=True )
    log_message("Auto-correlations flagged successfully")
    autocorr_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after flagging autocorr")
    report_flag(autocorr_flagging_summary, 'field')

    


@time_execution
def flagging():
    
     
    # if not use_aoflagger:

    log_message(f"Quacking 2s from each scan")
    flagdata(
        vis = vis, mode='quack', quackinterval=2.0, quackmode='beg',
        quackincrement=True,
        )
    flagdata(
        vis = vis, mode='quack', quackinterval=2.0, quackmode='endb',
        quackincrement=True,
        )
    log_message("Finished quacking")

    flagmanager(vis=vis, mode='save', versionname="after_quacking")

    quacking_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after quacking")
    report_flag(quacking_flagging_summary, 'field')

    if os.path.exists(manual_file):
        log_message(f"Flagging file {manual_file} exists")
        log_message(f"Flagging using {manual_file}")
        flagdata(vis = vis, mode='list',inpfile=manual_file)
        flagmanager(vis=vis, mode='save', versionname="after_manual_flagging")
        manual_flagging_summary = flagdata(vis=vis, mode='summary')
        log_message("======>>>REPORTING FLAGGING STATS after manual flagging")
        report_flag(manual_flagging_summary, 'field')

    else:
        log_message("Manual flagging file not supplied")


def calc_flagged_data(field):

    # Get the scan data
    tb = casatools.table()
    tb.open(vis + '/ANTENNA')
    antenna_names = tb.getcol('NAME')
    # antenna_names = [int(name) for name in antenna_names.tolist()]
    tb.close()

    try:
        for antenna in antenna_names:
            log_message(f"======>>>Calculating the flagging statistics for scans in antenna {antenna}")
            flagged_vis = flagdata(vis=vis,mode='summary',field=field,antenna=antenna)
            for key in sorted(flagged_vis['scan']):
                value = flagged_vis['scan'][key]
                flagged_scan = value['flagged']
                total_scan = value['total']
                ratio = flagged_scan / total_scan
                log_message(f"{ratio * 100:.2f}% of antenna {antenna} in scan {key} are flagged")
    except Exception as e:
        logging.warning(f"======>>>Exception exception {e}: Antenna {antenna} may not have data due to flagging")
        pass

@time_execution
def execute_aoflagger_strategy():

    """
    Flags using aoflagger
    """
    
    main_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(main_dir)
    # bright_source_strategy = os.path.join(project_root, 'data', 'flagging', '4C.rfis')
    # faint_source_strategy = os.path.join(project_root,'data','flagging','faint_sources.rfis')
    
    """ For AOflagger 3 """
    lua_strategy = os.path.join(project_root,'data','flagging','generic_strategy.lua')

    bright_source_strategy = faint_source_strategy = lua_strategy
    """ End """

    strategies = {
        "bright_source_strategy": bright_source_strategy,
        "faint_source_strategy": faint_source_strategy,
    }

    missing_files = [name for name, path in strategies.items() if not os.path.exists(path)]

    if missing_files:
        missing_list = ", ".join(missing_files)
        raise FileNotFoundError(f"Missing strategy file(s): {missing_list}")
    else:
        for name, path in strategies.items():
            log_message(f"Loaded strategy file: {path}")
    
    
    try:

        container = aoflagger_path
        log_message(f"Checking for container at: {container}")
        if os.path.exists(container):
            log_message(f"Found {container}")
            singularity_bind = os.path.join(os.path.dirname(os.path.dirname(aoflagger_path)))
            log_message(f"You are binding singularity to {singularity_bind}")
        else:
            log_message(f"{container} not found")
    except FileNotFoundError:
        log_message(f"Singularity container not found",level="ERROR")

    fields  = getfields()
    num_threads = get_number_of_threads()
    
    
    phase_calibrator_keys = [key for key, value in fields.items() if value in phase_calibrator]
    fringe_finder_keys = [key for key, value in fields.items() if value in fringe_finder]
    target_keys = [key for key, value in fields.items() if value in target]
    
    # log_message(phase_calibrator_keys,fringe_finder_keys,target_keys)

    bright_strategy_phasecal = ['aoflagger', '-v', '-indirect-read', '-fields', ','.join(map(str, phase_calibrator_keys)), '-strategy', bright_source_strategy, vis]
    faint_strategy = ['aoflagger', '-v', '-indirect-read', '-fields',','.join(map(str, target_keys)), '-strategy', faint_source_strategy, vis]
    bright_strategy_fringefinder = ['aoflagger', '-v', '-indirect-read', '-fields', ','.join(map(str, fringe_finder_keys)), '-strategy', bright_source_strategy, vis]
    

    # aoflagger_cmds = ['aoflagger', '-j', f'{num_threads}', '-indirect-read', '-strategy', flagging_strategy, vis]

    insert_position = 1  # Insert after 'aoflagger'
    # Insert the verbosity flag at the specified position if verbosity is enabled
    if verbosity==True:
        aoflagger_cmds.insert(insert_position, '-v')


    for field in fields.values():
        
        # Determine the appropriate strategy based on the type of field
        if field in phase_calibrator:
            flagging_strategy = bright_strategy_phasecal
        elif field in fringe_finder:
            flagging_strategy = bright_strategy_fringefinder
        elif field in target:
            flagging_strategy = faint_strategy
        else:
            log_message(f"No strategy defined for field {field}",level="ERROR")
            continue
                
        strategy_index = flagging_strategy.index('-strategy') + 1
        strategy_file = flagging_strategy[strategy_index]
        log_message(f"Flagging {field} using strategy: {strategy_file}")
        
        flagging_strategy.insert(1, '-j')
        flagging_strategy.insert(2, str(num_threads))

        command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + flagging_strategy 

        try:
            log_message(f"Executing: {' '.join(command_to_execute)}")
            process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            log_message(f"stdout:", stdout)
            log_message(f"stderr:", stderr)

            return_code = process.returncode
            if return_code == 0:
                log_message(f"Strategy executed successfully. Output:\n{stdout}")
            else:
                log_message(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}",level="ERROR")

            log_message(f"Finished flagging field {field}")

            flagmanager(vis=vis, mode='save', versionname=f"after_automatic_flagging_{field}")
            aoflagger_flagging_summary = flagdata(vis=vis, mode='summary',field=field)
            log_message("======>>>REPORTING FLAGGING STATS after automatic flagging")
            # calc_flagged_data(phase_calibrator)
            # calc_flagged_data(target)
            report_flag(aoflagger_flagging_summary, 'field')

        except Exception as e:
            logging.critical(f"An error occurred: {e}")





@time_execution
def flag_edge_channels():

    _,  nchan = get_msinfo()
    edge_channels = int(nchan[0]*(edge_channel_fraction))
    log_message(f"You are flagging edge channels {edge_channels}")
    start = str(edge_channels-1)
    end = str(nchan[0] - edge_channels)
    flagdata(vis=vis,mode='manual',spw=f"*:0~{start};{end}~{nchan[0]-1}",flagbackup=False)

    edge_channel_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after flagging the edge channels")
    report_flag(edge_channel_flagging_summary, 'field')
    
def gencal_tsys_gc():
    
    """
    This function generates the system temperatures and gaincurve calibration tables

    """
    

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    tsys_caltable = vis.replace('.ms','.tsys'); gcal_caltable = vis.replace('.ms','.gcal')
  
  
    if not os.path.exists(tsys_caltable):
        gencal(vis=vis, caltable=tsys_caltable, caltype='tsys', uniform = False)

    if not os.path.exists(gcal_caltable):
        gencal(vis =vis, caltable=gcal_caltable, caltype='gc', infile= f'{experiment}.gc')
    
    # Plot the caltable
    for m in ['time','frequency']:
        plotfile = os.path.join(plots_dir, f"{vis.replace('.ms', f'_tsys_{m}.png')}")
        if not os.path.exists(plotfile):
            plotms(
                vis=tsys_caltable, yaxis='tsys', xaxis=m, gridcols=3, gridrows=3, coloraxis='corr',
                iteraxis='antenna', highres=True, showgui=False, dpi=800, width=1500, height=750, plotfile=plotfile,
                overwrite=True,  
            ) 
    cal_tables_dict[tsys_caltable] = "nearest,nearest"
    cal_tables_dict[gcal_caltable] = "nearest"
    log_message(f"Cal tables {tsys_caltable} and {gcal_caltable} added to cal_tables_dict {cal_tables_dict}")
    


   
   
def run_accor():
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    plotfile = os.path.join(plots_dir, f"{vis.replace('.ms', f'_accor.png')}")
    accor_caltable = vis.replace('.ms','.accor')
    
    if not os.path.exists(accor_caltable):
        log_message("Running task accor ")
        accor(vis=vis, caltable=accor_caltable, solint='30s')
        
    plotms(vis=accor_caltable, xaxis='time', yaxis = 'amp', iteraxis='antenna',coloraxis='spw',
           highres=True, showgui=False, dpi=800, width=1500, height=750, plotfile=plotfile,
            overwrite=True)
    
    smoothcal_caltable =  vis.replace('.ms','_smmoth.accor')
    plotfile_smoothcal = os.path.join(plots_dir, f"{vis.replace('.ms', f'_smooth_accor.png')}")
    
    if not os.path.exists(smoothcal_caltable):
        log_message("Running smoothcal")
        smoothcal(vis=vis, tablein=accor_caltable, caltable=smoothcal_caltable, smoothtype='median', smoothtime=1800.0)
    
    plotms(vis=smoothcal_caltable, xaxis='time', yaxis = 'amp', iteraxis='antenna',coloraxis='spw',
           highres=True, showgui=False, dpi=800, width=1500, height=750, plotfile=plotfile_smoothcal,
            overwrite=True)
    
    cal_tables_dict[smoothcal_caltable] = 'nearest'
    
def applycal_tsys_gc():
    


    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    # log_message(table,interp)
    # log_message(os.getcwd())
    log_message(f"======>>>Applying {table} using interpolation {interp}")  
    applycal(vis = vis, field = '',gaintable=table,interp = interp, parang = True,
    )
    tsys_gc_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message(f"======>>>REPORTING FLAGGING STATS after applying {table}")
    report_flag(tsys_gc_flagging_summary, 'field')

def tec_corrections():

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    tec_caltable = vis.replace('.ms','.tec')

    log_message("Downloading tec files")
    from private import tec_maps

    tec_image, _, _ = tec_maps.create(vis=vis)
  
    log_message("Generating tec solutions")
  
    if not os.path.exists(tec_caltable):
        gencal(vis=vis, caltable=tec_caltable, caltype='tecim', infile=tec_image)

    cal_tables_dict[tec_caltable] = ""

    log_message(f"Applying tec calibration table {tec_caltable}")
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    log_message(f"======>>>Applying {table} using interpolation {interp}")  
    applycal(vis=vis, field="",gaintable=table,interp=interp,parang=True)

def plot_sbd(plotfile,timerange,datacolumn):
    
    try:
        plotms(
            vis=vis, xaxis='frequency', yaxis='phase', antenna = f"{refant.split(',')[0].strip()}&*", 
            timerange=timerange, correlation='LL',avgtime='1200',
            showgui=False, coloraxis='spw', overwrite=True, ydatacolumn = datacolumn,
            gridcols=3, gridrows=3, iteraxis='baseline', width=1920, height=1080, 
            plotfile=plotfile.replace('.png',timerange)+'.png',
        )
    except Exception as plotms_error:
        logging.critical(f"Error occurred during plotms: {plotms_error}")
        


@time_execution
def sbd_fringefit():

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    
    ##### Search a good timerange to use for sbd 
    sbd_search = os.path.join(working_directory).rstrip('/') + '/' + 'sbd_search_plots'
    if not os.path.exists(sbd_search):
        os.makedirs(sbd_search)
        
    global timerange_sbd
    timerange_sbd = search_sbd_fringefit_soln(vis,fringe_finder,refant=refant,minsnr=5.0,interval=60.0,sbd_search=sbd_search)
    
    
    sbd_plotfile_before = f"{plots_dir}/before_sbd_fringefit.png"
    sbd_table = vis.replace('.ms', '.sbd')

    
    if not os.path.exists(sbd_table):
        log_message(f"Fringefitting and writing caltable: {sbd_table}")
        fringefit(
            vis=vis, caltable=sbd_table, solint='inf',
            zerorates=True, timerange=timerange_sbd, refant=refant,
            minsnr=snr_sbd, parang=True
        )
        
        plot_sbd(sbd_plotfile_before.replace('.png',timerange_sbd)+'.png',timerange_sbd,'data')


    else:
        log_message(f"Caltable {sbd_table} exists. Will not write a new one")

    ### create an empty dict to hold the cal tables
    if not use_casa: 
        global cal_tables_dict
        cal_tables_dict = {}
        cal_tables_dict[sbd_table] = "nearest"
        log_message(f"Cal table {sbd_table} added to cal_tables_dict {cal_tables_dict}")
    else:
        cal_tables_dict[sbd_table] = "nearest"
        log_message(f"Cal table {sbd_table} added to cal_tables_dict {cal_tables_dict}")

@time_execution
def applycal_sbd_fringe():

    """
    Applying the sbd calibration table to the data and plots the cor rected scan

    """
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    sbd_plotfile_after =f"{plots_dir}/after_sbd_fringefit.png"

    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())

    log_message(f"======>>>Applying {table} using interpolation {interp}")  
    applycal(vis = vis, field = '',gaintable=table,interp = interp, parang = True,
    )

    
    plot_sbd(sbd_plotfile_after,timerange_sbd,'corrected')
    sbd_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after applying sbd corrections")
    report_flag(sbd_flagging_summary, 'field')

@time_execution
def mbd_fringefit():
    """
    Performs a global fringe fit on all the data and plots the delays, phases and rates

    """

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    mbd_table = vis.replace('.ms', '.mbd')
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    log_message(f"======>>>Running global fring")
    # log_message(f"Using solution interval {solint}")
    if not os.path.exists(mbd_table):
        log_message(f"Fringefitting and making {mbd_table}")
        fringefit(
            vis=vis, caltable=mbd_table, solint=solint,
            zerorates=False, field=phase_calibrator, refant=refant, minsnr=snr_mbd, combine='spw',
            corrdepflags=True,
            gaintable=table,
            interp=interp, parang=True,
        )

    for m in ['delay', 'phase', 'rate']:
        plotfile = f"{plots_dir}/{vis.replace('.ms', '')}_mbd_{m}.png"
        plotms(
            vis=mbd_table, yaxis=m, xaxis='time', gridcols=3, gridrows=3,
            coloraxis='corr', iteraxis='antenna', highres=True, showgui=False,  width=1920, height=1080,
            overwrite=True, plotfile=plotfile,
        )
    else:
        log_message("Multiband fringe successfully completed")
 
    cal_tables_dict[mbd_table] = "linear"


@time_execution
def applycal_mbd_fringe():
    """
    Applying all the global fringe fit solutions to the data and plotting the data 
    """
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    mbd_plotfile = f'{plots_dir}/applied_mbd.png'
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    fields = phase_calibrator + ',' + target
    log_message(f"======>>>Applying {table} using interpolation {interp}")    
    log_message(f"======>>>Applying to {fields} ")

    ## Fix this by checking the position of the dictionary
    ## also dont hardcode the num spw
    nspw,nchan = get_msinfo()
    log_message(f"Applying to spws: {nspw}")
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0]]
        log_message(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0]]
        log_message(f"spw mapping is {spwmap}")

    applycal(
            vis = vis, field = fields, gaintable=table,
            interp = interp, spwmap = spwmap , parang = True,
        )

    plotms(
            vis=vis, xaxis='frequency', yaxis='phase', antenna = f"{refant.split(',')[0].strip()}&*", ydatacolumn='corrected',
            correlation='LL', gridcols=3, gridrows=3,showgui=False, coloraxis='spw', iteraxis='baseline', avgchannel=nchan,
            plotfile=mbd_plotfile,overwrite=True, width=1920, height=1080, avgtime='9999',
        ) 
    
    mbd_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after applying mbd corrections")

    calc_flagged_data(phase_calibrator)
    calc_flagged_data(target)
    report_flag(mbd_flagging_summary, 'field')

@time_execution
def bpass():
    
    """
    Calculate the bandpass corrections and plots the solutions
    """


    bpass_table = vis.replace(".ms",".bpass")

    log_message("Calculating bandpass solutions")
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())

    nspw,_ = get_msinfo()
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0]]
        log_message(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0]]
        log_message(f"spw mapping is {spwmap}")

    if not os.path.exists(bpass_table):
        # log_message(f"{bpass_table} exists. Will not create a new one")
        bandpass(
            vis = vis, bandtype = 'B', solint= 'inf', minsnr=3.0, solnorm = True, 
            field = phase_calibrator,
            refant=refant, caltable = bpass_table,gaintable = table, 
            interp = interp,spwmap = spwmap, parang=True 
            )
    
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'

    for m in ['amp','phase']:
        plotfile = f"{plots_dir}/{vis.replace('.ms', '')}_bpass_{m}.png"
        plotms(
                vis=bpass_table, yaxis=m, xaxis='frequency', gridcols=3, gridrows=3, 
                coloraxis='spw',iteraxis='antenna', highres=True, showgui=False, width=1920, height=1080,
                overwrite=True, plotfile=plotfile,
            )  
    
    cal_tables_dict[bpass_table] = "nearest,nearest"

@time_execution
def applycal_bpass():

    """
    Applying the derived bandpass corrections
    """
    log_message("Applying bandpass solutions")

    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    log_message(f"======>>>Applying {table} using interpolation {interp}")   

    nspw,_ = get_msinfo()
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0],[]]
        log_message(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0],[]]
        log_message(f"spw mapping is {spwmap}")

    applycal(vis = vis, field = '', gaintable = table,interp = interp,
            spwmap = spwmap,parang = True,
        )
    
    bpass_flagging_summary = flagdata(vis=vis, mode='summary')
    log_message("======>>>REPORTING FLAGGING STATS after applying bpass corrections")
    report_flag(bpass_flagging_summary, 'field')


@time_execution
def after_cal_plots():

    """
    Make plots to check the calibration

    """

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
   
    sources = [phase_calibrator, target]
    yaxis = ['amp', 'phase']

    for source in sources:
        for y_value in yaxis:                
            plotfile = f"{plots_dir}/{vis.replace('.ms', '')}_{source}_{y_value}.png"
            plotms(vis=vis, xaxis='frequency', yaxis=y_value,antenna = f"{refant.split(',')[0].strip()}&*", ydatacolumn='corrected',
                correlation='LL', showgui=False, coloraxis='spw', avgtime='9999', field=source,
                gridcols=3, gridrows=3, iteraxis='baseline', plotfile=plotfile, overwrite=True, width=1920, height=1080)


def getimaging_params():

    try:
        ms = casatools.ms()
        tb = casatools.table()
        ms.open(vis)
        max_uv = ms.getdata('uvdist')['uvdist'].max()
        ms.close()

        tb.open(vis + '/SPECTRAL_WINDOW')
        chan_freq = tb.getcol('CHAN_FREQ')
        highest_freq = chan_freq.max()
        tb.close()

        # 3.6e6 converts the degrees to mas
        # 5 is the sampling
        cell_size = ((c.value / highest_freq) / max_uv) * (180. / np.pi) * (3.6e6 / 5)
        cell_size = np.round(cell_size)
        log_message("You are using a cell size of:", cell_size)

    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")

    return cell_size



@time_execution
def make_map(vis=vis,source=phase_calibrator):

                
    if use_tclean == True:
        pass
    
    if use_wsclean == True:

        log_message("Using WSCLEAN for imaging")

        """
        Find all Measurement Sets in 'split_ms' and make images for each.
        """
        split_dir = os.path.join(working_directory, 'split_ms')

        if not os.path.exists(split_dir):
            log_message("'split_ms' directory not found. Run split_calibrated_ms() first.")
            return

        ms_files = sorted([f for f in os.listdir(split_dir) if f.endswith('.ms')])
        if not ms_files:
            log_message("No .ms files found in 'split_ms' — nothing to image.",level="ERROR")
            return

        log_message(f"Found {len(ms_files)} measurement sets to image: {ms_files}")

        # Create the maps directory if needed
        maps_dir = os.path.join(working_directory, "images")
        os.makedirs(maps_dir, exist_ok=True)

        # Loop through each split MS file
        for msfile in ms_files:
            vis_path = os.path.join(split_dir, msfile)
            imagename = os.path.join(maps_dir, msfile.replace('.ms', ''))
            log_message(f"Starting imaging for {vis_path}")

            try:
                cell = get_imaging_cellsize(vis_path)
                imsize = 640
                if not os.path.exists(imagename + '-image.fits'):
                    log_message(f"Making {imagename}")

                    wsclean_cmd = [
                        'wsclean', '-log-time',
                        '-size', str(imsize), str(imsize),
                        '-name', imagename,
                        '-scale', str(cell),
                        '-mgain', '0.8',
                        '-niter', '1',
                    ]

                    if verbosity:
                        wsclean_cmd.insert(2, '-quiet')

                    num_threads = get_number_of_threads()
                    threads_to_use = str(int(num_threads / 4))
                    wsclean_cmd += ["-j", threads_to_use]

                    wsclean_cmd = wsclean_cmd + [vis_path]
                    run_wsclean(wsclean_sif, wsclean_cmd)

                    wsclean_fitsfile = imagename + '-image.fits'
                    get_im_stats(wsclean_fitsfile)
                    plot_fits(wsclean_fitsfile)

            except Exception as e:
                log_message(f"Error imaging {msfile}: {e}",level="ERROR")
                # continue

        log_message("All split Measurement Sets imaged successfully.")

    

        
def split_calibrated_ms():


    msmd = msmdtool()
    tb = tbtool()

    tb.open(vis)
    datacolumn = "corrected" if "CORRECTED_DATA" in tb.colnames() else "data"
    tb.close()

    targets = [t.strip() for t in target.split(',') if t.strip()]
    sources_to_split = targets + [phase_calibrator.strip(), fringe_finder.strip()]

    log_message(f"Sources to split: {sources_to_split}")

    # --- Create output directory ---
    split_dir = "split_ms"
    if not os.path.exists(split_dir):
        os.makedirs(split_dir)
        log_message(f" Created directory: {split_dir}")

    original_dir = os.getcwd()

    msmd.open(vis)
    available_sources = msmd.namesforfields()
    msmd.close()

    os.chdir(split_dir)
    log_message(f" Changed to directory: {os.getcwd()}")

    # --- Split each source ---
    for src in sources_to_split:
        if src not in available_sources:
            log_message(f"Source '{src}' not found in {vis}, skipping.")
            continue

        outputvis = f"{src}.ms"
        log_message(f"Splitting source '{src}' ---> {outputvis}")

        if not os.path.exists(outputvis):
            split(
                vis=os.path.abspath(os.path.join("..", vis)),  # ensure full path
                outputvis=outputvis,
                field=src,
                datacolumn=datacolumn,
            )

            log_message(f"Saved source {src} as {outputvis}")
        else:
            log_message(f"Split source {outputvis} exists.")

    # --- Return to original directory ---
    os.chdir(original_dir)
    log_message(f"Returned to directory: {original_dir}")

    log_message("All specified sources processed successfully.")
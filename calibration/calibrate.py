
import logging
from utils.helper_functions import *

def set_working_dir():

    """
    Creates a working dir if one does not exist
    """

    os.makedirs(working_directory)

    try:
        os.chdir(working_directory)
        logging.info(f"Changed working directory to {working_directory}")
    except Exception as e:
        logging.error(f"An error occurred while changing directory: {e}")
    
    logging.info(f"Setting logfile in working dir")

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)


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

    helper_scripts = 'casa-vlbi-master.zip'
    helper_scripts_dir = 'casa-vlbi'
    repo_url = 'https://github.com/jive-vlbi/casa-vlbi/archive/refs/heads/master.zip'

    # Download and extract the helper scripts if not already done
    if not os.path.exists(helper_scripts_dir):
        if not os.path.exists(helper_scripts):
            subprocess.run(['wget', '-c', repo_url, '-O', helper_scripts], check=True)

        with zipfile.ZipFile(helper_scripts, 'r') as zip_ref:
            zip_ref.extractall()
            logging.info("Zipped file extracted")

        # Move the extracted folder to the desired directory
        extracted_dir = helper_scripts.strip('.zip')
        if os.path.exists(extracted_dir):
            shutil.move(extracted_dir, helper_scripts_dir)
            logging.info("======>>> renamed to casa-vlbi")
        else:
            logging.info("======>>> Extraction failed or the directory was not found.")
    else:
        logging.info("======>>> JIVE helper scripts already downloaded: ", helper_scripts_dir)

    # Append the helper scripts directory to sys.path
    sys.path.append(helper_scripts_dir)

    # Import necessary functions from the casa-vlbi package
    try:
        from casavlbitools.fitsidi import append_tsys, append_gc, convert_flags
        from casavlbitools.casa import convert_gaincurve
        logging.info("======>>> imported casa-vlbi tools.")
    except ImportError as e:
        logging.error(f"======>>> importing casa-vlbi tools: {e}")
        sys.exit(1)

    # Check for UVFLG file and print
    if os.path.exists(uvflg_file):
        logging.info(f"======>>> UVFLG File: {uvflg_file}")
    else:
        logging.info("======>>> No uvflg file found")


    fitsfiles = glob.glob(os.path.join(idifitsfiles, f'{experiment}_1_1.IDI*'))
    fitsidifiles = natsorted(fitsfiles)

    # Convert flags
    try:
        convert_flags(infile=uvflg_file, idifiles=fitsidifiles, outfp=sys.stdout, outfile='{}_apriori.flag'.format(experiment))
        logging.info("======>>> Flag conversion completed.")
    except Exception as e:
        logging.info(f"======>>> Error during flag conversion: {e}")


    filename = fitsidifiles[0]
    hdul = fits.open(filename)

    extension_name = 'SYSTEM_TEMPERATURE'

    if any(extension_name == ext.header.get('EXTNAME') for ext in hdul):
        print(f"'{extension_name}' exists in the FITS file.")
        hdul.close()
    else:
        print(f"Extension '{extension_name}' does not exist in the FITS file.")
    
        hdul.close()

        print("Attaching TSYS table")
        for i in fitsidifiles:
            append_tsys(antab_file,idifiles=i)
        
        print("Finished attaching TSYS table")

    with fits.open(filename, mode='update') as hdul:
            gain_curve_exists = False
            try:
                hdu = hdul['GAIN_CURVE']
                gain_curve_exists = True
                print("GAIN_CURVE table found.")
            except KeyError:
                print("GAIN_CURVE table does not exist")

            if not gain_curve_exists:
                print("Attaching GAIN_CURVE table")
                append_gc(antab_file, fitsidifiles[0])  # Gain curve requires only one of the fits-idi files
                print("Finished attaching GAIN_CURVE table")

    # Convert gain curves and flag files
    logging.info("Converting gain curves")
    gc_table = f'{experiment}.gc'
    if not os.path.exists(gc_table):
        convert_gaincurve(antab_file, gc_table, min_elevation=0.0, max_elevation=90.0)



    """
    loops over all the FITS-IDI files and checks if the system temperatures
    have successfully been appended
    """
    # extension_name = 'SYSTEM_TEMPERATURE'
    # missing_extensions = []
    # for filename in fitsidifiles:
    #     hdul = fits.open(filename)

    #     # checks the extension through the ext.header.get ... loop ext through entire hdul
    #     if any(extension_name == ext.header.get('EXTNAME') for ext in hdul):
    #         logging.info(f"======>>>{extension_name}' exists in the FITS file.")
    #     else:
    #         logging.info(f"======>>>Extension '{extension_name}' does not exist in the {filename} file.")
    #         missing_extensions.append(filename)

    #     # Close the FITS file
    #     hdul.close()
        
    # # Print the filenames that do not contain the 'SYSTEM_TEMPERATURE' extension
    # logging.info("======>>> Filenames with missing 'SYSTEM_TEMPERATURE' extension:", missing_extensions)

    
@time_execution
def makems(vis,splitvis=None):

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(working_directory).rstrip('/') + '/' + 'calibration_dir'

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    fitsfiles = glob.glob(os.path.join(idifitsfiles, f'{experiment}_1_1.IDI*'))
    fitsidifiles = natsorted(fitsfiles)

    if use_casa == True:
        logging.info("======>>> use CASA has been requested")
        logging.info("======>>> Assuming TSYS and GC already attached to fitsidifiles")
        if not os.path.exists(vis):
            print(f"Making {vis}")
            casatasks.importfitsidi(
                vis= vis, fitsidifile=fitsidifiles,scanreindexgap_s=15.0,constobsid=True)
            listfile = vis.replace(".ms","_listobs.list")
            casatasks.listobs(vis = vis, listfile = listfile, overwrite=True)

    else:
        logging.info("======>>> Using UVFITS from AIPS")
        logging.info(f"======>>>Using {uvfits_file}")
        if not os.path.exists(vis):
            logging.info(f"======>>>Making {vis}")
            casatasks.importuvfits(vis=vis, fitsfile=uvfits_file)
            listfile = vis.replace(".ms","_listobs.list")
            casatasks.listobs(vis = vis, listfile = listfile, overwrite=True)

    if splitvis and not os.path.exists(splitvis):
        if not os.path.exists(splitvis):
            logging.info(f"Averaging to {timebin} and {width} channels")
            casatasks.split(
                vis = vis, outputvis = splitvis, timebin=timebin, width=width,
                datacolumn='data'
            )

            listfile = splitvis.replace(".ms","_split_listobs.list")
            casatasks.listobs(vis = splitvis, listfile = listfile, overwrite=True)


            
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


    logging.info("======>>> Plot visibilities to check bad data")

    sources = [phase_calibrator,target]

    for spw in range(0,nspw):
        for source in sources:
            plotfile = f"{flags_dir}/spw_{spw}.png" if save_as is None else f"{flags_dir}/{save_as}_{source}_spw_{spw}.png"
            plotms(vis=vis, xaxis='channel', yaxis='amp', field=source, iteraxis='antenna', gridcols=3, 
                spw=str(spw),gridrows=3, plotfile=plotfile, width=1500, height=750, dpi=300, showgui=False, 
                overwrite=True)

    logging.info("======>>> Finished plotting the visibilities")

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



@time_execution
def flagging():

     
    if not use_aoflagger:
        logging.info("Flagging the auto-correlations")
        casatasks.flagdata(
                vis = vis, autocorr=True )
        logging.info("Auto-correlations flagged successfully")

        autocorr_flagging_summary = flagdata(vis=vis, mode='summary')
        logging.info("======>>>REPORTING FLAGGING STATS after flagging autocorr")
        report_flag(autocorr_flagging_summary, 'field')


    # logging.info(f"Quacking every {integration_time}s from each scan")
    # casatasks.flagdata(
    #     vis = vis, mode='quack', quackinterval=integration_time, quackmode='beg',
    #     quackincrement=True,
    #     )
    # casatasks.flagdata(
    #     vis = vis, mode='quack', quackinterval=integration_time, quackmode='endb',
    #     quackincrement=True,
    #     )
    # logging.info("Finished quacking")

    # flagmanager(vis=vis, mode='save', versionname="after_quacking")

    # quacking_flagging_summary = flagdata(vis=vis, mode='summary')
    # logging.info("======>>>REPORTING FLAGGING STATS after quacking")
    # report_flag(quacking_flagging_summary, 'field')

    if os.path.exists(manual_file):
        logging.info(f"Flagging file {manual_file} exists")
        logging.info(f"Flagging using {manual_file}")
        casatasks.flagdata(vis = vis, mode='list',inpfile=manual_file)
        flagmanager(vis=vis, mode='save', versionname="after_manual_flagging")
        manual_flagging_summary = flagdata(vis=vis, mode='summary')
        logging.info("======>>>REPORTING FLAGGING STATS after manual flagging")
        report_flag(manual_flagging_summary, 'field')

    else:
        logging.info("Manual flagging file not supplied")

def antenna_flag(antenna):
    """
    Use this if you wish to flag some antennas
    """
    logging.info(f"You are flagging antennas {antenna}")
    flagdata(vis=vis,mode='manual',antenna=antenna)

    antenna_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info(f"======>>>REPORTING FLAGGING STATS after flagging {antenna}")
    report_flag(antenna_flagging_summary, 'field')

   
@time_execution
def execute_aoflagger_strategy():

    """
    Flags using aoflagger
    """

    try:

        container = aoflagger_path
        print(f"Checking for container at: {container}")
        if os.path.exists(container):
            logging.info(f"Found {container}")
            singularity_bind = os.path.join(os.path.dirname(os.path.dirname(aoflagger_path)))
            logging.info(f"You are binding singularity to {singularity_bind}")
        else:
            print(f"{container} not found")
    except FileNotFoundError:
        logging.critical(f"Singularity container not found")

    fields  = getfields()
    num_threads = get_number_of_threads()


    # phase_calibrator_keys = [key for key, value in fields.items() if value in phase_calibrator]
    # fringe_finder_keys = [key for key, value in fields.items() if value in fringe_finder]
    # target_keys = [key for key, value in fields.items() if value in target]
    # bright_strategy_phasecal = ['aoflagger', '-v', '-indirect-read', '-fields', ','.join(map(str, phase_calibrator_keys)), '-strategy', bright_source_strategy, vis]
    # faint_strategy = ['aoflagger', '-v', '-indirect-read', '-fields',','.join(map(str, target_keys)), '-strategy', faint_source_strategy, vis]
    # bright_strategy_fringefinder = ['aoflagger', '-v', '-indirect-read', '-fields', ','.join(map(str, fringe_finder_keys)), '-strategy', bright_source_strategy, vis]
    

    aoflagger_cmds = ['aoflagger', '-j', f'{num_threads}', '-indirect-read', '-strategy', flagging_strategy, vis]

    insert_position = 1  # Insert after 'aoflagger'
    # Insert the verbosity flag at the specified position if verbosity is enabled
    if verbosity==True:
        aoflagger_cmds.insert(insert_position, '-v')


    # for field in fields.values():
        
    #     # Determine the appropriate strategy based on the type of field
    #     if field in phase_calibrator:
    #         strategy = bright_strategy_phasecal
    #     elif field in fringe_finder:
    #         strategy = bright_strategy_fringefinder
    #     elif field in target:
    #         strategy = faint_strategy
    #     else:
    #         logging.critical(f"No strategy defined for field {field}")
            

        # logging.info(f"Flagging {field}")
    logging.info(f"Using strategy {flagging_strategy}")
    command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + aoflagger_cmds

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

        # logging.info(f"Finished flagging field {field}")

    except Exception as e:
        logging.critical(f"An error occurred: {e}")

    flagmanager(vis=vis, mode='save', versionname="after_automatic_flagging")

    aoflagger_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after automatic flagging")
    report_flag(aoflagger_flagging_summary, 'field')


def calc_flagged_data(field):

    # Get the scan data
    tb = casatools.table()
    tb.open(vis + '/ANTENNA')
    antenna_names = tb.getcol('NAME')
    # antenna_names = [int(name) for name in antenna_names.tolist()]
    tb.close()

    try:
        for antenna in antenna_names:
            print(f"======>>>Calculating the flagging statistics for scans in antenna {antenna}")
            flagged_vis = flagdata(vis=vis,mode='summary',field=field,antenna=antenna)
            for key in sorted(flagged_vis['scan']):
                value = flagged_vis['scan'][key]
                flagged_scan = value['flagged']
                total_scan = value['total']
                ratio = flagged_scan / total_scan
                print(f"{ratio * 100:.2f}% of antenna {antenna} in scan {key} are flagged")
    except Exception as e:
        print(f"======>>>Exception exception {e}: Antenna {antenna} may not have data due to flagging")
        pass



@time_execution
def flag_edge_channels():

    _,  nchan = get_msinfo()
    edge_channels = int(nchan[0]*(edge_channel_fraction))
    logging.info(f"You are flagging edge channels {edge_channels}")
    start = str(edge_channels-1)
    end = str(nchan[0] - edge_channels)
    flagdata(vis=vis,mode='manual',spw=f"*:0~{start};{end}~{nchan[0]-1}",flagbackup=False)

    edge_channel_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after flagging the edge channels")
    report_flag(edge_channel_flagging_summary, 'field')
    
def gencal_tsys_gc():
    
    """
    This function generates the system temperatures and gaincurve calibration tables

    """

    global cal_tables_dict
    cal_tables_dict = {}

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(working_directory).rstrip('/') + '/' + 'calibration_dir'

    tsys_caltable = vis.replace('.ms','.tsys'); gcal_caltable = vis.replace('.ms','.gcal')
    
    if not os.path.exists(tsys_caltable):
        gencal(vis=vis, caltable=tsys_caltable, caltype='tsys', uniform = False)

    if not os.path.exists(gcal_caltable):
        gencal(vis =vis, caltable=gcal_caltable, caltype='gc', infile= f'{experiment}.gc')

    # Plot the caltable
    for m in ['frequency','time']:
        plotfile = os.path.join(calibration_dir, f"{vis.replace('.ms', '')}_tsys_{m}.png")
        if not os.path.exists(plotfile):
            plotms(
                vis=f'{experiment}.tsys', yaxis='tsys', xaxis=m, gridcols=3, gridrows=3, coloraxis='corr',
                iteraxis='antenna', highres=True, showgui=False, dpi=800, width=1500, height=750, plotfile=plotfile,
                overwrite=True,  
            ) 
    cal_tables_dict[tsys_caltable] = "nearest,nearest"
    cal_tables_dict[gcal_caltable] = "nearest"
    logging.info(f"Cal tables {tsys_caltable} and {gcal_caltable} added to cal_tables_dict {cal_tables_dict}")
    
def applycal_tsys_gc():


    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())

    logging.info(f"======>>>Applying {table} using interpolation {interp}")  
    casatasks.applycal(vis = vis, field = '',gaintable=table,interp = interp, parang = True,
    )
    tsys_gc_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after applying tsys and gc")
    report_flag(tsys_gc_flagging_summary, 'field')

@time_execution
def sbd_fringefit():

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(working_directory).rstrip('/') + '/' + 'calibration_dir'


    sbd_plotfile_before = f"{calibration_dir}/before_sbd_fringefit.png"

    sbd_table = vis.replace('.ms', '_sbd.gcal')
    try:
        casaplotms.plotms(
            vis=vis, xaxis='frequency', yaxis='phase', antenna='EF&*', 
            timerange=timerange, correlation='LL',avgtime='1200',
            showgui=False, plotfile= sbd_plotfile_before, coloraxis='spw', overwrite=True,
            gridcols=3, gridrows=3, iteraxis='baseline', width=1920, height=1080
        )
    except Exception as plotms_error:
        logging.critical(f"Error occurred during plotms: {plotms_error}")

    
    sbd_table = vis.replace('.ms', '_sbd.gcal')
   
    if not os.path.exists(sbd_table):
        # logging.info(f"{sbd_table} exists. Will not create a new one")
        casatasks.fringefit(
            vis=vis, caltable=sbd_table, solint='inf',
            zerorates=True, timerange=timerange, refant=refant,
            minsnr=snr_sbd, parang=True
        )


    ### create an empty dict to hold the cal tables
    if not use_casa: 
        global cal_tables_dict
        cal_tables_dict = {}
        cal_tables_dict[sbd_table] = "nearest"
        logging.info(f"Cal table {sbd_table} added to cal_tables_dict {cal_tables_dict}")
    else:
        cal_tables_dict[sbd_table] = "nearest"
        logging.info(f"Cal table {sbd_table} added to cal_tables_dict {cal_tables_dict}")

@time_execution
def applycal_sbd_fringe():

    """
    Applying the sbd calibration table to the data and plots the cor rected scan

    """
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(plots_dir,'calibration_dir')

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)

    sbd_plotfile_after =f"{calibration_dir}/after_sbd_fringefit.png"

    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())

    logging.info(f"======>>>Applying {table} using interpolation {interp}")  
    casatasks.applycal(vis = vis, field = '',gaintable=table,interp = interp, parang = True,
    )
    
    casaplotms.plotms(vis=vis, xaxis='frequency', yaxis='phase', antenna='EF&*', ydatacolumn='corrected',
        timerange=timerange, correlation='LL',showgui=False, coloraxis='spw',avgtime='1200', width=1920, height=1080,
        gridcols=3, gridrows=3, iteraxis='baseline',plotfile=sbd_plotfile_after,overwrite=True
        )
    
    sbd_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after applying sbd corrections")
    report_flag(sbd_flagging_summary, 'field')

@time_execution
def mbd_fringefit():
    """
    Performs a global fringe fit on all the data and plots the delays, phases and rates

    """

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(plots_dir,'calibration_dir')

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)


    mbd_table = vis.replace('.ms', '_mbd.gcal')
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    logging.info(f"======>>>Running global fring")
    if not os.path.exists(mbd_table):
        logging.info(f"Fringefitting and making {mbd_table}")
        casatasks.fringefit(
            vis=vis, caltable=mbd_table, solint=solint,
            zerorates=False, field=phase_calibrator, refant=refant, minsnr=snr_mbd, combine='spw',
            corrdepflags=True,
            gaintable=table,
            interp=interp, parang=True,
        )

    for m in ['delay', 'phase', 'rate']:
        plotfile = f"{calibration_dir}/{vis.replace('.ms', '')}_mbd_{m}.png"
        casaplotms.plotms(
            vis=mbd_table, yaxis=m, xaxis='time', gridcols=3, gridrows=3,
            coloraxis='corr', iteraxis='antenna', highres=True, showgui=False,  width=1920, height=1080,
            overwrite=True, plotfile=plotfile,
        )
    else:
        logging.info("Multiband fringe successfully completed")
 
    cal_tables_dict[mbd_table] = "linear"


@time_execution
def applycal_mbd_fringe():
    """
    Applying all the global fringe fit solutions to the data and plotting the data 
    """
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(plots_dir,'calibration_dir')

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)


    mbd_plotfile = f'{calibration_dir}/applied_mbd.png'
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    fields = phase_calibrator + ',' + target
    logging.info(f"======>>>Applying {table} using interpolation {interp}")    
    logging.info(f"======>>>Applying to {fields} ")

    ## Fix this by checking the position of the dictionary
    ## also dont hardcode the num spw
    nspw,_ = get_msinfo
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")

    casatasks.applycal(
            vis = vis, field = fields, gaintable=table,
            interp = interp, spwmap = spwmap , parang = True,
        )

    # casaplotms.plotms(
    #         vis=vis, xaxis='frequency', yaxis='phase', antenna='EF&*', ydatacolumn='corrected',
    #         correlation='LL', gridcols=3, gridrows=3,showgui=False, coloraxis='spw',
    #         plotfile=mbd_plotfile,overwrite=True, width=1920, height=1080
    #     ) 
    
    mbd_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after applying mbd corrections")
    report_flag(mbd_flagging_summary, 'field')

@time_execution
def bpass():
    
    """
    Calculate the bandpass corrections and plots the solutions
    """


    bpass_table = vis.replace(".ms","_gcal.bpass")

    logging.info("Calculating bandpass solutions")
    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())

    nspw,_ = get_msinfo()
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")

    if not os.path.exists(bpass_table):
        # logging.info(f"{bpass_table} exists. Will not create a new one")
        casatasks.bandpass(
            vis = vis, bandtype = 'B', solint= 'inf', minsnr=3.0, solnorm = True, 
            # field = phase_calibrator + ',' + fringe_finder, 
            field = phase_calibrator,
            refant=refant, caltable = bpass_table,gaintable = table, 
            interp = interp,spwmap = spwmap, parang=True 
            )
    
    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(plots_dir,'calibration_dir')

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)

    for m in ['amp','phase']:
        plotfile = f"{calibration_dir}/{vis.replace('.ms', '')}_bpass_{m}.png"
        casaplotms.plotms(
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
    logging.info("Applying bandpass solutions")

    table = list(cal_tables_dict.keys())
    interp = list(cal_tables_dict.values())
    logging.info(f"======>>>Applying {table} using interpolation {interp}")   

    nspw,_ = get_msinfo()
    
    if use_casa == True:
        spwmap = [[],[],[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")
    else:
        spwmap = [[], nspw*[0]]
        logging.info(f"spw mapping is {spwmap}")

    casatasks.applycal(vis = vis, field = '', gaintable = table,interp = interp,
            spwmap = spwmap,parang = True,
        )
    
    bpass_flagging_summary = flagdata(vis=vis, mode='summary')
    logging.info("======>>>REPORTING FLAGGING STATS after applying bpass corrections")
    report_flag(bpass_flagging_summary, 'field')


@time_execution
def after_cal_plots():

    """
    Make plots to check the calibration

    """

    plots_dir = os.path.join(working_directory).rstrip('/') + '/' + 'plots'
    calibration_dir = os.path.join(plots_dir,'calibration_dir')

    if not os.path.exists(calibration_dir):
        os.makedirs(calibration_dir)


    sources = [phase_calibrator, target]
    yaxis = ['amp', 'phase']

    for source in sources:
        for y_value in yaxis:                
            plotfile = f"{calibration_dir}/{vis.replace('.ms', '')}_{source}_{y_value}.png"
            plotms(vis=vis, xaxis='frequency', yaxis=y_value, antenna='EF&*', ydatacolumn='corrected',
                correlation='LL', showgui=False, coloraxis='spw', avgtime='9999', field=source,
                gridcols=3, gridrows=3, iteraxis='baseline', plotfile=plotfile, overwrite=True, width=1920, height=1080)


# def getimaging_params():

    # try:
    #     ms = casatools.ms()
    #     tb = casatools.table()
    #     ms.open(vis)
    #     max_uv = ms.getdata('uvdist')['uvdist'].max()
    #     ms.close()

    #     tb.open(vis + '/SPECTRAL_WINDOW')
    #     chan_freq = tb.getcol('CHAN_FREQ')
    #     highest_freq = chan_freq.max()
    #     tb.close()

    #     # 3.6e6 converts the degrees to mas
    #     # 5 is the sampling
    #     cell_size = ((c.value / highest_freq) / max_uv) * (180. / np.pi) * (3.6e6 / 5)
    #     cell_size = np.round(cell_size)
    #     logging.info("You are using a cell size of:", cell_size)

    # except Exception as e:
    #     logging.critical(f"An unexpected error occurred: {e}")

    # return cell_size


"""
Yu need to remove this function and use wsclean here
"""

@time_execution
def dirty_map(source):

    
    
    msmd.open(vis)
    field_id = msmd.fieldsforname(source)[0]
    msmd.close()


    
    dirty_maps_dir = os.path.join(working_directory,'dirty_maps')
    if not os.path.exists(dirty_maps_dir):
        os.makedirs(dirty_maps_dir)

    imagename = f"{dirty_maps_dir}/{source}_dirty_map"

    if use_tclean == True:
        if not os.path.exists(imagename):
            tclean(vis= phasecal_ms, imagename=imagename,imsize=imsize, cell=cell,
                gridder='standard',weighting='briggs',robust=robust,niter=0, field = str(field_id)
                )
            fitsname = imagename+'.fits'
            exportfits(imagename=imagename+'.image',fitsimage=fitsname,overwrite=True)
            get_im_stats(fitsname)
            plot_fits(fitsname)


    if use_wsclean == True:
        if not os.path.exists(imagename+'-image.fits'):
                
            logging.info(f"Making {imagename}")
            # Insert the verbosity flag at the specified position if verbosity is enabled
            
            wsclean_cmd = ['wsclean', '-log-time','-size', f'{imsize[0]}', f'{imsize[1]}','-name',f'{imagename}','-scale', f'{cell}',\
                                '-mgain', '0.8', '-niter', '0' , '-field',f'{field_id}',f'{vis}']
            insert_position = 2
            if verbosity==True:
                wsclean_cmd.insert(insert_position, '-quiet')

            
            run_wsclean(wsclean_sif,wsclean_cmd)

        wsclean_fitsfile = imagename+'-image.fits'
        get_im_stats(wsclean_fitsfile)
        plot_fits(wsclean_fitsfile)




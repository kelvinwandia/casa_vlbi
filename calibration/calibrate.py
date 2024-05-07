


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

 



def makems(vis,splitvis=None):

    if not os.path.exists(vis):
        logging.info(f"Making {vis}")
        casatasks.importuvfits(
            vis=vis, fitsfile=uvfits_file
        )


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

        try:
            if phase_calibrator in fields.values() and target in fields.values():
                logging.info(f"{phase_calibrator} and {target} found in {vis}")

        except Exception as e:
            logging.critical(f"{fields} not present in {vis}")

        
        return fields




def flagging():


    logging.info("Flagging the auto-correlations")
    casatasks.flagdata(
            vis = vis, autocorr=True )
    logging.info("Auto-correlations flagged successfully")

    
    logging.info(f"Quacking every {integration_time}s from each scan")
    casatasks.flagdata(
        vis = vis, mode='quack', quackinterval=integration_time, quackmode='beg',
        quackincrement=True,
        )
    casatasks.flagdata(
        vis = vis, mode='quack', quackinterval=integration_time, quackmode='endb',
        quackincrement=True,
        )
    logging.info("Finished quacking")


    if os.path.exists(manual_file):
        logging.info(f"Flagging file {manual_file} exists")
        logging.info(f"Flagging using {manual_file}")
        casatasks.flagdata(
            vis = vis, mode='list',inpfile=manual_file
        )
    else:
        print("Manual flagging file not supplied")



def execute_aoflagger_strategy():

    try:
        # container = os.path.join(aoflagger_path.rstrip('/'), aoflagger_sif)
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

    phase_calibrator_keys = [key for key, value in fields.items() if value in phase_calibrator]
    fringe_finder_keys = [key for key, value in fields.items() if value in fringe_finder]
    target_keys = [key for key, value in fields.items() if value in target]
    bright_strategy = ['aoflagger', '-v', '-indirect-read', '-fields', ','.join(map(str, phase_calibrator_keys)), '-strategy', bright_source_strategy, vis]
    faint_strategy = ['aoflagger', '-v', '-indirect-read', '-fields',','.join(map(str, target_keys)), '-strategy', faint_source_strategy, vis]


    for field in fields.values():
        
        # Determine the appropriate strategy based on the type of field
        if field in phase_calibrator:
            strategy = bright_strategy
        elif field in fringe_finder:
            strategy = bright_strategy
        elif field in target:
            strategy = faint_strategy
        else:
            logging.critical(f"No strategy defined for field {field}")
            

        logging.info(f"Flagging {field}")
        logging.info(f"Using strategy {strategy}")
        command_to_execute = ['singularity', 'exec', '-B', singularity_bind, container] + strategy

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

            logging.info(f"Finished flagging field {value}")

        except Exception as e:
            logging.critical(f"An error occurred: {e}")





def sbd_fringefit():

    try:
        casaplotms.plotms(
            vis=vis, xaxis='frequency', yaxis='phase', antenna='EF&*', 
            timerange=timerange, correlation='LL', avgtime='1200',
            showgui=False, plotfile=vis.replace('.ms','.png'), coloraxis='spw', overwrite=True,
            gridcols=3, gridrows=3, iteraxis='baseline'
        )
    except Exception as plotms_error:
        logging.critical(f"Error occurred during plotms: {plotms_error}")

    try:
        sbd_table = vis.replace('.ms', '.sbd')
        os.system(f'rm -rf {sbd_table}.*')
    except Exception as rm_error:
        logging.critical(f"Error occurred during file removal: {rm_error}")

    try:
        casatasks.fringefit(
            vis=vis, caltable=sbd_table, solint='inf',
            zerorates=True, timerange=timerange, refant=refant,
            minsnr=snr_sbd, parang=True
        )
    except Exception as fringefit_error:
        logging.critical(f"Error occurred during fringefit: {fringefit_error}")

def applycal_sbd_fringe():

    """
    Applying the sbd calibration table to the data and plots the cor rected scan

    """

    sbd_table = vis.replace(".ms",".sbd")
    sbd_plotfile_corr = sbd_table.replace("time_jump.png","_corrected_time_jump.png")


    casatasks.applycal(
        vis = vis, field = '',gaintable=sbd_table,
        interp = ['nearest'], parang = True,
        )
    
    casaplotms.plotms(
        vis=vis, xaxis='frequency', yaxis='phase', antenna='EF&*', ydatacolumn='corrected',
        timerange=timerange, 
        correlation='LL',showgui=False, coloraxis='spw',avgtime='1200',
        gridcols=3, gridrows=3, iteraxis='baseline',
        plotfile=sbd_plotfile_corr+'_corrected_time_jump.png',overwrite=True
        )
    


def mbd_fringefit():
    """
    Performs a global fringe fit on all the data and plots the delays, phases and rates

    """


    sbd_table = vis.replace('.ms', '.sbd')
    mbd_table = vis.replace('.ms', '.mbd')

    try:
        if not os.path.exists(sbd_table):
            raise FileNotFoundError(f"The table file '{sbd_table}' does not exist. You need to run instr fring")

        os.system('rm -rf {}*'.format(mbd_table))

        casatasks.fringefit(
            vis=vis, caltable=mbd_table, solint=solint,
            zerorates=False, field=phase_calibrator, refant=refant, minsnr=snr_mbd, combine='spw',
            corrdepflags=True,
            gaintable=[sbd_table],
            interp=['nearest'], parang=True,
        )

        for m in ['delay', 'phase', 'rate']:
            plotfile = '{}_mbd_{}'.format(vis.replace(".ms", ""), m) + ".png"
            casaplotms.plotms(
                vis=mbd_table, yaxis=m, xaxis='time', gridcols=3, gridrows=3,
                coloraxis='corr', iteraxis='antenna', highres=True, showgui=False, dpi=800, width=1500,
                height=750, overwrite=True, plotfile=plotfile,
            )
        else:
            logging.info("Multiband fringe successfully completed")
    except FileNotFoundError as e:
        print(f"Error: {e}")

def applycal_mbd_fringe():
    """
    Applying all the global fringe fit solutions to the data and plotting the data 
    """
    # logging.info("Applying mbd corrections")


    mbd_table = vis.replace(".ms",".mbd")
    sbd_table = vis.replace(".ms",".sbd")

    casatasks.applycal(
            vis = vis, field = phase_calibrator + ',' + target, 
            gaintable=[sbd_table,mbd_table],
            interp = ['nearest','linear'], spwmap = [[], 8*[0]], parang = True,
        )
    mbd_plotfile_corr = vis.replace(".ms","_corrected_phases.png")

    casaplotms.plotms(
            vis=sbd_table, xaxis='frequency', yaxis='phase', antenna='EF&*', ydatacolumn='corrected',
            correlation='LL', gridcols=3, gridrows=3,showgui=False, coloraxis='spw',
            plotfile=mbd_plotfile_corr,overwrite=True
        ) 
    
    

def bpass():
    
    """
    Calculate the bandpass corrections and plots the solutions
    """


    mbd_table = vis.replace(".ms",".mbd")
    sbd_table = vis.replace(".ms",".sbd")
    bpass_table = vis.replace(".ms",".bpass")


    os.system(f'rm -rf {bpass_table}*')

    logging.info("Calculating bandpass solutions")

    casatasks.bandpass(
        vis = vis, bandtype = 'B', solint= 'inf', minsnr=3.0, solnorm = True, field = phase_calibrator, 
        refant=refant, caltable = bpass_table,gaintable = [sbd_table,mbd_table], 
        interp = ['nearest','linear'],spwmap = [[], 8*[0]], parang=True 
        )
    
    for m in ['amp','phase']:
        casaplotms.plotms(
                vis=bpass_table, yaxis=m, xaxis='frequency', gridcols=3, gridrows=3, 
                coloraxis='spw',iteraxis='antenna', highres=True, showgui=False, dpi=800, width=1500, 
                height=750, overwrite=True, plotfile='{}_{}'.format(vis.replace(".ms",".bpass"),m)+'.png',
            )  

def applycal_bpass():

    """
    Applying the derived bandpass corrections
    """
    logging.info("Applying bandpass solutions")


    mbd_table = vis.replace(".ms",".mbd")
    sbd_table = vis.replace(".ms",".sbd")
    bpass_table = vis.replace(".ms",".bpass")


    casatasks.applycal(
            vis = vis, field = '', gaintable = [sbd_table,mbd_table,bpass_table],
            interp = ['nearest','linear','nearest,nearest'],
            spwmap = [[], 8*[0],[]],
            parang = True,
        )


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

def mytclean(source,niter):

    if niter == 0:
        # cell_size = getimaging_params()
        cell_size='1mas'
        logging.info(f"Making dirty map for {source}")
        
    imagename = source+f"_num_{niter}_iterations"
    logging.info(f"Making {imagename}")
    
    os.system(f"rm -r {imagename}.*")
    

   
    im_size = [int(x) for x in imsize.split(',')]
    logging.info("Running tclean")
    casatasks.tclean(
        vis=vis,imsize=im_size,imagename=imagename,cell=cell_size,
        niter=niter, deconvolver='clark',interactive=False, gridder='standard',
        field=source,
        )  



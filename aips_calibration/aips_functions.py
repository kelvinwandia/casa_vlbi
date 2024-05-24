
def logfile():

    """
    Creates a logfile for all the tasks in the cwd
    """

    timestamp = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    try:
        logging(f"Creating AIPS logfile {AIPS.log.name}")
        AIPS.log = open(os.getcwd()+'/'+f'AIPS_{timestamp}.log','w')
    except Exception as e:
        logging.error(f"An error occurred while creating AIPS logfile: {e}")
        



def set_working_dir(working_dir):

    """
    Creates a working dir if one does not exist
    """

    os.makedirs(working_dir)

    try:
        os.chdir(working_dir)
        logging.info(f"Changed working directory to {working_dir}")
    except Exception as e:
        logging.error(f"An error occurred while changing directory: {e}")
    
    logging.info(f"Setting logfile in working dir")

    logfile()


def TV(condition):
    
    """
    Controls tv TV based on the condition
    Args:
        condition (int): If True, the TV will be started; if False, it will be closed.
    Returns:
        str: A message indicating whether the TV was started or closed.
    """
    tv = AIPSTV()

    if condition==True:
        tv.start()
        logging.info("AIPTV started")
        return "TV started."
    else:
        # Close the TV
        tv.kill()
        logging.info("AIPTV killed")
        return "TV closed."
    
def get_table(data,table):

    """
    Gets the table in uvdata

    Args:
        table(str): the table name to get

    Returns:
        table version
    """

    ver = 0
    for i in range(len(data.tables)):
        try:
            if table in data.tables[i][1]:
                ver = data.tables[i][0]
                logging.info(f"Highest table of type {table} is {ver}")
        except Exception as e:
            logging.error(f"Error {e} occured when getting {table}")

    return ver

sorted_indata = None

def sort_to_tb():
    global sorted_indata  # 
    if sorted_indata is not None:
        # If the sorted data is already available, return it directly
        return sorted_indata

    loaded_indata = AIPSUVData(experiment,inclass,inseq,indisk)
    
    if loaded_indata.header['sortord'] =='TB':
        print('Sort order is TB')
        indata = loaded_indata
    else: 
        print("Sorting data to TB")
        sorted_indata = AIPSUVData(experiment,'UVSRT',inseq,indisk)
        uvsrt = AIPSTask('UVSRT')
        uvsrt.indata = loaded_indata
        uvsrt.sort = 'TB'
        if sorted_indata.exists():
            logging.info(f"Found prior {sorted_indata}, zapping it")
            sorted_indata.zap(force=True)
        uvsrt.go()
        uvsrt.outdata = sorted_indata
        indata = sorted_indata
    
    
    # set_indata.inclass = inclass
    # set_indata.indisk = indisk
    # set_indata.inseq = inseq

    return indata
    

def set_indata():
    """
    Sets the default values for indisk and inseq
    Checks if the UVDATA is sorted and if not sorts to TB

    """
    indata = sort_to_tb()
    tasav_indata = AIPSUVData(experiment,'TASAV',inseq,indisk)
    set_indata.indata = indata
    set_indata.tasav_indata = tasav_indata


    # logging.info(f"Setting indata as {indata} and tasav as {tasav_indata}")
    # set_indata.inclass = inclass
    # set_indata.indisk = indisk
    # set_indata.inseq = inseq




def cleanup():

    """
    Sets the default values for indisk and inseq
    """
    indata = set_indata.indata
    tasav_indata = set_indata.tasav_indata

    if indata.exists():
        try:
            logging.info(f"Zapping {indata} and clearing state")
            indata.clrstat()
            indata.zap(force=True)
        except Exception as e:
            logging.error(f"Unable to zap {indata}") 
    else:
        logging.info(f"No {indata} to zap")

    if tasav_indata.exists():
        try:
            logging.info(f"Zapping {tasav_indata} and clearing state")
            tasav_indata.clrstat()
            tasav_indata.zap(force=True)
        except Exception as e:
            logging.error(f"Unable to zap {tasav_indata}") 

    else:
        print(f"No {tasav_indata} to zap")



def load_fitsfiles(file_extension):

    logging.info(f"The file extension is: {file_extension}")
    """
    This function loads either idi/fits files or the pipeline calibration in form of TASAV

    Args:
        zap_old_data (bool): deletes old uvdata/tasav data
    """


    pattern = f"{experiment.strip()}_{pointing.strip()}_1.{file_extension.strip()}*"
    logging.info(f"Searching for fitsfiles with extension {pattern}")

    try:
        fitsfiles = glob.glob(fitsfiles_dir.rstrip('/')+'/'+pattern)
        # logging.info(f"Found fitsfiles in fitsfiles_dir: {fitsfiles}")

        if not fitsfiles:
            fitsfiles = glob.glob(working_dir.rstrip('/')+'/'+pattern)
            logging.info(f"Found fitsfiles in working_dir: {fitsfiles}")

        if not fitsfiles:
            raise FileNotFoundError("No fitsfiles were found")

        if len(fitsfiles)>1:
            # fitsfiles = sorted(fitsfiles,key=lambda x: int(re.findall(r'\d+$', x)[0]))
            fitsfiles = natsorted(fitsfiles)
            if len(fitsfiles)>11:
                for i in range(0,len(fitsfiles),10):
                    chunk = fitsfiles[i:i+10]
                    logging.info("Chunk %d: %s", (i // 10) + 1, ", ".join(chunk))
        else:
            logging.info(f"Found fitsfiles: {fitsfiles[0]}")

    except Exception as e:
        logging.error(f"An error occurred while getting the fitsfiles: {e}")
        logging.error(traceback.format_exc())  



    indata = AIPSUVData(experiment,inclass,inseq,indisk)

    if zap_data==True:
        if indata.exists():
            logging.info(f"Zapping requested")
            indata.clrstat()
            indata.zap(force=True)
        else:
            logging.info("No existing UVDATA found to zap")


    if file_extension == file_extension_for_cal:
        if indata.exists():
            logging.info("Zapping existing UVDATA")
            indata.clrstat()
            indata.zap(force=True)
        else:
            logging.info("No existing UVDATA found to zap")


    try:
        logging.info("Executing task FITLD")
        fitld = AIPSTask('FITLD')
        fitld.digicor = -1
        fitld.douvcomp = -1 
        fitld.ncount = len(fitsfiles) 
        fitld.outname, fitld.outclass, fitld.outseq,fitld.outdisk = \
            indata.name,indata.klass,indata.seq,indata.disk
        fitld.doconcat = 1
        fitld.clint = 0.25 # set aips CL interval to 15 seconds
        fitld.datain = fitsfiles[0]

        if indata.exists():
            logging.info("UVDATA exists, will not write a new one")
            
        else:
            logging.info(f"Loading {fitsfiles}")
            fitld.go()
    except Exception as e:
        logging.error(f"Exception {e} occured while executing FITLD")


def load_tasav():


    tasav_file_path = os.path.join(fitsfiles_dir, tasav_file)
    if os.path.exists(tasav_file_path):
        logging.info(f"Found tasav file in fitsfiles_dir: {tasav_file_path}")
    else:
        raise FileNotFoundError("No tasav file found in either fitsfiles_dir or working_dir")
    # indata = AIPSUVData(experiment,inclass,inseq,indisk)
    tasav_indata = AIPSUVData(experiment,'TASAV',inseq,indisk)
    # tasav_indata = set_indata.tasav_indata

    # if zap_data==True:
    #     try:
    #         logging.info(f"Zapping requested")
    #         tasav_indata.clrstat()
    #         tasav_indata.zap(force=True)
    #     except Exception as e:
    #         logging.error(f"{e}")


    try:
        logging.info(f"Loading {tasav_file}")
        fitld = AIPSTask('FITLD')
        fitld.digicor = -1
        fitld.douvcomp = -1 
        fitld.ncount = 1
        fitld.outname, fitld.outclass, fitld.outseq,fitld.outdisk = \
                tasav_indata.name,tasav_indata.klass,tasav_indata.seq,tasav_indata.disk    
        fitld.doconcat = 1
        fitld.clint = 0.25
        fitld.datain = tasav_file_path
        if tasav_indata.exists():
            logging.info(f"TASAV file already exists. Will not write a new one")
        else:
            fitld.go()

    except Exception as e:
        logging.error(f"Exception {e} occured while attempting to load {tasav_file}")



def get_obs_params():

    """
    Gets the observation date and the antennas used in the observation
    
    Returns:
        obs_date (int; year,month,day): used in downloading ionex files for TECOR
        antennas (list): observing antennas, useful for pb correction
        refant (int): the reference antenna
        refant_indices (list): list of antennas that could be used as refants
            refant_indices is used in FRING


    
    FITS problem when reading the AN table  -- IDIFITS work
    UnicodeDecodeError: 'utf-8' codec can't decode byte 0xb2 in position 1: invalid start byte

    The above exception was the direct cause of the following exception:

    SystemError: <class 'UnicodeDecodeError'> returned a result with an error set

    The above exception was the direct cause of the following exception:

    """
    set_indata()
    indata = set_indata.indata
    tasav_indata = set_indata.tasav_indata

    obs_date = indata.header.date_obs
    date_obj = datetime.strptime(obs_date,"%Y-%m-%d")
    year = date_obj.year
    month = date_obj.month
    day = date_obj.day
    # logging.info(f"These data were observed on {year}_{month}_{day}")
    day_of_year = date_obj.timetuple().tm_yday

    # NB: Use AN table from TASAV - FITS imported from CASA have AN that is not understood
    # Something is off with epoch 1
    # Antennas that participated
    
    # if params['file_extension'] == "gv020a":
    #     obs_antennas = ['JB', 'WB', 'EF', 'ON', 'MC', 'TR', 'NT', 'AR', 'GB']
    # else:
    obs_antennas = []
    for antenna in tasav_indata.antennas:
        obs_antennas.append(antenna)

    logging.info(f"{obs_antennas} participated in the observation")
    logging.info(f"{searchants} will be used as substitute if {refant} not found")
    # Set reference antenna
   
    indices = [obs_antennas.index(ant) for ant in searchants]
    refant_indices =[i+1 for i in indices]
    

    try:
        find_refant = next((ant for ant in searchants if ant in obs_antennas),None)
        if find_refant is not None:
            refant_index = obs_antennas.index(find_refant)+1 # AIPS begins indexing at 1 
            indices = [obs_antennas.index(ant) for ant in searchants]
            refant_indices =[i+1 for i in indices]
        else:
            raise ValueError("Reference antenna not found")
    except ValueError as v:
        logging.error(f"Caught value error {v}\nSupply a reference antenna")


    logging.info(f"The refant indices are {refant_indices} and the refant is {refant} with index {refant_index}")


    # TODO: Get the integration time from the data
    # Get the solution interval
    # data_keys = indata.header.keys()
    # for key in data_keys:
    #     value = indata.header[key]
    #     print(f"{key}:{value}")

    return year,day_of_year,obs_antennas,refant_index,refant_indices

def zap_cal_tables(tab_name):
    """
    Zaps existing calibration tables
    Checks whether the input is a list and zap otherwise checks if its a single string
    (AIPS tablename have 2 char) and converts it to a list

    Args:
        tab_name: str (for single table) or list of table name(s) to delete
    """
    set_indata()
    indata = set_indata.indata

    if isinstance(tab_name,list):
        logging.info(f"You are zapping tables {tab_name}")
        for i in tab_name:
            table = get_table(indata,i)
            if table != 0:
                logging.info(f"Extension table {i} will be zapped")
                indata.zap_table(i,-1)

    elif isinstance(tab_name,str):
        if len(tab_name)!=2:
            logging.info("Invalid table name: enter tables as lists")
        else:
            table_name = [tab_name]
            for i in tab_name:
                table = get_table(indata,i)
                if table != 0:
                    logging.warning(f"Extension table {i} will be zapped")
                    indata.zap_table(i,-1)



def runindxr():
    """
    Checks if the index table exists, then deletes it and builds a new one
    
    This is important since we are going to zap all CL tables before doing tacop
    from TASAV file - we can build the pristine CL using INDXR
    """
    set_indata()
    indata = set_indata.indata
    nx_ver = get_table(indata,'NX')

    
    if nx_ver != 0:
        logging.info("Zapping index table")
        indata.zap_table('NX',-1)

    # build a new table - should also build a new CL1
    
    indxr = AIPSTask('INDXR')
    indxr.indata = indata
    # indxr.cparm[3] = integration_time
    logging.info("Building new index and pristine CL using INDXR ")
    indxr.go()


   
def runtacop():

    """
    Copies CL2 which contains the amplitude and parallactic angle corrections
    Also copies BP1 which contains the bandpass corrections and FG1

    """

    set_indata()
    indata = set_indata.indata
    tasav_indata = set_indata.tasav_indata


    # Zap all tables
    
    tables_to_zap = ['BP','FG','CL','SN','TE']

    zap_cal_tables(tables_to_zap)
    # for inext,invers in tables_to_copy.items():
    #     indata.zap_table(inext,-1)

    runindxr()

    tables_to_copy = {'CL':2, 'BP':1, 'FG':1}
    logging.info(f"You are copying tables: {tables_to_copy}")

    tacop = AIPSTask('TACOP')
    tacop.indata = tasav_indata
    tacop.outdata = indata
    tacop.ncount = 1

    # Copy the tables 
    for inext, invers in tables_to_copy.items():
        tacop.inext = inext
        tacop.invers = invers
        tacop.go()
    
    logging.info(f"Tables {tables_to_copy} copied")


def download_ionex_files():
    
    """
    Downloads ionex files. Able to download multiple files if observation spans days
    The TEC file format: YYYY/DDD/jplgDDD0.YYi.Z
    Uses Zcat to decompress the file
    Commented code downloads from newer website - if NASA changes TEC files location
    Returns:
        TEC file downloaded
    """
    
    year,day_of_year,_,_,_ = get_obs_params()

    ionex_files =[]
    day_of_year = [day_of_year]
    for day in day_of_year:

        # Get the correct file based on day, use zfill since you cant turn str to int and keep preceding zero
        num_day_digits = len(str(day))
        if num_day_digits == 1:
            day = str(day).zfill(3)
        elif num_day_digits == 2:
            day = str(day).zfill(3)

        # tec_file = f'jplg{day}0.{str(year)[-2:]}i'+'.Z'
        # tec_file_url = f"https://cddis.nasa.gov/archive/gnss/products/ionex/{year}/{day}/"
        # tec_file_url = tec_file_url+tec_file

        # os.system(f"rm -rf {tec_file}")
        # try:
        #     # subprocess.run(["wget","-c",tec_file_url],check=True)
        #     ionex_files.append(tec_file)
        #     print(f"{tec_file} downloaded successfully.")
        # except subprocess.CalledProcessError:
        #     print(f"Failed to download {tec_file}.")


        tec_file = f'jplg{day}0.{str(year)[-2:]}i'+'.Z'

        tec_file_path = os.path.join(working_dir,tec_file)
        ### Good idea to force download of a new tec file
        os.system(f"rm -r {tec_file_path}")
        if os.path.exists(tec_file_path):
            logging.info(f"Searching for {tec_file} in {tec_file_path}")
            logging.info(f" Found downloaded tecfile {tec_file_path} ")
            logging.info(f"IONEX file {tec_file} exists. Will not download another")
            input_file = tec_file
            output_file = tec_file.strip('.Z')
            # decompress the file
            logging.info("Decompressing the file")
            try:
                with open(output_file, "w") as output_file_obj:
                    subprocess.run(["zcat", input_file], stdout=output_file_obj, check=True, text=True)
                    ionex_files.append(output_file)
                logging.info(f"File {input_file} has been decompressed and saved as {output_file}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running zcat for {input_file}: {e}")
            except Exception as e:
                logging.error(f"An unexpected error occurred: {e}")

        else:
            ftp_url = f"ftp://gdc.cddis.eosdis.nasa.gov/gnss/products/ionex/{year}/{day}/{tec_file}"
            logging.info(f"Downloading {tec_file}")
            curl_command = f'curl -u "anonymous:daip@nrao.edu" --ftp-ssl-reqd "{ftp_url}" -o "{tec_file}"'
            logging.info(f"Executing {curl_command}")
            try:
                timeout = 3600
                logging.warning(f"FTP will break if no response is received after {timeout}s")
                subprocess.run(curl_command, shell=True, check=True,timeout=timeout)
                logging.info(f"{tec_file} downloaded successfully.")

                input_file = tec_file
                output_file = tec_file.strip('.Z')
                # decompress the file
                logging.info("Decompressing the file")
            
                with open(output_file, "w") as output_file_obj:
                    subprocess.run(["zcat", input_file], stdout=output_file_obj, check=True, text=True)
                    ionex_files.append(output_file)
                    logging.info(f"File {input_file} has been decompressed and saved as {output_file}")
   
            except subprocess.TimeoutExpired:
                logging.error("FTP operation time out")

            except subprocess.CalledProcessError as e:
                logging.error(f"Error downloading file: {e}")


    return ionex_files

def runtecor():

    """
    Calibrates ionospheric delays and Faraday rotation and writes a new CL table
    Using the FIRST file as the INFILE and setting NFILE to len(ionex_files) will find all the files
    Args:
        Uses ionex_files downloaded from NASA
    """

    set_indata()
    indata = set_indata.indata
    ionex_files = download_ionex_files()

    tecor = AIPSTask('TECOR')
    tecor.indata = indata
    # if len(ionex_files) > 1:
    #     for i in range(len(ionex_files)):
    #         tecor.infile = ionex_files[i]
    # else:
    tecor.infile= ionex_files[0]
    tecor.nfiles = len(ionex_files)
    tecor.aparm[1] = 1 # calculate the dispersive delay corrections like VLBATECR
    tecor.gainver = get_table(indata,'CL')
    tecor.gainuse = get_table(indata,'CL')+1
    logging.info(f"Running TECOR and making CL table {get_table(indata,'CL')+1}")
    tecor.go()


def fring_instr(calsour):

    """
    Corrects the phases by calibrating for delays and rates
    Args:
        timerange (int list): The timerange to search for solutions
        type (str): for controlling FRING; instrumental or global fringefitting
    """

    set_indata()
    indata = set_indata.indata
    
    _,_,_,refant_index,refant_indices = get_obs_params()

    if isinstance(calsour, str):
        calsour = [calsour]

    timerange_str = fring_timerange
    timerange = [int(value) for value in timerange_str.split(',')]

    fring = AIPSTask('FRING')
    fring.indata = indata
    fring.docal = 1
    fring.gainuse = get_table(indata,'CL')
    fring.weightit = 1
    fring.refant = refant_index
    fring.solint = 0 # 0 means 10 min
    fring.dparm[9] = 1 # do not fit rates
    fring.timerang[1:] = timerange
    fring.snver = get_table(indata,'SN')+1
    fring.calsour[1:] = calsour
    fring.aparm[6] = 3 # print in detail
    fring.aparm[7] = 7 # snr
    fring.cmethod = 'dft'
    fring.search[1:] = refant_indices


    if len(timerange) !=8:
        logging.error(f"Invalid format for AIPS timerang")
        logging.error("Timerange should consist of 8 entries")
        logging.error("First four entries specify the start day, hour, minute and second")
        logging.error("and the last four give the end day, hour, minute and second")
        logging.critical("Please supply a valid timerange")

    else:
        logging.info(f"Running FRING using timerange: {timerange} and  refant indexed as: {refant_index}")
        logging.info(f"FRING corrections will be derived using timerange: {timerange}")

        fring.go()


 
def apply_solutions(calsour):

    """
    Takes the SN table and makes a new CL table

    Args:
        calsour (str): source used to derive the calibrations 
        sources (str): sources to which the calibrations are to be applied
        interpol (str): the interpolation scheme to be used
        opcode (str): either CALI or CALP; CALP does not flag the CL table
    """


    # NB: This is important since for instrumental fring you generally dont provide
    # the calibrator

    # if calsour.lower() == "phase_calibrator" or calsour.lower() == "calibrator":
    #     calsour = phase_calibrator
    # else:
    #         calsour = ""

    # print(f"The calsour is : {calsour} and is of type {type(calsour)}")

    if isinstance(calsour, str):
        calsour = [calsour]

    if not target or not calsour:
        logging.error("Target or phase calibrator is empty")
    
    else:
        sources = target+calsour
        logging.info(f"Sources: {sources}")

    set_indata()
    indata = set_indata.indata

    _,_,_,refant_index,_ = get_obs_params()

    clcal = AIPSTask('CLCAL')
    clcal.indata = indata
    clcal.calsour[1:]=calsour
    clcal.interpol = ""
    clcal.sources[1:] = sources 
    clcal.snver = get_table(indata,'SN')
    clcal.gainver = get_table(indata,'CL')
    clcal.gainuse = get_table(indata,'CL')+1
    clcal.opcode = "CALI"
    clcal.refant = refant_index
    logging.info(f"Calsour is {calsour}")
    logging.info(f"Applying calibration solutions to {sources} and writing CL table {get_table(indata,'CL')+1}")
    clcal.go()


def global_fring(calsour):

    """
    Corrects the phases by calibrating for delays and rates
    Args:
        timerange (int list): The timerange to search for solutions
        type (str): for controlling FRING; instrumental or global fringefitting
    """

    set_indata()
    indata = set_indata.indata
    
    _,_,_,refant_index,refant_indices = get_obs_params()

    if isinstance(calsour, str):
        calsour = [calsour]

    fring = AIPSTask('FRING')
    fring.indata = indata
    fring.docal = 1
    fring.gainuse = get_table(indata,'CL')
    fring.weightit = 1
    fring.refant = refant_index
    fring.solint = fring_solint
    fring.aparm[9] = 1 # turn on search window
    fring.dparm[1:] = [1,2000,400,integration_time]
    fring.calsour[1:] = calsour
    fring.aparm[5] = 1 # combine all IF
    fring.aparm[6] = 3 # print in detail
    fring.snver = get_table(indata,'SN')+1
    fring.aparm[6] = 3 # print in detail
    fring.aparm[7] = fring_snr # snr
    fring.cmethod = 'dft'
    fring.search[1:] = refant_indices


    logging.info(f"Deriving fring solutions from {calsour}")
    logging.info(f"Using solint {fring_solint} and {refant_index}")
    logging.info(f"Solutions with snr < {fring_snr} will be rejected")
    logging.info(f"Will search {refant_indices} if solutions not found for {refant_index}")

    # print(f"The source used for deriving corrections is {calsour}")
    fring.go()




def runbpass(calsour):

    """
    Calculate bandpass corrections
    """

    set_indata()
    indata = set_indata.indata

    _,_,_,refant_index,_ = get_obs_params()

    if isinstance(calsour, str):
        calsour = [calsour]

    bpass = AIPSTask('BPASS')
    bpass.indata = indata
    bpass.calsour[1:] = calsour
    bpass.docal = 1
    bpass.refant = refant_index
    bpass.solint = -1 # use whole time range
    bpass.weightit = 1
    bpass.cmethod = 'dft'
    bpass.soltype = 'l1r'
    bpass.bpassprm[1] = 0
    bpass.bpassprm[2] = 1
    bpass.bpassprm[10] = 1
    logging.info(f"Deriving solutions from {calsour}")
    bpass.go()


def runsplat_init(target,phase_calibrator,fringe_finder):

    if isinstance(target, str):
        target = [target]
    if isinstance(phase_calibrator, str):
        phase_calibrator = [phase_calibrator]
    if isinstance(fringe_finder, str):
        fringe_finder = [fringe_finder]

    sources = target+phase_calibrator+fringe_finder

    logging.info(f"Sources: {sources} will be SPLAT")

    set_indata()
    indata = set_indata.indata

    splat_file = AIPSUVData(experiment,'SPLAT',1,1)
    if splat_file.exists():
        try:
            logging.info(f"SPLAT file {splat_file} exists. Zapping {splat_file}")
            splat_file.clrstat()
            splat_file.zap(force=True)
        except Exception as e:
            logging.error(f"Unable to zap {splat_file}") 

    doband = -1; docal =1; doflag = 1; flg_table = get_table(indata,'FG'); cal_table = get_table(indata,'CL')
    logging.info("You are applying apriori flags and calibration tables")
    logging.info(f"Applying caltable: {cal_table}, apriori flags: {flg_table}")
    logging.info(f"Bandpass calibration is disabled: doband={doband}")

    splat = AIPSTask('SPLAT')
    splat.indata = indata
    splat.source[1:] = sources
    splat.docal = docal
    splat.gainuse = cal_table
    splat.doband = doband
    splat.bpver = -1 
    splat.flagv = flg_table
    splat.go()

def runsplat_final(target,phase_calibrator,fringe_finder):

    if isinstance(target, str):
        target = [target]
    if isinstance(phase_calibrator, str):
        phase_calibrator = [phase_calibrator]
    if isinstance(fringe_finder, str):
        fringe_finder = [fringe_finder]

    sources = target+phase_calibrator+fringe_finder

    logging.info(f"Sources: {sources} will be SPLAT")

    set_indata()
    indata = set_indata.indata

    splat_file = AIPSUVData(experiment,'SPLAT',1,1)
    if splat_file.exists():
        try:
            logging.info(f"SPLAT file {splat_file} exists. Zapping {splat_file}")
            splat_file.clrstat()
            splat_file.zap(force=True)
        except Exception as e:
            logging.error(f"Unable to zap {splat_file}") 

    
    doband = 1; docal = 1; cal_table = get_table(indata,'CL'); bpver=get_table(indata,'BP');  doflg = -1
        
    logging.info("You are applying calibration and bandpass tables")
    logging.info("Note the table numbering has changed as previous tables were destroyed during conversion to ms")
    logging.info(f"Applying  final caltable: {cal_table} and bandpass table: {bpver}")
    logging.info(f"Bandpass calibration is disabled: doband={doband}")

    splat = AIPSTask('SPLAT')
    splat.indata = indata
    splat.source[1:] = sources
    splat.docal = docal 
    splat.gainuse = cal_table
    splat.doband = doband
    splat.bpver = bpver
    splat.flagv = doflg
    splat.go()


def runfittp(file_extension):

    """
    Write either image of UVFITS to the working directory
    """

    set_indata()

    inseq=indisk = 1
    indata = AIPSUVData(experiment,'SPLAT',inseq,indisk)
    fittp_output =  working_dir +'/'+experiment+'.'+file_extension

    os.system(f"rm -r {fittp_output}")

    fittp = AIPSTask('FITTP')
    fittp.indata = indata
    fittp.dataout = fittp_output
    if not os.path.exists(fittp_output):
        logging.info(f"Writing {fittp_output} to disk")
        fittp.go()
    else:
        logging.info(f"{fittp_output} exists, will not write a new one")

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



def load_fitsfiles(fitsfiles_dir,file_extension):

    logging.info(f"The file extension is: {file_extension}")
    """
    This function loads either idi/fits files or the pipeline calibration in form of TASAV

    Args:
        zap_old_data (bool): deletes old uvdata/tasav data
    """

    fitsfiles = [os.path.join(fitsfiles_dir, f) for f in os.listdir(fitsfiles_dir) 
             if f.endswith(f".{file_extension}")]
    fitsfiles = natsorted(fitsfiles)

    logging.info(f"Found FITS files: {fitsfiles}")

    if not fitsfiles:
        logging.error(f"No FITS files with extension {file_extension} found in {fitsfiles_dir}")
        raise FileNotFoundError("No FITS files to load!")

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



    indata = AIPSUVData(experiment,inclass,inseq,indisk)

    if zap_data==True:
        if indata.exists():
            logging.info(f"Zapping requested")
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


def create_antab_file(workin_dir):

    set_indata()
    indata = set_indata.indata

    if not indata.exists():
        raise FileNotFoundError(f"UVDATA {experiment}.{inclass}.{seq} not found on disk {indisk}")
    
    outtext = os.path.join(workin_dir,f"{experiment}.antab")

    try:
        ty_versions = indata.table('TY',0)  
        if not ty_versions:
            logging.warning(f"No TY table found in {indata.name}.{indata.klass}.{indata.seq}. Skipping ANTAB generation.")
            return
        

        iantb = AIPSTask('IANTB')
        iantb.inname = indata.name
        iantb.inclass = indata.klass
        iantb.inseq = indata.seq
        iantb.indisk = indata.disk
        iantb.antennas = []
        iantb.subarray = 0
        iantb.freqid = -1
        iantb.tyver = 0
        iantb.gcver = 0
        iantb.outtext = outtext

        logging.info(f"Generating ANTAB file: {outtext} from {indata.name}.{indata.klass}.{indata.seq}")
        iantb.go()
        logging.info("ANTAB file generated successfully")

    except Exception as e:
        logging.error(f"Exception while running IANTB: {e}")



def runfittp():

    """
    Write either image of UVFITS to the working directory
    """

    set_indata()
    indata = set_indata.indata

    inseq=indisk = 1
    fittp_output =  working_dir +'/'+experiment+'_exported.UVFITS'

    os.system(f"rm -r {fittp_output}")

    fittp = AIPSTask('FITTP')
    fittp.indata = indata
    fittp.dataout = fittp_output
    if not os.path.exists(fittp_output):
        logging.info(f"Writing {fittp_output} to disk")
        fittp.go()
    else:
        logging.info(f"{fittp_output} exists, will not write a new one")




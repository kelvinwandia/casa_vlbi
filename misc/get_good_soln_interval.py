
from casatasks import *
import casatools

vis = '/raid1/scratch/kelvinw/rsg12_1/rsg12_1.ms'

def query_ms_correct_ant_tables():

    """
    This function opens the ms using the table tools, corrects the antenna
    tables, gets the number of spectral windows which are important during mapping
    in fringefitting and gets the telescope name
    
    """


    tb = casatools.table()

    print("Adding antenna dish diameters to measurement set")

    evn_antennas = ['JB','WB','EF','ON','NT','SR','MC','O6','SH','T6','UR','KM','TR','YS','IR','AR','HH','SV','ZC','BD','GB']
    diams = [76.,25.,100.,25.,32.,65.,32.,20.,25.,65.,25.,40.,32.,40.,32.,305.,26.,32.,32.,32.,100.]
    
    # Antennas that participated in the observation
    tb.open(vis+'/ANTENNA')
    ant_name = tb.getcol('NAME')
    print(ant_name)
    tb.close()

    ### This portion needs fixing - wrong antenna diameters are being placed into the measurement set

    matching_antenna_names = []
    matching_antenna_diams = []
    added_antennas = set() # ensures there is no duplicate
    diameter_map = dict(zip(evn_antennas, diams))

    # Find matching antennas and their diameters
    for name in ant_name:
        if name in diameter_map:
            matching_antenna_names.append(name)
            matching_antenna_diams.append(diameter_map[name])

    print("Matching antenna names:", matching_antenna_names)
    print("Matching antenna diameters:", matching_antenna_diams)

    # # print(matching_antenna_diams)
    
    print("The antennas used in the observations are: ",dict(zip(matching_antenna_names,matching_antenna_diams)))

    # Put the diams in the measurement set
    tb.open(vis+'/ANTENNA',nomodify=False)
    tb.putcol('DISH_DIAMETER',matching_antenna_diams)
    tb.close()

    print("Finished adding antenna dish diameters")


query_ms_correct_ant_tables()
import os,re,sys,json,time,subprocess,glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg') 
from astropy.io import fits
import casatasks
import casatools
from astropy.coordinates import SkyCoord
from astropy.modeling import models
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from scipy.constants import c
from scipy.special import j1 # bessel func of order 1


import logging
from utils.helper_functions import *


vis = '/raid1/scratch/kelvinw/gv020_working_dir/gv020a_working_dir/M15.ms'
### TODO: phasecenters should be supplied by the gaia querying script
offset_sources_coords=['21h29m58.246512s +12d10m01.2339s','21h30m01.203493s +12d10m38.1592s',
            '21h29m58.312403s +12d10m02.6740s','21h29m51.9034555s +12d10m17.13240s',
            '21h30m02.085700s +12d09m04.2203s']

# offset_sources_coords=['21h29m58.246512s +12d10m01.2339s']

pointing_centre = ['21h29m58.350000s +12d10m01.50000s']

def load_primary_beams(pb_file):

    """
    Loads the primary beam json file and gets the station names in the measurement set

    Returns:
        stations: names of antennas that participated in the observation
        antenna_parameters: antenna names, primary beam model and diameters
    """
    pb_json_file = pb_file
    try:
        with open(pb_json_file, 'r') as file:
            primary_beams = json.load(file)
            logging.info(f"======>>> {pb_json_file} found")

    except FileNotFoundError:
        logging.info(f"======>>> {pb_json_file} not found")
        return {}
   

    # Read the stations from the measurement set
    tb = casatools.table()
    tb.open(f'{vis}/ANTENNA',nomodify=False)
    stations = tb.getcol('STATION')
    tb.close()

    antenna_parameters = {}
    antennas_json_file =[]
    pb_params_to_access = ['antenna','diameter','pb_model','pb_freq']
    for antenna, pb_params in primary_beams.items():
        if antenna in stations:
            antennas_json_file.append(antenna)
            antenna_values = {}
            for param in pb_params_to_access:
                if param in pb_params['L']:
                    value = pb_params['L'][param]
                    antenna_values[param] = value
                antenna_parameters[antenna] = antenna_values

    # sort stations to match antennas in primary beams json file
    # sift and return only pb params for antennas that observed
    # reorganise antennas also so correct tables are made - mismatch between ordering of
    # antennas in json file and stations
    sorted_antennas = sorted(stations.tolist(),key=lambda x:antennas_json_file.index(x))
    stations = sorted_antennas
    
    output_file ="stations_pb_params.txt"
    os.system(f'rm -r {output_file}')
    with open(output_file,'w') as file:
        for antenna,values in antenna_parameters.items():
            file.write(f"Antenna: {antenna}\n")
            for param, value in values.items():
                file.write(f"{param}: {value}\n")
            file.write("\n")
    
    return stations, antenna_parameters


def angsep():

    """
    Calculates the angular separation between provided coordinates and the pointing centre and
    converts the coordinates from equatorial to rad using astropy
    # follow this link for the reference eqn
    # https://link.springer.com/referenceworkentry/10.1007/978-0-387-35973-1_761

    Args:
        offset_source: the equatorial coordinates (J2000) of the source

    Returns:
        separation_angle: the angular separation
    """

    separation_angle =  []


    ra_pointing_centre, dec_pointing_centre = pointing_centre[0].split(' ')
    
    
    pointing_centre_skycoord_obj = SkyCoord(ra_pointing_centre,dec_pointing_centre, unit=(u.hourangle,u.deg),frame='icrs')
    
    for source_coords in offset_sources_coords:
        
        source_coords_skycoord_obj = SkyCoord(source_coords.split()[0],source_coords.split()[1],unit=(u.hourangle,u.deg),frame='icrs')
        
        ## This formular needs fixing -- has been fixed
        vector_dot_product =  np.cos(pointing_centre_skycoord_obj.dec.radian)*np.cos(source_coords_skycoord_obj.dec.radian)*\
            np.cos( pointing_centre_skycoord_obj.ra.radian- source_coords_skycoord_obj.ra.radian) + \
            np.sin(pointing_centre_skycoord_obj.dec.radian)*np.sin(source_coords_skycoord_obj.dec.radian)
       
        angle = np.arccos(vector_dot_product)
        # print(f"The offset angle is {angle*180*60/np.pi} arcmin")

        offset_angles_dict = {
            'coordinates':source_coords_skycoord_obj.to_string('hmsdms'),
            'angular_separation_arcmin': angle*180*60/np.pi, # arcmin
            'angular_separation_radians': angle
        }

        separation_angle.append(offset_angles_dict)
    # print(separation_angle)
    return separation_angle


def calculate_pb_attenuations():
    
    """
    Calculates the attenuations using two primary beam models: A Gaussian for EVN antennas
    and Airy function for GBT adn Arecibo

    Calls: load_primary_beam and angsep functions which return
        antenna_parameters: the name, diameter and primary_beam model (Gaussian or Airy)\
        separation_angle: offset angle from the phase centre

    
    """

    # Calling load_primary_beam and angsep

    _, antenna_parameters = load_primary_beams(pb_file) 
    separation_angle = angsep()
    attenuations = {}
    separation_angle_radians = []

    for result in separation_angle:
        angsep_radians = result['angular_separation_radians']
        angsep_coordinates = result['coordinates']
        separation_angle_radians.append(angsep_radians)
    separation_dict = dict(zip([result['coordinates'] for result in separation_angle], separation_angle_radians))    # print(separation_angle_radians)
    # print(separation_dict)

    for antenna, values in antenna_parameters.items():
        # print(antenna,values)
        if 'diameter' in values and 'pb_model' in values and 'pb_freq' in values:
            diameter = values['diameter']
            model = values['pb_model']
            freq = values['pb_freq']
            attenuations[antenna] = {}

            for coordinates, offset_rad in separation_dict.items():
                # print(offset_rad)
                wavelength = c/freq
                P = I = 1

                if model == 'G':
                    factor_gauss = 4 * np.log(2) * diameter**2 * offset_rad**2
                    attenuation = P * np.exp(-factor_gauss / wavelength**2)
                elif model == 'B':
                    factor_bessel = (np.pi / wavelength)*diameter*np.sin(offset_rad)
                    attenuation = I*(2*j1(factor_bessel)/factor_bessel)
                else:
                    attenuation=1.0
                
                attenuations[antenna][f'{coordinates}'] = np.sqrt(attenuation) 

    # for antenna, attenuation in attenuations.items():
    #     print(antenna,attenuation)


    attenuations_file = "stations_pb_attenuation.txt"
    os.system(f'rm -r {attenuations_file}')
    with open(attenuations_file,'a') as file:
        file.write(f"Attenuations: \n")
        for antenna, attenuation in attenuations.items():
            file.write(f"Antenna: {antenna}, Attenuation: {attenuation}\n")
        file.write("\n")

    return attenuations

@time_execution
def gencal_pb_table():
    
    """
    Corrects the loss of flux due to attenuations of the primary beam for each offset source
    and generates a calibration table with the corrections

    Calls:
        calculate_pb_attenuations

    Returns:
        caltables (dict): dict of keys:coordinates (str) and values:caltables (list)  of all the generated calibration tables
    """
    

    # stations,_ =  load_primary_beams()
    attenuations = calculate_pb_attenuations()
    
    # Re-organise the dictionary to get the coordinates as the keys and for get the attenuation for 
    # each antenna
    offset_data = {}
    for antenna, offset_attenuations in attenuations.items():
        # print(antenna,offset_attenuations)
        for offset, attenuation in offset_attenuations.items():
            # print(offset,attenuation)
            offset_data.setdefault(offset, {})[antenna] = attenuation

    attenuations_file = "stations_pb_reoganised_attenuation.txt"
    os.system(f'rm -r {attenuations_file}')
    with open(attenuations_file,'a') as file:
        # file.write(f"Attenuations: \n")
        for coordinate, attenuation in offset_data.items():
            file.write(f"Coordinate: {coordinate}, Attenuation: {attenuation}\n")
        file.write("\n")
    caltables = {}

    for coord, antennas_atten in offset_data.items():
        caltable = coord.replace(' ','')+'.pbcorr'
        # os.system(f'rm -r {caltable}')
        if not os.path.exists(caltable):
        # Get the station name and attenuation and generate a caltable
            for station, antenna_atten in antennas_atten.items():
                logging.info(f"======>>> Generating caltable {caltable} for antenna {station}")
                casatasks.gencal(
                    vis = vis, parameter=antenna_atten, antenna=station,caltype='amp',
                    caltable = caltable
                )
        else:
            # logging.info(f"======>>> Caltable {caltable} exists. Will not generate a new one")
            logging.info(f"======>>> Caltables exists. Will not generate a new one")
            

        if coord not in caltables:
            caltables[coord] = []
        caltables[coord].append(caltable)


        """
        Applies the pbcorr calibration tables and maps the sources
        Will phaseshift the data to the offset position and applycalibrations using mstransform

        Calls:
            pbcorr
        """

        
        for coordinate, table in caltables.items():
            # Write cal_table as txt for docallib in mstransform
            # cal_table = 'caltable_'+coordinate.replace(" ","")+'.txt'
            # os.system(f"rm -r {cal_table}")
            # logging.info(f"======>>> Writing calfile {cal_table}")
            # if not os.path.exists(cal_table):
            #     with open(cal_table,'w') as file:
            #         cal_file = "caltable='"+''.join(table[0])+"'"
            #         file.write(cal_file+'\n')
          
            phaseshifted_ms = coordinate.replace(' ','')+'_phaseshifted.ms'
            
            # subprocess.run(['rm','-r',phaseshifted_ms])
            # os.system(f"rm -r {phaseshifted_ms}*")
            phasecenter = 'J2000'+ ' '+ coordinate
            # if not os.path.exists(phaseshifted_ms):
            #     logging.info(f"======>>> Phaseshifting {vis} to {phasecenter}")
            #     casatasks.phaseshift(
            #         vis = vis, outputvis = phaseshifted_ms, datacolumn='corrected',
            #         phasecenter = phasecenter
            #     )
            
            transformed_ms = coordinate.replace(' ','')+'_transformed'+'.ms'
            # os.system(f"rm -r {transformed_ms}*")
            # if not os.path.exists(transformed_ms):
            #     # logging.info(f"======>>>Transforming {transformed_ms} and applying {cal_table}")
            #     # casatasks.mstransform(
            #     #     vis = phaseshifted_ms,outputvis = transformed_ms,
            #     #     timeaverage=True, timebin='20s',datacolumn='data',
            #     #     chanaverage=True, chanbin=512, docallib = True,
            #     #     callib = cal_file )
            #     casatasks.split(vis=phaseshifted_ms,outputvis=transformed_ms, datacolumn='data',timebin='16s',width=8)
            #     os.system(f"rm -r {phaseshifted_ms}*")
            #     casatasks.applycal(vis=transformed_ms,gaintable=table)

                
            # helper_functions.run_singularity_container
            print(f"Imaging {transformed_ms}")
            imagename = "map_"+coordinate.replace(' ','')
            os.system(f"rm -r {imagename}")
            casatasks.tclean(
                vis = transformed_ms,imagename= imagename, datacolumn='corrected',
                cell = cell, imsize=imsize, deconvolver='clark',
                niter=0
            )
            fitsimage = imagename+'.fits'
            casatasks.exportfits(imagename=imagename+'.image',fitsimage=fitsimage,overwrite=True)
            get_im_stats(imagename+'.image')
            # helper_functions.plot_fits(fitsimage)
            


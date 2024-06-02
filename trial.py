
# import subprocess, os

# jive_plotter = '/home/kelvin/Desktop/singularity/jiveplot_latest.sif'
# vis = '/home/kelvin/Desktop/gv020_working_dir/dont_delete/trial.ms'

# def run_jplotter(command):

#     """
#     Runs wsclean commands 
#     """

#     container = jive_plotter
#     if os.path.exists(container):
#         singularity_bind = os.path.join(os.path.dirname(os.path.dirname(jive_plotter)))

#     command_to_execute = ['singularity', 'run', '-B', singularity_bind, container] + command
#     try:
#         print("Executing: %s", ' '.join(command_to_execute))
#         process = subprocess.Popen(command_to_execute, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#         stdout, stderr = process.communicate()
#         print("stdout: %s", stdout)
#         print("stderr: %s", stderr)

#         return_code = process.returncode
#         if return_code == 0:
#             print(f"Strategy executed successfully. Output:\n{stdout}")
#         else:
#             print(f"Error executing strategy. Return code: {return_code}\nError message: {stderr}")  

#     except Exception as e:
#         print(f"An error occurred: {e}")




# jplotter_cmd = ['jplotter','pt','wt','ms',vis]

# run_jplotter(jplotter_cmd)

# from casatasks import *
# import casatools
# import matplotlib.pyplot as plt
# import numpy as np

# caltable = '/home/kelvin/Desktop/gv020_working_dir/gv020b/trial.sbd'

# tb = casatools.table()

# tb.open(caltable)
# snr = tb.getcol('SNR')
# scanno = tb.getcol('SCAN_NUMBER')
# snr_shape = tb.getcol('SNR').shape 
# scanno_shape = tb.getcol('SCAN_NUMBER').shape
# tb.close()


# # # plt.hist(snr, bins=50, histtype='step', label='sbd_snr' )
# # # plt.show()
# # # print(scanno)
# print("SNR shape:", snr_shape)
# print("SCAN_NUMBER shape:", scanno_shape)

# # Reshape SNR array to flatten the antenna dimension
# snr_flat = snr.reshape(8, 4680).ravel()
# # print(snr_flat)

# scanno_repeated = np.tile(scanno,8)
# # print(scanno_repeated)
# # print(len(scanno_repeated))

# print("Length of flattened SNR:", len(snr_flat))
# print("Length of repeated SCAN_NUMBER:", len(scanno_repeated))


# snr_by_scan = {}

# for scan,snr_value in zip(scanno_repeated,snr_flat):
#     if scan in snr_by_scan:
#         snr_by_scan[scan].append(snr_value)
#     else:
#         snr_by_scan[scan] = [snr_value]

# scan_number_to_check = 1
# snr_values_for_scan = snr_by_scan.get(scan_number_to_check, [])
# print(f"SNR values for scan number {scan_number_to_check}: {snr_values_for_scan}")



"""
Different bit of code
"""
# # Reshape SNR to separate the antenna dimension and extract the first antenna's SNR values
# print("SNR shape:", snr_shape)
# print("Type of snr:", type(snr))
# print("Shape of snr:", snr.shape)
# snr_first_antenna = snr[0, 0, :]

# # Verify the length of the first antenna's SNR array matches the SCAN_NUMBER array length
# print("Length of first antenna's SNR array:", len(snr_first_antenna))  # Should match scanno length

# # Initialize a dictionary to store SNR values by scan number for the first antenna
# snr_by_scan_first_antenna = {}

# # Populate the dictionary with SNR values for each scan number
# for scan_index in range(len(scanno)):
#     scan_number = scanno[scan_index]
#     snr_value = snr_first_antenna[scan_index]
    
#     if scan_number not in snr_by_scan_first_antenna:
#         snr_by_scan_first_antenna[scan_number] = []
    
#     snr_by_scan_first_antenna[scan_number].append(snr_value)

# # Example usage: Access and print SNR values for a specific scan number
# scan_number_to_check = 1
# snr_values_for_scan = snr_by_scan_first_antenna.get(scan_number_to_check, [])
# print(f"SNR values for scan number {scan_number_to_check} for the first antenna: {snr_values_for_scan}")

# # Optional: Save the SNR values to a file for further analysis
# import json

# output_filename = 'snr_first_antenna.json'
# with open(output_filename, 'w') as outfile:
#     json.dump(snr_by_scan_first_antenna, outfile)

# print(f"SNR values for the first antenna have been saved to {output_filename}")



# import glob
# from scipy import stats
# import numpy as np
# msname = '/home/kelvin/Desktop/gv020_working_dir/gv020d/gv020d.ms'
# msmd.open(msname)
# scans = msmd.scansforfield(field='J2139+1423')
# nscans = len(scans)
# scans = scans.astype(int).tolist()

# for i in scans:
#     print(f"Using scan {str(i)}")
#     caltable_name = f'trial_{str(i)}.sbd'
#     if not os.path.exists(caltable_name):
#         fringefit(vis=msname,caltable=caltable_name,zerorates=True,refant='EF,JB,WB',parang=True,minsnr=7,scan=str(i),field='J2139+1423')



# caltables = glob.glob('*.sbd')

# for caltable in caltables:

#     tb = casatools.table()

#     tb.open(caltable)
#     snr = tb.getcol('SNR').ravel()
#     snr_shape = tb.getcol('SNR').shape 
#     tb.close()


#     # plt.hist(snr, bins=50, histtype='step', label=f'sbd_snr_{caltable}' )
#     # plt.show()
#     print(f"SNR for caltable {caltable}\n is :{snr}")

#     # print( 'median = {0}'.format( np.median( snr ) ) )
#     print( 'P(<=7) = {0}'.format( stats.percentileofscore( snr, 7) ) )


# import numpy as np
# from scipy import stats
# import casatools
# import glob

# caltables = glob.glob('*.sbd')
# max_percentage = 0
# best_caltable = None
# percentiles = {}

# for caltable in caltables:
#     tb = casatools.table()
#     tb.open(caltable)
#     snr = tb.getcol('SNR').ravel()
#     tb.close()

#     # Remove NaN values from the snr array
#     snr = snr[~np.isnan(snr)]

#     # Calculate percentage of meaningful (non-zero) data
#     total_data_points = len(snr)
#     meaningful_data_points = np.sum(snr != 0)
#     meaningful_percentage = (meaningful_data_points / total_data_points) * 100

#     if meaningful_percentage > max_percentage:
#         max_percentage = meaningful_percentage
#         best_caltable = caltable

#     # Calculate the percentage of scores greater than 7
#     percentile = 100 - stats.percentileofscore(snr, 50)
#     percentiles[caltable] = percentile

#     print(f'{caltable} - P(>7) = {percentile}, Meaningful Data Percentage = {meaningful_percentage:.2f}%')

# print(f'\nThe caltable with the highest percentage of meaningful data is: {best_caltable} with {max_percentage:.2f}% meaningful data.')

# # Optional: print the percentile score of the best caltable
# if best_caltable in percentiles:
#     print(f'The percentage of scores greater than 7 for the best caltable {best_caltable} is {percentiles[best_caltable]}')


# caltables = glob.glob('*.sbd')

# for caltable in caltables:
#     tb = casatools.table()
#     tb.open(caltable)
#     snr = tb.getcol('SNR')  # Get SNR array
#     snr_shape = tb.getcol('SNR').shape 
#     antennas = tb.getcol('ANTENNA1') # Get antenna array
#     antennas_shape = tb.getcol('ANTENNA1').shape  
#     tb.close()

#     unique_antennas = np.unique(antennas)  # Get unique antennas
#     print(len(unique_antennas))
#     print(snr_shape)
#     print(antennas_shape)

#     for antenna in unique_antennas:
#         # Filter SNR values for the current antenna

#         snr_antenna = snr[antennas == antenna]
#         print(snr_antenna)

#         # Calculate the percentage of SNR for the current antenna
#         total_points = len(snr_antenna)
#         non_zero_points = np.count_nonzero(snr_antenna)
#         percentage = (non_zero_points / total_points) * 100

#         print(f"Caltable: {caltable}, Antenna: {antenna}, Percentage of SNR: {percentage:.2f}%")




# max_percentage = 0
# best_caltable = None
# percentiles = {}

# for tablename in caltables:
#     tb.open(tablename + '/ANTENNA')
#     antenna_names = tb.getcol('NAME')
#     tb.close()
#     tb.open(tablename)
#     antenna_ids = tb.getcol('ANTENNA1')
#     # times  = tb.getcol('TIME')
#     flags = tb.getcol('FLAG')
#     delays = tb.getcol('FPARAM')
#     snrs = tb.getcol('SNR')
#     tb.close()
#     # Analyse number of good solutions:
#     good_frac = []
#     good_snrs = []
#     for i, ant_id in enumerate(np.unique(antenna_ids)):
#         cond = antenna_ids == ant_id
#         # t = times[cond]
#         f = flags[0, 0, :][cond]
#         p = delays[0, 0, :][cond]
#         snr = snrs[0, 0, :][cond]
#         frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
#         snr_mean = np.nanmean(snr[~f])
#         good_frac.append(frac)
#         good_snrs.append(snr_mean)
#     sort_idx = np.argsort(good_frac)[::-1]
#     print('Antennas sorted by % of good solutions:')
#     for i in sort_idx:
#         print('{0:3}: {1:4.1f}, <SNR> = {2:4.1f}'.format(antenna_names[i],
#                                                             good_frac[i],
#                                                             good_snrs[i]))
#     if good_frac[sort_idx[0]] < 90:
#         print('Small fraction of good solutions with selected refant!')
#         print('Please inspect antennas to select optimal refant')
#         print('You may want to use refantmode= flex" in default_params')


# max_good_antennas = 0
# least_flagged_percentage = 100
# best_caltable = None
# best_antenna_flags = None

# for tablename in caltables:
#     tb.open(tablename + '/ANTENNA')
#     antenna_names = tb.getcol('STATION')
#     tb.close()
#     tb.open(tablename)
#     antenna_ids = tb.getcol('ANTENNA1')
#     flags = tb.getcol('FLAG')
#     delays = tb.getcol('FPARAM')
#     snrs = tb.getcol('SNR')
#     tb.close()
    
#     # Analyze number of good solutions for each antenna
#     good_antennas = 0
#     total_unflagged_percentage = 0
#     antenna_flags = {}

#     for i, ant_id in enumerate(np.unique(antenna_ids)):
#         cond = antenna_ids == ant_id
#         f = flags[0, 0, :][cond]
#         snr = snrs[0, 0, :][cond]
#         unflagged_frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
#         if unflagged_frac == 100:
#             good_antennas += 1
#         total_unflagged_percentage += unflagged_frac
#         antenna_flags[antenna_names[i]] = unflagged_frac
    
#     # Calculate the average unflagged percentage across all antennas
#     avg_unflagged_percentage = total_unflagged_percentage / len(np.unique(antenna_ids))
    
#     # Update best_caltable if necessary
#     if good_antennas > max_good_antennas or (good_antennas == max_good_antennas and avg_unflagged_percentage < least_flagged_percentage):
#         max_good_antennas = good_antennas
#         least_flagged_percentage = avg_unflagged_percentage
#         best_caltable = tablename
#         best_antenna_flags = antenna_flags

# # Print the antennas on source and the percentage of unflagged data
# print(f"The calibration table with the most antennas on source and the least flagged data is '{best_caltable}':")
# for antenna, unflagged_percentage in best_antenna_flags.items():
#     print(f"- Antenna {antenna}: {unflagged_percentage:.2f}% unflagged data")




# import glob
# from scipy import stats
# import numpy as np
# msname = '/home/kelvin/Desktop/gv020_working_dir/gv020d/gv020d.ms'
# msmd.open(msname)
# scans = msmd.scansforfield(field='J2139+1423')
# nscans = len(scans)
# scans = scans.astype(int).tolist()

# max_good_antennas = 0
# least_flagged_percentage = 100
# smallest_zero_snr_percentage = 100
# best_caltable = None
# best_antenna_flags = None


# for i in scans[0:2]:
#     print(f"Using scan {str(i)}")
#     tablename = f'trial_{str(i)}.sbd'
#     if not os.path.exists(tablename):
#         fringefit(vis=msname,caltable=tablename,zerorates=True,refant='EF,JB,WB',parang=True,minsnr=20,scan=str(i),field='J2139+1423')

#     tb.open(tablename + '/ANTENNA')
#     antenna_names = tb.getcol('STATION')
#     tb.close()
#     tb.open(tablename)
#     antenna_ids = tb.getcol('ANTENNA1')
#     flags = tb.getcol('FLAG')
#     delays = tb.getcol('FPARAM')
#     snrs = tb.getcol('SNR')
#     tb.close()
    
#     # Analyze number of good solutions for each antenna
#     good_antennas = 0
#     total_unflagged_percentage = 0
#     antenna_flags = {}

#     for i, ant_id in enumerate(np.unique(antenna_ids)):
#         cond = antenna_ids == ant_id
#         f = flags[0, 0, :][cond]
#         snr = snrs[0, 0, :][cond]
#         unflagged_frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
#         if unflagged_frac == 100:
#             good_antennas += 1
#         total_unflagged_percentage += unflagged_frac
#         antenna_flags[antenna_names[i]] = unflagged_frac
    
#     # Calculate the average unflagged percentage across all antennas
#     avg_unflagged_percentage = total_unflagged_percentage / len(np.unique(antenna_ids))
    
#     # Calculate the percentage of zero SNR values
#     total_snrs = np.prod(snrs.shape)
#     zero_snrs_count = np.count_nonzero(snrs == 0)
#     zero_snr_percentage = (zero_snrs_count / total_snrs) * 100
    
#     # Update best_caltable based on the criteria
#     if good_antennas > max_good_antennas \
#         or (good_antennas == max_good_antennas and avg_unflagged_percentage < least_flagged_percentage) \
#         or (good_antennas == max_good_antennas and avg_unflagged_percentage == least_flagged_percentage and zero_snr_percentage < smallest_zero_snr_percentage):
#         max_good_antennas = good_antennas
#         least_flagged_percentage = avg_unflagged_percentage
#         smallest_zero_snr_percentage = zero_snr_percentage
#         best_caltable = tablename
#         best_antenna_flags = antenna_flags

# # Print the best calibration table with all criteria
# if best_caltable:
#     print(f"The calibration table with the most antennas on source, the least flagged data, and the smallest percentage of zero SNR values is '{best_caltable}':")
#     for antenna, unflagged_percentage in best_antenna_flags.items():
#         print(f"- Antenna {antenna}: {unflagged_percentage:.2f}% unflagged data")
#     print(f"Smallest percentage of zero SNR values: {smallest_zero_snr_percentage:.2f}%")
# else:
#     print("No calibration table found meeting all criteria.")




# max_good_antennas = 0
# least_flagged_percentage = 100
# smallest_zero_snr_percentage = 100
# best_caltable = None
# best_antenna_flags = None

# for tablename in caltables:
#     tb.open(tablename + '/ANTENNA')
#     antenna_names = tb.getcol('STATION')
#     tb.close()
#     tb.open(tablename)
#     antenna_ids = tb.getcol('ANTENNA1')
#     flags = tb.getcol('FLAG')
#     delays = tb.getcol('FPARAM')
#     snrs = tb.getcol('SNR')
#     tb.close()
    
#     # Analyze number of good solutions for each antenna
#     good_antennas = 0
#     total_unflagged_percentage = 0
#     antenna_flags = {}

#     for i, ant_id in enumerate(np.unique(antenna_ids)):
#         cond = antenna_ids == ant_id
#         f = flags[0, 0, :][cond]
#         snr = snrs[0, 0, :][cond]
#         unflagged_frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
#         if unflagged_frac == 100:
#             good_antennas += 1
#         total_unflagged_percentage += unflagged_frac
#         antenna_flags[antenna_names[i]] = unflagged_frac
    
#     # Calculate the average unflagged percentage across all antennas
#     avg_unflagged_percentage = total_unflagged_percentage / len(np.unique(antenna_ids))
    
#     # Calculate the percentage of zero SNR values
#     total_snrs = np.prod(snrs.shape)
#     zero_snrs_count = np.count_nonzero(snrs == 0)
#     zero_snr_percentage = (zero_snrs_count / total_snrs) * 100
    
#     # Update best_caltable based on the criteria
#     if good_antennas > max_good_antennas \
#         or (good_antennas == max_good_antennas and avg_unflagged_percentage < least_flagged_percentage) \
#         or (good_antennas == max_good_antennas and avg_unflagged_percentage == least_flagged_percentage and zero_snr_percentage < smallest_zero_snr_percentage):
#         max_good_antennas = good_antennas
#         least_flagged_percentage = avg_unflagged_percentage
#         smallest_zero_snr_percentage = zero_snr_percentage
#         best_caltable = tablename
#         best_antenna_flags = antenna_flags

# # Print the best calibration table with all criteria
# if best_caltable:
#     print(f"The calibration table with the most antennas on source, the least flagged data, and the smallest percentage of zero SNR values is '{best_caltable}':")
#     for antenna, unflagged_percentage in best_antenna_flags.items():
#         print(f"- Antenna {antenna}: {unflagged_percentage:.2f}% unflagged data")
#     print(f"Smallest percentage of zero SNR values: {smallest_zero_snr_percentage:.2f}%")
# else:
#     print("No calibration table found meeting all criteria.")



# import glob, os, re
# from scipy import stats
# import numpy as np

# from casatasks import *
# import casatools

# msmd = casatools.msmetadata()
# tb = casatools.table()

# msname = '/raid1/scratch/kelvinw/gv020_working_dir/gv020b/gv020b.ms'

# # Determine the directory containing the .ms file
# ms_dir = os.path.dirname(msname)

# new_working_dir = os.path.join(ms_dir, 'sbd_files')
# os.makedirs(new_working_dir, exist_ok=True)
# os.chdir(new_working_dir)

# # Print the new current working directory to verify
# print(f"Current working directory: {os.getcwd()}")


# msmd.open(msname)
# scans = msmd.scansforfield(field='J2139+1423')
# nscans = len(scans)
# scans = scans.astype(int).tolist()

# best_scan = None
# least_flagged_percentage = 100
# refant = "EF,JB,WB,GB,MC"
# for scan in scans:
#     print(f"Calculating single band delay solutions for {str(scan)}")
#     tablename = f'trial_{str(scan)}.sbd'
#     if not os.path.exists(tablename):
#         fringefit(vis=msname,caltable=tablename,zerorates=True,refant=refant,parang=True,minsnr=20,scan=str(scan),field='J2139+1423')

#     tb.open(tablename + '/ANTENNA')
#     antenna_names = tb.getcol('STATION')
#     tb.close()
#     tb.open(tablename)
#     antenna_ids = tb.getcol('ANTENNA1')
#     flags = tb.getcol('FLAG')
#     snrs = tb.getcol('SNR')
#     tb.close()
    
#     good_frac = []
#     good_snrs = []

#     total_flagged_percentage = 0

#     for i, ant_id in enumerate(np.unique(antenna_ids)):
#         cond = antenna_ids == ant_id
#         # t = times[cond]
#         f = flags[0, 0, :][cond]
#         snr = snrs[0, 0, :][cond]
#         frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
#         flagged_frac = 100 - frac
#         snr_mean = np.nanmean(snr[~f])
#         good_frac.append(frac)
#         good_snrs.append(snr_mean)
#         total_flagged_percentage += flagged_frac
    
#     avg_flagged_percentage = total_flagged_percentage / len(np.unique(antenna_ids))

#     if avg_flagged_percentage < least_flagged_percentage:
#         least_flagged_percentage = avg_flagged_percentage
#         best_scan = scan

#     sort_idx = np.argsort(good_frac)[::-1]
#     print(f"Antennas sorted by % of good solutions for scan: {scan}")

#     for i in sort_idx:
#         print(f"{antenna_names[i]:<3}: {good_frac[i]:4.1f}, <SNR> = {good_snrs[i]:4.1f}")

# if best_scan is not None:
#     print(f"The best scan is {best_scan} with the least flagged data percentage of {least_flagged_percentage:.2f}%")
# else:
#     print("No scans found.")



"""
This code will sort the scans and write them to a txt file
"""


import glob, os, re
from scipy import stats
import numpy as np

from casatasks import *
import casatools

msmd = casatools.msmetadata()
tb = casatools.table()

msname = '/raid1/scratch/kelvinw/gv020_working_dir/gv020b/gv020b.ms'

# Determine the directory containing the .ms file
ms_dir = os.path.dirname(msname)

new_working_dir = os.path.join(ms_dir, 'sbd_files')
os.makedirs(new_working_dir, exist_ok=True)
os.chdir(new_working_dir)

# Print the new current working directory to verify
print(f"Current working directory: {os.getcwd()}")


msmd.open(msname)
scans = msmd.scansforfield(field='J2139+1423')
nscans = len(scans)
scans = scans.astype(int).tolist()


results = []
refant = "EF,JB,WB,GB,MC"

for scan in scans:
    print(f"Calculating single band delay solutions for {str(scan)}")
    tablename = f'trial_{str(scan)}.sbd'
    if not os.path.exists(tablename):
        fringefit(vis=msname, caltable=tablename, zerorates=True, refant=refant, parang=True, minsnr=20, scan=str(scan), field='J2139+1423')

    tb.open(tablename + '/ANTENNA')
    antenna_names = tb.getcol('STATION')
    tb.close()
    tb.open(tablename)
    antenna_ids = tb.getcol('ANTENNA1')
    flags = tb.getcol('FLAG')
    snrs = tb.getcol('SNR')
    tb.close()
    
    good_frac = []
    good_snrs = []

    total_flagged_percentage = 0

    for i, ant_id in enumerate(np.unique(antenna_ids)):
        cond = antenna_ids == ant_id
        f = flags[0, 0, :][cond]
        snr = snrs[0, 0, :][cond]
        frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
        flagged_frac = 100 - frac
        snr_mean = np.nanmean(snr[~f])
        good_frac.append(frac)
        good_snrs.append(snr_mean)
        total_flagged_percentage += flagged_frac
    
    avg_flagged_percentage = total_flagged_percentage / len(np.unique(antenna_ids))

    results.append((scan, avg_flagged_percentage, good_frac, good_snrs, antenna_names))

# Sort results based on the least flagged data percentage
results.sort(key=lambda x: x[1])

best_scan = results[0][0]
least_flagged_percentage = results[0][1]


output_file = "sorted_scans_info.txt"
os.system(f"rm -r {output_file}")
with open(output_file, 'w') as f:
    for result in results:
        scan, avg_flagged_percentage, good_frac, good_snrs, antenna_names = result
        f.write(f"Antennas sorted by % of good solutions for scan: {scan}\n")

        sort_idx = np.argsort(good_frac)[::-1]
        for i in sort_idx:
            f.write(f"    {antenna_names[i]:<3}: {good_frac[i]:4.1f}%, <SNR> = {good_snrs[i]:4.1f}\n")

    if results:
        best_scan = results[0][0]
        least_flagged_percentage = results[0][1]
        f.write(f"\nThe best scan is {best_scan} with the least flagged data percentage of {least_flagged_percentage:.2f}%\n")
    else:
        f.write("No scans found.\n")

for result in results:
    scan, avg_flagged_percentage, good_frac, good_snrs, antenna_names = result
    print(f"Antennas sorted by % of good solutions for scan: {scan}")

    sort_idx = np.argsort(good_frac)[::-1]
    for i in sort_idx:
        print(f"{antenna_names[i]:<3}: {good_frac[i]:4.1f}, <SNR> = {good_snrs[i]:4.1f}")

if best_scan is not None:
    print(f"The best scan is {best_scan} with the least flagged data percentage of {least_flagged_percentage:.2f}%")
else:
    print("No scans found.")

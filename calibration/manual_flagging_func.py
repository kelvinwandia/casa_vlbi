import os, re
import casatools
import casaplotms
from casatasks import *

"""
This function reads the casalogfile and writes the selected points from viewer for flagging
The assumption here is that the autocorrelations have been flagged
"""
working_dir = '/raid1/scratch/kelvinw/rsg12_1'
log = 'casa-20240704-084616.log'
outputfile = 'rsg12_1.flag'
log_file = os.path.join(working_dir,log)
flagging_file = os.path.join(working_dir,outputfile)
def read_logfile():


    # Initialize a list to store lines containing the keyword 'locate'
    locate_lines = []
    

    # Open and read the log file
    with open(log_file, 'r') as file:
        for line in file:
            if 'locate' in line:
                scan_index = line.find('Scan')
                corr_index = line.find('X')
                # Extract the substring starting from "Scan"
                modified_line = line[scan_index:corr_index]
                modified_line = re.sub(r"(Field=\S+)\s\[\d+\]", r"\1", modified_line)
                modified_line = re.sub(r'Time=\S+\s', '', modified_line)
                modified_line = re.sub(r'BL=\S+\s\&\s\S+\s\[(.*?)\]', r'BL=\1', modified_line)
                # Remove any remaining square brackets
                modified_line = modified_line.replace('[', '').replace(']', '')
                modified_line = re.sub(r'Spw=(\d+)\sChan=(\d+)', r"spw='\1:\2~\2'", modified_line)
                modified_line = re.sub(r'Freq=\S+\s', '', modified_line)
                # Change Scan to scan
                modified_line = modified_line.replace('Scan=', 'scan=')
                # Change Field to field and add single quotes to the value
                modified_line = re.sub(r'Field=(\S+)', r"field='\1'", modified_line)
                # Change Corr to correlation and add single quotes to the value
                modified_line = re.sub(r'Corr=(\S+)', r"correlation='\1'", modified_line)
                modified_line = re.sub(r'BL=(\S+)', r"antenna='\1'", modified_line)
                # Add the modified line to the list
                # remove the word scan -- flagging file doesnt recognise that
                modified_line = re.sub(r'scan=\S+\s', '', modified_line)
                locate_lines.append(modified_line)

    # Print all the lines with the keyword 'locate'
    for line in locate_lines:
        print(line)

    with open(flagging_file, 'w') as outfile:
        for line in locate_lines:
            outfile.write(line + '\n')

read_logfile()
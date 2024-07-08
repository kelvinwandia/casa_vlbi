import os, re
import casatools
import casaplotms
from casatasks import *

"""
This function reads the casalogfile and writes the selected points from viewer for flagging
The assumption here is that the autocorrelations have been flagged

Seems like the field is not working -- specify field when calling flagdata in CASA
not the issue ... try and figure out whats happening

"""
working_dir = '/raid1/scratch/kelvinw/gv020_working_dir/gv020a_flagging_working_dir'
casalog = 'casa-20240708-043801.log'
outputfile = 'gv020a_1s_avg.flag'
log_file = os.path.join(working_dir,casalog)
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
                # Remove any remaining square brackets
                # modified_line = modified_line.replace('[', '').replace(']', '')
                modified_line = re.sub(r'Spw=(\d+)\sChan=(\d+)', r"spw='\1:\2~\2'", modified_line)
                # modified_line = re.sub(r'Spw=(\d+)\sChan=(\d+)', r"spw='\1:\2'", modified_line)

                modified_line = re.sub(r'Freq=\S+\s', '', modified_line)
                # Change Scan to scan
                modified_line = modified_line.replace('Scan=', 'scan=')
                # Change Field to field and add single quotes to the value
                modified_line = re.sub(r'Field=(\S+)', r"field='\1'", modified_line)
                # Change Corr to correlation and add single quotes to the value
                modified_line = re.sub(r'Corr=(\S+)', r"correlation='\1'", modified_line)
                # Find and replace the BL field value
                modified_line = re.sub(r'BL=\d+@(.*?)\s&\s\d+@(\w+)\s\[.*?(?=\s)', r'BL=\1&\2', modified_line)

                modified_line = re.sub(r'BL=(\S+)', r"antenna='\1'", modified_line)
                # Add the modified line to the list
                # remove the word scan -- flagging file doesnt recognise that
                modified_line = re.sub(r'scan=\S+\s', '', modified_line)
                locate_lines.append(modified_line)


    os.system(f"rm -r {flagging_file}")

    with open(flagging_file, 'w') as outfile:
        for line in locate_lines:
            if line.strip() == "":
                continue
            modified_line = f"mode='manual' {line}"
            print(modified_line)
            outfile.write(modified_line.strip() + '\n')


def remove_duplicates(input_file, output_file):
    seen_lines = set()  # To keep track of seen lines

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip empty lines
            if line.strip() == "":
                continue

            # Remove duplicates
            if line not in seen_lines:
                seen_lines.add(line)
                outfile.write(line)
read_logfile()
input_file = flagging_file
cleaned_flagging_file = 'gv020a_1s_final.flag'
output_file = os.path.join(working_dir,cleaned_flagging_file)
remove_duplicates(input_file, output_file)


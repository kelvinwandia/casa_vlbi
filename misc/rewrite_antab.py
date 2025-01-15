def extract_and_modify_columns(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Check if the line starts with '!' (comment line)
            if line.strip().startswith('!'):
                # Write the comment line as is to the output file
                outfile.write(line)
            else:
                # Split the line by spaces
                columns = line.strip().split()
                # Ensure there are at least 4 columns
                if len(columns) >= 4:
                    # Convert the third and fourth columns to floats, divide by 2, and format them back as strings
                    columns[2] = str(float(columns[2]) / 2)
                    columns[3] = str(float(columns[3]) / 2)
                    # Write the first four modified columns to the output file, joined by a space
                    outfile.write(' '.join(columns[:4]) + '\n')
                else:
                    # Optionally, log or print lines with fewer than 4 columns (e.g., blank lines)
                    print(f"Skipping line (not enough columns): {line.strip()}")


# Specify the input and output file paths
input_file = '/raid1/scratch/kelvinw/casa_vlbi/data/gb.antab'
output_file = '/raid1/scratch/kelvinw/casa_vlbi/data/gb_corrected.antab'

extract_and_modify_columns(input_file,output_file)

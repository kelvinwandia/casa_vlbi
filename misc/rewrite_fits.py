
# import astropy.io.fits as fits
# import os, glob

# fits_path = '/mirror/scratch/kelvinw/gv020_fitsfiles/gv020d'
# fitsfiles = glob.glob(os.path.join(fits_path,'*.IDI*'))
# output_fits_file = 'rewritten_fits.idi'
# # print(fitsfiles[0])

# for fitsfile in fitsfiles:
#     with fits.open(fitsfile,mode='update') as hdul:
#         uv_data = hdul['UV_DATA'].data
#         # print(uv_data.columns)
#         print("Baseline values:", uv_data['BASELINE'])
#         valid_rows = uv_data['BASELINE'] != 0
#         # print(valid_rows)
#         uv_data_clean = uv_data[valid_rows]
#         hdul['UV_DATA'].data = uv_data_clean
#         hdul.writeto(output_fits_file, overwrite=True)
#         print(f"Processed and cleaned FITS file: {fitsfile}")


import astropy.io.fits as fits
import os
import glob

# Define input and output directories
fits_path = '/mirror/scratch/kelvinw/gv020_fitsfiles/gv020d'
fitsfiles = glob.glob(os.path.join(fits_path, '*.IDI*'))

output_dir = fits_path
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Iterate over each FITS file
for fitsfile in fitsfiles:
    print(f"Processing file: {fitsfile}")

    # Open the FITS file in read-only mode to avoid in-place modification
    with fits.open(fitsfile, mode='readonly') as hdul:
        if 'UV_DATA' not in hdul:
            print(f"'UV_DATA' extension not found in {fitsfile}. Skipping.")
            continue

        uv_data = hdul['UV_DATA'].data

        # Check if 'BASELINE' column exists
        if 'BASELINE' not in uv_data.columns.names:
            print(f"'BASELINE' column not found in {fitsfile}. Skipping.")
            continue

        # Identify rows where BASELINE is not zero
        valid_rows = uv_data['BASELINE'] != 0

        # If all rows are zero, trigger the else block
        if not valid_rows.any():
            print(f"No cleaning required for {fitsfile} (all BASELINE values are zero).")
        else:
            print(f"Cleaning required for {fitsfile}. Removing zero BASELINE rows.")
            
            # Keep only valid rows
            uv_data_clean = uv_data[valid_rows]

            # Copy the entire HDUList to preserve the file structure
            hdul_new = fits.HDUList([hdu.copy() for hdu in hdul])

            # Replace the UV_DATA with the cleaned version
            hdul_new['UV_DATA'].data = uv_data_clean

            # Define the output file path
            output_fits_file = os.path.join(output_dir, os.path.basename(fitsfile))

            # Write the cleaned FITS file
            hdul_new.writeto(output_fits_file, overwrite=True)
            print(f"Cleaned file saved to {output_fits_file}")

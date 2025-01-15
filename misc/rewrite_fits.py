import astropy.io.fits as fits
import os, glob

fits_path = '/mirror/scratch/kelvinw/gv020_fitsfiles/gv020d'
fitsfiles = glob.glob(os.path.join(fits_path,'*.IDI*'))
output_fits_file = 'rewritten_fits.idi'
# print(fitsfiles[0])

for fitsfile in fitsfiles:
    with fits.open(fitsfile,mode='update') as hdul:
        uv_data = hdul['UV_DATA'].data
        # print(uv_data.columns)
        print("Baseline values:", uv_data['BASELINE'])
        valid_rows = uv_data['BASELINE'] != 0
        # print(valid_rows)
        uv_data_clean = uv_data[valid_rows]
        hdul['UV_DATA'].data = uv_data_clean
        # hdul.writeto(output_fits_file, overwrite=True)
        print(f"Processed and cleaned FITS file: {fitsfile}")
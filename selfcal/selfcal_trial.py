
import casatools
from casatasks import * 
import os,glob,subprocess
import numpy as np
from casaplotms import plotms

from typing import Tuple, List


working_directory = '/raid1/scratch/kelvinw/k2_18b/working_dir/selfcal'
vis='//raid1/scratch/kelvinw/k2_18b/working_dir/selfcal/K2-18_split_2s_4.ms'
# phase_calibrator = '1848+283'
target = 'K2-18'
refant = 'ea18, ea22, ea07'

cell = '3.07arcsec'
overwrite = True

os.chdir(working_directory)

def get_msinfo():

    tb = casatools.table()
    nchan = []
    msmd = casatools.msmetadata()
    msmd.open(vis)
    bandwidth = msmd.bandwidths()
    nspw = len(bandwidth)
    for spw in range(nspw):
        nchan.append(msmd.nchan(spw))
    msmd.close()
    return nspw,nchan




def do_gaincal(solint: str,refant: str,minsnr: int) ->str:
    
    caltable = f'selfcal_solint_{solint}.tb'
    if not os.path.exists(caltable):
        gaincal(vis = vis, 
                caltable= caltable, 
                solint = solint, 
                refant = refant, 
                calmode ='p',
                gaintype='G',
                minsnr = minsnr
                )
        

    casa_plotting(caltable)

    return caltable

def casa_plotting(caltable: str, antennas_per_plot: int = 9, coloraxis: list = ['spw','scan'] ) -> None:

    msmd = casatools.msmetadata()
    msmd.open(vis)
    num_ants = len(msmd.antennaids())
    msmd.close()

    for i in range(0, num_ants, antennas_per_plot):
        for color in coloraxis:
            plotms(
                vis=caltable,
                xaxis='time',
                yaxis='phase',
                iteraxis='antenna',
                gridcols=3,
                gridrows=3,
                coloraxis=color,
                antenna=f"{i}~{min(i + antennas_per_plot - 1, num_ants - 1)}",
                plotfile=caltable.replace('.tb', f'_{color}_antennas_{i}_{i + antennas_per_plot - 1}.png'),
                width = 1500,
                height = 750,
                dpi = 300,
                showgui = False,
                overwrite =  True,
            )

solution_intervals =['inf','3s','6s','12s','24s','48s','96s','192s','384s']
for solint in solution_intervals:
    do_gaincal(solint=solint,refant='ea18',minsnr=1)

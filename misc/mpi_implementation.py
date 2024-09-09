import os, glob
from mpi4py import MPI
from casatasks import *
import casatools



def makems(vis,fitsidifile):
    print(f"Making {vis} from {fitsidifile}")


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    filepath = '/home/kelvin/Desktop/EVN_continuum_pt1_casa6/fitsfiles/'
    fitsfiles = glob.glob(os.path.join(filepath,'*.IDI*'))

    num_files = len(fitsfiles)
    # print(num_files)
    files_per_process = num_files // size
    extra_files = num_files % size

    start_index = rank*files_per_process + min(rank,extra_files)
    end_index = start_index + files_per_process + (1 if rank<extra_files else 0)
    files_to_process = fitsfiles[start_index:end_index]
    # print(files_to_process)

    # for fitsfile in files_to_process:
    #     # print(fitsfile)
    importfitsidi(vis=fitsfile[0].replace('.IDI','.ms'),fitsidifile=fitsfiles)

    partition()


main()




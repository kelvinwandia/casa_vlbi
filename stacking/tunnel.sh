#!/bin/bash
#SBATCH --time=504:00:00
#SBATCH --mem=
#SBATCH -w compute-0-100
#SBATCH --job-name=jupyter_lab
#SBATCH --output=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.out
#SBATCH --error=/share/nas/kelvinw/sbatch/logs/out-slurm_%j.err


echo ">>> activate conda"
source /home/kelvinw/.bashrc
conda activate meti
cat /etc/hosts


## get tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -s)

## print tunneling instructions to jupyter-log-{jobid}.txt
echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport galahad
    -----------------------------------------------------------------
    "
## start an ipcluster instance and launch jupyter server
jupyter-lab --no-browser --port=$ipnport --ip=$ipnip

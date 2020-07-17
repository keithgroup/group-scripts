#!/bin/bash
#SBATCH --job-name=example-calculation-ts
#SBATCH --output=out-example-calculation-ts.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=0-00:10:00
#SBATCH --cluster=smp

cd $SBATCH_O_WORKDIR
module purge
module load openmpi/3.1.4
module load orca/4.2.0

cp $SLURM_SUBMIT_DIR/*.inp $SLURM_SCRATCH

export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
# Suppresses OpenFabrics (openib) absence nonfatal error.
export OMPI_MCA_btl_base_warn_component_unused=0
# Makes all error messages print.
export OMPI_MCA_orte_base_help_aggregate=0

cd $SLURM_SCRATCH
$(which orca) *.inp


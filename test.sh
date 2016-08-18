#!/bin/bash
#MSUB -l nodes=1:ppn=8
#MSUB -l walltime=00:05:00
#MSUB -e ./error.log
#      if keyword omitted : default is submitting directory
#MSUB -o ./output.log
#      if keyword omitted : default is submitting directory

### start of jobscript
    cd $PBS_O_WORKDIR
    echo "workdir: $PBS_O_WORKDIR"
    NSLOTS=8
    echo "running on $NSLOTS cpus ..."
    
#MAIN SCRIPT
date
/opt/parastation/bin/mpiexec -np $NSLOTS disStrange_mpi_test.exe  disStrange_mpi.ini
date
exit

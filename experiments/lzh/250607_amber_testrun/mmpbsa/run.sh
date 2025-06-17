#!/bin/bash

cd /workspace/MD/experiments/250607_amber_testrun1/mmpbsa/output
mpirun -np 10 $AMBERHOME/bin/MMPBSA.py.MPI -O -i ../mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat \
    -sp ../../processed/com_solv.top -cp ../../processed/com.top -rp ../../processed/rec.top -lp ../../processed/lig.top \
    -y ../../output/sim.mdcrd \
    -eo FINAL_RESULTS_MMPBSA.csv > progress.log 2>&1
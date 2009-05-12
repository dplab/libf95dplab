#! /bin/bash 
#
#  Script to run test problem.
#
if [ $# -lt 2 ]; then
    echo 'need arguments:  data-file num_procs'
    exit 1
fi
#
DATAFILE=$1
NUM_PROC=$2
#
MAINDIR=`pwd`
BINARY=./linear-solver.x
MPIRUN=/opt/openmpi-1.2.7/bin/mpirun 
#
cp $DATAFILE input.txt
$MPIRUN -np $NUM_PROC $BINARY 

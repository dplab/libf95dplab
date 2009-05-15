#! /bin/bash 
#
#  Script to run tests for linear-solver.x
#
INPUTS=(\
   tridiag-5-hbc-Arhs.inp\
   tridiag-5-hbc-nArhs.inp\
   tridiag-5-nbc-Arhs.inp\
   tridiag-5-nbc-nArhs.inp\
)
#
for inp in ${INPUTS[*]}; do
    bn=`basename $inp .inp`
    echo "------------------------$inp"
    echo -n "testing $inp ... "
    ./Run.sh $inp 1 > $bn-1.out
    ./Run.sh $inp 2 > $bn-2.out
    ./Run.sh $inp 3 > $bn-3.out
    ./Run.sh $inp 4 > $bn-4.out
    echo "done"
done

#!/bin/bash

###############################################################################
#                                                                             #
# This script collects the results of VBHF calculations for quantum sampling  #
# inputs.                                                                     #
# The equilibrium structure is supposed to be calculated, see below to change #
# the filename as needed.                                                     #
#                                                                             #
###############################################################################

if [[ -f results_test ]]; then
	rm results_test
fi

MODE=1

while [[ MODE -le 144 ]];
do
	echo "mode" $MODE

	abs=`ls all_0/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk -F "_" '{ print $4 }'`

	ALL_0=`grep "CONVERGED TOTAL VB/HF ENERGY"  all_0/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`
	ALL_1=`grep "CONVERGED TOTAL VB/HF ENERGY"  all_1/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`
	ALL_m1=`grep "CONVERGED TOTAL VB/HF ENERGY"  all_-1/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`
	MONO_0=`grep "CONVERGED TOTAL VB/HF ENERGY"  mono_0/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`
	MONO_1=`grep "CONVERGED TOTAL VB/HF ENERGY"  mono_1/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`
	MONO_m1=`grep "CONVERGED TOTAL VB/HF ENERGY"  mono_-1/result_"$MODE"_*.out | sort -g -t "_" -k 4 | awk '{ print $6 }'`

	echo $abs | tr " " "\n" > abs
	
	echo $ALL_0 | tr " " "\n" > ALL_0
	echo $ALL_1 | tr " " "\n" > ALL_1
	echo $ALL_m1 | tr " " "\n" > ALL_m1
	echo $MONO_0 | tr " " "\n" > MONO_0
	echo $MONO_1 | tr " " "\n" > MONO_1
	echo $MONO_m1 | tr " " "\n" > MONO_m1

	echo "Mode" $MODE >> results_test
	paste abs ALL_0 ALL_1 ALL_m1 MONO_0 MONO_1 MONO_m1 >> results_test

	((MODE++))

done

rm abs ALL_0 ALL_1 ALL_m1 MONO_0 MONO_1 MONO_m1

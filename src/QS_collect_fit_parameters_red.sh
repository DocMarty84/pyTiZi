#!/bin/bash

###############################################################################
#                                                                             #
# This script collects the linear fit parameters of EA, IP and PE for quantum #
# sampling.                                                                   #
#                                                                             #
###############################################################################

if [[ -f data_fit.* ]]; then
	rm data_fit.*
fi

for x in `ls *.plt | sort -g -t "_" -k 2`
do
	NAME=`echo $x | cut -d "." -f1`
	NUMBER=`echo $NAME | cut -d "_" -f2 | cut -d "-" -f1`
	TYPE=`echo $NAME | cut -d "." -f1 | cut -d "-" -f2`
	gnuplot $x
	mv fit.log $NAME.fit

	a=`grep "+/-" $NAME.fit | grep "a" | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | awk '{print $3}'`
	
	echo $NUMBER $b $a >> data_fit.$TYPE

done

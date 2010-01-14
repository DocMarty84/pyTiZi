#!/bin/bash

###############################################################################
#                                                                             #
# This script collects the linear fit parameters of EA, IP and PE for quantum #
# sampling.                                                                   #
#                                                                             #
###############################################################################

if [[ -f *.data_fit ]]; then
	rm data_fit.*
fi

for x in `ls *.plt | sort -g -t "_" -k 2`
do
	NAME=`echo $x | cut -d "." -f1`
	NUMBER=`echo $NAME | cut -d "_" -f2`
	gnuplot $x
	mv fit.log $NAME.fit

	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 1p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 1p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.IP_cluster

	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 2p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 2p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.EA_cluster
	
	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 3p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 3p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.IP_alone

	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 4p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 4p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.EA_alone
	
	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 5p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 5p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.P_plus

	a=`grep "+/-" $NAME.fit | grep "a" | sed -n 6p | awk '{print $3}'`
	b=`grep "+/-" $NAME.fit | grep "b" | sed -n 6p | awk '{print $3}'`
	echo $NUMBER $b $a >> data_fit.P_minus
done

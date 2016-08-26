#!/bin/bash

if [[ $# == 2 ]]; then
	echo "use inputfile $1.mech. additional options are $2"
	$FlameManLib/ScanMan -i $1.mech -3rsS -$2 &> $1.out
	echo "mv" $1.pre "$myData"
	mv $1.pre $myData
elif [[ $# == 1 ]]; then
	echo "use inputfile $1.mech."
	$FlameManLib/ScanMan -i $1.mech -3srS &> $1.out
	echo "mv" $1.pre "$myData"
	mv $1.pre $myData
else
	echo usage: ScanMan.script '<inputfile>.mech <options>'
	echo "       moves outputfile to $myData"
fi

#!/bin/bash

version=${1}
outfile=${2:-./task_to_resubmit.txt}

rm $outfile
touch $outfile

for dir in ~/nobackup/crab/*'_'$version'_'*/; do
	echo "Checking job status: " $dir
	crab status -d "$dir" > tmp
	if !(grep -q COMPLETED tmp); then
		echo $dir | cat >> "${outfile}"
	fi
	rm tmp
done

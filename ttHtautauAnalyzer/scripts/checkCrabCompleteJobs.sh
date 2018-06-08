#!/bin/bash

version=${1}
outfile=${2:-./task_to_purge.txt}
crabdir=${3:-/afs/cern.ch/work/z/ztao/public/workspace/crab}

rm $outfile
touch $outfile

for dir in $crabdir/*'_'$version'_'*/; do
	echo "Checking job status: " $dir
	crab status -d "$dir" > tmp_p
	if (grep -q COMPLETED tmp_p); then
		echo $dir | cat >> "${outfile}"
	fi
	rm tmp_p
done

#!/bin/bash

version=${1}
crabdir=${2:-/afs/cern.ch/work/z/ztao/public/workspace/crab}

for dir in $crabdir/*'_'$version'_'*/; do
	crab status -d "$dir"
done

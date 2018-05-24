#!/bin/bash

version=${1}

for dir in ~/nobackup/crab/*'_'$version'_'*/; do
	crab status -d "$dir"
done

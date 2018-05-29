#!/bin/bash

file=${1:-./task_to_resubmit.txt}

while read dir; do
	crab status -d "$dir"
done <${file}

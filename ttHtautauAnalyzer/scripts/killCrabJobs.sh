#!/bin/bash
list_to_kill=$1

while read dir; do
	echo "Killing crab job: " $dir
	crab kill -d "$dir"
done <${list_to_kill}

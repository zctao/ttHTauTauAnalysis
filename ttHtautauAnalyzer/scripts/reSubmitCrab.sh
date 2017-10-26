#!/bin/bash
infile=${1:-./task_to_resubmit.txt}

while read line; 
do
	#echo "$line"
	crab resubmit -d "$line" #--sitewhitelist=T2_EE_Estonia
done < $infile

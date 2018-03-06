#!/bin/bash
infile=${1:-./task_to_resubmit.txt}

while read line; 
do
	#echo "$line"
	crab resubmit -d "$line" #--siteblacklist=T3_US_UCR #--sitewhitelist=T2_EE_Estonia
done < $infile

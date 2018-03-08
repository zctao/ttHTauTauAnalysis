#!/bin/bash

cleanlist=$1

while read line;
do
	eos root://cmseos.fnal.gov rm -r $line
done < $cleanlist

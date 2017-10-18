#!/bin/bash

# hadd root files in eos storage given a list of directories
# usage: ./haddEOSRoot.sh <target root file name> <list of directories in eos>

target=$1
inlist=$2

touch tmp

while read line;
do
	#xrdfsls -u $line | grep '.root' >> $tmp
	xrdfs root://cmseos.fnal.gov ls -u $line | grep '.root' >> tmp
done < $inlist

hadd $target @tmp
rm tmp

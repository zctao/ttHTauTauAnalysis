#!/bin/bash

# hadd crab output ntuples given analysis type and samples names
# usage: ./haddCrabOutputs.sh <anatype> <samplelist.txt> <out.root> <sample1> <sample2> ...
# example: ./haddCrabOutputs.sh 1l2tau ../data/SampleList_FastSim2016.txt test.root ttH_p1 ttH_p2 ttH_p3

analysis_type=$1
samplelist=$2
outname=$3
shift
shift
shift
samples=$@

#echo $analysis_type
#echo $samplelist
#echo $samples

# make list of directories containing root files to be addes
dumpCrabOutputDirList.py $analysis_type $samples -l $samplelist -s >> list.tmp

# hadd root files
haddEOSRoot.sh $outname list.tmp

rm list.tmp

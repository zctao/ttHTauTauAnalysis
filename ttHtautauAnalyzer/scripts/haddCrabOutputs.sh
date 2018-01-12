#!/bin/bash

# hadd crab output ntuples given analysis type and samples names
# usage: ./haddCrabOutputs.sh <anatype> <samplelist.txt> <eosdirectory> <version> <prefix> <log> <out> <sample1> <sample2> ...
# example: ./haddCrabOutputs.sh 1l2tau ../data/SampleList_FastSim2016.txt /store/user/ztao/ttHtautau_80X/ train_nov2017 crab_m17_ nlist.log test ttH_p1 ttH_p2 ttH_p3

analysis_type=$1
samplelist=$2
eosdir=$3
version=$4
prefix=$5
log=$6
outname=$7
shift
shift
shift
shift
shift
shift
shift
samples=$@

# make list of directories containing root files to be addes
dumpCrabOutputDirList.py $analysis_type $samples -l $samplelist -p $prefix -r $eosdir -s >> $outname'.list'

# hadd root files
haddEOSRoot.sh $outname'.root' $outname'.list'

if [ ! -f $outname'.root' ]; then
	echo $outname'.root is not created. exit.'
	rm $outname'.list'
	exit
fi

# sample full name
fn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" $samples)
#echo $fn

# make directory
directory=$eosdir$fn'/'$version'/'
eos root://cmseos.fnal.gov mkdir $directory

# copy ntuple to eos directory
xrdcp -f $outname'.root' root://cmseos.fnal.gov/$directory
xrdcp -f $outname'.list' root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
# TODO check if directory exists
mkdir ~/scratch/$fn
mkdir ~/scratch/$fn'/'$version
mv $outname'.root' ~/scratch/$fn'/'$version/.
mv $outname'.list' ~/scratch/$fn'/'$version/.

# print file name
echo $directory$outname'.root' >> $log

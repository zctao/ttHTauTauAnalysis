#!/bin/bash

sample=$1
index=$2
version=$3
inputfile=$4

workdir=$(pwd)
echo Work directory $workdir

cmssw=CMSSW_9_4_7

echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630

xrdcp -s root://cmseos.fnal.gov//store/user/ztao/Condor/$cmssw.tgz .
tar -xf $cmssw.tgz
rm $cmssw.tgz
cd $cmssw/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh` 2>&1 > /dev/null
scramexit=$?
if [[ $scramexit -ne 0 ]]; then
    echo scram exited with code $scramexit
    exit $scramexit
fi

echo "CMSSW: "$CMSSW_BASE
echo Arguments passed to this script are: $samples $index $version $inputfile

cd $workdir
# copy input file
echo $inputfiles
xrdcp -s root://cmsxrootd.fnal.gov/$inputfile input.root

cp -r $cmssw/src/ttHTauTauAnalysis/ttHtautauAnalyzer/test/qg/ .

cmsRun $cmssw/src/ttHTauTauAnalysis/ttHtautauAnalyzer/test/analyzer2017_cfg.py SampleName=$sample doCutFlow=True isData=False inputFiles=file:input.root

echo "*******************************************"
# output file default name: output_$sample.root
file=output_$sample.root
echo "xrdcp -f output_$sample.root root://cmseos.fnal.gov//store/user/ztao/Condor/eventNtuples/condor_${version}_${sample}_incl/output_${sample}_${index}.root"
xrdcp -f ${file} root://cmseos.fnal.gov//store/user/ztao/Condor/eventNtuples/condor_${version}_${sample}_incl/output_${sample}_${index}.root 2>&1 > /dev/null
xrdexit=$?
if [[ $xrdexit -ne 0 ]]; then
	rm *.root
	echo xrdcp exited with code $xrdexit
	exit $xrdexit
fi
rm ${file}

cd $workdir
rm -rf $cmssw
rm -r qg/
rm input.root

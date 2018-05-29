#!/bin/bash

version=$1
logname=$2
shift
shift
samples=$@

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
echo Arguments passed to this script are: $version $samples

cd $workdir
produceMVANtuplesv2.py $samples --datasetlist $cmssw/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/DatasetList_2017reMiniAODv2.csv --version $version -o ./ --transfer_inputs -c NA JESUp JESDown TESUp TESDown --log $logname

echo "*******************************************"
# copy output to eos
# all output root files should be under $workdir/$version/
#eos root://cmseos.fnal.gov mkdir -p /store/user/ztao/Condor/mvaNtuples/$version
outdir=root://cmseos.fnal.gov//store/user/ztao/Condor/mvaNtuples/$version
cd $workdir/$version/
for file in *.root
do
	echo "xrdcp -f ${file} ${outdir}/${file}"
	xrdcp -f ${file} ${outdir}/${file} 2>&1 > /dev/null
	xrdexit=$?
	if [[ $xrdexit -ne 0 ]]; then
		rm *.root
		echo xrdcp exited with code $xrdexit
		exit $xrdexit
	fi
	rm ${file}	
done
cd $workdir
rm -rf $cmssw
rm -rf $version/

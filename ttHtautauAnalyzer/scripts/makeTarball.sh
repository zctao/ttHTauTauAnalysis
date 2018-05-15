#!/bin/bash

cmssw=${1:-CMSSW_9_4_7}
log=${2:-savedfiles.log}

tar -zcvf $cmssw.tgz $cmssw \
--exclude=src/ttH \
--exclude=src/RecoEgamma \
--exclude=src/PhysicsTools \
--exclude=src/HadTop \
--exclude=src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles \
--exclude=external \
--exclude=tmp \
--exclude=doc \
--exclude=test \
--exclude-vcs \
> $log

# copy to eos
xrdcp -f CMSSW_9_4_7.tgz root://cmseos.fnal.gov//store/user/ztao/Condor/

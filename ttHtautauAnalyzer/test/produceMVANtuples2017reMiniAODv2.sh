#!/bin/bash

redirector='root://cmsxrootd.fnal.gov/'
eostopdir='/store/user/ztao/ttHtaus_94X/'
datasetlist='../dataFiles/DatasetList_2017reMiniAODv2.csv'
outdir='/uscms/home/ztao/nobackup/mvaNtuples/'
version='may2018'

samples='ttHJetToNonbb data_e'
#corrections='jesup jesdown tesup tesdown'

mkdir -p $outdir$version

produceMVANtuplesv2.py $samples --datasetlist $datasetlist --version $version -o $outdir -r $redirector -t $eostopdir #-c $corrections

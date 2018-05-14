#!/bin/bash

redirector='root://cmsxrootd.fnal.gov/'
eostopdir='/store/user/ztao/ttHtaus_94X/'
datasetlist='../dataFiles/DatasetList_2017reMiniAODv2.csv'
outdir='/uscms/home/ztao/nobackup/mvaNtuples/'
version='may2018'

samples='ttHJetToNonbb data_e'
#samples='ttHJetToNonbb TTW TTW_psw TTZ TTZ_M1to10 TTGJets WZ ZZ ZZ_ext WW WWW WWZ WZZ ZZZ TTWW TTTT tZq WWds TTToDiLep TTToDiLep_psw TTToSemiLep TTToSemiLep_psw TTToHad TTToHad_psw ST_sLep ST_sLep_psw ST_tT ST_tTbar ST_tWT ST_tWT_psw ST_tWTbar ST_tWTbar_psw DYJets_M50 DYJets_M50_ext DYJets_M10to50 data_e data_mu data_dieg data_dimu data_mueg'
#corrections='jesup jesdown tesup tesdown'

mkdir -p $outdir$version

produceMVANtuplesv2.py $samples --datasetlist $datasetlist --version $version -o $outdir -r $redirector -t $eostopdir #-c $corrections

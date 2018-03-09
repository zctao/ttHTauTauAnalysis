#!/bin/bash

version=2017

outdir='/uscms/home/ztao/nobackup/mvaNtuples/M17/jan2018/'
#outdir='./'
# timestamp
timestamp=`date +%Y%m%d`

mkdir -p $outdir

samples='ttH TTW TTZ TTGJets TGJets WG ZG WZ ZZ WW WWds WpWp WZZ WWZ WWW ZZZ tZq TTTT fakes_data data_obs'
#samples='fakes_data'
corrections='jesup jesdown tesup tesdown'
datasets='SingleMuon SingleElectron DoubleMuon DoubleEG MuonEG'
ntuplelist='../test/ntuplelist_m17_20180307.log'

produceMVANtuples.py 1l2tau $samples -d $datasets -l $ntuplelist -o $outdir -v $version #-s -c $corrections 

# mv mva ntuple list
mv $outdir'mvaNtuples.txt' $outdir'mvaNtuples_'$timestamp'.txt'

# 2lss1tau

# 3l1tau

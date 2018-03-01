#!/bin/bash

# my eos directory
eosrootdir='/store/user/ztao/ttHtautau_80X/'

# sample list
samplelist='../data/SampleList_FastSim2016.txt'

nlist=genlist_fs16.log

# add time stamp to log file
timestamp=`date +%Y%m%d`
echo $timestamp | cat >> $nlist

# ttH
haddCrabOutputs.sh gen $samplelist $eosrootdir Gen crab_fs16_ $nlist genNtuple_ttH_fs16 ttH_p1 ttH_p2 ttH_p3

# TTZ
haddCrabOutputs.sh gen $samplelist $eosrootdir Gen crab_fs16_ $nlist genNtuple_TTZ_fs16 TTZ

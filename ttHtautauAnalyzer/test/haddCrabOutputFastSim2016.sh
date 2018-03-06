#!/bin/bash

# my eos directory
eosrootdir='/store/user/ztao/ttHtautau_80X/'

# sample list
samplelist='../data/SampleList_FastSim2016.txt'

# version
version='train_feb2018'

# crab job prefix
prefix='crab_fs16_'

# timestamp
timestamp=`date +%Y%m%d`
#echo $timestamp
nlist=ntuplelist_fs16_$timestamp.log
touch $nlist

###############
# hadd ntuples

# ttH
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist $eosrootdir $version $prefix $nlist ntuple_ttH_1l2tau_loose ttH_p1 ttH_p2 ttH_p3
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist $eosrootdir $version $prefix $nlist ntuple_ttH_2lss1tau_loose ttH_p1 ttH_p2 ttH_p3

# TT_DiLep
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TT_DiLep_1l2tau_loose TT_DiLep_p1 TT_DiLep_p2 TT_DiLep_p3
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TT_DiLep_2lss1tau_loose TT_DiLep_p1 TT_DiLep_p2 TT_DiLep_p3

# TT_SemiLep
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TT_SemiLep_1l2tau_loose TT_SemiLep_p1 TT_SemiLep_p2 TT_SemiLep_p3
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TT_SemiLep_2lss1tau_loose TT_SemiLep_p1 TT_SemiLep_p2 TT_SemiLep_p3

# TTW
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TTW_1l2tau_loose TTW
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TTW_2lss1tau_loose TTW

# TTZ
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TTZ_1l2tau_loose TTZ
haddCrabOutputs.sh 2lss1tau $samplelist $eosrootdir $version $prefix $nlist ntuple_TTZ_2lss1tau_loose TTZ

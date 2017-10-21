#!/bin/bash

# my eos directory
eosrootdir='/store/user/ztao/ttH_80X/'

# sample list
samplelist='../data/SampleList_FastSim2016.txt'

# sample full names
tthfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" ttH_p1)
ttdilepfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_DiLep_p1)
ttsemilepfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_SemiLep_p1)
ttwfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTW)
ttzfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTZ)


# ttH
# make directory
directory=$eosrootdir$tthfn'/train_oct2017/'
eos root://cmseos.fnal.gov mkdir $directory
# hadd ntuples
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist ntuple_ttH_1l2tau_loose.root ttH_p1 ttH_p2 ttH_p3
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist ntuple_ttH_2lss1tau_loose.root ttH_p1 ttH_p2 ttH_p3
# copy ntuple to eos directory
xrdcp -f ntuple_ttH_1l2tau_loose.root root://cmseos.fnal.gov/$directory
xrdcp -f ntuple_ttH_2lss1tau_loose.root root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
read something
mv ntuple_ttH_1l2tau_loose.root ~/scratch/.
mv ntuple_ttH_2lss1tau_loose.root ~/scratch/.

# TT_DiLep
# make directory
directory=$eosrootdir$ttdilepfn'/train_oct2017/'
eos root://cmseos.fnal.gov mkdir $directory
# hadd ntuples
# 1l2tau
haddCrabOutputs.sh 1l2tau $samplelist ntuple_TT_DiLep_1l2tau_loose.root TT_DiLep_p1 TT_DiLep_p2 TT_DiLep_p3
# 2lss1tau
haddCrabOutputs.sh 2lss1tau $samplelist ntuple_TT_DiLep_2lss1tau_loose.root TT_DiLep_p1 TT_DiLep_p2 TT_DiLep_p3
# copy ntuple to eos directory
xrdcp -f ntuple_TT_DiLep_1l2tau_loose.root root://cmseos.fnal.gov/$directory
xrdcp -f ntuple_TT_DiLep_2lss1tau_loose.root root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
read something
mv ntuple_TT_DiLep_1l2tau_loose.root ~/scratch/.
mv ntuple_TT_DiLep_2lss1tau_loose.root ~/scratch/.

# TT_SemiLep
# make directory
directory=$eosrootdir$ttsemilepfn'/train_oct2017/'
eos root://cmseos.fnal.gov mkdir $directory
# hadd ntuples
haddCrabOutputs.sh 1l2tau $samplelist ntuple_TT_SemiLep_1l2tau_loose.root TT_SemiLep_p1 TT_SemiLep_p2 TT_SemiLep_p3
haddCrabOutputs.sh 2lss1tau $samplelist ntuple_TT_SemiLep_2lss1tau_loose.root TT_SemiLep_p1 TT_SemiLep_p2 TT_SemiLep_p3
# copy ntuple to eos directory
xrdcp -f ntuple_TT_SemiLep_1l2tau_loose.root root://cmseos.fnal.gov/$directory
xrdcp -f ntuple_TT_SemiLep_2lss1tau_loose.root root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
read something
mv ntuple_TT_SemiLep_1l2tau_loose.root ~/scratch/.
mv ntuple_TT_SemiLep_2lss1tau_loose.root ~/scratch/.

# TTW
# make directory
directory=$eosrootdir$ttwfn'/train_oct2017/'
eos root://cmseos.fnal.gov mkdir $directory
# hadd ntuples
haddCrabOutputs.sh 1l2tau $samplelist ntuple_TTW_1l2tau_loose.root TTW
haddCrabOutputs.sh 2lss1tau $samplelist ntuple_TTW_2lss1tau_loose.root TTW
# copy ntuple to eos directory
xrdcp -f ntuple_TTW_1l2tau_loose.root root://cmseos.fnal.gov/$directory
xrdcp -f ntuple_TTW_2lss1tau_loose.root root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
read something
mv ntuple_TTW_1l2tau_loose.root ~/scratch/.
mv ntuple_TTW_2lss1tau_loose.root ~/scratch/.

# TTZ
# make directory
directory=$eosrootdir$ttzfn'/train_oct2017/'
eos root://cmseos.fnal.gov mkdir $directory
# hadd ntuples
haddCrabOutputs.sh 1l2tau $samplelist ntuple_TTZ_1l2tau_loose.root TTZ
haddCrabOutputs.sh 2lss1tau $samplelist ntuple_TTZ_2lss1tau_loose.root TTZ
# copy ntuple to eos directory
xrdcp -f ntuple_TTZ_1l2tau_loose.root root://cmseos.fnal.gov/$directory
xrdcp -f ntuple_TTZ_2lss1tau_loose.root root://cmseos.fnal.gov/$directory

echo 'Finished copying ntuples to eos. Move local files to ~/scratch/'
read something
mv ntuple_TTZ_1l2tau_loose.root ~/scratch/.
mv ntuple_TTZ_2lss1tau_loose.root ~/scratch/.



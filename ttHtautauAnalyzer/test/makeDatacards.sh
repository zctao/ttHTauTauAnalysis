#!/bin/bash

date=${1:-2018jun26}
mvaNtupleVersion=${2:-jun2018v4}
outdir=${3:-/uscms/home/ztao/nobackup/datacards/$date}

# 1l2tau
echo "Making datacards for 1l2tau"
makeHistograms.py signal_1l2tau mvaOutput 0 1 --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_1l2tau_41p53invfb_noRebin_$date.root --systematics | tee yields_1l2tau_$date.txt
# rebin
binDatacards.py 1l2tau datacards_1l2tau_41p53invfb_noRebin_$date.root -o datacards_1l2tau_41p53invfb_Binned_$date.root -p -c ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs

# 2lss1tau
echo "Making datacards for 2lss1tau"
makeHistograms.py signal_2lss1tau mvaOutput 0 1 --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data flips_data data_obs --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --genTauSplit --version $mvaNtupleVersion -o datacards_2lss1tau_41p53invfb_noRebin_$date.root --closureFile ../dataFiles/datacard/Closure_FR_syst/Closure_FR_lepton_syst_2lss1tau_41p5invfb_jun2018v4.root --systematics | tee yields_2lss1tau_$date.txt
# rebin
binDatacards.py 2lss1tau datacards_2lss1tau_41p53invfb_noRebin_$date.root -o datacards_2lss1tau_41p53invfb_Binned_$date.root -p -c ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data flips_data data_obs

# 3l1tau
echo "Making datacards for 3l1tau"
makeHistograms.py signal_3l1tau mvaOutput 0 1 --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --genTauSplit --version $mvaNtupleVersion -o datacards_3l1tau_41p53invfb_noRebin_$date.root --closureFile ../dataFiles/datacard/Closure_FR_syst/Closure_FR_lepton_syst_3l1tau_41p5invfb_jun2018v4.root --systematics | tee yields_3l1tau_$date.txt
# rebin
binDatacards.py 3l1tau datacards_3l1tau_41p53invfb_noRebin_$date.root -o datacards_3l1tau_41p53invfb_Binned_$date.root -p -c ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs

# 2l2tau
echo "Making datacards for 2l2tau"
makeHistograms.py signal_2l2tau mvaOutput 0 1 --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_2l2tau_41p53invfb_noRebin_$date.root --systematics | tee yields_2l2tau_$date.txt
# rebin
binDatacards.py 2l2tau datacards_2l2tau_41p53invfb_noRebin_$date.root -o datacards_2l2tau_41p53invfb_Binned_$date.root -p -c ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs

mkdir -p $outdir
echo "Moving outputs to "$outdir
mv datacards_*.root $outdir
mv yields_*.txt $outdir

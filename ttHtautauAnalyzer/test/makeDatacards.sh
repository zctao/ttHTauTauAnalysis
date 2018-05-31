#!/bin/bash

date=${1:-2018may31}
mvaNtupleVersion=${2:-may2018v2_1}

# 1l2tau
echo "Making datacards for 1l2tau"
makeHistograms.py signal_1l2tau mvaOutput 0 1 --channels ttH TTW TTZ EWK Rares tH VH Conversion ggH fakes_data --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_1l2tau_41p53invfb_noRebin_$date.root --systematics | tee datacards_1l2tau_$date.txt

# 2lss1tau
echo "Making datacards for 2lss1tau"
makeHistograms.py signal_2lss1tau mvaOutput 0 1 --channels ttH TTW TTZ EWK Rares tH VH Conversion ggH fakes_data flips_data --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_2lss1tau_41p53invfb_noRebin_$date.root --systematics | tee datacards_2lss1tau_$date.txt

# 3l1tau
echo "Making datacards for 3l1tau"
makeHistograms.py signal_3l1tau mvaOutput 0 1 --channels ttH TTW TTZ EWK Rares tH VH Conversion ggH fakes_data --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_3l1tau_41p53invfb_noRebin_$date.root --systematics | tee datacards_3l1tau_$date.txt

# 2l2tau
echo "Making datacards for 2l2tau"
makeHistograms.py signal_2l2tau mvaOutput 0 1 --channels ttH TTW TTZ EWK Rares tH VH Conversion ggH fakes_data --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --hDecaySplit --version $mvaNtupleVersion -o datacards_2l2tau_41p53invfb_noRebin_$date.root --systematics | tee datacards_2l2tau_$date.txt

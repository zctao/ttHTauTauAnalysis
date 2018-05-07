#!/bin/bash

# my eos directory
eosrootdir='/store/user/ztao/ttHtaus_94X/'

# sample list
samplelist='../dataFiles/DatasetList_2017reMiniAODv2.csv'

# version
version='may2018'

# crab job prefix
#prefix='crab_2018may_'
prefix='crab_2018may'

# timestamp
timestamp=`date +%Y%m%d`
#echo $timestamp
nlist=ntuplelist_2017remaodv2_$timestamp.log
touch $nlist

###############
# hadd ntuples
analysis='incl'

#MC
samples='ttHJetToNonbb TTW TTW_psw TTZ TTZ_psw WZ TTToDiLep TTToDiLep_psw TTToSemiLep TTToSemiLep_psw TTToHad TTToHad_psw'
# TTGJets TGJets WG ZG ZZ WW WWds WpWp WZZ WWZ WWW ZZZ tZq TTTT'
#corrections='jesup jesdown tesup tesdown'
corrections=''
skipNominal=false

# Data
datasets='mu e dimu dieg mueg'
dataset_era='2017b 2017c 2017d 2017e 2017f'

# MC
for sample in $samples
do
	start=`date +%s`
	
	if [ "$skipNominal" != true ]; then
		haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$analysis $sample
	fi
	for cor in $corrections
	do
		haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$cor'_'$analysis $sample'_'$cor
	done
	
	end=`date +%s`
	echo 'hadd sample' $sample 'takes' $((end-start))'s'
done

# Data	
for dataset in $datasets
do
	collection=''
	for era in $dataset_era
	do
		collection=$collection'data_'$dataset'_'$era' '
	done
	
	start=`date +%s`
	haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix'_' $nlist ntuple_'data_'$analysis $collection
	end=`date +%s`
	echo 'hadd' $dataset 'dataset takes' $((end-start))'s'
done

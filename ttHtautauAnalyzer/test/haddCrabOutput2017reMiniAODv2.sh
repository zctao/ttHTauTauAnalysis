#!/bin/bash

# sample list
samplelist=${1:-'../dataFiles/DatasetList_2017reMiniAODv2.csv'}
#samplelist='../dataFiles/DatasetList_2017reMiniAODv2.csv'

# my eos directory
eosrootdir='/store/user/ztao/ttHtaus_94X/'

# version
version='may2018v2'

# crab job prefix
prefix='crab_2018mayv2_'

# timestamp
timestamp=`date +%Y%m%d`
#echo $timestamp
nlist=ntuplelist_2017remaodv2_$version_$timestamp.log
touch $nlist

###############
# hadd ntuples
analysis='incl'

#MC
samples=''
#samples='ttHJetToNonbb TTW TTW_psw TTZ TTZ_psw WZ WZTo3LNu TTGJets TTZ_M1to10 ZZ ZZ_ext WW WWTo2L2Nuds WZZ WWZ WWW ZZZ tZq TTTT TTWW ggHZZ4l WWds' #
#samples='ttHToNonbb TTToDiLep TTToDiLep_psw TTToSemiLep TTToSemiLep_psw TTToHad TTToHad_ps DYJets_M50 DYJets_M50_ext DYJets_M10to50 ST_sLep ST_sLep ST_tT ST_tTbar ST_tWT ST_tWT ST_tWTbar ST_tWTbar'
# TGJets WG ZG WJets
# ST_tWll WpWpJJ VHToNonbb tHW tHq

#corrections='jesup jesdown tesup tesdown'
corrections=''
skipNominal=false

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
datasets=''
#datasets='mu e dimu dieg mueg'
dataset_era=''
#dataset_era='2017b 2017c 2017d 2017e 2017f'

#for dataset in $datasets
#do
#	collection=''
#	for era in $dataset_era
#	do
#		collection=$collection'data_'$dataset'_'$era' '
#	done
#	
#	start=`date +%s`
#	haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection
#	end=`date +%s`
#	echo 'hadd' $dataset 'dataset takes' $((end-start))'s'
#done

collection_e='data_e_2017b data_e_2017c data_e_2017c_missingLumis data_e_2017d data_e_2017e data_e_2017f'
#haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection_e

collection_mu='data_mu_2017b data_mu_2017c data_mu_2017d data_mu_2017d_missingLumis data_mu_2017e data_mu_2017f data_mu_2017f_missingLumis'
#haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection_mu

collection_mueg='data_mueg_2017b data_mueg_2017c data_mueg_2017c_missingLumis data_mueg_2017d data_mueg_2017e data_mueg_2017f'
#haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection_mueg

collection_dimu='data_dimu_2017b data_dimu_2017c data_dimu_2017d data_dimu_2017e data_dimu_2017e_missingLumis data_dimu_2017f'
#haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection_dimu

collection_dieg='data_dieg_2017b data_dieg_2017c data_dieg_2017d data_dieg_2017e data_dieg_2017f data_dieg_2017f_missingLumis'
#haddCrabOutputs.sh $analysis $samplelist $eosrootdir $version $prefix $nlist ntuple_'data_'$analysis $collection_dieg

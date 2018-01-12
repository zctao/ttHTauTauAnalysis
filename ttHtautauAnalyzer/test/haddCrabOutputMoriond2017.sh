#!/bin/bash

# my eos directory
eosrootdir='/store/user/ztao/ttHtautau_80X/'

# sample list
samplelist='../data/SampleList_Moriond17.txt'


# version
version='dec2017'

# crab job prefix
prefix='crab_m17_'

# timestamp
timestamp=`date +%Y%m%d`
#echo $timestamp
nlist=ntuplelist_m17_$timestamp.log
touch $nlist

###############
# hadd ntuples
analysis='1l2tau 2lss1tau 3l1tau'
#analysis=''

#MC
samples='ttH TTW TTZ TTGJets TGJets WG ZG WZ ZZ WW WWds WpWp WZZ WWZ WWW ZZZ tZq TTTT'
#samples=''
corrections='jesup jesdown tesup tesdown'
#corrections=''
skipNominal=false

# Data
channels='fakes_data flips_data data_obs'
#channels=''
datasets='mu e dimu dieg mueg'
dataset_era='2016b 2016c 2016d 2016e 2016f 2016g 2016h_v2 2016h_v3'

for anatype in $analysis
do
	echo 'process' $anatype
	# MC
	for sample in $samples
	do
		start=`date +%s`
		if [[ $sample != "TTW" && $sample != "TTGJets" ]]; then
			if [ "$skipNominal" != true ]; then
				haddCrabOutputs.sh $anatype $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$anatype $sample
			fi
			for cor in $corrections
			do
				haddCrabOutputs.sh $anatype $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$cor'_'$anatype $sample'_'$cor
			done
		else
			if [ "$skipNominal" != true ]; then
				haddCrabOutputs.sh $anatype $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$anatype $sample $sample'_ext'
			fi
			for cor in $corrections
			do
				haddCrabOutputs.sh $anatype $samplelist $eosrootdir $version $prefix $nlist ntuple_$sample'_'$cor'_'$anatype $sample'_'$cor $sample'_ext_'$cor
			done
		fi
		end=`date +%s`
		echo 'hadd sample' $sample 'takes' $((end-start))'s'
	done

	# Data	
	for channel in $channels
	do
		if [[ $channel == 'flips_data' && $anatype != "2lss1tau" ]]; then
			continue
		fi
		
		for dataset in $datasets
		do
			collection=''
			for era in $dataset_era
			do
				collection=$collection$channel'_'$dataset'_'$era' '
			done
			
			start=`date +%s`
			haddCrabOutputs.sh $anatype $samplelist $eosrootdir $version $prefix $nlist ntuple_$channel'_'$anatype $collection
			end=`date +%s`
			echo 'hadd' $channel $dataset 'dataset takes' $((end-start))'s'
		done
	done
	
done

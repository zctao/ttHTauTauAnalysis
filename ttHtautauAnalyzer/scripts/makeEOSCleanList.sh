#!/bin/bash

oldestdate=$1
eostopdir=${2:-/store/user/ztao/ttHtautau_80X/}
outfile=${3:-eoscleanlist_$oldestdate.txt}

rm $outfile
touch $outfile

eosls='eos root://cmseos.fnal.gov ls'
samples=$($eosls $eostopdir)

for sample in $samples; do
	echo "checking" $sample
	for subdir in $eostopdir$sample/crab_*/; do
		crabdirs=$($eosls "$subdir")
		#echo $crabdirs
		for crabdir in $crabdirs; do
			timestamps=$($eosls $eostopdir$sample/$crabdir/)
			#echo $eostopdir$sample/$crabdir/
			hasNewerOnes=0
			for timestamp in $timestamps; do
				#echo $timestamp
				date=$(cut -d'_' -f1 <<<"$timestamp" )
				if [ "$date" -lt "$oldestdate" ]; then
					echo $eostopdir$sample/$crabdir/$timestamp/ | cat >> "${outfile}"
				else
					hasNewerOnes=1
				fi
				
			done
			if [ $hasNewerOnes -eq 0 ]; then
				echo "WARNING!" $eostopdir$sample/$crabdir/ "does not have ntuples produced after" $oldestdate
			fi
		done
	done
done

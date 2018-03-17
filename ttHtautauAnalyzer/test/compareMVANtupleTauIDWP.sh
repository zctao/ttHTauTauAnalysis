#!/bin/bash

samples="TT_DiLep TT_SemiLep ttH TTZ TTW"
ntupleDir="/uscms/home/ztao/nobackup/mvaNtuples/TestMar2018/"
taus="taup taum"
variables="easym pt ldgtrkpt ldgtrkE decaymode"
decaymodes="0 1 10"

for sample in $samples; do
	for tau in $taus; do
		for var in $variables; do
			#echo $ntuple $var
		    compareTreeBranches.py $ntupleDir"mvaVars_"$sample"_1l2tau.root" $tau"_"$var $tau"_"$var --cut1 $tau"_tightWP==0" --cut2 $tau"_tightWP==1" -o /uscms/home/ztao/public_html/BDTFeb2018_noGenMatch/plots/$sample/compare_$tau"_"$var.pdf

			for mode in $decaymodes; do
				compareTreeBranches.py $ntupleDir"mvaVars_"$sample"_1l2tau.root" $tau"_"$var $tau"_"$var --cut1 $tau"_tightWP==0&&"$tau"_decaymode=="$mode --cut2 $tau"_tightWP==1&&"$tau"_decaymode=="$mode -o /uscms/home/ztao/public_html/BDTFeb2018_noGenMatch/plots/$sample/compare_$tau"_"$var"_"$mode.pdf
			done
		done
	done
done

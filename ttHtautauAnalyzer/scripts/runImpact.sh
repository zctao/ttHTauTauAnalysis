#!/bin/bash

# Run under CombineHarvester/ttH_htt/ after data cards are made with CombineHarvester
label=${1}
inputcard=${2}

combineTool.py -M T2W -i $inputcard.txt

workdir=$(pwd)
mkdir impact_$label
cd impact_$label

combineTool.py -M Impacts -d $workdir/$inputcard.root -m 125 --expectSignal 1 --allPars --parallel 8 -t -1 --doInitialFit --robustFit 1

combineTool.py -M Impacts -d $workdir/$inputcard.root -m 125 --expectSignal 1 --allPars --parallel 8 -t -1 --robustFit 1 --doFits

combineTool.py -M Impacts -m 125 -d $workdir/$inputcard.root -o impacts.json

plotImpacts.py -i impacts.json -o impacts_$label

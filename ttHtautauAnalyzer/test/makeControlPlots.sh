#!/bin/bash

analysis=${1}
variable=${2}
mvaNtupleVersion=${3}
xlabel=${4}
xmin=${5}
xmax=${6}
nbins=${7:-100}

outdir="/uscms/home/ztao/public_html/Datacards/CR_plots/$mvaNtupleVersion/$analysis"
mkdir -p $outdir

#channels="ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs"
channels="ttH TTW TTWW TTZ EWK Rares conversions fakes_data data_obs"
if [ "$analysis" = "2lss1tau" ]; then
	channels="$channels flips_data"
fi

makeHistograms.py control_$analysis $variable $xmin $xmax -b $nbins --channels $channels --luminosity 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --version $mvaNtupleVersion -o $outdir/${variable}_${analysis}_ctrl_41p53invfb.root

if [ "$variable" = "mvaOutput" ]; then
	binDatacards.py $analysis $outdir/${variable}_${analysis}_ctrl_41p53invfb.root -o $outdir/${variable}_${analysis}_ctrl_41p53invfb_binned.root -p -c $channels --nosystematics
	mv $outdir/${variable}_${analysis}_ctrl_41p53invfb_binned.root $outdir/${variable}_${analysis}_ctrl_41p53invfb.root
fi

# plot
label="$analysis CR"
plotvar=$variable
if [ "$variable" = "mvaOutput" ]; then
	plotvar="x"
fi

echo plotDatacards.py $plotvar $outdir/${variable}_${analysis}_ctrl_41p53invfb.root -t "\"${label}"\" -c $channels -o $outdir/${variable}_${analysis}_fakeCR.pdf -x $xlabel

plotDatacards.py $plotvar $outdir/${variable}_${analysis}_ctrl_41p53invfb.root -t "\"${label}"\" -c $channels -o $outdir/${variable}_${analysis}_fakeCR.pdf -x $xlabel



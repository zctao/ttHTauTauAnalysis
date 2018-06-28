#!/bin/bash

version=${1:-jun2018v4}
outdir=${2:-/uscms/home/ztao/public_html/Datacards/ClosureTest/jun2018v4}

mkdir -p $outdir

makeClosureShapesv2.py 1l2tau -l 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv --version $version -o $outdir/ClosureTest_FR_1l2tau_41p53invfb.root -p

makeClosureShapesv2.py 2lss1tau -l 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv --version $version -o $outdir/ClosureTest_FR_2lss1tau_41p53invfb.root -p

makeClosureShapesv2.py 3l1tau -l 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv --version $version -o $outdir/ClosureTest_FR_3l1tau_41p53invfb.root -p

makeClosureShapesv2.py 2l2tau -l 41.53 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv --version $version -o $outdir/ClosureTest_FR_2l2tau_41p53invfb.root -p

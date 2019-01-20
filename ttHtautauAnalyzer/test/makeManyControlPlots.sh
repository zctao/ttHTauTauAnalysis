#!/bin/bash

./makeControlPlots.sh 1l2tau mvaOutput jun2018v4 BDT 0 1
./makeControlPlots.sh 1l2tau mTauTauVis jun2018v4 "m_{#tau#tau}^{vis} [GeV]" 0 200 20
./makeControlPlots.sh 1l2tau nJet jun2018v4 nJet 1.5 12.5 11

./makeControlPlots.sh 2lss1tau mvaOutput jun2018v4 BDT 0 1
./makeControlPlots.sh 2lss1tau nJet jun2018v4 nJet 1.5 12.5 11

./makeControlPlots.sh 3l1tau mvaOutput jun2018v4 BDT 0 1
./makeControlPlots.sh 3l1tau nJet jun2018v4 nJet 1.5 12.5 11

./makeControlPlots.sh 2l2tau mvaOutput jun2018v4 BDT 0 1
./makeControlPlots.sh 2l2tau nJet jun2018v4 nJet 1.5 12.5 11

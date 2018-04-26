#!/bin/bash

### trigger scale factors
if [ ! -d "triggerSF" ]; then
	mkdir triggerSF
fi
## 1l2tau
# single lepton trigger
wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Muon/Run2016BtoH/Muon_Mu22OR_eta2p1_eff.root?raw=true -O triggerSF/Muon_Mu22OR_eta2p1_eff.root
wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Electron/Run2016BtoH/Electron_Ele25WPTight_eff.root?raw=true -O triggerSF/Electron_Ele25WPTight_eff.root
# lepton leg of lepton+tau cross trigger
wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Electron/Run2016BtoH/Electron_Ele24_eff.root?raw=true -O triggerSF/Electron_Ele24_eff.root
wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Muon/Run2016BtoH/Muon_Mu19leg_2016BtoH_eff.root?raw=true -O triggerSF/Muon_Mu19leg_2016BtoH_eff.root
# tau leg of lepton+tau cross trigger
cp /afs/cern.ch/work/t/tstreble/public/triggerSF_weights/trigger_sf_*.root triggerSF/.

# TODO
### lepton ID efficiencies

### lepton fake rate

### electron charge flip rate

### tau ID efficiency and fake rate

### b-tagging
#wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/CSVv2_94XSF_V1_B_F.csv
#wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/DeepCSV_94XSF_V1_B_F.csv # not working

### pile up

wget https://raw.githubusercontent.com/cms-jet/QGLDatabase/master/SQLiteFiles/QGL_cmssw8020_v2.db

### event BDT weights
if [ ! -d "evtMVAWeights" ]; then
	mkdir evtMVAWeights
fi
cd evtMVAWeights

# 1l2tau
if [ ! -d "1l2tau" ]; then
	mkdir 1l2tau
fi
wget https://raw.githubusercontent.com/cms-ttH/ttH-TauRoast/master/data/bdt.vvt/sklearn_tt.xml -O 1l2tau/sklearn_tt.xml
wget https://raw.githubusercontent.com/cms-ttH/ttH-TauRoast/master/data/bdt.vvt/sklearn_ttZ.xml -O 1l2tau/sklearn_ttZ.xml

# 2lss1tau
#if [ ! -d "2lss1tau" ]; then
#	mkdir 2lss1tau
#fi

# 3l1tau
if [ ! -d "3l1tau" ]; then
	mkdir 3l1tau
fi
wget https://raw.githubusercontent.com/CERN-PH-CMG/cmgtools-lite/80X/TTHAnalysis/data/kinMVA/tth/3l_ttbar_BDTG.weights.xml -O 3l1tau/3l_ttbar_BDTG.weights.xml
wget https://raw.githubusercontent.com/CERN-PH-CMG/cmgtools-lite/80X/TTHAnalysis/data/kinMVA/tth/3l_ttV_BDTG.weights.xml -O 3l1tau/3l_ttV_BDTG.weights.xml

cd -

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

### pile up

### event BDT weights

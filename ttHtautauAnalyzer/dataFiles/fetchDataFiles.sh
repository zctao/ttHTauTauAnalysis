#!/bin/bash

### Dataset list
#wget https://docs.google.com/spreadsheets/d/1p4jsl1u7DqJ6e6MgiiSCFrOLrReJnvKvkW54noOwxTY/export?format=csv -O DatasetList_2017reMiniAODv2.csv

### QGTagger database
#wget https://raw.githubusercontent.com/cms-jet/QGLDatabase/master/SQLiteFiles/QGL_cmssw8020_v2.db -O ../test/qg/QGL_cmssw8020_v2.db

##########

### mva
dir_mva=mva

mkdir -p $dir_mva
cd $dir_mva
# HadTop Tagger
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml

# Hj Tagger
#scp ztao@lxplus.cern.ch:/afs/cern.ch/user/b/binghuan/public/TTHLep/weights_IHEP/*.xml .
cd -

mkdir -p $dir_mva/1l2tau/
cd $dir_mva/1l2tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/1l_2tau_XGB_plainKin_evtLevelTT_TTH_13Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/1l_2tau_XGB_HTT_evtLevelSUM_TTH_VT_17Var.xml
cd -

mkdir -p $dir_mva/2lss1tau/
cd $dir_mva/2lss1tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelTTV_TTH_15Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelTT_TTH_16Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelSUM_TTH_M_18Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_HTT_evtLevelSUM_TTH_M_19Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_HTTMEM_evtLevelSUM_TTH_M_20Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_JointBDT_plainKin_1B_M.xml
cd -

mkdir -p $dir_mva/3l1tau/
cd $dir_mva/3l1tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelTTV_TTH_13Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelTT_TTH_15Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelSUM_TTH_M_12Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_JointBDT_plainKin_1B_M.xml
cd -

mkdir -p $dir_mva/2l2tau/
cd $dir_mva/2l2tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelTTV_TTH_14Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelTT_TTH_11Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelSUM_TTH_VT_13Var.xml
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_JointBDT_plainKin_1B_VT.xml
cd -

### fake rate
# UPDATE ME

### scale factors
dir_sf=trigger_sf
mkdir -p $dir_sf
# trigger sf UPDATE ME
## 1l2tau
# single lepton trigger
#wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Muon/Run2016BtoH/Muon_Mu22OR_eta2p1_eff.root?raw=true -O $dir_sf/Muon_Mu22OR_eta2p1_eff.root
#wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Electron/Run2016BtoH/Electron_Ele25WPTight_eff.root?raw=true -O $dir_sf/Electron_Ele25WPTight_eff.root
# lepton leg of lepton+tau cross trigger
#wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Electron/Run2016BtoH/Electron_Ele24_eff.root?raw=true -O $dir_sf/Electron_Ele24_eff.root
#wget https://github.com/CMS-HTT/LeptonEfficiencies/blob/master/Muon/Run2016BtoH/Muon_Mu19leg_2016BtoH_eff.root?raw=true -O $dir_sf/Muon_Mu19leg_2016BtoH_eff.root
# tau leg of lepton+tau cross trigger
#cp /afs/cern.ch/work/t/tstreble/public/triggerSF_weights/trigger_sf_*.root $dir_sf/.

# b-tagging
#wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/DeepCSV_94XSF_V1_B_F.csv # not working

### lepton ID scale factor and fake rate
mkdir -p lepton_sf
# UPDATE ME

### Data cards
mkdir -p datacard
# UPDATEME

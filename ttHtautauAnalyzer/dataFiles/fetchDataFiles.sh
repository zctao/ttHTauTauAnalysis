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
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/master/data/evtLevel_2018March/1l_2tau_XGB_HTT_evtLevelSUM_TTH_VT_17Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/1l_2tau_XGB_plainKin_evtLevelTT_TTH_13Var.xml
cd -

mkdir -p $dir_mva/2lss1tau/
cd $dir_mva/2lss1tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/master/data/evtLevel_2018March/2lss_1tau_XGB_HTT_evtLevelSUM_TTH_M_18Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelTTV_TTH_15Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelTT_TTH_16Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_plainKin_evtLevelSUM_TTH_M_18Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_HTT_evtLevelSUM_TTH_M_19Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_HTTMEM_evtLevelSUM_TTH_M_20Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2lss_1tau_XGB_JointBDT_plainKin_1B_M.xml
cd -

mkdir -p $dir_mva/3l1tau/
cd $dir_mva/3l1tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/master/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelSUM_TTH_M_12Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelTTV_TTH_13Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_plainKin_evtLevelTT_TTH_15Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/3l_1tau_XGB_JointBDT_plainKin_1B_M.xml
cd -

mkdir -p $dir_mva/2l2tau/
cd $dir_mva/2l2tau/
wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/master/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelSUM_TTH_VT_13Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelTTV_TTH_14Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_plainKin_evtLevelTT_TTH_11Var.xml
#wget https://raw.githubusercontent.com/HEP-KBFI/tth-htt/7b35df2ea61eac0c75cfd036b9ae89e363be57cb/data/evtLevel_2018March/2l_2tau_XGB_JointBDT_plainKin_1B_VT.xml
cd -

### fake rate
# lepton fake rate
wget https://github.com/HEP-KBFI/tth-htt/blob/master/data/FR_lep_ttH_mva090_2017_CERN_2018May29.root?raw=true -O FR_lep_ttH_mva090_2017_CERN_2018May29.root
#OUTDATED scp ztao@lxplus.cern.ch:/afs/cern.ch/user/g/gpetrucc/public_html/drop/plots/ttH/94X/ttH/lepMVA/v1.0.1/fr-comb/ttH_fr-v1_0_1.root .
# tau fake rate
#scp ztao@lxplus.cern.ch:/afs/cern.ch/user/v/veelken/public/ttHAnalysis/jetToTauFakeRate/FR_tau_2017_v1.root .

# electron charge mis-id rate
#scp ztao@lxplus.cern.ch:/afs/cern.ch/work/s/ssawant/public/ttHAnalysis/eChargeMisId/ElectronChargeMisIdRates_2017.root .

### scale factors
dir_sf=trigger_sf
mkdir -p $dir_sf
# loose vs reco
## muon
#wget http://www.hep.uniovi.es/sscruz/ttH/may28/scaleFactors.root
## 1l2tau
# single lepton trigger and lepton leg of the cross triggers
#scp ztao@lxplus.cern.ch:/afs/cern.ch/user/v/veelken/public/triggerSFs2017/*.root $dir_sf/.
# tau leg of lepton+tau cross trigger
# retrieve from https://github.com/truggles/TauTriggerSFs2017/tree/tauTriggers2017_MCv2_PreReMiniaod

# b-tagging
#wget https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation94X/DeepCSV_94XSF_V1_B_F.csv # not working

### lepton ID scale factor and fake rate
mkdir -p lepton_sf
# UPDATE ME
#scp ztao@lxplus.cern.ch:/afs/cern.ch/work/s/sesanche/public/forTTH/SFs_may17/lepMVAEffSF_e_2lss.root lepton_sf/.
#scp ztao@lxplus.cern.ch:/afs/cern.ch/work/s/sesanche/public/forTTH/SFs_may17/lepMVAEffSF_m_2lss.root lepton_sf/.
#scp ztao@lxplus.cern.ch:/afs/cern.ch/work/s/sesanche/public/forTTH/SFs_may17/lepMVAEffSF_e_3l.root lepton_sf/.
#scp ztao@lxplus.cern.ch:/afs/cern.ch/work/s/sesanche/public/forTTH/SFs_may17/lepMVAEffSF_m_3l.root lepton_sf/.


### Data cards
mkdir -p datacard
# UPDATEME

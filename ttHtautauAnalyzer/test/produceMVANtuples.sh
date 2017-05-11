#!/bin/bash
echo 'ttH'
root -b '../macro/makeMVANtuple.cc+("/eos/uscms/store/user/ztao/ttH_80X/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_ttH.root","mvaVars_ttH_loose.root",1)'

echo 'TT'
root -b '../macro/makeMVANtuple.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_TT_DiLep.root","mvaVars_TTDilep_loose.root",1)'
root -b '../macro/makeMVANtuple.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_TT_SemiLep.root","mvaVars_TTSemilep_loose.root",1)'

echo 'TTW'
root -b '../macro/makeMVANtuple.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/train_may2017/loose/output_TTW.root","mvaVars_TTW_loose.root",1)'

echo 'TTZ'
root -b '../macro/makeMVANtuple.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train_may2017/loose/output_TTZ.root","mvaVars_TTZ_loose.root",1)'

#echo 'WZ'
#'/eos/uscms/store/user/ztao/ttH_80X/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/train_may2017/loose/output_WZ.root'

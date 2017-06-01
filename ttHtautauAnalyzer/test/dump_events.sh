#!/bin/bash
echo 'dump TTH events'
root -b '../macro/dumpEvents.cc+("/eos/uscms/store/user/ztao/ttH_80X/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_ttH.root","../python/EventDump_ttH.py",1)'

echo 'dump TTbar events'
root -b '../macro/dumpEvents.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_TT_DiLep.root","../python/EventDump_TTDilep.py",1)'
root -b '../macro/dumpEvents.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_TT_SemiLep.root","../python/EventDump_TTSemilep.py",1)'

echo 'dump TTV events'
root -b '../macro/dumpEvents.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/train_may2017/loose/output_TTW.root","../python/EventDump_TTW.py",1)'
root -b '../macro/dumpEvents.cc+("/eos/uscms/store/user/ztao/ttH_80X/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/train_may2017/loose/output_TTZ.root","../python/EventDump_TTZ.py",1)'

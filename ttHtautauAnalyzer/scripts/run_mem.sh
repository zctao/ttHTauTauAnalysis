#!/bin/bash

config=${1:-'ttHTauTauAnalysis/ttHtautauAnalyzer/mem_cfg.py'}

input=${2:-'/eos/uscms/store/user/ztao/ttH_80X/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/train_may2017/loose/output_ttH.root'}

maxevent=${3}
startevent=${4}
output=${5}

runMEM $config $input $maxevents $startevent $output

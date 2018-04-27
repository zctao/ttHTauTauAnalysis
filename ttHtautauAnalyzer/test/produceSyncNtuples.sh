#!/bin/bash

outdir=/uscms/home/ztao/nobackup/ttHTT_syncNtuple/94X/26apr2018/
mkdir -p $outdir

echo 'producing sync ntuple with event selection turned off'
cmsRun analyzer2017_cfg.py doSync=True TurnOffEvtSel=True
mv output_.root $outdir"output_sync.root"

makeSyncNtuple -d $outdir -o syncNtuple_object.root --makeObjectNtuple true
mv syncNtuple_object.root $outdir

echo "producing inclusive sync ntuple"
cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_incl doCutFlow=True AnalysisType=inclusive
mv output_sync_event_incl.root $outdir

makeSyncNtuple -d $outdir -o syncNtuple_event.root  --make1l2tau true --make2lss1tau true --make3l1tau true --make2l2tau true
mv syncNtuple_event.root $outdir

#echo "producing sync ntuple for inclusive 1l2tau"
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_1l2tau_incl doCutFlow=True AnalysisType=1l2tau SelectionRegion=inclusive_1l2tau
#mv output_sync_event_1l2tau_incl.root $outdir

#echo "producing sync ntuple for inclusive 2lss1tau"
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2lss1tau_incl doCutFlow=True AnalysisType=2lss1tau SelectionRegion=inclusive_2lss1tau
#mv output_sync_event_2lss1tau_incl.root $outdir

#echo "producing sync ntuple for inclusive 3l1tau"
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_3l1tau_incl doCutFlow=True AnalysisType=3l1tau SelectionRegion=inclusive_3l1tau
#mv output_sync_event_3l1tau_incl.root $outdir

#echo "producing sync ntuple for inclusive 2l2tau"
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2l2tau_incl doCutFlow=True AnalysisType=2l2tau SelectionRegion=inclusive_2l2tau
#mv output_sync_event_2l2tau_incl.root $outdir

#makeSyncNtuple -d $outdir -o syncNtuple_event.root  --make1l2tau true --make2lss1tau true --make3l1tau true --make2l2tau true
#mv syncNtuple_event.root $outdir

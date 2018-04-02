#!/bin/bash

outdir=/uscms/home/ztao/nobackup/ttHTT_syncNtuple/94X/1apr2018/

echo 'producing sync ntuple with event selection turned off'
cmsRun analyzer2017_cfg.py doSync=True TurnOffEvtSel=True
mv output_.root $outdir"output_sync.root"

makeSyncNtuple -d $outdir -o syncNtuple_object.root --makeObjectNtuple true
mv syncNtuple_object.root $outdir

echo "producing sync ntuple for inclusive 1l2tau"
cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_1l2tau_incl doCutFlow=True AnalysisType=1l2tau SelectionRegion=inclusive_1l2tau
mv output_sync_event_1l2tau_incl.root $outdir

echo "producing sync ntuple for inclusive 2lss1tau"
cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2lss1tau_incl doCutFlow=True AnalysisType=2lss1tau SelectionRegion=inclusive_2lss1tau
mv output_sync_event_2lss1tau_incl.root $outdir

echo "producing sync ntuple for inclusive 3l1tau"
cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_3l1tau_incl doCutFlow=True AnalysisType=3l1tau SelectionRegion=inclusive_3l1tau
mv output_sync_event_3l1tau_incl.root $outdir

makeSyncNtuple -d $outdir -o syncNtuple_event.root  --make1l2tau true --make2lss1tau true --make3l1tau true
mv syncNtuple_event.root $outdir




#echo 'producing sync ntuple for 2lss1tau signal region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2lss1tau_sr doCutFlow=True
#mv output_sync_event_2lss1tau_sr.root $outdir

#echo 'producing sync ntuple for 2lss1tau fake extrapolation region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2lss1tau_fake doCutFlow=True doSystematics=False SelectionRegion=control_fake_2lss1tau
#mv output_sync_event_2lss1tau_fake.root $outdir

#echo 'producing sync ntuple for 2lss1tau charge flip control region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_2lss1tau_flip doCutFlow=True doSystematics=False SelectionRegion=control_2los1tau
#mv output_sync_event_2lss1tau_flip.root $outdir

#echo 'producing sync ntuple for 1l2tau signal region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_1l2tau_sr doCutFlow=True doSystematics=True AnalysisType=1l2tau SelectionRegion=signal_1l2tau
#mv output_sync_event_1l2tau_sr.root $outdir

#echo 'producing sync ntuple for 1l2tau fake extrapolation region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_1l2tau_fake doCutFlow=True doSystematics=False AnalysisType=1l2tau SelectionRegion=control_fake_1l2tau
#mv output_sync_event_1l2tau_fake.root $outdir

#echo 'producing sync ntuple for 3l1tau signal region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_3l1tau_sr doCutFlow=True doSystematics=True AnalysisType=3l1tau SelectionRegion=signal_3l1tau
#mv output_sync_event_3l1tau_sr.root $outdir

#echo 'producing sync ntuple for 3l1tau fake extrapolation region'
#cmsRun analyzer2017_cfg.py doSync=True SampleName=sync_event_3l1tau_fake doCutFlow=True doSystematics=False AnalysisType=3l1tau SelectionRegion=control_fake_3l1tau
#mv output_sync_event_3l1tau_fake.root $outdir

#makeSyncNtuple -d $outdir -o syncNtuple_event.root --makeObjectNtuple false
#mv syncNtuple_event.root $outdir

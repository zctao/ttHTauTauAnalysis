#!/bin/bash
echo 'producing sync ntuple with event selection turned off'
#cmsRun analyzer2016_cfg.py doSync=True TurnOffEvtSel=True doSystematics=False
#root -b '../macro/makeFlatNtuple.cc+("output_.root","syncNtuple.root")'
#mv output_.root ~/nobackup/ttHTT_syncNtuple/80X/Test/output_sync.root
#mv syncNtuple.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 2lss1tau signal region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_2lss1tau_sr doCutFlow=True
mv output_sync_event_2lss1tau_sr.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 2lss1tau fake extrapolation region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_2lss1tau_fake doCutFlow=True doSystematics=False SelectionRegion=control_fake_2lss1tau
mv output_sync_event_2lss1tau_fake.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 2lss1tau charge flip control region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_2lss1tau_flip doCutFlow=True doSystematics=False SelectionRegion=control_2los1tau
mv output_sync_event_2lss1tau_flip.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 1l2tau signal region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_1l2tau_sr doCutFlow=True doSystematics=True AnalysisType=1l2tau SelectionRegion=signal_1l2tau
mv output_sync_event_1l2tau_sr.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 1l2tau fake extrapolation region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_1l2tau_fake doCutFlow=True doSystematics=False AnalysisType=1l2tau SelectionRegion=control_fake_1l2tau
mv output_sync_event_1l2tau_fake.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 3l1tau signal region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_3l1tau_sr doCutFlow=True doSystematics=True AnalysisType=3l1tau SelectionRegion=signal_3l1tau
mv output_sync_event_3l1tau_sr.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for 3l1tau fake extrapolation region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_3l1tau_fake doCutFlow=True doSystematics=False AnalysisType=3l1tau SelectionRegion=control_fake_3l1tau
mv output_sync_event_3l1tau_fake.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

cd ../macro/
root -b <<EOF
.x makeSyncNtuples.cc+("~/nobackup/ttHTT_syncNtuple/80X/Test/","syncNtuple_event.root")
.q
EOF
mv syncNtuple_event.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.
cd -

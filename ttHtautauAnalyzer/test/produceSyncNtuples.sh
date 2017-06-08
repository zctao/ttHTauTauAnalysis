#!/bin/bash
#echo 'producing sync ntuple with event selection turned off'
#cmsRun analyzer2016_cfg.py doSync=True TurnOffEvtSel=True doSystematics=False
#root -b '../macro/makeFlatNtuple.cc+("output_.root","syncNtuple.root")'
#mv output_.root ~/nobackup/ttHTT_syncNtuple/80X/Test/output_sync.root
#mv syncNtuple.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for signal region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_sr doCutFlow=True
mv output_sync_event_sr.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for fake extrapolation region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_fake doCutFlow=True doSystematics=False SelectionRegion=control_1lfakeable
mv output_sync_event_fake.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

echo 'producing sync ntuple for charge flip control region'
cmsRun analyzer2016_cfg.py doSync=True SampleName=sync_event_flip doCutFlow=True doSystematics=False SelectionRegion=control_2los1tau
mv output_sync_event_flip.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.

cd ../macro/
root -b 'makeSyncNtuples.cc+("~/nobackup/ttHTT_syncNtuple/80X/Test/","syncNtuple_event.root")'
mv syncNtuple_event.root ~/nobackup/ttHTT_syncNtuple/80X/Test/.
cd -

#!/bin/bash
echo 'Comparing object selections in ntuples...'
root -b '../macro/syncNtupleComparer.cc(
"", 
"~/nobackup/ttHTT_syncNtuple/80X/Test/syncNtuple.root","eventTree", 
"~/nobackup/ttHTT_syncNtuple/80X/Summer16/syncNtuple.root","syncTree", 
"","", 
"","", 
"~/public_html/Googoo/syncObjSel/"
)'

echo 'Comparing event selection ntuples...'

echo 'Signal region : '
root -b '../macro/syncNtupleComparer.cc(
"2lss1tau_SR/",
"~/nobackup/ttHTT_syncNtuple/80X/Test/syncNtuple_sr.root","eventTree",
"~/nobackup/ttHTT_syncNtuple/80X/Summer16/syncNtuple_event.root", "syncTree_2lSS1tau_SR",
"","",
"","",
"~/public_html/Googoo/syncEvtSel/"
)'

echo 'Fake extrapolation region : '
root -b '../macro/syncNtupleComparer.cc(
"2lss1tau_Fake/",
"~/nobackup/ttHTT_syncNtuple/80X/Test/syncNtuple_fake.root","eventTree",
"~/nobackup/ttHTT_syncNtuple/80X/Summer16/syncNtuple_event.root", "syncTree_2lSS1tau_Fake",
"","",
"","",
"~/public_html/Googoo/syncEvtSel/"
)'

echo 'Charge flip region : '
root -b '../macro/syncNtupleComparer.cc(
"2lss1tau_Flip/",
"~/nobackup/ttHTT_syncNtuple/80X/Test/syncNtuple_flip.root","eventTree",
"~/nobackup/ttHTT_syncNtuple/80X/Summer16/syncNtuple_event.root", "syncTree_2lSS1tau_Flip",
"","",
"","",
"~/public_html/Googoo/syncEvtSel/"
)'

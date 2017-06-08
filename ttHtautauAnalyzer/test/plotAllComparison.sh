#! /bin/sh

# Signal Region
echo 'Signal Region:'
root -b '../macro/compareNtuples.C(\
"2lSS1tau_SR/",
"syncNtuple_event.root","syncTree_2lSS1tau_SR", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_Summer16_BDTs.root ","syncTree_2lSS1tau_SR", \
"","", \
"","", \
"~/public_html/2017Jun07/syncEvtSel/"
)'

# Fake
echo 'Fake Control Region:'
root -b '../macro/compareNtuples.C(\
"2lSS1tau_Fake/",
"syncNtuple_event.root","syncTree_2lSS1tau_Fake", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_Summer16_BDTs.root ","syncTree_2lSS1tau_Fake", \
"","", \
"","", \
"~/public_html/2017Jun07/syncEvtSel/"
)'

# Charge Flip
echo 'Charge Flip Control Region:'
root -b '../macro/compareNtuples.C(\
"2lSS1tau_Flip/",
"syncNtuple_event.root","syncTree_2lSS1tau_Flip", \
"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_Summer16_BDTs.root ","syncTree_2lSS1tau_Flip", \
"","", \
"","", \
"~/public_html/2017Jun07/syncEvtSel/"
)'

#! /bin/sh

cornell_ntuple_obj="/afs/cern.ch/user/z/ztao/public/ttHTT_Sync/2017/syncNtuple_object.root"
llr_ntuple_obj="/afs/cern.ch/user/c/cmartinp/public/sync_2017/syncNtuple_ttH_object_v3.root"
tallinn_ntuple_obj="/afs/cern.ch/user/k/kaehatah/public/sync_2017/sync_Tallinn_v9.root"

scp ztao@lxplus.cern.ch:$cornell_ntuple_obj cornell_obj.root
scp ztao@lxplus.cern.ch:$llr_ntuple_obj llr_obj.root
scp ztao@lxplus.cern.ch:$tallinn_ntuple_obj tallinn_obj.root

echo 'Object Selection'
root -b <<EOF
.x ../macro/compareNtuples.C(\
"cornell_obj.root","Cornell",\
"llr_obj.root","LLR",\
"tallinn_obj.root","Tallinn",\
"","",\
"syncTree","~/public_html/syncNtuples/2018Apr10/syncObjSel/")
EOF

# Signal Region
#echo 'Signal Region:'
#root -b <<EOF
#.x ../macro/compareNtuples.C(\
#"2lSS1tau_SR/",\
#"syncNtuple_event.root","syncTree_2lSS1tau_SR",\
#"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_#Summer16_MEM_BDTs.root","syncTree_2lSS1tau_SR",\
#"","",\
#"","",\
#"~/public_html/2017Jun07/syncEvtSel/"\
#)
#.q
#EOF

# Fake
#echo 'Fake Control Region:'
#root -b <<EOF
#.x ../macro/compareNtuples.C(\
#"2lSS1tau_Fake/",\
#"syncNtuple_event.root","syncTree_2lSS1tau_Fake",\
#"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_#Summer16_BDTs.root","syncTree_2lSS1tau_Fake",\
#"","",\
#"","",\
#"~/public_html/2017Jun07/syncEvtSel/"\
#)
#.q
#EOF

# Charge Flip
#echo 'Charge Flip Control Region:'
#root -b <<EOF
#.x ../macro/compareNtuples.C(\
#"2lSS1tau_Flip/",\
#"syncNtuple_event.root","syncTree_2lSS1tau_Flip",\
#"/afs/cern.ch/work/t/tstreble/public/syncNtuple_ttH_Htautau/syncNtuple_event_ttH_80X_#Summer16_BDTs.root","syncTree_2lSS1tau_Flip",\
#"","",\
#"","",\
#"~/public_html/2017Jun07/syncEvtSel/"\
#)
#.q
#EOF


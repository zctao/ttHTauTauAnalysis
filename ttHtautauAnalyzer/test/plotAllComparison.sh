#! /bin/sh

#cornell_ntuple_obj="/afs/cern.ch/work/z/ztao/public/ttHTT_Sync/2017/syncNtuple_object_cornell_v11.root"
#llr_ntuple_obj="/afs/cern.ch/user/c/cmartinp/public/sync_2017/syncNtuple_ttH_object_v11.root"
#tallinn_ntuple_obj="/afs/cern.ch/user/k/kaehatah/public/sync_2017/sync_Tallinn_v28.root"
#ihep_ntuple_obj="/afs/cern.ch/user/b/binghuan/public/TTHLep/2018Sync/IHEP_ttHsyncV11.root"

#scp ztao@lxplus.cern.ch:$cornell_ntuple_obj cornell_obj.root
#scp ztao@lxplus.cern.ch:$llr_ntuple_obj llr_obj.root
#scp ztao@lxplus.cern.ch:$tallinn_ntuple_obj tallinn_obj.root
#scp ztao@lxplus.cern.ch:$ihep_ntuple_obj ihep_obj.root

#echo 'Object Selection'
#root -b <<EOF
#.x ../macro/compareNtuples.C(\
#"cornell_obj.root","Cornell",\
#"llr_obj.root","LLR",\
#"tallinn_obj.root","Tallinn",\
#"ihep_obj.root","IHEP",\
#"syncTree","~/public_html/syncNtuples/2018Apr20/syncObjSel/",\
#true)
#EOF

#cornell_ntuple_evt="/afs/cern.ch/work/z/ztao/public/ttHTT_Sync/2017/syncNtuple_event_cornell_v12.root"
#llr_ntuple_evt="/afs/cern.ch/user/c/cmartinp/public/sync_2017/syncNtuple_event_ttH_tautau_v3.root"
#tallinn_ntuple_evt="/afs/cern.ch/user/k/kaehatah/public/sync_2017/sync_Tallinn_v29.root"
#ihep_ntuple_evt="/afs/cern.ch/user/b/binghuan/public/TTHLep/2018Sync/IHEP_EvtSync_V8.root"

#scp ztao@lxplus.cern.ch:$cornell_ntuple_evt cornell_evt.root
#scp ztao@lxplus.cern.ch:$llr_ntuple_evt llr_evt.root
#scp ztao@lxplus.cern.ch:$tallinn_ntuple_evt tallinn_evt.root
#scp ztao@lxplus.cern.ch:$ihep_ntuple_evt ihep_evt.root

# 2lSS1tau SR
mkdir -p ~/public_html/syncNtuples/2018Jun06/syncEvtSel/2lSS1tau_SR
echo '2lSS1tau Signal Region:'
root -b <<EOF
.x ../macro/compareNtuples.C(\
"tallinn_evt.root","Tallinn",\
"cornell_evt.root","Cornell",\
"ihep_evt.root","IHEP",\
"","",\
"syncTree_2lSS1tau_SR","~/public_html/syncNtuples/2018Jun06/syncEvtSel/2lSS1tau_SR/",\
false)
EOF

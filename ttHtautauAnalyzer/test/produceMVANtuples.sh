#!/bin/bash

redirector='root://cmsxrootd.fnal.gov/'
eostopdir='/store/user/ztao/ttH_80X/'
era='train_oct2017'
samplelist='../data/SampleList_FastSim2016.txt'
outdir='/uscms/home/ztao/nobackup/mvaNtuples/'

echo 'ttH'
tthfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" ttH_p1)
ttHntuple_1l2tau=$redirector$eostopdir$tthfn'/'$era'/ntuple_ttH_1l2tau_loose.root'
ttHntuple_2lss1tau=$redirector$eostopdir$tthfn'/'$era'/ntuple_ttH_2lss1tau_loose.root'

root -b <<EOF
.x ../macro/makeMVANtuple.cc+("$ttHntuple_1l2tau","ttH","1l2tau",1,1,"$outdir")
.x ../macro/makeMVANtuple.cc+("$ttHntuple_2lss1tau","ttH","2lss1tau",1,1,"$outdir")
.q
EOF

echo 'TT'
ttdlfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_DiLep_p1)
TTdlntuple_1l2tau=$redirector$eostopdir$ttdlfn'/'$era'/ntuple_TT_DiLep_1l2tau_loose.root'
TTdlntuple_2lss1tau=$redirector$eostopdir$ttdlfn'/'$era'/ntuple_TT_DiLep_2lss1tau_loose.root'

ttslfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_SemiLep_p1)
TTslntuple_1l2tau=$redirector$eostopdir$ttslfn'/'$era'/ntuple_TT_SemiLep_1l2tau_loose.root'
TTslntuple_2lss1tau=$redirector$eostopdir$ttslfn'/'$era'/ntuple_TT_SemiLep_2lss1tau_loose.root'

root -b <<EOF
.x ../macro/makeMVANtuple.cc+("$TTdlntuple_1l2tau","TT_DiLep","1l2tau",1,0,"$outdir")
.x ../macro/makeMVANtuple.cc+("$TTdlntuple_2lss1tau","TT_DiLep","2lss1tau",1,0,"$outdir")
.x ../macro/makeMVANtuple.cc+("$TTslntuple_1l2tau","TT_SemiLep","1l2tau",1,0,"$outdir")
.x ../macro/makeMVANtuple.cc+("$TTslntuple_2lss1tau","TT_SemiLep","2lss1tau",1,0,"$outdir")
.q
EOF

echo 'TTW'
ttwfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTW)
ttwntuple_1l2tau=$redirector$eostopdir$ttwfn'/'$era'/ntuple_TTW_1l2tau_loose.root'
ttwntuple_2lss1tau=$redirector$eostopdir$ttwfn'/'$era'/ntuple_TTW_2lss1tau_loose.root'

root -b <<EOF
.x ../macro/makeMVANtuple.cc+("$ttwntuple_1l2tau","TTW","1l2tau",1,1,"$outdir")
.x ../macro/makeMVANtuple.cc+("$ttwntuple_2lss1tau","TTW","2lss1tau",1,1,"$outdir")
.q
EOF

echo 'TTZ'
ttzfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTZ)
ttzntuple_1l2tau=$redirector$eostopdir$ttzfn'/'$era'/ntuple_TTZ_1l2tau_loose.root'
ttzntuple_2lss1tau=$redirector$eostopdir$ttzfn'/'$era'/ntuple_TTZ_2lss1tau_loose.root'

root -b <<EOF
.x ../macro/makeMVANtuple.cc+("$ttzntuple_1l2tau","TTZ","1l2tau",1,1,"$outdir")
.x ../macro/makeMVANtuple.cc+("$ttzntuple_2lss1tau","TTZ","2lss1tau",1,1,"$outdir")
.q
EOF

#!/bin/bash

redirector='root://cmsxrootd.fnal.gov/'
eostopdir='/store/user/ztao/ttHtautau_80X/'
era='train_feb2018'
samplelist='../data/SampleList_FastSim2016.txt'
#outdir='/uscms/home/ztao/nobackup/mvaNtuples/Feb2018_noGenMatch/'
outdir='/uscms/home/ztao/nobackup/mvaNtuples/Feb2018/'

echo 'ttH'
tthfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" ttH_p1)
ttHntuple_1l2tau=$redirector$eostopdir$tthfn'/'$era'/ntuple_ttH_1l2tau_loose.root'
ttHntuple_2lss1tau=$redirector$eostopdir$tthfn'/'$era'/ntuple_ttH_2lss1tau_loose.root'

makeMVANtuple -i $ttHntuple_1l2tau -o $outdir'mvaVars_ttH_1l2tau.root' --anatype 1l2tau --seltype loose_1l2tau --xsection 0.215 --systematics false #--mc_matching false #--loose_selection true
#makeMVANtuple -i $ttHntuple_2lss1tau -o $outdir'mvaVars_ttH_2lss1tau.root' --anatype 2lss1tau --seltype loose_2lss1tau --systematics false --loose_selection true

echo 'TT'
ttdlfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_DiLep_p1)
TTdlntuple_1l2tau=$redirector$eostopdir$ttdlfn'/'$era'/ntuple_TT_DiLep_1l2tau_loose.root'
TTdlntuple_2lss1tau=$redirector$eostopdir$ttdlfn'/'$era'/ntuple_TT_DiLep_2lss1tau_loose.root'

makeMVANtuple -i $TTdlntuple_1l2tau -o $outdir'mvaVars_TT_DiLep_1l2tau.root' --anatype 1l2tau --seltype loose_1l2tau --xsection 87.3 --systematics false --mc_matching false #--loose_selection true 
#makeMVANtuple -i $TTdlntuple_2lss1tau -o $outdir'mvaVars_TT_DiLep_2lss1tau.root' --anatype 2lss1tau --seltype loose_2lss1tau --systematics false --loose_selection true

ttslfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TT_SemiLep_p1)
TTslntuple_1l2tau=$redirector$eostopdir$ttslfn'/'$era'/ntuple_TT_SemiLep_1l2tau_loose.root'
TTslntuple_2lss1tau=$redirector$eostopdir$ttslfn'/'$era'/ntuple_TT_SemiLep_2lss1tau_loose.root'

makeMVANtuple -i $TTslntuple_1l2tau -o $outdir'mvaVars_TT_SemiLep_1l2tau.root' --anatype 1l2tau --seltype loose_1l2tau --xsection 182. --systematics false --mc_matching false #--loose_selection true 
#makeMVANtuple -i $TTslntuple_2lss1tau -o $outdir'mvaVars_TT_SemiLep_2lss1tau.root' --anatype 2lss1tau --seltype loose_2lss1tau --systematics false --loose_selection true

echo 'TTW'
ttwfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTW)
ttwntuple_1l2tau=$redirector$eostopdir$ttwfn'/'$era'/ntuple_TTW_1l2tau_loose.root'
ttwntuple_2lss1tau=$redirector$eostopdir$ttwfn'/'$era'/ntuple_TTW_2lss1tau_loose.root'

makeMVANtuple -i $ttwntuple_1l2tau -o $outdir'mvaVars_TTW_1l2tau.root' --anatype 1l2tau --seltype loose_1l2tau --xsection 0.204 --systematics false #--mc_matching false #--loose_selection true
#makeMVANtuple -i $ttwntuple_2lss1tau -o $outdir'mvaVars_TTW_2lss1tau.root' --anatype 2lss1tau --seltype loose_2lss1tau --systematics false --loose_selection true

echo 'TTZ'
ttzfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTZ)
ttzntuple_1l2tau=$redirector$eostopdir$ttzfn'/'$era'/ntuple_TTZ_1l2tau_loose.root'
ttzntuple_2lss1tau=$redirector$eostopdir$ttzfn'/'$era'/ntuple_TTZ_2lss1tau_loose.root'

makeMVANtuple -i $ttzntuple_1l2tau -o $outdir'mvaVars_TTZ_1l2tau.root' --anatype 1l2tau --seltype loose_1l2tau --xsection 0.253 --systematics false #--mc_matching false #--loose_selection true
#makeMVANtuple -i $ttzntuple_2lss1tau -o $outdir'mvaVars_TTZ_2lss1tau.root' --anatype 2lss1tau --seltype loose_2lss1tau --systematics false --loose_selection true

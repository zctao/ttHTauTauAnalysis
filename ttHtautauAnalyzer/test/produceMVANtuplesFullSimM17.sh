#!/bin/bash

redirector='root://cmsxrootd.fnal.gov/'
eostopdir='/store/user/ztao/ttHtautau_80X/'
version='jan2018'
samplelist='../data/SampleList_Moriond17.txt'
outdir='/uscms/home/ztao/nobackup/mvaNtuples/Feb2018_FullSim_noGenMatch/'
#outdir='/uscms/home/ztao/nobackup/mvaNtuples/Feb2018_FullSim/'

echo 'ttH'
tthfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" ttH)
ttHntuple_1l2tau=$redirector$eostopdir$tthfn'/'$version'/ntuple_ttH_1l2tau.root'

makeMVANtuple -i $ttHntuple_1l2tau -o $outdir'mvaVars_ttH_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 0.215 --systematics false --mc_matching false

echo 'TT'
ttdlfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTJets_DiLep)
TTdlntuple_1l2tau=$redirector$eostopdir$ttdlfn'/'$version'/ntuple_TTJets_DiLep_1l2tau.root'

makeMVANtuple -i $TTdlntuple_1l2tau -o $outdir'mvaVars_TTJets_DiLep_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 87.3 --systematics false --mc_matching false 

ttltfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTJets_LepT)
TTltntuple_1l2tau=$redirector$eostopdir$ttltfn'/'$version'/ntuple_TTJets_LepT_1l2tau.root'

makeMVANtuple -i $TTltntuple_1l2tau -o $outdir'mvaVars_TTJets_LepT_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 182. --systematics false --mc_matching false

ttltbarfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTJets_LepTbar)
TTltbarntuple_1l2tau=$redirector$eostopdir$ttltbarfn'/'$version'/ntuple_TTJets_LepTbar_1l2tau.root'

makeMVANtuple -i $TTltbarntuple_1l2tau -o $outdir'mvaVars_TTJets_LepTbar_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 182. --systematics false --mc_matching false

echo 'TTW'
ttwfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTW)
ttwntuple_1l2tau=$redirector$eostopdir$ttwfn'/'$version'/ntuple_TTW_1l2tau.root'

makeMVANtuple -i $ttwntuple_1l2tau -o $outdir'mvaVars_TTW_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 0.204 --systematics false --mc_matching false

echo 'TTZ'
ttzfn=$(python -c "import sys; from ttHTauTauAnalysis.ttHtautauAnalyzer.crab_utils import getSampleFullname; print getSampleFullname(sys.argv[1],'$samplelist')" TTZ)
ttzntuple_1l2tau=$redirector$eostopdir$ttzfn'/'$version'/ntuple_TTZ_1l2tau.root'

makeMVANtuple -i $ttzntuple_1l2tau -o $outdir'mvaVars_TTZ_1l2tau.root' --anatype 1l2tau --seltype signal_1l2tau --xsection 0.253 --systematics false --mc_matching false

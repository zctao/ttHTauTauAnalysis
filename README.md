# ttHTauTauAnalysis

## Installation

Setup CMSSW environment:

	cmsrel CMSSW_8_0_26_patch1
	cd CMSSW_8_0_26_patch1/src/
	cmsenv
	git cms-init

For MET Correction:

	git cms-merge-topic cms-met:METRecipe_8020 -u
	git cms-merge-topic cms-met:METRecipe_80X_part2 -u

Electron MVA ID:

	git cms-merge-topic ikrav:egm_id_80X_v2

Get Analyzer:

	git clone https://github.com/zctao/ttHTauTauAnalysis.git

MiniAODHelper:

	git clone https://github.com/cms-ttH/MiniAOD.git
	(currently use branch CMSSW_8_0_24_v1_sync)

LeptonID package shared with ND ttH-Multilepton group:

	git clone https://github.com/cms-ttH/ttH-LeptonID.git ttH/LeptonID

Compile:

	scram b -j 32

Add the area containing the Electron MVA weights:

	cd $CMSSW_BASE/external
	cd slc6_amd64_gcc530/
	git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
	cd data/RecoEgamma/ElectronIdentification/data
	git checkout egm_id_80X_v1
	cd $CMSSW_BASE/src

Get data files before running analyzer:

	cd ttHTauTauAnalysis/ttHtautauAnalyzer/data
	./fetchDataFiles.sh
	cd -

When running with CRAB one needs to add the following option to the crab config file: config.JobType.sendExternalFolder = True This is needed until the PR including this ID will be integrated in CMSSW/cms-data.

## Usage

Produce sync ntuples:

	./ttHtautauAnalyzer/test/produceSyncNtuples.sh

Submit CRAB jobs:

	python crabSubmitter.py --samples='<SampleLists>' --channels='<channels>'

	python crabSubmitter.py -h for help
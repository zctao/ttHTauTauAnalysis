# ttHTauTauAnalysis

## Installation

Setup CMSSW environment:

	cmsrel CMSSW_9_4_7
	cd CMSSW_9_4_7/src/
	cmsenv
	git cms-init

Electron Scale & Smearing + MVA ID:

	git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
	git cms-merge-topic cms-egamma:Egamma80XMiniAODV2_946

Get Analyzer:

	git clone https://github.com/zctao/ttHTauTauAnalysis.git
	cd ttHTauTauAnalysis
	git checkout cmssw_9_4_x
	cd -

Add MVA package if necessary:

	git clone https://github.com/zctao/ttHTauTauMVA.git

Add MEM interface if necessary:

	git clone https://github.com/zctao/ttHTauTau_MEM_Interface.git

Add package for kinematic hadronic top fit:

	git clone https://github.com/zctao/HTT_kinfit.git HadTop/HTT_kinfit

LeptonID package:

	git clone https://github.com/zctao/ttH-LeptonID.git ttH/LeptonID
	
For lepton-tau cross trigger scale factor:

	git clone -b tauTriggers2017_MCv2_PreReMiniaod https://github.com/truggles/TauTriggerSFs2017 TriggerSF/TauTriggerSFs2017

Compile:

	scram b -j 32

Add the area containing the Electron MVA weights:

	cd $CMSSW_BASE/external
	cd slc6_amd64_gcc630/
	git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
	cd data/RecoEgamma/PhotonIdentification/data
	git checkout CMSSW_9_4_0_pre3_TnP
	cd $CMSSW_BASE/external
	cd slc6_amd64_gcc630/
	git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
	cd data/RecoEgamma/ElectronIdentification/data
	git checkout CMSSW_9_4_0_pre3_TnP
	cd $CMSSW_BASE/src

Get data files before running analyzer:

	cd ttHTauTauAnalysis/ttHtautauAnalyzer/data
	./fetchDataFiles.sh
	cd -

When running with CRAB one needs to add the following option to the crab config file: config.JobType.sendExternalFolder = True This is needed until the PR including this ID will be integrated in CMSSW/cms-data.
It'd also be necessary to remove old MVA weights in RecoEgamma/PhotonIdentification/data/ and RecoEgamma/ElectronIdentification/data/ to reduce crab sandbox tarball size (max 100 MB).

## Usage

* General structures:

For 2017 analysis

    cmsRun analyzer2017_cfg.py

Runs ttHtautauAnalyzer implemented in plugins/ and produces an event ntuple with the ntuple format defined in interface/eventNtuple.h.
It applies a loose selection on the input sample (pass\_ttH\_ltau\_inclusive\_selection() implemented in src/EventSelector.cc), hopefully covering all event categories for the interest of this analysis.

The inclusive event ntuple needs to be further processed by binaries to make synchronization ntuples:

    makeSyncNtuple

or to make ntuples for training MVA or for producing datacards:

    makeMVANtuple

The ntuple formats are defined in interface/syncNtuple.h and interface/mvaNtuple.h respectively.

Event selections for a region of interest (e.g. signal region selection for 2lss1tau category) are applied in either of the binaries. The type of the event selections can be specified by the command line arguments.

Scale factors, event weights, systematics, etc. are also handled during this step.

* Produce sync ntuples:

In the test/ directory

	./produceSyncNtuples.sh <output directory> <input root file>

* CRAB jobs

All the following scripts assume CRAB work area: ~/nobackup/crab/

Submit CRAB jobs to produce event ntuples:
	   
	submitCrabJobsv2.py <datasetlist.csv> --sample_list <sample_list> --prefix <job_label>
	
	submitCrabJobsv2.py -h for help

Check CRAB job status:

    checkCrabJobStatus.sh <job_label>

It generates a list of incomplete jobs, named task\_to\_resubmit.txt by default.

    checkCrabJobStatusFromFile.sh <list_of_jobs>

Same as the above, except it only loops over the jobs in the list. It generates an updated list of incomplete jobs, named task\_to\_resubmit\_updated.txt by default.

Resubmit CRAB jobs:

    reSubmitCrab.sh <list_of_jobs>

Resubmit CRAB jobs in the list. By default it reads task\_to\_resubmit.txt, unless specified otherwise.

Collect and hadd CRAB outputs on EOS:

    haddCrabOutputs.py <datasetList.csv> <crab_job_label> <samples>

	haddCrabOutputs.py -h to show other optional arguments


* mvaNtuples:

Make flat mva ntuple from event ntuple: 

     produceMVANtuplesv2.py <crab_job_label> <list of samples> --datasetlist <datasetlist.csv> 

	 produceMVANtuplesv2.py -h to show other optional arguments

This Python script is a wrapper of binary 'makeMVANtuple' in /bin that does the actual production of the mva ntuples.
The script picks the ntuple files on EOS produced by CRAB jobs, given the <crab job label> (e.g. 2019jan) and the samples names (e.g. TTW).

By default, it runs all analysis categories, including 1l2tau, 2lss1tau, 3l1tau, 2l2tau, control\_ttW and control\_ttZ. The optional argument

    --analysis <list_of_analysis_categories>

can be used to specify a subset of the analysis categories.

For each category, by default event selections "signal" and "control" are applied to the MC samples, "signal", "control", "application\_fake" and "control\_fakeAR" are applied to the collision data. The optional argument

    --selection <selection_type>

can be used to specify a particular selection type. E.g. --selection application_fake on TTbar samples for closure studies.

The selection types are grouped in two categories: "datacard" and "control". For MC samples, "datacard" category includes only "signal" selection, "control" includes only "control" selections; For collision data samples, "datacard" category includes "signal" and "application_fake" selections, "control" includes "control" and "control_fakeAR" selcetions.
The optional argument

    --jobtypes <jtype>

can be used to pick either or both of the categories. By default, both categories are included. NOTE: this option is overwritten by --selection <selection_type> if <selection_type> is provided (i.e. not None).

To use LPC Condor to produce mva ntuples:

    prepareCondorJobs.py <crab_job_label> <list of samples> --datasetlist <datasetlist.csv> --fname\_jdl <condor.jdl> --fname\_sh <makentuple.sh>

This script has almost the same arguments as produceMVANtuplesv2.py. It generates one condor config file <condor.jdl>, and a bash script <makentuple.sh> that sets up the CMSSW work environment on the worker node and then runs produceMVANtuplesv2.py

* Datacards:

Make datacards using flat mvaNtuples:

    makeHistograms.py <selection_type> <variable_name> <xmin> <xmax> --channels <list_of_channels> --luminosity <lumi> --version <mvaNtuple_label>

    makeHistograms.py -h to show other optional arguments

To make datacards <selection_type> would be signal_<analysis_category>.
And <variable_name> would be the "mvaOutput" in the mva ntuple, with xmin = 0 and xmax = 1.

An example of running the above scripts to make datacards can be found in /test/makeDatacards.sh
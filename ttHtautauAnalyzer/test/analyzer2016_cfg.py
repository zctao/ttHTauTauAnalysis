### CMSSW_8_0_26_patch1

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ttH")

### initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### Run in unscheduled mode
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True)

###
options = VarParsing.VarParsing('analysis')

options.register('Debug', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Debug mode")

options.register('isData', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on collision data or not")

options.register('SampleName','',  #'ttH'
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Sample name")

options.register('TauESType', 'NA',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Tau energy scale: NA, tauESUp, tauESDown")

options.register('JECType', 'NA',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "JEC type: NA, JESUp, JESDown, JERUp, JERDown")

options.register("doJERSmearing", False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "apply jet energy smearing for MC or not")

options.register('TurnOffEvtSel', False, 
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Turn off event selection")

options.register('TurnOffHLTCut', True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Turn off HLT path check in event selection")

options.register('AnalysisType', '2lss1tau',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Analysis type currently supported: 1l2tau, 2lss1tau, 3l1tau")

options.register('SelectionRegion', 'signal_2lss1tau',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Which selection region to apply: signal_2lss1tau, control_2los1tau, control_fake_2lss1tau, signal_1l2tau, control_fake_1l2tau, signal_3l1tau, control_fake_3l1au, loose_2lss1tau, loose_1l2tau")

options.register('doSystematics', True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Include systematics or not")

options.register('doCutFlow', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Fill Cutflow histogram or not")

options.register('doSync', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "For Synchronization")

options.register('Is2016H', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "If the processed dataset is 2016H or not")

options.register("HIPSafeMediumMuon", True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "A switch for normal or HIP safe medium muon definions")

###
options.maxEvents = -1
options.inputFiles='file:/uscms/home/ztao/nobackup/datasample/ttH_80X/ttHnonbb.root'
#options.inputFiles='/store/mc/RunIISummer16MiniAODv2/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00D10AF2-76BE-E611-8EFB-001E67457DFA.root'
#options.inputFiles='/store/mc/RunIISummer16MiniAODv2/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/14466B11-62E3-E611-85A2-02163E011900.root'

# get and parse the command line arguments
options.parseArguments()

### Global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

MCGT = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
DATAGT = '80X_dataRun2_2016SeptRepro_v7'
if options.Is2016H:
    DATAGT = '80X_dataRun2_Prompt_v16'

process.GlobalTag.globaltag = DATAGT if options.isData else MCGT

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
)

### Filters for running on data
if options.isData:
    # primary vertex filter
    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                               vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                               minimumNDOF = cms.uint32(4) ,
                                               maxAbsZ = cms.double(24), 
                                               maxd0 = cms.double(2) 
    )

### Inputs
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
)

### JEC
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

JECLevel = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
if options.isData:
    JECLevel.append('L2L3Residual')

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', JECLevel, 'None')
)

### MET Correction and Uncertainty
from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams_Data
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# re-miniAOD dataset slimmedMETs already included mu correction. Only update JEC
if options.isData: 
    runMetCorAndUncFromMiniAOD(process, isData=options.isData,)
else:
    # re-correct MET based on bad muons on the fly for MC
    from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation

    muonRecoMitigation(
        process = process,
        pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
        runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
        selection="",
        muonCollection="", 
        cleanCollName="cleanMuonsPFCandidates",
        cleaningScheme="computeAllApplyClone", 
        postfix=""
    )

    process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
    process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

    runMetCorAndUncFromMiniAOD(process,
                               isData=False,
                               pfCandColl="cleanMuonsPFCandidates",
                               recoMetFromPFCs=True,
                               postfix=""#"MuClean"
    )

    process.mucorMET = cms.Sequence(
        process.badGlobalMuonTaggerMAOD *
        process.cloneGlobalMuonTaggerMAOD *
        #process.badMuons * # If you are using cleaning mode "all", uncomment this line
        process.cleanMuonsPFCandidates *
        process.fullPatMetSequence#MuClean
    )

### load the analysis
# electron MVA developed by the EGamma POG
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# LeptonID producer from ttH Multi-lepton group
process.load("ttH.LeptonID.ttHLeptons_cfi")
# Analyzer
process.load("ttHTauTauAnalysis.ttHtautauAnalyzer.ttHtaus_cfi")

### update parameter sets
process.ttHLeptons.rhoParam = "fixedGridRhoFastjetCentralNeutral"
process.ttHLeptons.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.ttHLeptons.LooseCSVWP = cms.double(0.5426)
process.ttHLeptons.MediumCSVWP = cms.double(0.8484)
process.ttHLeptons.IsHIPSafe = cms.bool(options.HIPSafeMediumMuon)

process.ttHtaus.electrons = cms.InputTag("ttHLeptons")
process.ttHtaus.muons = cms.InputTag("ttHLeptons")
process.ttHtaus.taus = cms.InputTag("ttHLeptons")
process.ttHtaus.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.ttHtaus.mets = cms.InputTag("slimmedMETs","","ttH")
process.ttHtaus.rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral")
process.ttHtaus.do_systematics = cms.bool(options.doSystematics)
process.ttHtaus.turn_off_event_sel = cms.bool(options.TurnOffEvtSel)
process.ttHtaus.sampleName = cms.string(options.SampleName)
process.ttHtaus.TauESType = cms.string(options.TauESType)
process.ttHtaus.JECType = cms.string(options.JECType)
process.ttHtaus.using_real_data = cms.bool(options.isData)
process.ttHtaus.analysis_type = cms.string(options.AnalysisType)
process.ttHtaus.selection_region = cms.string(options.SelectionRegion)
process.ttHtaus.turn_off_HLT_cut = cms.bool(options.TurnOffHLTCut)
process.ttHtaus.debug_mode = cms.bool(options.Debug)
process.ttHtaus.do_sync = cms.bool(options.doSync)
process.ttHtaus.doCutFlow = cms.bool(options.doCutFlow)
process.ttHtaus.doJERsmear = cms.bool(options.doJERSmearing)
# CSV WPs for Summer16 80X MC + ReReco Data
process.ttHtaus.csv_loose_wp = cms.double(0.5426)
process.ttHtaus.csv_medium_wp = cms.double(0.8484)
process.ttHtaus.csv_tight_wp = cms.double(0.9535)
if options.isData:
    print 'skip'
#    process.ttHtaus.MET_filters.append("Flag_badMuons")
#    process.ttHtaus.MET_filters.append("Flag_duplicateMuons") 
else:
    process.ttHtaus.badmu = cms.InputTag("badGlobalMuonTaggerMAOD","bad","ttH")
    process.ttHtaus.clonemu = cms.InputTag("cloneGlobalMuonTaggerMAOD","bad","ttH")

    ######## FIXME ########
# for 2016H data only
#if options.Is2016H:
    #process.ttHtaus.HLT_electron_muon_triggers = cms.vstring(
    #    ['HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v',
    #     'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v'])


### Output
out_file = 'output_' + options.SampleName + '.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(out_file)
)

### Performance
#process.Timing = cms.Service("Timing",
#                             summaryOnly = cms.untracked.bool(True),
#                             useJobReport = cms.untracked.bool(True)
#)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        oncePerEventMode=cms.untracked.bool(True)
#                                        )
#process.Tracer = cms.Service('Tracer')

#Path
if options.isData:
    process.p = cms.Path(
        process.primaryVertexFilter *
        process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
        process.fullPatMetSequence *
        process.electronMVAValueMapProducer *
        process.ttHLeptons *
        process.ttHtaus
    )
else:
    process.p = cms.Path(
        process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
        #process.fullPatMetSequence *
        process.mucorMET  *
        process.electronMVAValueMapProducer *
        process.ttHLeptons *
        process.ttHtaus
    )

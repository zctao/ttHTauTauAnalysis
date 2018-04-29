### CMSSW_9_4_4

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

options.register('AnalysisType', 'NA',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Analysis type currently supported: inclusive, 1l2tau, 2lss1tau, 3l1tau, 2l2tau")

options.register('SelectionRegion', 'NA',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Which selection region to apply: NA, inclusive_1l2tau, inclusive_2lss1tau, inclusive_3l1tau, inclusive_2l2tau")

#options.register('doSystematics', True,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.bool,
#                 "Include systematics or not")

options.register('doCutFlow', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Fill Cutflow histogram or not")

options.register('doSync', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "For Synchronization")

###
options.maxEvents = -1
options.inputFiles='file:/uscms/home/ztao/nobackup/datasample/ttH_94X/ttHJetToNonbb.root'
#/store/mc/RunIIFall17MiniAOD/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/0CF65340-0200-E811-ABB7-0025905C53F0.root

# get and parse the command line arguments
options.parseArguments()

### Global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

MCGT = '94X_mc2017_realistic_v13'
DATAGT = '94X_dataRun2_v6'

process.GlobalTag.globaltag = DATAGT if options.isData else MCGT

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
)

############################
### Filters for running on data
if options.isData:
    # primary vertex filter
    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                               vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                               minimumNDOF = cms.uint32(4) ,
                                               maxAbsZ = cms.double(24), 
                                               maxd0 = cms.double(2) 
    )
##############################

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

## Quark-gluon likelihood
# add the database object
qgDatabaseVersion = 'cmssw8020_v2'

from CondCore.CondDB.CondDB_cfi import *
#CondDB.connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000')
CondDB.connect = cms.string('sqlite:qg/QGL_'+qgDatabaseVersion+'.db')
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource", CondDB,
                                toGet = cms.VPSet(),
)

for type in ['AK4PFchs','AK4PFchs_antib']:
    process.QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
        record = cms.string('QGLikelihoodRcd'),
        tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
        label  = cms.untracked.string('QGL_'+type)
    )))

process.es_prefer_QGL = cms.ESPrefer("PoolDBESSource","QGPoolDBESSource")
  
# load QGTagger
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.QGTagger.srcVertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices")
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

### Electron MVA VID-based receipe
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

# ##################################################
# # rerun tau ID
from ttHTauTauAnalysis.ttHtautauAnalyzer.runTauIdMVA import *
na = TauIDEmbedder(process, cms, 
                   debug=True,
                   toKeep = ["dR0p32017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()

### MET correction 
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=options.isData,
                           )

### load the analysis
# LeptonID producer from ttH Multi-lepton group
process.load("ttH.LeptonID.ttHLeptons_cfi")
# Analyzer
process.load("ttHTauTauAnalysis.ttHtautauAnalyzer.ttHtaus_cfi")

### update parameter sets
process.ttHLeptons.rhoParam = "fixedGridRhoFastjetAll"
# if rerun tauID
process.ttHLeptons.taus = cms.InputTag("NewTauIDsEmbedded")
#
process.ttHLeptons.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.ttHLeptons.LooseCSVWP = cms.double(0.1522)  # DeepCSV WP
process.ttHLeptons.MediumCSVWP = cms.double(0.4941) # DeepCSV WP
process.ttHLeptons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values")
process.ttHLeptons.mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories")
#process.ttHLeptons.IsHIPSafe = cms.bool(options.HIPSafeMediumMuon)

process.ttHtaus.electrons = cms.InputTag("ttHLeptons")
process.ttHtaus.muons = cms.InputTag("ttHLeptons")
process.ttHtaus.taus = cms.InputTag("ttHLeptons")
process.ttHtaus.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
process.ttHtaus.mets = cms.InputTag("slimmedMETs","","ttH")
#process.ttHtaus.rho = cms.InputTag("fixedGridRhoFastjetAll")
#process.ttHtaus.do_systematics = cms.bool(options.doSystematics)
process.ttHtaus.turn_off_event_sel = cms.bool(options.TurnOffEvtSel)
process.ttHtaus.sample_name = cms.string(options.SampleName)
process.ttHtaus.TauESType = cms.string(options.TauESType)
process.ttHtaus.JECType = cms.string(options.JECType)
process.ttHtaus.using_collision_data = cms.bool(options.isData)
process.ttHtaus.analysis_type = cms.string(options.AnalysisType)
process.ttHtaus.selection_region = cms.string(options.SelectionRegion)
process.ttHtaus.turn_off_HLT_cut = cms.bool(options.TurnOffHLTCut)
process.ttHtaus.debug_mode = cms.bool(options.Debug)
process.ttHtaus.do_sync = cms.bool(options.doSync)
process.ttHtaus.doCutFlow = cms.bool(options.doCutFlow)
process.ttHtaus.doJERsmear = cms.bool(options.doJERSmearing)
process.ttHtaus.verbosity = cms.int32(1)
# DeepCSV WPs 
process.ttHtaus.csv_loose_wp = cms.double(0.1522)
process.ttHtaus.csv_medium_wp = cms.double(0.4941)
process.ttHtaus.csv_tight_wp = cms.double(0.8001)


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
        process.egmGsfElectronIDSequence *
        process.rerunMvaIsolationSequence * process.NewTauIDsEmbedded * # *getattr(process, "NewTauIDsEmbedded")
        process.fullPatMetSequence *
        process.ttHLeptons *
        process.QGTagger *
        process.ttHtaus
    )
else:
    process.p = cms.Path(
        process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
        process.egmGsfElectronIDSequence *
        process.rerunMvaIsolationSequence * process.NewTauIDsEmbedded * # *getattr(process, "NewTauIDsEmbedded")
        process.fullPatMetSequence *
        process.ttHLeptons *
        process.QGTagger *
        process.ttHtaus
    )

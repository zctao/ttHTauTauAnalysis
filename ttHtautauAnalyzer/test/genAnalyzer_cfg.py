import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ttH")

options = VarParsing.VarParsing('analysis')

options.register('pdgid', 25,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "pdgID of the target boson")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

options.maxEvents = -1
#options.inputFiles = 'file:/uscms/home/ztao/nobackup/datasample/ttH_80X/ttZ.root'
options.inputFiles = 'file:/uscms/home/ztao/nobackup/datasample/ttH_80X/ttHnonbb.root'
# ttHnonbb
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00D10AF2-76BE-E611-8EFB-001E67457DFA.root'
# ttZ
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/0603A111-42B5-E611-9369-0CC47A7EED28.root'

options.outputFile = 'GenParticles.root'

options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
)

process.GenAna = cms.EDAnalyzer('ttHGenAnalyzer',
                                #packed = cms.InputTag("packedGenParticles"),
                                pruned = cms.InputTag("prunedGenParticles"),
                                boson = cms.int32(options.pdgid)
)

process.p = cms.Path(process.GenAna)

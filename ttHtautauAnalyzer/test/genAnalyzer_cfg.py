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
# ttHnonbb
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00D10AF2-76BE-E611-8EFB-001E67457DFA.root'
#options.inputFiles = '/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_100897.root'
# ttZ
#'root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/0603A111-42B5-E611-9369-0CC47A7EED28.root'
#options.inputFiles = '/store/user/matze/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/faster_v9_ttZ_maod_54aa74f75231422e9f4d3766cb92a64a/ttZ_maod_13294.root'
options.inputFiles = ['/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_102533.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_104121.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_105619.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_107191.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_110805.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_112797.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_115969.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_117451.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_119576.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_120913.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_122224.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_123727.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_124419.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_126292.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_128763.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_132163.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_133992.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_135717.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_139897.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_140042.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_142815.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_143601.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_147490.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_151458.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_154249.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_157031.root',
'/store/user/matze/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/faster_v8_ttH_maod_p1_3a2fa29ab1d54ae0995b28f27b405be9/ttH_maod_161398.root']


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
                                boson = cms.int32(options.pdgid),
                                tauptcut = cms.double(0.), # 20.
                                tauetacut = cms.double(1000000.) # 2.3
)

process.p = cms.Path(process.GenAna)

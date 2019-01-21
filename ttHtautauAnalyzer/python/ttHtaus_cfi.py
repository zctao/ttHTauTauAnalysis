# Configurations for ttH tautau analysis
#
# Zhengcheng Tao
#

import FWCore.ParameterSet.Config as cms

ttHtaus = cms.EDAnalyzer('ttHtautauAnalyzer',
                         # loose selection
                         looseSelection = cms.bool(False),
                         # Generic parameters
                         sample_name = cms.string(''),
                         verbosity = cms.int32(0),
                         using_collision_data = cms.bool(False),
                         debug_mode = cms.bool(False),
                         do_sync = cms.bool(False),
                         turn_off_event_sel = cms.bool(False),
                         doCutFlow = cms.bool(True),
                         # triggers
                         print_HLT_event_path = cms.bool(False),
                         HLT_config_tag = cms.string('HLT'),
                         filter_config_tag = cms.string('PAT'),
                         # tauES
                         TauESUnc = cms.double(0.03),
                         # CSV WP
                         csv_loose_wp = cms.double(0.460),
                         csv_medium_wp = cms.double(0.800),
                         csv_tight_wp = cms.double(0.935),
                         # InputTags 
                         pv = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         pileup = cms.InputTag("slimmedAddPileupInfo"),
                         rho = cms.InputTag("fixedGridRhoFastjetAll"),
                         electrons = cms.InputTag("slimmedElectrons"),
                         muons = cms.InputTag("slimmedMuons"),
                         taus = cms.InputTag("slimmedTaus"),
                         jets = cms.InputTag("slimmedJets"),
                         mets = cms.InputTag("slimmedMETs"),
                         pfcand = cms.InputTag("packedPFCandidates"),
                         beamspot = cms.InputTag("offlineBeamSpot"),
                         packedgen = cms.InputTag("packedGenParticles"),
                         prunedgen = cms.InputTag("prunedGenParticles"),
                         badmu = cms.InputTag(""),
                         clonemu = cms.InputTag("")
)

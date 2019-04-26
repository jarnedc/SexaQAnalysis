import FWCore.ParameterSet.Config as cms

process = cms.Process("MCCand")

#process.load("Configuration.StandardSequences.FakeConditions_cff")
#process.load("Configuration.StandardSequences.GeometryIdeal_cff")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
#process.load("SimTracker.TrackHistory.GenTrackMatcher_cfi")
process.load("SexaQAnalysis.Analyzer_V0Fitter.producer_trackAssociatior_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///user/jdeclerc/CMSSW_8_0_30/src/STEP2_Sexaq/SUS-RunIISummer16DR80Premix-00068_1_374.root'),
    duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('standAloneAssociatorProducer.root')
)

process.genParticles.src = cms.InputTag("generator")


process.generalGenTrackMatcher.trackAssociator = cms.untracked.InputTag("trackAssociatorByHitsblqblq")

process.p = cms.Path(process.genTrackMatcher)
process.ep = cms.EndPath(process.out)
process.genParticles.saveBarCodes = True

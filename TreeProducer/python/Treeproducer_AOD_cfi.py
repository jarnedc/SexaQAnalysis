import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_AOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    lambdaKshortCollection = cms.InputTag("sMassFilter"),
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    isData = cms.untracked.bool(True),
    triggerName = cms.vstring('HLT_DiCentralPFJet170_v1'),
    genCollection =  cms.InputTag("genParticles","","HLT"),
)

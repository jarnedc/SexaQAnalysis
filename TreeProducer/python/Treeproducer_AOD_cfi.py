import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_AOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    prescales        = cms.InputTag("patTrigger", "", "RECO"), #not found
    prescales2        = cms.InputTag("patTrigger", "", "PAT"), #not found
    genjetCollection  = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    isData = cms.untracked.bool(True),
    triggerName = cms.vstring('HLT_DiCentralPFJet170_v1'),
    Partons_Source =  cms.InputTag("genParticles","","HLT"),
)

import FWCore.ParameterSet.Config as cms

analyzer = cms.EDAnalyzer('Analyzer',
    isData = cms.untracked.bool(True),
    genCollection =  cms.InputTag("genParticles","","HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    trackCollection    = cms.InputTag("generalTracks", "", "RECO"),    
)

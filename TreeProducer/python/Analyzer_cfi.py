import FWCore.ParameterSet.Config as cms

analyzer = cms.EDAnalyzer('Analyzer',
    isData = cms.untracked.bool(True),
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    resonCandidates = cms.InputTag("rMassFilter"),
    sexaqCandidates = cms.InputTag("sMassFilter"),
)

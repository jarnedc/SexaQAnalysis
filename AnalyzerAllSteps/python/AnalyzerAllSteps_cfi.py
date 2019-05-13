import FWCore.ParameterSet.Config as cms

AnalyzerAllSteps = cms.EDAnalyzer('AnalyzerAllSteps',
    isData = cms.untracked.bool(True),
    beamspot = cms.InputTag("offlineBeamSpot"),
    offlinePV = cms.InputTag("offlinePrimaryVertices","","RECO"),
    genCollection_GEN =  cms.InputTag("genParticles","","GEN"),
    genCollection_SIM_GEANT =  cms.InputTag("genParticlesPlusGEANT","","SIM"),
    generalTracksCollection =  cms.InputTag("generalTracks","","RECO"),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles","SEXAQ"),
    V0KsCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    V0LCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    trackAssociators = cms.InputTag("trackingParticleRecoTrackAsssociation"),
    TrackingParticles = cms.InputTag("mix","MergedTrackTruth")
)

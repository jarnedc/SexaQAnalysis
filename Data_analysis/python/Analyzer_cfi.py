import FWCore.ParameterSet.Config as cms

analyzer = cms.EDAnalyzer('Analyzer',
    isData = cms.untracked.bool(True),
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    resonCandidates = cms.InputTag("rMassFilter", "sVertexCompositePtrCandidate","SEXAQ"),
    sexaqCandidates = cms.InputTag("sMassFilter", "sVertexCompositePtrCandidate","SEXAQ"),
   
    #triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    #vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    #trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    #lambdaKshortCollection = cms.InputTag("sMassFilter","sVertexCompositePtrCandidate","SEXAQ"),
    #lambdaCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    #kshortCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    #sCollection = cms.InputTag("lambdaKshortVertexFilter","Sparticles","SEXAQ"),
    #sTrackCollection = cms.InputTag("lambdaKshortVertexFilter","sParticlesTracks","SEXAQ"),
    #triggerName = cms.vstring('HLT_DiCentralPFJet170_v1'),
    #genCollection =  cms.InputTag("genParticles","","HLT"),
)

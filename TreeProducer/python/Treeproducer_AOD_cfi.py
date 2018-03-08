import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_AOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    lambdaKshortCollection = cms.InputTag("sMassFilter","sVertexCompositePtrCandidate","SEXAQ"),
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    sCollection = cms.InputTag("lambdaKshortVertexFilter","Sparticles","SEXAQ"),
    sTrackCollection = cms.InputTag("lambdaKshortVertexFilter","sParticlesTracks","SEXAQ"),
    isData = cms.untracked.bool(True),
    triggerName = cms.vstring('HLT_DiCentralPFJet170_v1'),
    genCollection =  cms.InputTag("genParticles","","HLT"),
)

import FWCore.ParameterSet.Config as cms

Analyzer_V0_angular_correlation = cms.EDAnalyzer('Analyzer_V0_angular_correlation',
    isData = cms.untracked.bool(True),
    #beamspot = cms.InputTag("offlineBeamSpot"),
    #vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    #resonCandidates = cms.InputTag("rMassFilter", "sVertexCompositePtrCandidate","SEXAQ"),
    #sexaqCandidates = cms.InputTag("sMassFilter", "sVertexCompositePtrCandidate","SEXAQ"),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles","SEXAQ"),
   
    #triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    #trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    #lambdaKshortCollection = cms.InputTag("sMassFilter","sVertexCompositePtrCandidate","SEXAQ"),
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda","RECO"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort","RECO"),
    kshortCollectionlambdaKshortFilter = cms.InputTag("lambdaKshortFilter","kshort","SEXAQ"),
    lambdaCollectionlambdaKshortFilter = cms.InputTag("lambdaKshortFilter","lambda","SEXAQ"),
    sCollectionMassFilter = cms.InputTag("sMassFilter","sVertexCompositePtrCandidate","SEXAQ"),
    rCollectionMassFilter = cms.InputTag("rMassFilter","sVertexCompositePtrCandidate","SEXAQ"),
    
    nPVsCollection = cms.InputTag("InitialProducer","nPVs","SEXAQ"),
    nelectronsCollection = cms.InputTag("InitialProducer","nelectrons","SEXAQ"),
    njetsCollection = cms.InputTag("InitialProducer","njets","SEXAQ"),
    nkshortsCollection = cms.InputTag("InitialProducer","nkshorts","SEXAQ"),
    nlambdasCollection = cms.InputTag("InitialProducer","nlambdas","SEXAQ"),
    nmuonsCollection = cms.InputTag("InitialProducer","nmuons","SEXAQ"),
    ntracksCollection = cms.InputTag("InitialProducer","ntracks","SEXAQ"),
    HTCollection = cms.InputTag("InitialProducer","HT","SEXAQ"),
    TKHTCollection = cms.InputTag("InitialProducer","TKHT","SEXAQ"),
    TwoTopJetsCollection = cms.InputTag("InitialProducer","TwoTopJets","SEXAQ"),
    METCollection = cms.InputTag("InitialProducer","MET","SEXAQ"),
    TKMETCollection = cms.InputTag("InitialProducer","TKMET","SEXAQ"),
    

    #sCollection = cms.InputTag("lambdaKshortVertexFilter","Sparticles","SEXAQ"),
    #sTrackCollection = cms.InputTag("lambdaKshortVertexFilter","sParticlesTracks","SEXAQ"),
    #triggerName = cms.vstring('HLT_DiCentralPFJet170_v1'),
    #genCollection =  cms.InputTag("genParticles","","HLT"),
)

import FWCore.ParameterSet.Config as cms

OverlapFilter = cms.EDFilter(
    'OverlapFilter',
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    #genCollection    = cms.InputTag("genCollection"),
    genCollection = cms.InputTag("genParticles","","HLT"),
    isData = cms.bool(True),
    minNrLambda = cms.uint32(1),
    minNrKshort = cms.uint32(1),
    prescaleFalse = cms.uint32(0) # 0 means no prescale, reject all
)

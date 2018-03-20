import FWCore.ParameterSet.Config as cms

lambdaKshortVertexFilter = cms.EDFilter(
    'LambdaKshortVertexFilter',
    lambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    kshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    maxchi2ndofVertexFit = cms.double(5.),
    isData = cms.bool(True)
)

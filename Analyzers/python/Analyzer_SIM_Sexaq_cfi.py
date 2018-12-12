import FWCore.ParameterSet.Config as cms

Analyzer_SIM_Sexaq = cms.EDAnalyzer('Analyzer_SIM_Sexaq',
    isData = cms.untracked.bool(True),
    genCollection_GEN =  cms.InputTag("genParticles","","HLT"),
    genCollection_SIM_GEANT =  cms.InputTag("genParticlesPlusGEANT","","SIM"),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles","SEXAQ"),
)

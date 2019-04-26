import FWCore.ParameterSet.Config as cms

Analyzer_GEN_Sexaq = cms.EDAnalyzer('Analyzer_GEN_Sexaq',
    isData = cms.untracked.bool(True),

    genCollection_GEN = cms.InputTag("genParticles","","GEN"),
    genCollection_SIM = cms.InputTag("genParticles","","SIM"),
    genCollection_PLUSGEANT = cms.InputTag("genParticlesPlusGEANT","","SIM"),
    genCollection_HLT = cms.InputTag("genParticles","","HLT"),
    sexaqCandidates = cms.InputTag("lambdaKshortVertexFilter", "sParticles","SEXAQ")
)

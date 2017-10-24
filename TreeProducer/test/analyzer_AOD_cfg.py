import FWCore.ParameterSet.Config as cms

process = cms.Process("HEXAQTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'file:///user/gflouris/Analysis/SIMPS/SignalProduction/testPG/CMSSW_8_0_21/src/NoPile/SUS-RunIISummer16DR80Premix-00068.root'
    )
)


# Analyzer
process.load("HexaAnalysis.TreeProducer.Analyzer_cfi")
process.analyzer.isData = cms.untracked.bool(False)
process.p = cms.Path(process.analyzer)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('MC_analysis_tree.root')
)

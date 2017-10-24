import FWCore.ParameterSet.Config as cms


isData = False

process = cms.Process("HEXAQTREE")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

if(isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'file:///user/gflouris/Analysis/SIMPS/SignalProduction/testPG/CMSSW_8_0_21/src/NoPile/SUS-RunIISummer16DR80Premix-00068.root'
    )
)


# Adaptive vertex finder
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexFinder_cfi')

# Tree producer
process.load("HexaAnalysis.TreeProducer.Treeproducer_AOD_cfi")
process.tree.isData = cms.untracked.bool(isData)
process.p = cms.Path(process.inclusiveVertexFinder + process.tree)
#process.p = cms.Path(process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('MC_tree.root')
)



##Keep edm output file only during debugging
process.out = cms.OutputModule("PoolOutputModule",

    outputCommands = cms.untracked.vstring(
     #'drop *',
     'keep *'

),
   fileName = cms.untracked.string("outputfile_debug.root")
)

process.output_step = cms.EndPath(process.out)

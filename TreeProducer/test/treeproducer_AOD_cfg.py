import FWCore.ParameterSet.Config as cms


### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## data or MC options
options.register(
	'isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',-1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
	'flag to indicate max events to process')


process = cms.Process("HEXAQTREE")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

if(options.isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_2017Repro_v4', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v7', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     #'file://./E4C41A5A-1976-E711-B680-02163E01439D.root'
     #'file:///user/gflouris/Analysis/SIMPS/SignalProduction/testPG/CMSSW_8_0_21/src/NoPile/SUS-RunIISummer16DR80Premix-00068.root'
     'file://./0679350A-05A3-E711-A1F8-0025905B85CA.root'

    )
)


# Adaptive vertex finder
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexFinder_cfi')

# Tree producer
process.load("HexaAnalysis.TreeProducer.Treeproducer_AOD_cfi")
process.tree.isData = cms.untracked.bool(options.isData)
process.p = cms.Path(process.inclusiveVertexFinder + process.tree)
#process.p = cms.Path(process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('tree.root')
)

# #Keep edm output file only during debugging
# process.out = cms.OutputModule("PoolOutputModule",
#     outputCommands = cms.untracked.vstring(
#         #'drop *',
#         'keep *'
# ),
#    fileName = cms.untracked.string("outputfile_debug.root")
# )

# process.output_step = cms.EndPath(process.out)

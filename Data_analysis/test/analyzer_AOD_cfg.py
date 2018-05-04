import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## data or MC options
options.register(
	'isData',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',-1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
	'flag to indicate max events to process')


process = cms.Process("SEXAQDATAANA")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

if(options.isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(2000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#/pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_Run2016G/180321_194524/0000/
#/pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/

#single muon
#rangeOfStrings=[str(i) for i in range(1,148)] #1 to 148
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_Run2016G/180504_024915/0000/events_skimmed_' + x + '.root' for x in rangeOfStrings])
#zero bias
rangeOfStrings=[str(i) for i in range(1,500)] #1 to 500
tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/ZeroBias/ZeroBias_280318/180504_025020/0000/events_skimmed_' + x + '.root' for x in rangeOfStrings])

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		*tupleOfFiles #the * is to unpack the tuple of filenames
#      "dccp:///pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_100.root"

#      'file:///user/lowette/SexaQ/CMSSW_8_0_30/src/SexaQAnalysis/TreeProducer/test/events_skimmed.root'

		
    )
)


# Analyzer
process.load("SexaQAnalysis.Data_analysis.Analyzer_cfi")
process.analyzer.isData = cms.untracked.bool(True)

process.p = cms.Path(
  process.analyzer
)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analysis_SingleMuon_Run2016G.root')
)

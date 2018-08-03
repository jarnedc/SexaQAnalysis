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
	
	
#options.isData==True


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

#Jarne

#SingleMuon
rangeOfStringsC=[str(i) for i in range(1,130)]
rangeOfStringsD=[str(i) for i in range(1,194)]
rangeOfStringsE=[str(i) for i in range(1,175)]
rangeOfStringsF=[str(i) for i in range(1,130)]
rangeOfStringsG=[str(i) for i in range(1,295)]

#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134251/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsC], ['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134343/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsD], ['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134430/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsE] )

#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134251/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsC, 'file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134343/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsD, 'file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134430/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsE, 'file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134500/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsF, 'file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134540/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsG] )

tupleOfFiles1=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134251/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsC] )
tupleOfFiles2=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134343/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsC] )
tupleOfFiles3=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_multicrab/180731_134430/0000/events_skimmed_' + x + '.root' for x in rangeOfStringsC] )

#zero bias
#rangeOfStrings=[str(i) for i in range(1,500)] #1 to 500
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/ZeroBias/ZeroBias_280318/180505_032510/0000/events_skimmed_' + x + '.root' for x in rangeOfStrings])

#MET
#rangeOfStrings=[str(i) for i in range(1,513)] #1 to 503
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/MET/MET_Run2016G/180505_234637/0000/events_skimmed_' + x + '.root' for x in rangeOfStrings])


#rangeOfStrings=[str(i) for i in range(1,145)] #1 to 145 
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_Run2016G/180629_143246/0000/MC_events_skimmed_' + x + '.root' for x in rangeOfStrings])

#Florian

#/pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_Run2016G/180505_032403/0000   range to 148 
#/pnfs/iihe/cms/store/user/jdeclerc/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Wjets_30032018/180505_035823/0000  range to 483 
#/pnfs/iihe/cms/store/user/jdeclerc/ZeroBias/ZeroBias_280318/180505_032510/0000   range to 500 
#/pnfs/iihe/cms/store/user/jdeclerc/MinBias_TuneCUETP8M1_13TeV-pythia8/MinBias_260318/180505_035725/0000    range up to 201





process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		*tupleOfFiles1 #the * is to unpack the tuple of filenames
#      "dccp:///pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_100.root"

#      'file:///user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/MC_events_skimmed.root',
#      'file:///user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/V0_correlation_ZeroBias_single_file.root',

		
    )
)


# Analyzer
process.load("SexaQAnalysis.V0_angular_correlation_analysis.Analyzer_cfi")
process.Analyzer_V0_angular_correlation.isData = cms.untracked.bool(True) ##############SET BACK TO TRUE################

process.p = cms.Path(
  process.Analyzer_V0_angular_correlation
)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analysed_SingleMuon_LambdaKshort_correlation_allFilters_2016data_multicrab_runC.root')
)

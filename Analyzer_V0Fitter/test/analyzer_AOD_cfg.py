import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')



## data or MC options
options.register(
	'isData',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',20,VarParsing.multiplicity.singleton,VarParsing.varType.int,
	'flag to indicate max events to process')
	
	
#options.isData==True


process = cms.Process("SEXAQDATAANA")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')
#process.load('SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi')
#process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag

if(options.isData==True):
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#combined files of the RECO = Step2 with the smeared vertex (~130k antiS produced)
	#'file:///pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step2/Sexaquark_13TeV_trial5_Vtx_smeared/combined/crmc_Sexaq_combined_trial5_Vtx_Smeared.root',
	#'file:///pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step2/Sexaquark_13TeV_trial6_Vtx_smeared/combined/crmc_Sexaq_combined_trial6_Vtx_Smeared_Step2.root',
	#'file:///pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Step2/Sexaquark_13TeV_trial7_Vtx_smeared/combined/crmc_Sexaq_combined_trial7_Vtx_Smeared_Step2.root',




	#test file for the FEVT with the tracking particles:
	'file:///user/jdeclerc/CMSSW_8_0_30/src/STEP2_Sexaq/SUS-RunIISummer16DR80Premix-00068_1_374.root'
	#the skimmed files of teh run with the smeared vertex
	#'file:///pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/Sexaquark_13TeV_trial6_Vtx_smeared/combined/crmc_Sexaq_combined_trial6_Vtx_Smeared_skimmed.root'
	#centrally produced SingleMuon MC:
	#'file:///user/lowette/SexaQ/RunIISummer16DR80_MinBias_TuneCUETP8M1_13TeV-pythia8_GEN-SIM-RECO_NoPU_RECO_80X_mcRun2_asymptotic_v14-v1_100000_00150044-D075-E611-AAE8-001E67505A2D.root'
    ),
    duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)


# Analyzer
process.load("SexaQAnalysis.Analyzer_V0Fitter.Analyzer_V0Fitter_cfi") 
process.Analyzer_V0Fitter.isData = cms.untracked.bool(True) ##############SET BACK TO TRUE################

process.p = cms.Path(
  #process.simHitTPAssocProducer*process.trackAssociatorByHits*process.Analyzer_V0Fitter
  process.Analyzer_V0Fitter
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('analyzed_sexaq_Step2_Vtx_smeared_edm.root')
)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analyzed_sexaq_Step2_Vtx_smeared.root')
#    fileName = cms.string('analyzed_sexaq_Skimmed_Vtx_smeared.root')
#    fileName = cms.string('analyzed_Single_Muon_Step2.root')
)

process.ep = cms.EndPath(process.out)

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




process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/jdeclerc/crmc_Sexaq/Skimmed/Sexaquark_13TeV_trial4_neutron_mass_0/events_skimmed_MC_Sexaq_IIDD.root',
  			
    ),
    duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)


# Analyzer
#process.load("SexaQAnalysis.V0_angular_correlation_analysis.Analyzer_cfi")
#process.Analyzer_V0_angular_correlation.isData = cms.untracked.bool(True) ##############SET BACK TO TRUE################

#process.p = cms.Path(
#  process.Analyzer_V0_angular_correlation
#)

process.load("SexaQAnalysis.Analyzers.Analyzer_SIM_Sexaq_cfi")
process.Analyzer_SIM_Sexaq.isData = cms.untracked.bool(True) ##############SET BACK TO TRUE################

process.p = cms.Path(
  process.Analyzer_SIM_Sexaq
)

# Output
process.TFileService = cms.Service('TFileService',
#    fileName = cms.string('analyzed_events_skimmed_sexaq_crmc.root')
    fileName = cms.string('events_analyzed_MC_Sexaq_IIDD.root')
#    fileName = cms.string('MinBias_no_pile_up.root')
#    fileName = cms.string('WJets_with_pile_up.root')
#    fileName = cms.string('XiGun_0_displacement.root')
)

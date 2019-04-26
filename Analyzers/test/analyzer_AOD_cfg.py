import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')



## data or MC options
options.register(
	'isData',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC')

options.register(
	'maxEvts',10000,VarParsing.multiplicity.singleton,VarParsing.varType.int,
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
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound'))

#Jarne

#SingleMuon
rangeOfStrings=[str(i) for i in range(1,74)] 
tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/SingleMuon/SingleMuon_Run2016G-07Aug17-v1/181021_130627/0000/events_skimmed_2016_trialD_' + x + '.root' for x in rangeOfStrings])

#Min Bias MC (no Pile-up)
#rangeOfStrings=[str(i) for i in (1,2,3,4,5,6,7,8,9,10,11,13,14,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,53,54,55,56,57,58,59,61,64,65,66,67,68,69,70,72,73,74,75,76,77,83,84,85,86,89,90,91,92,93,95,96,97)] #starting from 14, because something is wrong with 13
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/MinBias_TuneCUETP8M1_13TeV-pythia8/MinBias_asData_2_010818/181019_181540/0000/Xi1820_0cm_displacement_' + x + '.root' for x in rangeOfStrings])

#W+jets (with pile-up)
#rangeOfStrings=[str(i) for i in range(1,483)] 
#tupleOfFiles=tuple(['file:///pnfs/iihe/cms/store/user/jdeclerc/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Wjets_01082018/181019_083349/0000/Xi1820_0cm_displacement_' + x + '.root' for x in rangeOfStrings])


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
	#	*tupleOfFiles #the * is to unpack the tuple of filenames
#      "dccp:///pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_100.root"      

	#low cross section (40mb), but no granddaughter cuts
	#'file:///user/jdeclerc/CMSSW_7_1_20_patch2/src/runGENSIM_Sexaq/sexaq_crmc_Sexaq_SIM_combined_xs40_withoutGranddaughterCuts.root'

	#high cross section (1200mb), but no granddaughter cuts
	#'file:///user/jdeclerc/CMSSW_7_1_20_patch2/src/runGENSIM_Sexaq/crmc_Sexaq_SIM_combined_xs1200_noGranddaughterCuts.root',	
	#'file:///user/jdeclerc/CMSSW_7_1_20_patch2/src/runGENSIM_Sexaq/crmc_Sexaq_SIM_combined_xs1200_noGranddaughterCuts001.root'	

	#high cross section (1200mb), with granddauhgter cuts
	#'file:///user/jdeclerc/CMSSW_7_1_20_patch2/src/runGENSIM_Sexaq_withCuts/sexaq_crmc_Sexaq_SIM_combined_xs1200_withGranddaughterCuts.root'

	#high cross section(1200mb), with the granddauhgter cuts (adapted the pt cut for the granddaughter to all be pt > 0.2) and added a deltaR(Ks,Lambda) cut
        'file:///user/jdeclerc/CMSSW_7_1_20_patch2/src/runGENSIM_Sexaq_withCuts/sexaq_crmc_Sexaq_SIM_combined_xs1200_withGranddaughterCuts.root'
    ),
    duplicateCheckMode = cms.untracked.string ("noDuplicateCheck")
)


# Analyzer
process.load("SexaQAnalysis.Analyzers.Analyzer_SIM_Sexaq_cfi") 
process.Analyzer_SIM_Sexaq.isData = cms.untracked.bool(True) ##############SET BACK TO TRUE################

process.p = cms.Path(
  process.Analyzer_SIM_Sexaq
)


# Output
process.TFileService = cms.Service('TFileService',
#    fileName = cms.string('analyzed_Xi1820.root')
#    fileName = cms.string('analyzed_sexaq_NO_PU.root')
#    fileName = cms.string('analyzed_sexaq_NO_PU_neutron_mass_0.root')
#    fileName = cms.string('analyzed_sexaq_GENSIM.root')
#    fileName = cms.string('analyzed_sexaq_GENSIM_single_file_4.root')
#    fileName = cms.string('analyzed_sexaq_Step2.root')
    fileName = cms.string('analyzed_sexaq_crmc_Sexaq_SIM_combined_xs1200_withGranddaughterCutsAndDeltaRCut.root')
#     fileName = cms.string('analyzed_sexaq_Skimmed_Vtx_smeared_trial_7.root')
#    fileName = cms.string('analyzed_sexaq_Step2_Vtx_smeared_skimmed.root')
#    fileName = cms.string('analyzed_sexaq_Step2_Vtx_smeared_status1_skimmed.root')
#    fileName = cms.string('analyzed_sexaq_minBias.root')
#    fileName = cms.string('analyzed_sexaq_Step2_file17.root')
)

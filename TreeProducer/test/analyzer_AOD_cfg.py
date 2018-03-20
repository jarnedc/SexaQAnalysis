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


process = cms.Process("SEXAQANA")

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
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      "dccp:///pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_100.root"
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_100.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_101.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_103.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_104.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_105.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_106.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_107.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_108.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_109.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_10.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_110.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_111.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_112.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_113.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_114.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_115.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_116.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_117.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_118.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_119.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_11.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_120.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_121.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_122.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_123.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_124.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_125.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_126.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_127.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_128.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_129.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_12.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_130.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_131.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_132.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_133.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_134.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_135.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_136.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_137.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_138.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_139.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_13.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_140.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_141.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_143.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_145.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_146.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_147.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_148.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_14.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_15.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_16.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_17.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_18.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_19.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_1.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_20.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_21.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_22.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_23.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_24.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_25.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_26.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_27.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_28.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_29.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_2.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_30.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_31.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_32.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_33.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_34.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_35.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_36.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_37.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_38.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_39.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_3.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_40.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_41.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_42.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_43.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_44.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_45.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_46.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_47.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_48.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_49.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_4.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_50.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_51.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_52.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_53.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_54.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_55.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_56.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_57.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_58.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_59.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_5.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_60.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_61.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_62.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_63.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_64.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_65.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_66.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_67.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_68.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_69.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_6.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_70.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_71.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_72.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_73.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_74.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_75.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_76.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_77.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_78.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_79.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_7.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_80.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_81.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_82.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_83.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_84.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_85.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_86.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_87.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_88.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_89.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_8.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_90.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_91.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_92.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_93.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_94.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_95.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_96.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_97.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_98.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_99.root",
      "dcap://maite.iihe.ac.be//pnfs/iihe/cms/store/user/lowette/SingleMuon/SingleMuon_Run2016G/180205_152747/0000/events_skimmed_9.root",
#      'file:///user/lowette/SexaQ/CMSSW_8_0_30/src/SexaQAnalysis/TreeProducer/test/events_skimmed.root'
    )
)


# Analyzer
process.load("SexaQAnalysis.TreeProducer.Analyzer_cfi")
process.analyzer.isData = cms.untracked.bool(True)

process.p = cms.Path(
  process.analyzer
)


# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analysis.root')
)

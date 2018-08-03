from CRABClient.UserUtilities import config

config = config()

pyCfgParams = ['isData=True']

config.General.requestName = 'SexaQ_MinBias_asData'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_data_cfg.py'

config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISummer16DR80-NoPU_RECO_80X_mcRun2_asymptotic_v14-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
config.Data.lumiMask = ''
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/jdeclerc/'
config.Data.publication = False
config.Data.outputDatasetTag = 'MinBias_asData_2_010818'

config.Site.storageSite = 'T2_BE_IIHE'
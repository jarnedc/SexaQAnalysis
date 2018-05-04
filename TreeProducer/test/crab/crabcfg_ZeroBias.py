from CRABClient.UserUtilities import config

config = config()

pyCfgParams = ['isData=True']

config.General.requestName = 'SexaQ_ZeroBias'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_data_cfg.py'

config.Data.inputDataset = '/ZeroBias/Run2016G-07Aug17-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/jdeclerc/'
config.Data.publication = False
config.Data.outputDatasetTag = 'ZeroBias_280318'

config.Site.storageSite = 'T2_BE_IIHE'

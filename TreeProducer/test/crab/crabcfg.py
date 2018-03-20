from CRABClient.UserUtilities import config

config = config()

pyCfgParams = ['isData=True']

config.General.requestName = 'SexaQ_MinBias'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_AOD_cfg.py'

config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISummer16DR80-NoPU_RECO_80X_mcRun2_asymptotic_v14-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
config.Data.lumiMask = '' #'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/lowette/'
config.Data.publication = False
config.Data.outputDatasetTag = 'MinBias'

config.Site.storageSite = 'T2_BE_IIHE'

from CRABClient.UserUtilities import config

config = config()

pyCfgParams = ['isData=True']

#config.General.requestName = 'SexaQ_SingleMuon'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_data_cfg.py'

#config.Data.inputDataset = '/SingleMuon/Run2016G-23Sep2016-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 500000
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = ''
config.Data.outLFNDirBase = '/store/user/jdeclerc/'
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuon_multicrab'
#config.Data.splitting = 'Automatic'

config.Site.storageSite = 'T2_BE_IIHE'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################


    config.General.requestName = 'SingleMuon_Multicrab_runB'
    config.Data.inputDataset = '/SingleMuon/Run2016B-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)


    config.General.requestName = 'SingleMuon_Multicrab_runC'
    config.Data.inputDataset = '/SingleMuon/Run2016C-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)


    config.General.requestName = 'SingleMuon_Multicrab_runD'
    config.Data.inputDataset = '/SingleMuon/Run2016D-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)


    config.General.requestName = 'SingleMuon_Multicrab_runE'
    config.Data.inputDataset = '/SingleMuon/Run2016E-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)

    
    config.General.requestName = 'SingleMuon_Multicrab_runF'
    config.Data.inputDataset = '/SingleMuon/Run2016F-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'SingleMuon_Multicrab_runG'
    config.Data.inputDataset = '/SingleMuon/Run2016G-23Sep2016-v1/AOD'
    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 4
    submit(config)

#    config.General.requestName = 'runC'
#    config.Data.inputDataset = '/DoubleMuParked/Run2012C-22Jan2013-v1/AOD'
#    config.Data.unitsPerJob = 500000
#    config.Data.totalUnits = 8
#    submit(config)

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'skimEtaMuMuGamma_pythiaGun'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'Run3_standardData.py'

#config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
config.Data.inputDataset = '/store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 5
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'skimEtaMuMuGamma_pythiaGun'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'

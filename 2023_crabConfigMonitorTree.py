from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'monitorSkim_2023C-v4'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'monitorTree.py'

config.Data.inputDataset = '/ScoutingPFMonitor/Run2023C-PromptReco-v4/MINIAOD'
config.Data.secondaryInputDataset  = '/ScoutingPFMonitor/Run2023C-v1/RAW'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'monitorSkim_2023C-v4'

# These values only make sense for processing data
#    Select input data based on a lumi mask
# Golden JSON Era C: https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_eraC_367095_368823_Golden.json
config.Data.lumiMask = 'Cert_Collisions2023_eraC_367095_368823_Golden.json'
# Golden JSON Era D: https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_eraD_369803_370790_Golden.json
#config.Data.lumiMask = 'Cert_Collisions2023_eraD_369803_370790_Golden.json'


# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'

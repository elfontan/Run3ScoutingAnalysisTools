from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'skimEtaMuMuGamma_pythiaGun_MINIAOD_3'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'Run3_standardData.py'

#config.Data.inputDataset = '/ScoutingPFRun3/Run2022F-v1/RAW'
config.Data.userInputFiles = [ 
            #'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700) 
            #'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_2/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700) 
            'root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(3500) 
        ]
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 100000 #10000 
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.ignoreLocality = True

# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'skimEtaMuMuGamma_pythiaGun_MINIAOD_3'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

config.Site.whitelist = ['T2_US*','T2_CH*']
# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'

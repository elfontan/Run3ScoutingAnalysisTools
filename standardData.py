import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing # ADDED                                                                                                          
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ADDED                                                                                                                                           
options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments           
options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('outFile',
                 'analysisEtaMuMuGamma.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file')

options.parseArguments()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source                                                                                                                                  
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(options.skipEvents), #added                                                          
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)
# end ADDED  

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(options.outFile)
)

#process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
#)
#
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#        "file:/eos/user/j/jfriesen/QCD_EtaMuMuGamma_TuneCP5_13p6TeV-pythia8_Run3Summer22EE_MiniAODv3/MiniAODv3_1.root",
#       "root://cmseos.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_987.root",
# )
#)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

# ---------------------------------
# NO for datasets in miniAOD format
# ---------------------------------
#process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
#process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

#process.gtStage2Digis.InputLabel = cms.InputTag( "hltStage2Digis" )
#process.gtStage2Digis.InputLabel = cms.InputTag( "AlgoBlkInputTag" )

#process.TFileService = cms.Service("TFileService", 
#    fileName = cms.string("analysisEtaMuMuGamma.root")
#)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun2_v2', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

#L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2"]
#L1Info = ["L1_SingleMu22"]
L1Info = ["L1_SingleMu22", "L1_DoubleMu_15_7", "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", "L1_DoubleMu4_SQ_OS_dR_Max1p2", "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"]

# Full list RUN 2:
# "L1_DoubleMu_12_5", "L1_DoubleMu_15_7", "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", "L1_DoubleMu4_SQ_OS_dR_Max1p2", "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", "L1_DoubleMu0_SQ", "L1_DoubleMu0er1p5_SQ_dR_Max1p4", "L1_DoubleMu0er2p0_SQ_dR_Max1p4", "L1_DoubleMu4p5_SQ_OS", "L1_TripleMu_5_3_3", "L1_TripleMu_5_5_3", "L1_QuadMu0", "L1_SingleMu22"

# RUN 3
#L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2"]

process.standardTree = cms.EDAnalyzer('gen_standardRecoTreeMakerRun3',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1              = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      #genP              = cms.untracked.InputTag("genParticles"),
                                      genP              = cms.untracked.InputTag("prunedGenParticles"),
                                      packedGenP        = cms.untracked.InputTag("packedGenParticles"),
                                      offlineMuons      = cms.untracked.InputTag("slimmedMuons"),
                                      offlinePhotons    = cms.untracked.InputTag("slimmedPhotons"),
                                      pfLabel           = cms.untracked.InputTag("packedPFCandidates"),
                                  )

process.p = cms.Path(process.standardTree)

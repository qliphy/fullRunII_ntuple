from WMCore.Configuration import Configuration
name = 'WWW'
steam_dir = 'xulyu'

config = Configuration()
config.section_("General")
config.General.requestName   = 'case1_3000'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Autumn18_V8_MC_L1FastJet_AK4PFchs.txt','Autumn18_V8_MC_L2Relative_AK4PFchs.txt','Autumn18_V8_MC_L3Absolute_AK4PFchs.txt','Autumn18_V8_MC_L1FastJet_AK8PFchs.txt','Autumn18_V8_MC_L2Relative_AK8PFchs.txt','Autumn18_V8_MC_L3Absolute_AK8PFchs.txt','Autumn18_V8_MC_L1FastJet_AK8PFPuppi.txt','Autumn18_V8_MC_L2Relative_AK8PFPuppi.txt','Autumn18_V8_MC_L3Absolute_AK8PFPuppi.txt','Autumn18_V8_MC_L1FastJet_AK4PFPuppi.txt','Autumn18_V8_MC_L2Relative_AK4PFPuppi.txt','Autumn18_V8_MC_L3Absolute_AK4PFPuppi.txt']
#config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/WWW-v2-3000/xulyu-crab_WWW-MA-3000-5-28028af67189b3de7224b79195bd0e1d/USER'
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =5
config.Data.totalUnits = -1
config.Data.publication = False
##config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/' + steam_dir + '/' + name + '/'
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'case1_3000'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

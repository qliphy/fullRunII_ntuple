from WMCore.Configuration import Configuration
name = 'WWW'
steam_dir = 'xulyu'

config = Configuration()
config.section_("General")
config.General.requestName   = 'multi_triboson'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Fall17_17Nov2017_V8_MC_L1FastJet_AK4PFchs.txt','Fall17_17Nov2017_V8_MC_L2Relative_AK4PFchs.txt','Fall17_17Nov2017_V8_MC_L3Absolute_AK4PFchs.txt','Fall17_17Nov2017_V8_MC_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017_V8_MC_L2Relative_AK8PFchs.txt','Fall17_17Nov2017_V8_MC_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017_V8_MC_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017_V8_MC_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017_V8_MC_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017_V8_MC_L1FastJet_AK4PFPuppi.txt','Fall17_17Nov2017_V8_MC_L2Relative_AK4PFPuppi.txt','Fall17_17Nov2017_V8_MC_L3Absolute_AK4PFPuppi.txt']
#config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis_multi.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.sendExternalFolder = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/multi_sim_5_22/xulyu-crab_multi_MA_5_24-05265a4e675f881e1e171fd06785811b/USER'
#config.Data.inputDBS = 'global'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =5
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/' + steam_dir + '/' + name + '/'
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'multi_sig'


config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

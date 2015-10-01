from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'LLGDV_Analysis_ZJetsNuNuHT400-600_noLeptonTrigger' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_11092015WithLeptons.py'
config.JobType.outputFiles = ['RecoOutput.root']

config.Data.inputDataset = '/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'LLGDV_Analysis_AnalysisNTuple_ZJetsNuNuHT400-600_Spring15_v03'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

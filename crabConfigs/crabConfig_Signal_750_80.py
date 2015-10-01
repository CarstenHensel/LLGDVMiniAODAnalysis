from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'LLGDV_Analysis_Signal_750_80' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_11092015WithLeptons.py'
config.JobType.outputFiles = ['RecoOutput.root']

config.Data.inputDataset = '/CRAB_PrivateMC/mhamer-LLG_750_80_Spring15_stepPAT-fb89f44b0d6970d718ed21d513cd1c9d/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'LLGDV_Analysis_AnalysisNTuple_Signal_750_80_Spring15_v03'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'LLGDV_Analysis_QCD_Pt120to170' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_27112015.py'
config.JobType.outputFiles = ['RecoOutput.root']
config.Data.inputDataset = '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.outputDatasetTag = 'LLGDV_Analysis_AnalysisNTuple_QCD_Pt120to170_Fall15'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

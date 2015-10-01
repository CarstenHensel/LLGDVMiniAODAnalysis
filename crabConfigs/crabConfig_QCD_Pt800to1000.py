from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'LLGDV_Analysis_QCD_Pt800to1000_withLeptonTrigger' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_11092015WithLeptons.py'
config.JobType.outputFiles = ['RecoOutput.root']
config.Data.inputDataset = '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25nsRecodebug_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'LLGDV_Analysis_AnalysisNTuple_QCD_Pt800to1000_Spring15_v03'
config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

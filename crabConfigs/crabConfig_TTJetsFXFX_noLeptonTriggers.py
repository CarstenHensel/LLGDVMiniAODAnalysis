from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'LLGDV_Analysis_TTBarFXFX_noLeptonTrigger' 
config.General.workArea = 'crabJobs'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'JobConfig_11092015.py'
config.JobType.outputFiles = ['RecoOutput.root']

config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/mhamer' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'LLGDV_Analysis_AnalysisNTuple_TTBarFXFX_Spring15_v03NLT'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'

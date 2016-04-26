import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('inputfile',
          '',
          VarParsing.multiplicity.list,
          VarParsing.varType.string,
          "Input File")

options.parseArguments()


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:./BE7E8874-DFB8-E511-B9DA-001E67DFF5D7.root',
    )
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process,dataFormat)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

for idmod in my_id_modules:
  setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)


process.demo = cms.EDAnalyzer('LLGDVMiniAODAnalysis',
    IsData = cms.bool(False),
    StoreGenData = cms.bool(False),
    RunLeptonTriggers = cms.bool(True),
    ignoreTriggers = cms.bool(True),
    bits = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects = cms.InputTag("selectedPatTrigger"),
    METFilters = cms.InputTag("TriggerResults", "","PAT"),
    pupInfo = cms.InputTag("slimmedAddPileupInfo"),
    conversions = cms.InputTag('allConversions'),
    GenEventInfo = cms.InputTag("generator"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secVertices = cms.InputTag("slimmedSecondaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    genJets = cms.InputTag("slimmedGenJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    pfCands = cms.InputTag("packedPFCandidates" ),
)


process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo )

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

process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag= 'MCRUN2_74_V9'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:./testfile.root',
    )
)

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),
   reverseDecision = cms.bool(False)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process,dataFormat)

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

for idmod in my_id_modules:
  setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)


process.demo = cms.EDAnalyzer('LLGDVMiniAODAnalysis',
    RunLeptonTriggers = cms.bool(False),
    bits = cms.InputTag("TriggerResults","","HLT"),
    METFilters = cms.InputTag("TriggerResults", "","PAT"),
    conversions = cms.InputTag('allConversions'),
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
    packed = cms.InputTag("packedGenParticles"),
    pruned = cms.InputTag("prunedGenParticles"),
    pfCands = cms.InputTag("packedPFCandidates" ),
)


process.p = cms.Path(process.HBHENoiseFilterResultProducer * process.ApplyBaselineHBHENoiseFilter * process.egmGsfElectronIDSequence * process.demo )

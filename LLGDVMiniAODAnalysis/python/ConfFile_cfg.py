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

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:./miniAOD-prod_PAT_1.root',
    )
)

process.demo = cms.EDAnalyzer('LLGDVMiniAODAnalysis',
    bits = cms.InputTag("TriggerResults","","HLT"),
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


process.p = cms.Path(process.demo)

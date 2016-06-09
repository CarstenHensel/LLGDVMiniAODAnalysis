// -*- C++ -*-
//
// Package:    LLGDVMiniAODAnalysis/LLGDVMiniAODAnalysis
// Class:      LLGDVMiniAODAnalysis
// 
/**\class LLGDVMiniAODAnalysis LLGDVMiniAODAnalysis.cc LLGDVMiniAODAnalysis/LLGDVMiniAODAnalysis/plugins/LLGDVMiniAODAnalysis.cc

 Description: 
 Simple analysis class to dump output in a plain rootfile

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Hamer
//         Created:  Tue, 13 Jan 2015 17:58:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TFile.h"
#include "TTree.h"

std::vector<double> CalculateVertex( std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> weight, std::vector<int> charge, std::vector<double> distance, int &nConsidered, double &weightednConsidered, std::vector<double> &error, double &maxScore );
//
//
// class declaration
//

class LLGDVMiniAODAnalysis : public edm::EDAnalyzer {
   public:
      explicit LLGDVMiniAODAnalysis(const edm::ParameterSet&);
      ~LLGDVMiniAODAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<reco::PFJetCollection> jetTokennoCHS_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secVtxToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<edm::TriggerResults> METFilterBits_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PupInfoToken_;
      //edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetToken electronToken_;
      edm::EDGetToken tauToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

      // the b-tagging algorithms for which we want the numbers
      std::vector<std::string> btagAlgorithms;


      // some flags:
      bool isData = false;
      bool storeGenData = false;
      bool ignoreTriggers = false;
      bool useCHSJets = true;

      // the output file and tree
      TFile *fOutput = new TFile("RecoOutput.root", "RECREATE");
      TTree *tOutput = new TTree("RecoData", "RecoData");
      TTree *tMetaData = new TTree("MetaData", "MetaData" );

      //setup the variables for the metadata tree first:
      int nEventsProcessed = 0;
      int nEventsAccepted = 0;


      
      // the variables used for output
      // pileup related variables;
      float NumberOfTrueInteractions;
      int NumberOfObservedInteractions;
      double GenLevel_HT;;

      // missing transverse energy
      std::vector<double> *met = new std::vector<double>;
      std::vector<double> *met_x = new std::vector<double>;
      std::vector<double> *met_y = new std::vector<double>;

      // jets passing tight jet id for MHT calculation
      std::vector<double> *tightJet_eta = new std::vector<double>;
      std::vector<double> *tightJet_phi = new std::vector<double>;
      std::vector<double> *tightJet_pt = new std::vector<double>;

      // the jet variables
      std::vector<double> *jet_eta = new std::vector<double>;
      std::vector<double> *jet_phi = new std::vector<double>;
      std::vector<std::vector<double> > *jet_pt = new std::vector<std::vector<double> >;
      std::vector<std::vector<double>* > *jet_btagInfo = new std::vector<std::vector<double>* >;
      std::vector<double> *jet_vertex_x = new std::vector<double>;
      std::vector<double> *jet_vertex_y = new std::vector<double>;
      std::vector<double> *jet_vertex_z = new std::vector<double>;
      std::vector<double> *jet_vertex_score = new std::vector<double>;
      std::vector<int> *jet_nCons = new std::vector<int>;
      std::vector<double> *jet_averageDistance = new std::vector<double>;
      std::vector<double> *jet_rmsDistance = new std::vector<double>;

      // the muon variables
      std::vector<double> *muon_px = new std::vector<double>;
      std::vector<double> *muon_py = new std::vector<double>;
      std::vector<double> *muon_pz = new std::vector<double>;
      std::vector<double> *muon_e = new std::vector<double>;
      std::vector<double> *muon_eta = new std::vector<double>;
      std::vector<double> *muon_phi = new std::vector<double>;
      std::vector<double> *muon_iso = new std::vector<double>;
      std::vector<double> *muon_charge = new std::vector<double>;
      std::vector<bool> *muon_isTightMuon = new std::vector<bool>;
      std::vector<bool> *muon_isMediumMuon = new std::vector<bool>;
      std::vector<bool> *muon_isLooseMuon = new std::vector<bool>;

      // the electron variables
      std::vector<double> *electron_px = new std::vector<double>;
      std::vector<double> *electron_py = new std::vector<double>;
      std::vector<double> *electron_pz = new std::vector<double>;
      std::vector<double> *electron_e = new std::vector<double>;
      std::vector<double> *electron_eta = new std::vector<double>;
      std::vector<double> *electron_phi = new std::vector<double>;
      std::vector<double> *electron_iso = new std::vector<double>;
      std::vector<double> *electron_charge = new std::vector<double>;
      std::vector<bool> *electron_isVeto = new std::vector<bool>;
      std::vector<bool> *electron_isLoose = new std::vector<bool>;
      std::vector<bool> *electron_isMedium = new std::vector<bool>;
      std::vector<bool> *electron_isTight = new std::vector<bool>;
      std::vector<bool> *electron_isHEEP = new std::vector<bool>;

      // the tau variables
      std::vector<double> *tau_px = new std::vector<double>;
      std::vector<double> *tau_py = new std::vector<double>;
      std::vector<double> *tau_pz = new std::vector<double>;
      std::vector<double> *tau_e = new std::vector<double>;


      // the trigger bits and names
      std::vector<int> *triggerBits = new std::vector<int>;
      std::vector<std::string> *triggerNames = new std::vector<std::string>;
      std::vector<std::string> *triggerNamesTree = new std::vector<std::string>;
      std::vector<int> *METFilterBits = new std::vector<int>;
      std::vector<std::string> *METFilterNames = new std::vector<std::string>;

      // the jet constituents 
      /*
      std::vector<std::vector<double> >* jet_constVertex_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pt      = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_eta     = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_phi     = new std::vector<std::vector<double> >;
      std::vector<std::vector<int> >*    jet_const_charge  = new std::vector<std::vector<int> >;
      std::vector<std::vector<double> >* jet_const_px = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_py = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pz = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pca0_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pca0_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pca0_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dxy = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dz = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_d = new std::vector<std::vector<double> >;
      */

      // the primary vertex information
      std::vector<double> *vertex_x = new std::vector<double>;
      std::vector<double> *vertex_y = new std::vector<double>;
      std::vector<double> *vertex_z = new std::vector<double>;
      std::vector<double> *vertex_nTracks = new std::vector<double>;
      std::vector<double> *vertex_pt = new std::vector<double>;
      std::vector<double> *vertex_ndof = new std::vector<double>;
      std::vector<double> *vertex_d0 = new std::vector<double>;
      std::vector<double> *vertex_dx = new std::vector<double>;
      std::vector<double> *vertex_dy = new std::vector<double>;
      std::vector<double> *vertex_dz = new std::vector<double>;


      // the secondary vertices
      std::vector<double> *secVertex_x = new std::vector<double>;
      std::vector<double> *secVertex_y = new std::vector<double>;
      std::vector<double> *secVertex_z = new std::vector<double>;
      std::vector<double> *secVertex_pt = new std::vector<double>;
      std::vector<double> *secVertex_ndof = new std::vector<double>;
      std::vector<double> *secVertex_chi2 = new std::vector<double>;
      std::vector<double> *secVertex_dx = new std::vector<double>;
      std::vector<double> *secVertex_dy = new std::vector<double>;
      std::vector<double> *secVertex_dz = new std::vector<double>;


      // trigger objects
      std::vector<std::string> *to_TriggerNames = new std::vector<std::string>;
      std::vector<std::vector<double> > *to_pt = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> > *to_eta = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> > *to_phi = new std::vector<std::vector<double> >;


      // MC truth information
      std::vector<double> *mct_px = new std::vector<double>;
      std::vector<double> *mct_py = new std::vector<double>;
      std::vector<double> *mct_pz = new std::vector<double>;
      std::vector<double> *mct_E = new std::vector<double>;
      std::vector<int> *mct_id = new std::vector<int>;
      std::vector<int> *mct_status = new std::vector<int>;
      std::vector<std::vector<int> > *mct_parentId = new std::vector<std::vector<int> >;
      std::vector<std::vector<int> > *mct_parentStatus = new std::vector<std::vector<int> >;
      std::vector<std::vector<double> > *mct_parentE = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> > *mct_parentPx = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> > *mct_parentPy = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> > *mct_parentPz = new std::vector<std::vector<double> >;


      // event metadata
      unsigned int RunNumber = 0;
      unsigned long long EventNumber = 0;
      unsigned int LuminosityBlock = 0;

      double sumOfWeights = 0.;
      double generatorWeight = 0.;
      std::vector<double> *generatorWeights = new std::vector<double>; 
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
LLGDVMiniAODAnalysis::LLGDVMiniAODAnalysis(const edm::ParameterSet& iConfig):
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetTokennoCHS_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsnochs"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  secVtxToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjects"))),
  METFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METFilters"))),
  genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  genEvtInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenEventInfo") )),
  PupInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pupInfo"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
  conversionsToken_(consumes<reco::ConversionCollection > (iConfig.getParameter<edm::InputTag>("conversions"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
  {
   //now do what ever initialization is needed
   
   // define the triggers we (might) want to use
   // these seem to be interesting for us
   /* triggerNames->push_back( "HLT_PFJet260_v1" );
    * triggerNames->push_back( "HLT_JetE30_NoBPTX_v1" );
    * triggerNames->push_back( "HLT_JetE30_NoBPTX3BX_NoHalo_v1" );
    * triggerNames->push_back( "HLT_JetE50_NoBPTX3BX_NoHalo_v1" );
    * triggerNames->push_back( "HLT_JetE70_NoBPTX3BX_NoHalo_v1" );
   */
   bool RunLeptonTriggers = iConfig.getParameter<bool>("RunLeptonTriggers");
   isData = iConfig.getParameter<bool>("IsData");
   storeGenData = iConfig.getParameter<bool>("StoreGenData");
   ignoreTriggers = iConfig.getParameter<bool>("IgnoreTriggers"); 
   useCHSJets = iConfig.getParameter<bool>("UseCHSJets");

   std::string trigVersion = (isData) ? "_v2" : "_v1";

   triggerNames->push_back("HLT_PFMET170_NoiseCleaned" );
   triggerNames->push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight" );
   triggerNames->push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight" );
   triggerNames->push_back("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight" );
   triggerNames->push_back("HLT_PFMET90_JetIdCleaned_PFMHT90_IDTight" );
   triggerNames->push_back("HLT_PFMET90_NoiseCleaned_PFMHT90_IDTight" );
   triggerNames->push_back("HLT_PFMET90_PFMHT90_IDTight" );

   if( RunLeptonTriggers ) {
    triggerNames->push_back("HLT_Ele22_eta2p1_WP75_Gsf");
    triggerNames->push_back("HLT_Ele27_WP85_Gsf");
    triggerNames->push_back("HLT_Ele27_eta2p1_WP75_Gsf_LooseIsoPFTau20");
    triggerNames->push_back("HLT_Ele27_eta2p1_WP75_Gsf");
    triggerNames->push_back("HLT_Ele32_eta2p1_WP75_Gsf");
    triggerNames->push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT");
    triggerNames->push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL");
    triggerNames->push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL");
    triggerNames->push_back("HLT_Mu50");
    triggerNames->push_back("HLT_Mu45_eta2p1");
   
   }
   //for( unsigned int itrig = 0; itrig < triggerNames->size(); ++itrig ) {
   //   triggerNames->at(itrig) += trigVersion;
   //}
   // the btag algorithms for which we want infos
   btagAlgorithms.push_back("pfJetBProbabilityBJetTags");
   btagAlgorithms.push_back("pfJetProbabilityBJetTags");
   btagAlgorithms.push_back("pfTrackCountingHighPurBJetTags");
   btagAlgorithms.push_back("pfTrackCountingHighEffBJetTags");
   btagAlgorithms.push_back("pfCombinedInclusiveSecondaryVertexV2BJetTags");
   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
      jet_btagInfo->push_back( new std::vector<double> );
   }

   tOutput -> Branch("TightJet_eta", &tightJet_eta );
   tOutput -> Branch("TightJet_phi", &tightJet_phi );
   tOutput -> Branch("TightJet_pt", &tightJet_pt );

   // set the output branches for the tree
   tOutput -> Branch("RecoJet_eta", &jet_eta );
   tOutput -> Branch("RecoJet_phi", &jet_phi );
   tOutput -> Branch("RecoJet_pt", &jet_pt );
   tOutput -> Branch("RecoJet_Vertex_x", &jet_vertex_x );
   tOutput -> Branch("RecoJet_Vertex_y", &jet_vertex_y );
   tOutput -> Branch("RecoJet_Vertex_z", &jet_vertex_z );
   tOutput -> Branch("RecoJet_Vertex_Score", &jet_vertex_score );
   tOutput -> Branch("RecoJet_nConsidered", &jet_nCons );
   tOutput -> Branch("RecoJet_AverageDistanceToVertex", &jet_averageDistance );
   tOutput -> Branch("RecoJet_RMSDistanceToVertex", &jet_rmsDistance );
   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
     std::string branchname = "RecoJet_btag_" + btagAlgorithms.at(iBtagAlgo);
     tOutput -> Branch( branchname.c_str(), &(jet_btagInfo->at(iBtagAlgo)) ); 
   }
   /*
   tOutput -> Branch("RecoJet_constVertex_x", &jet_constVertex_x );
   tOutput -> Branch("RecoJet_constVertex_y", &jet_constVertex_y );
   tOutput -> Branch("RecoJet_constVertex_z", &jet_constVertex_z );
   tOutput -> Branch("RecoJet_const_pt", &jet_const_pt );
   tOutput -> Branch("RecoJet_const_charge", &jet_const_charge );
   tOutput -> Branch("RecoJet_const_px", &jet_const_px );
   tOutput -> Branch("RecoJet_const_py", &jet_const_py );
   tOutput -> Branch("RecoJet_const_pz", &jet_const_pz );
   tOutput -> Branch("RecoJet_const_pca0_x", &jet_const_pca0_x );
   tOutput -> Branch("RecoJet_const_pca0_y", &jet_const_pca0_y );
   tOutput -> Branch("RecoJet_const_pca0_z", &jet_const_pca0_z );
   tOutput -> Branch("RecoJet_const_closestVertex_dxy", &jet_const_closestVertex_dxy );
   tOutput -> Branch("RecoJet_const_closestVertex_dz", &jet_const_closestVertex_dz );
   tOutput -> Branch("RecoJet_const_closestVertex_d", &jet_const_closestVertex_d );
   tOutput -> Branch("RecoJet_const_eta", &jet_const_eta );
   tOutput -> Branch("RecoJet_const_phi", &jet_const_phi );
   */ 
   tOutput -> Branch("PUINFO_NumberOfTrueInteractions", &NumberOfTrueInteractions );
   tOutput -> Branch("PUINFO_NumberOfObservedInteractiosn", &NumberOfObservedInteractions );
   tOutput -> Branch("GenLevel_HT", &GenLevel_HT );
   tOutput -> Branch("RecoMuon_px", &muon_px );
   tOutput -> Branch("RecoMuon_py", &muon_pz );
   tOutput -> Branch("RecoMuon_pz", &muon_py );
   tOutput -> Branch("RecoMuon_E", &muon_e );
   tOutput -> Branch("RecoMuon_eta", &muon_eta );
   tOutput -> Branch("RecoMuon_phi", &muon_phi );
   tOutput -> Branch("RecoMuon_iso", &muon_iso);
   tOutput -> Branch("RecoMuon_charge", &muon_charge);
   tOutput -> Branch("RecoMuon_isTightMuon", &muon_isTightMuon );
   tOutput -> Branch("RecoMuon_isMediumMuon", &muon_isMediumMuon );
   tOutput -> Branch("RecoMuon_isLooseMuon", &muon_isLooseMuon );
   tOutput -> Branch("RecoElectron_px", &electron_px );
   tOutput -> Branch("RecoElectron_py", &electron_pz );
   tOutput -> Branch("RecoElectron_pz", &electron_py );
   tOutput -> Branch("RecoElectron_E", &electron_e );
   tOutput -> Branch("RecoElectron_eta", &electron_eta );
   tOutput -> Branch("RecoElectron_phi", &electron_phi );
   tOutput -> Branch("RecoElectron_iso", &electron_iso );
   tOutput -> Branch("RecoElectron_charge", &electron_charge );
   tOutput -> Branch("RecoElectron_isVeto", &electron_isVeto );
   tOutput -> Branch("RecoElectron_isLoose", &electron_isLoose );
   tOutput -> Branch("RecoElectron_isMedium", &electron_isMedium );
   tOutput -> Branch("RecoElectron_isTight", &electron_isTight );
   tOutput -> Branch("RecoElectron_isHEEP", &electron_isHEEP );
   tOutput -> Branch("RecoTau_px", &tau_px );
   tOutput -> Branch("RecoTau_py", &tau_py );
   tOutput -> Branch("RecoTau_pz", &tau_pz );
   tOutput -> Branch("RecoTau_e", &tau_e );
   tOutput -> Branch("TriggerBits", &triggerBits );
   tOutput -> Branch("TriggerNames", &triggerNamesTree );
   tOutput -> Branch("METFilterBits", &METFilterBits );
   tOutput -> Branch("METFilterNames", &METFilterNames );
   
   tOutput -> Branch("RecoVertex_x", &vertex_x );
   tOutput -> Branch("RecoVertex_y", &vertex_y );
   tOutput -> Branch("RecoVertex_z", &vertex_z );
   tOutput -> Branch("RecoVertex_ndof", &vertex_ndof ); 
   tOutput -> Branch("RecoVertex_d0", &vertex_d0 );
   tOutput -> Branch("RecoVertex_xError", &vertex_dx );
   tOutput -> Branch("RecoVertex_yError", &vertex_dy );
   tOutput -> Branch("RecoVertex_zError", &vertex_dz );
   tOutput -> Branch("RecoVertex_nTracks", &vertex_nTracks );
   tOutput -> Branch("RecoVertex_pt", &vertex_pt );
   
   tOutput -> Branch("RecoSecVertex_x", &secVertex_x );
   tOutput -> Branch("RecoSecVertex_y", &secVertex_y );
   tOutput -> Branch("RecoSecVertex_z", &secVertex_z );
   tOutput -> Branch("RecoSecVertex_ndof", &secVertex_ndof );
   tOutput -> Branch("RecoSecVertex_chi2", &secVertex_chi2 );
   tOutput -> Branch("RecoSecVertex_pt", &secVertex_pt );
   tOutput -> Branch("RecoSecVertex_xError", &secVertex_dx );
   tOutput -> Branch("RecoSecVertex_yError", &secVertex_dy );
   tOutput -> Branch("RecoSecVertex_zError", &secVertex_dz );
 
   tOutput -> Branch("TriggerObject_TriggerName", &to_TriggerNames );
   tOutput -> Branch("TriggerObject_pt", &to_pt );
   tOutput -> Branch("TriggerObject_eta", &to_eta );
   tOutput -> Branch("TriggerObject_phi", &to_phi );

   if( storeGenData && !isData ) {
     tOutput -> Branch("GenLevel_px", &mct_px );
     tOutput -> Branch("GenLevel_py", &mct_py );
     tOutput -> Branch("GenLevel_pz", &mct_pz );
     tOutput -> Branch("GenLevel_E", &mct_E );
     tOutput -> Branch("GenLevel_PDGID", &mct_id );
     tOutput -> Branch("GenLevel_status", &mct_status );
     tOutput -> Branch("GenLevel_ParentId", &mct_parentId );
     tOutput -> Branch("GenLevel_ParentPx", &mct_parentPx );
     tOutput -> Branch("GenLevel_ParentPy", &mct_parentPy );
     tOutput -> Branch("GenLevel_ParentPz", &mct_parentPz );
     tOutput -> Branch("GenLevel_ParentE", &mct_parentE );
     tOutput -> Branch("GenLevel_ParentStatus", &mct_parentStatus );
   }

   tOutput -> Branch("MET", &met );
   tOutput -> Branch("MET_x", &met_x );
   tOutput -> Branch("MET_y", &met_y );
   tOutput -> Branch("EventNumber", &EventNumber );
   tOutput -> Branch("RunNumber", &RunNumber );
   tOutput -> Branch("LuminosityBlock", &LuminosityBlock );
   tOutput -> Branch("GeneratorWeights", &generatorWeights );
   tOutput -> Branch("GeneratorWeight", &generatorWeight );

   // the metadata tree
   tMetaData -> Branch("nEventsProcessed", &nEventsProcessed );
   tMetaData -> Branch("SumOfWeights", &sumOfWeights );
   tMetaData -> Branch("nEventsAccepted", &nEventsAccepted );
}



LLGDVMiniAODAnalysis::~LLGDVMiniAODAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LLGDVMiniAODAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
   nEventsProcessed += 1;
   // clear all variables
   met->clear(); 
   met_x -> clear();
   met_y -> clear();
   muon_px->clear();
   muon_py->clear();
   muon_pz->clear();
   muon_e->clear();
   muon_eta->clear();
   muon_phi->clear();
   muon_iso->clear();
   muon_charge->clear();
   muon_isTightMuon->clear();
   muon_isMediumMuon->clear();
   muon_isLooseMuon->clear();
   electron_px->clear();
   electron_py->clear();
   electron_pz->clear();
   electron_e->clear();
   electron_eta->clear();
   electron_phi->clear();
   electron_iso->clear();
   electron_charge->clear();
   electron_isVeto->clear();
   electron_isLoose->clear(); 
   electron_isMedium->clear();
   electron_isTight->clear();
   electron_isHEEP->clear();
   tau_px->clear();
   tau_py->clear();
   tau_pz->clear();
   tau_e->clear();
   triggerBits->clear();
   triggerNamesTree->clear();
   METFilterBits->clear();
   tightJet_eta->clear();
   tightJet_phi->clear();
   tightJet_pt->clear();
   jet_eta->clear();
   jet_phi->clear();
   jet_pt->clear();
   jet_vertex_x->clear();
   jet_vertex_y->clear();
   jet_vertex_z->clear();
   jet_vertex_score->clear();
   jet_nCons->clear();
   jet_averageDistance->clear();
   jet_rmsDistance->clear();
   NumberOfTrueInteractions = -1.;
   NumberOfObservedInteractions = -1;
   GenLevel_HT = -1.;
   /*
   jet_constVertex_x->clear();
   jet_constVertex_y->clear();
   jet_constVertex_z->clear();
   jet_const_pt->clear();
   jet_const_eta->clear();
   jet_const_phi->clear();
   jet_const_charge->clear();
   jet_const_px->clear();
   jet_const_py->clear();
   jet_const_pz->clear();
   jet_const_pca0_x->clear();
   jet_const_pca0_y->clear();
   jet_const_pca0_z->clear();
   jet_const_closestVertex_dxy->clear();
   jet_const_closestVertex_dz->clear();
   jet_const_closestVertex_d->clear();
   */
   vertex_x -> clear();
   vertex_y -> clear();
   vertex_z -> clear();
   vertex_d0 -> clear();
   vertex_dx -> clear();
   vertex_dy -> clear();
   vertex_dz -> clear();
   vertex_nTracks -> clear();
   vertex_pt -> clear();
   vertex_ndof -> clear();
   secVertex_x -> clear();
   secVertex_y -> clear();
   secVertex_z -> clear();
   secVertex_ndof -> clear();
   secVertex_chi2 -> clear();
   secVertex_pt -> clear();
   secVertex_dx -> clear();
   secVertex_dy -> clear();
   secVertex_dz -> clear();
   to_TriggerNames -> clear();
   to_pt -> clear();
   to_eta -> clear();
   to_phi -> clear();
   generatorWeights -> clear();
   generatorWeight = 0.;
   if( storeGenData && !isData ) {
     mct_px -> clear();
     mct_py -> clear();
     mct_pz -> clear();
     mct_E -> clear();
     mct_id -> clear();
     mct_status -> clear();
     mct_parentId -> clear();
     mct_parentPx -> clear();
     mct_parentPy -> clear();
     mct_parentPz -> clear();
     mct_parentE -> clear();
     mct_parentStatus -> clear();
   }


   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
     jet_btagInfo->at(iBtagAlgo)->clear();
   }


   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

   // get the weight and PU info
   if( !isData ) {
     //generator weight:
     edm::Handle<GenEventInfoProduct> genEvtInfo; 
     iEvent.getByToken( genEvtInfoToken_, genEvtInfo );

     generatorWeight = genEvtInfo->weight();
     sumOfWeights += generatorWeight;
   
     const std::vector<double>& evtWeights = genEvtInfo->weights();
     for( unsigned int iWeight = 0; iWeight < evtWeights.size(); ++iWeight ) {
       generatorWeights->push_back( evtWeights.at(iWeight ) );
     }
  
     // pileup info:
     Handle<std::vector< PileupSummaryInfo > >  PupInfo;
     iEvent.getByToken(PupInfoToken_, PupInfo);

     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) { 
          NumberOfTrueInteractions = PVI->getTrueNumInteractions();
          NumberOfObservedInteractions = PVI->getPU_NumInteractions();
       }
     }
   }

   else {
     generatorWeight = 1.;
     sumOfWeights += 1.;
   }
   //edm::Handle<LHEEventProduct> lheEventProduct;
   //iEvent.getByToken( lheEventToken_, lheEventProduct );
   //std::cout << "compare " << generatorWeight << " with " << lheEventProduct->hepeup().XWGTUP << std::endl;


   // get all the physics objects from the miniAOD
   Handle<reco::PFJetCollection> jetsnoCHS;
   if( !useCHSJets ) {
     iEvent.getByToken( jetTokennoCHS_, jetsnoCHS );
   }

   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   
   Handle<pat::JetCollection> jets;
   iEvent.getByToken( jetToken_, jets );
   
   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken( vtxToken_, vertices );
   
   Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
   iEvent.getByToken( secVtxToken_, secVertices );

   Handle<pat::METCollection> mets;
   iEvent.getByToken( metToken_, mets );

   edm::Handle<edm::TriggerResults> evTriggerBits;
   iEvent.getByToken( triggerBits_, evTriggerBits );
  
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   edm::Handle<edm::TriggerResults> evMETFilterBits;
   iEvent.getByToken( METFilterBits_, evMETFilterBits );
  
   Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken( genJetToken_, genJets );
  
   Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   
   Handle<edm::View<reco::GenParticle> > pruned;
   Handle<edm::View<pat::PackedGenParticle> > packed;
   
   if( !isData ) {
     iEvent.getByToken(prunedGenToken_, pruned);
     if( storeGenData ) {
       iEvent.getByToken(packedGenToken_, packed);
     }
   }


   edm::EventAuxiliary aux = iEvent.eventAuxiliary();
   edm::EventID id = aux.id();
  
   EventNumber = id.event();
   RunNumber = id.run();
   LuminosityBlock = id.luminosityBlock();


   
   const edm::TriggerNames &filterNames = iEvent.triggerNames(*evMETFilterBits);
   bool passETMissFilter = true;
   for(unsigned int i = 0; i < evMETFilterBits->size(); ++i ) {
   // std::cout <<" testing met filter " << filterNames.triggerName(i) << std::endl;
    if(    filterNames.triggerName(i) == "Flag_HBHENoiseFilter" 
        || filterNames.triggerName(i) == "Flag_HBHENoiseIsoFilter"
        || filterNames.triggerName(i) == "Flag_CSCTightHalo2015Filter"
        || filterNames.triggerName(i) == "Flag_EcalDeadCellTriggerPrimitiveFilter"
        || filterNames.triggerName(i) == "Flag_goodVertices"
        || filterNames.triggerName(i) == "Flag_eeBadScFilter"
        || filterNames.triggerName(i) == "Flag_chargedHadronTrackResolutionFilter"
        || filterNames.triggerName(i) == "Flag_muonBadTrackFilter" ) {
          if( !evMETFilterBits->accept(i) ) passETMissFilter = false;    
        }
   } 
   if( !passETMissFilter && !ignoreTriggers ) return;

   const edm::TriggerNames &names = iEvent.triggerNames(*evTriggerBits);
   bool passTrigger = false;
   //for(unsigned int i = 0; i < evTriggerBits->size(); ++i ) {
   //  std::cout << "got trigger " << names.triggerName(i) << std::endl;
   //}
   for( unsigned int j = 0; j < triggerNames->size(); ++j ) {
    for(unsigned int i = 0; i < evTriggerBits->size(); ++i ) {
        //std::cout << "got trigger " << names.triggerName(i) << std::endl;
        if( names.triggerName(i).find( triggerNames->at(j) ) != std::string::npos ) {
          triggerNamesTree->push_back( names.triggerName(i) );
          triggerBits->push_back( ( (evTriggerBits->accept(i)) ? 1 : 0) );
          if( evTriggerBits->accept(i) ) passTrigger = true;
          continue;
        }
      }
   }
   // store only events in the ntuple which pass the trigger
   if( !passTrigger && !ignoreTriggers ) return;
  
   // now fill the primary vertex information
   int firstGoodVertexIdx = -1;
   int iVtx = 0;
   for( const reco::Vertex &v : *vertices ) {
      bool isFake = v.isFake();
      if( !isFake && v.ndof() >= 4. && v.position().Rho() <= 2.0 && fabs(v.position().z()) <= 24. ) {
        if( firstGoodVertexIdx == -1 ) firstGoodVertexIdx = iVtx;
        vertex_x -> push_back( v.x() );
        vertex_y -> push_back( v.y() );
        vertex_z -> push_back( v.z() );
        vertex_dx -> push_back( v.xError() );
        vertex_dy -> push_back( v.yError() );
        vertex_dz -> push_back( v.zError() );
        vertex_ndof -> push_back( v.ndof() );
        vertex_d0 -> push_back( v.position().rho() );
        vertex_nTracks -> push_back( v.nTracks() );
        vertex_pt -> push_back( v.p4().pt() );
      }
      iVtx ++;
   }

   // manual implementation of good vertex filter
   if( firstGoodVertexIdx == -1 && !ignoreTriggers ) return;

   // fill the trigger objects
   //std::string trigVersion = (isData) ? "_v2" : "_v1";
   //to_TriggerNames->push_back("HLT_Mu50" + trigVersion );
   to_TriggerNames->push_back("HLT_Mu50" );
  

   for( unsigned int k = 0; k < to_TriggerNames->size(); ++k ) {
      std::vector<double> pt, eta, phi;
      to_pt->push_back( pt );
      to_eta->push_back( eta );
      to_phi->push_back( phi );
   }

   for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
      obj.unpackPathNames(names);
      //std::cout << "\t   Collection: " << obj.collection() << std::endl;
      //std::cout << "\t   Type IDs:   ";
      //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
      //std::cout << std::endl;
      //std::cout << "\t   Filters:    ";
      //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
      //std::cout << std::endl;
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      std::vector<std::string> pathNamesLast = obj.pathNames(true);
      //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
         for( unsigned int k = 0; k < to_TriggerNames->size(); ++k ) {
           if( pathNamesAll[h].find(to_TriggerNames->at(k)) != std::string::npos ) {
              bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
              //bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
              bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
              //bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
              if( isLF || isBoth ) {
                 to_pt->at(k).push_back( obj.pt() );
                 to_eta->at(k).push_back( obj.eta() );
                 to_phi->at(k).push_back( obj.phi() );
              }
           }
         }
         //std::cout << "   " << pathNamesAll[h];
         //if (isBoth) std::cout << "(L,3)";
         //if (isL3 && !isBoth) std::cout << "(*,3)";
         //if (isLF && !isBoth) std::cout << "(L,*)";
         //if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
      }
      //std::cout << std::endl;
   }
   
   // muons
   // currently using the muon id taken from here:
   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
   
   for( const pat::Muon &m : *muons ) {
      if( !m.isPFMuon() ) continue;
      if( !(m.isGlobalMuon() || m.isTrackerMuon()) ) continue;
      if( !m.isLooseMuon() ) continue;
      if( m.pt() < 10. ) continue;
      if( fabs(m.eta()) > 2.5 ) continue;
      double pfRelIso = ( m.pfIsolationR04().sumChargedHadronPt 
                        + std::max(0., m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt ) ) 
                        / m.pt();
     
      muon_px->push_back( m.px() );
      muon_py->push_back( m.pz() );
      muon_pz->push_back( m.py() );
      muon_e->push_back( m.energy() );
      muon_phi->push_back( m.phi() );
      muon_eta->push_back( m.eta() );
      muon_iso->push_back( pfRelIso );
      muon_charge->push_back( m.charge() );
      muon_isLooseMuon->push_back( m.isLooseMuon() );
      muon_isMediumMuon->push_back( m.isMediumMuon() );
      muon_isTightMuon->push_back( m.isTightMuon(vertices->at(firstGoodVertexIdx)) );
   }
  


   // electrons
   // implemented the electron id taken from here:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
   // for isolation: https://www.sprace.org.br/twiki/bin/view/Main/ElectronIDRun2
   edm::Handle<reco::ConversionCollection> conversions;
   iEvent.getByToken( conversionsToken_, conversions );

   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
   iEvent.getByToken( eleVetoIdMapToken_, veto_id_decisions );
   iEvent.getByToken( eleLooseIdMapToken_, loose_id_decisions );
   iEvent.getByToken( eleMediumIdMapToken_, medium_id_decisions );
   iEvent.getByToken( eleTightIdMapToken_, tight_id_decisions );
   iEvent.getByToken( eleHEEPIdMapToken_, heep_id_decisions );


   for( size_t i = 0; i < electrons->size(); ++i ) {
      const auto e = electrons->ptrAt(i);
      if( e->pt() < 10 ) continue;
      if( fabs( e->eta()) > 2.5 ) continue;
      if( ! (*veto_id_decisions)[e] ) continue;


      reco::GsfElectron::PflowIsolationVariables pfIso = e->pfIsolationVariables();
      double absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
      double reliso = absiso/e->pt();

      electron_px->push_back( e->px() );
      electron_py->push_back( e->py() );
      electron_pz->push_back( e->pz() );
      electron_e->push_back( e->energy() );
      electron_phi->push_back( e->superCluster()->phi() );
      electron_eta->push_back( e->superCluster()->eta() );
      electron_charge->push_back( e->charge() );
      electron_isVeto->push_back( (*veto_id_decisions)[e]);
      electron_isLoose->push_back( (*loose_id_decisions)[e]);
      electron_isMedium->push_back( (*medium_id_decisions)[e]);
      electron_isTight->push_back( (*tight_id_decisions)[e]);
      electron_isHEEP->push_back( (*heep_id_decisions)[e]);
      electron_iso->push_back( reliso );
   }

   // and the taus:
   // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
   for( const pat::Tau &tau : *taus ) {
    if( tau.pt() < 20. ) continue;
    if( !( tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT") >= 0.5 && tau.tauID("againstMuonLoose3") >= 0.5 && tau.tauID("againstElectronVLooseMVA6") >= 0.5 ) ) continue;
    
    tau_px->push_back( tau.px() );
    tau_py->push_back( tau.py() );
    tau_pz->push_back( tau.pz() );
    tau_e->push_back( tau.energy() );
   }

   // jets
   // using the tight selection from:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
   int ctrJet = -1;
   if( useCHSJets ) {
   for( const pat::Jet &j : *jets ) {
    
     bool hasLargeMuonFraction = false;
     bool hasLargeEMFraction = false;
     bool hasSmallNeutralMultiplicity = false;

     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.99 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.9 ) hasLargeEMFraction = true; 
     if( j.numberOfDaughters() <= 1 ) continue;
     //if( j.muonEnergyFraction() >= 0.8 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) hasLargeMuonFraction = true;

     if( fabs(j.eta()) < 2.4 ) {
        if( j.chargedEmEnergyFraction() >= 0.99 ) continue;
        if( j.chargedEmEnergyFraction() >= 0.9 ) hasLargeEMFraction = true;
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
     }
     if( fabs(j.eta()) > 3.0 ) {
        if( j.neutralMultiplicity() <= 10 ) hasSmallNeutralMultiplicity = true;
     }
     if( j.pt() < 10. ) continue;
     
     tightJet_pt->push_back( j.pt() );
     tightJet_eta->push_back( j.eta() );
     tightJet_phi->push_back( j.phi() );


     if( hasSmallNeutralMultiplicity || hasLargeMuonFraction || hasLargeEMFraction ) continue;


     ctrJet += 1;
     
     std::vector<double> constVert_x;
     std::vector<double> constVert_y;
     std::vector<double> constVert_z;
     std::vector<double> const_pt;
     std::vector<double> const_eta;
     std::vector<double> const_phi;
     std::vector<int> const_charge;
     std::vector<double> const_px;
     std::vector<double> const_py;
     std::vector<double> const_pz;
     std::vector<double> const_pca0_x;
     std::vector<double> const_pca0_y;
     std::vector<double> const_pca0_z;
     std::vector<double> constVert_closestVertex_dxy;
     std::vector<double> constVert_closestVertex_dz;
     std::vector<double> constVert_closestVertex_d;
     std::vector<double> average_distance( vertices->size(), 0. );
     
     // loop over the jet constituents to find the vertex closest to each:
     for( unsigned int iD = 0; iD < j.numberOfDaughters(); ++iD ) {
        
        const pat::PackedCandidate &dau1 = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(iD));
        pat::PackedCandidate dau(dau1);
      
        const_px.push_back( dau.px() );
        const_py.push_back( dau.py() );
        const_pz.push_back( dau.pz() );
        const_pca0_x.push_back( dau.vertex().x() );
        const_pca0_y.push_back( dau.vertex().y() );
        const_pca0_z.push_back( dau.vertex().z() );

        // minimum distances (total, xy and z) for the jet constituent to any vertex
        double dMin = 100000.;
        double dxyMin = 10000.;
        double dzMin = 10000.;
        int ctr = -1;
       
        // get the unity vector pointing in the direction of the momentum and a reference point to build a self-made pseudo track
        // the 'track' is then (x(t), y(t), z(t) ) = (x,y,z) + t*(px,py,pz);
        // tmin, the time parameter for the point of closest approach is determined by minimising d = sqrt( (vx - x(t))^2 + (vy - y(t))^2 + (vz - z(t))^2);

        double jetVertex_x = -10000.;
        double jetVertex_y = -10000.;
        double jetVertex_z = -10000.;
        // first loop over the primary vertices
        for( const reco::Vertex &v : *vertices ) {
            ctr += 1;
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.position().x();
            double pv_y = v.position().y();
            double pv_z = v.position().z();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);

            // if the vertex is closer than the current reference vertex, set dMin, dxyMin, dzMin, and also change the vertex of reference
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              jetVertex_x = v.position().x();
              jetVertex_y = v.position().y();
              jetVertex_z = v.position().z();
             
            }
        }
  
        // now do the same for the secondary vertices
        // however, for some reason, I have to use additional variable here, jet constitutent won't accept a VertexCompositePtrCandidate as a new reference
        for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.vx();
            double pv_y = v.vy();
            double pv_z = v.vz();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);
            
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              ctr = -1;  
              jetVertex_x = v.vx();
              jetVertex_y = v.vy();
              jetVertex_z = v.vz();
          }
        }
        
        // now fill all the variables for the jet constituent
        constVert_closestVertex_dxy.push_back( dxyMin );
        constVert_closestVertex_dz.push_back( dzMin );
        constVert_closestVertex_d.push_back( dMin );
        constVert_x.push_back( jetVertex_x );
        constVert_y.push_back( jetVertex_y );
        constVert_z.push_back( jetVertex_z );
        const_pt.push_back( dau.pt() );
        const_eta.push_back( dau.eta() );
        const_phi.push_back( dau.phi() );
        const_charge.push_back( dau.charge() );
     }

     // now calculate the jet-vertex:
     int nCons = 0;
     double weightednCons = 0.;
     std::vector<double> error(3,0.);
     double theScore = 0.;
     std::vector<double> position = CalculateVertex( constVert_x, constVert_y, constVert_z, const_pt, const_charge, constVert_closestVertex_d, nCons, weightednCons, error, theScore );
     
     double averageDistance = 0.;
     double rmsDistance = 0.;
     double chargedConsts = 0.;
     // now I know the jet vertex - calculated average distance to the highest score vertex
     for( unsigned int iConst = 0; iConst < const_px.size(); ++iConst ) {
        if( const_charge.at(iConst) == 0 ) continue;
        double const_p = sqrt( const_px.at(iConst)*const_px.at(iConst) + const_py.at(iConst)*const_py.at(iConst) + const_pz.at(iConst)*const_pz.at(iConst) );
        double px = const_px.at(iConst)/const_p;
        double py = const_py.at(iConst)/const_p;
        double pz = const_pz.at(iConst)/const_p;
        double x = const_pca0_x.at(iConst);
        double y = const_pca0_y.at(iConst);
        double z = const_pca0_z.at(iConst);
        double target_x = position.at(0);
        double target_y = position.at(1);
        double target_z = position.at(2);
            
        double tmin = - ( px*(x-target_x) + py*(y-target_y) + pz*(z-target_z) );
        double dx_min = target_x - x - tmin*px;
        double dy_min = target_y - y - tmin*py;
        double dz_min = target_z - z - tmin*pz;
           
        double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
        double d = sqrt( dxy*dxy + dz_min*dz_min);
        averageDistance += d;
        chargedConsts += 1.;
     }
     for( unsigned int iConst = 0; iConst < const_px.size(); ++iConst ) {
        if( const_charge.at(iConst) == 0 ) continue;
        double const_p = sqrt( const_px.at(iConst)*const_px.at(iConst) + const_py.at(iConst)*const_py.at(iConst) + const_pz.at(iConst)*const_pz.at(iConst) );
        double px = const_px.at(iConst)/const_p;
        double py = const_py.at(iConst)/const_p;
        double pz = const_pz.at(iConst)/const_p;
        double x = const_pca0_x.at(iConst);
        double y = const_pca0_y.at(iConst);
        double z = const_pca0_z.at(iConst);
        double target_x = position.at(0);
        double target_y = position.at(1);
        double target_z = position.at(2);
            
        double tmin = - ( px*(x-target_x) + py*(y-target_y) + pz*(z-target_z) );
        double dx_min = target_x - x - tmin*px;
        double dy_min = target_y - y - tmin*py;
        double dz_min = target_z - z - tmin*pz;
           
        double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
        double d = sqrt( dxy*dxy + dz_min*dz_min);
        rmsDistance += (d - averageDistance)*(d - averageDistance);
     }
     rmsDistance /= chargedConsts;
     rmsDistance = sqrt( rmsDistance );
     averageDistance /= chargedConsts;

     // and fill allthe jet variables
     std::vector<double> jpts;

     jecUnc->setJetEta(j.eta() );
     jecUnc->setJetPt(j.pt()); // here you must use the CORRECTED jet pt
     double unc = jecUnc->getUncertainty(true);
     
     jpts.push_back( j.pt() );
     jpts.push_back( j.pt()*(1.+unc));
     jpts.push_back( j.pt()*(1.-unc));
     
     jet_pt->push_back( jpts );
     jet_eta->push_back( j.eta() );
     jet_phi->push_back( j.phi() );
     for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) { 
       jet_btagInfo->at(iBtagAlgo)->push_back( j.bDiscriminator( btagAlgorithms.at(iBtagAlgo) ) );
     }
     jet_vertex_x->push_back( position.at(0) );
     jet_vertex_y->push_back( position.at(1) );
     jet_vertex_z->push_back( position.at(2) );
     jet_vertex_score->push_back( theScore );
     jet_nCons->push_back( nCons );
     jet_averageDistance->push_back( averageDistance );
     jet_rmsDistance->push_back( rmsDistance );
    /* 
     jet_constVertex_x->push_back( constVert_x ); 
     jet_constVertex_y->push_back( constVert_y ); 
     jet_constVertex_z->push_back( constVert_z ); 
     jet_const_closestVertex_dxy->push_back(constVert_closestVertex_dxy);
     jet_const_closestVertex_dz->push_back(constVert_closestVertex_dz);
     jet_const_closestVertex_d->push_back(constVert_closestVertex_d);
     jet_const_pt->push_back( const_pt );
     jet_const_eta->push_back( const_eta );
     jet_const_phi->push_back( const_phi );
     jet_const_charge->push_back( const_charge );
     jet_const_px->push_back( const_px );
     jet_const_py->push_back( const_py );
     jet_const_pz->push_back( const_pz );
     jet_const_pca0_x->push_back( const_pca0_x );
     jet_const_pca0_y->push_back( const_pca0_y );
     jet_const_pca0_z->push_back( const_pca0_z );
     */ 
  }
  }
   else {
   for( const reco::Jet &jr : *jetsnoCHS ) {

     const pat::Jet j( jr );
     
     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.90 ) continue;
     if( j.numberOfDaughters() <= 1 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) continue;
     
     if( fabs(j.eta()) < 2.4 ) {
        if( j.chargedEmEnergyFraction() >= 0.9 ) continue;
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
     }
     if( fabs(j.eta()) > 3.0 ) {
        if( j.neutralMultiplicity() <= 10 ) continue;
     }
     if( j.pt() < 10. ) continue;
     
     ctrJet += 1;
     
     std::vector<double> constVert_x;
     std::vector<double> constVert_y;
     std::vector<double> constVert_z;
     std::vector<double> const_pt;
     std::vector<double> const_eta;
     std::vector<double> const_phi;
     std::vector<int> const_charge;
     std::vector<double> const_px;
     std::vector<double> const_py;
     std::vector<double> const_pz;
     std::vector<double> const_pca0_x;
     std::vector<double> const_pca0_y;
     std::vector<double> const_pca0_z;
     std::vector<double> constVert_closestVertex_dxy;
     std::vector<double> constVert_closestVertex_dz;
     std::vector<double> constVert_closestVertex_d;
     std::vector<double> average_distance( vertices->size(), 0. );
     
     // loop over the jet constituents to find the vertex closest to each:
     for( unsigned int iD = 0; iD < j.numberOfDaughters(); ++iD ) {
        
        const pat::PackedCandidate &dau1 = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(iD));
        pat::PackedCandidate dau(dau1);
      
        const_px.push_back( dau.px() );
        const_py.push_back( dau.py() );
        const_pz.push_back( dau.pz() );
        const_pca0_x.push_back( dau.vertex().x() );
        const_pca0_y.push_back( dau.vertex().y() );
        const_pca0_z.push_back( dau.vertex().z() );

        // minimum distances (total, xy and z) for the jet constituent to any vertex
        double dMin = 100000.;
        double dxyMin = 10000.;
        double dzMin = 10000.;
        int ctr = -1;
       
        // get the unity vector pointing in the direction of the momentum and a reference point to build a self-made pseudo track
        // the 'track' is then (x(t), y(t), z(t) ) = (x,y,z) + t*(px,py,pz);
        // tmin, the time parameter for the point of closest approach is determined by minimising d = sqrt( (vx - x(t))^2 + (vy - y(t))^2 + (vz - z(t))^2);

        double jetVertex_x = -10000.;
        double jetVertex_y = -10000.;
        double jetVertex_z = -10000.;
        // first loop over the primary vertices
        for( const reco::Vertex &v : *vertices ) {
            ctr += 1;
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.position().x();
            double pv_y = v.position().y();
            double pv_z = v.position().z();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);

            // if the vertex is closer than the current reference vertex, set dMin, dxyMin, dzMin, and also change the vertex of reference
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              jetVertex_x = v.position().x();
              jetVertex_y = v.position().y();
              jetVertex_z = v.position().z();
             
            }
        }
  
        // now do the same for the secondary vertices
        // however, for some reason, I have to use additional variable here, jet constitutent won't accept a VertexCompositePtrCandidate as a new reference
        for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
            double x = dau.vertex().x();
            double y = dau.vertex().y();
            double z = dau.vertex().z();
            double px = dau.px()/dau.p();
            double py = dau.py()/dau.p();
            double pz = dau.pz()/dau.p();
            double pv_x = v.vx();
            double pv_y = v.vy();
            double pv_z = v.vz();
            double tmin = - ( px*(x-pv_x) + py*(y-pv_y) + pz*(z-pv_z) );
            double dx_min = pv_x - x - tmin*px;
            double dy_min = pv_y - y - tmin*py;
            double dz_min = pv_z - z - tmin*pz;
           
            double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
            double d = sqrt( dxy*dxy + dz_min*dz_min);
            
            if( d < dMin ) {
              dMin = d;
              dxyMin = dxy;
              dzMin = dz_min;
              ctr = -1;  
              jetVertex_x = v.vx();
              jetVertex_y = v.vy();
              jetVertex_z = v.vz();
          }
        }
        
        // now fill all the variables for the jet constituent
        constVert_closestVertex_dxy.push_back( dxyMin );
        constVert_closestVertex_dz.push_back( dzMin );
        constVert_closestVertex_d.push_back( dMin );
        constVert_x.push_back( jetVertex_x );
        constVert_y.push_back( jetVertex_y );
        constVert_z.push_back( jetVertex_z );
        const_pt.push_back( dau.pt() );
        const_eta.push_back( dau.eta() );
        const_phi.push_back( dau.phi() );
        const_charge.push_back( dau.charge() );
     }

     // now calculate the jet-vertex:
     int nCons = 0;
     double weightednCons = 0.;
     std::vector<double> error(3,0.);
     double theScore = 0.;
     std::vector<double> position = CalculateVertex( constVert_x, constVert_y, constVert_z, const_pt, const_charge, constVert_closestVertex_d, nCons, weightednCons, error, theScore );
     
     double averageDistance = 0.;
     double rmsDistance = 0.;
     double chargedConsts = 0.;
     // now I know the jet vertex - calculated average distance to the highest score vertex
     for( unsigned int iConst = 0; iConst < const_px.size(); ++iConst ) {
        if( const_charge.at(iConst) == 0 ) continue;
        double const_p = sqrt( const_px.at(iConst)*const_px.at(iConst) + const_py.at(iConst)*const_py.at(iConst) + const_pz.at(iConst)*const_pz.at(iConst) );
        double px = const_px.at(iConst)/const_p;
        double py = const_py.at(iConst)/const_p;
        double pz = const_pz.at(iConst)/const_p;
        double x = const_pca0_x.at(iConst);
        double y = const_pca0_y.at(iConst);
        double z = const_pca0_z.at(iConst);
        double target_x = position.at(0);
        double target_y = position.at(1);
        double target_z = position.at(2);
            
        double tmin = - ( px*(x-target_x) + py*(y-target_y) + pz*(z-target_z) );
        double dx_min = target_x - x - tmin*px;
        double dy_min = target_y - y - tmin*py;
        double dz_min = target_z - z - tmin*pz;
           
        double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
        double d = sqrt( dxy*dxy + dz_min*dz_min);
        averageDistance += d;
        chargedConsts += 1.;
     }
     for( unsigned int iConst = 0; iConst < const_px.size(); ++iConst ) {
        if( const_charge.at(iConst) == 0 ) continue;
        double const_p = sqrt( const_px.at(iConst)*const_px.at(iConst) + const_py.at(iConst)*const_py.at(iConst) + const_pz.at(iConst)*const_pz.at(iConst) );
        double px = const_px.at(iConst)/const_p;
        double py = const_py.at(iConst)/const_p;
        double pz = const_pz.at(iConst)/const_p;
        double x = const_pca0_x.at(iConst);
        double y = const_pca0_y.at(iConst);
        double z = const_pca0_z.at(iConst);
        double target_x = position.at(0);
        double target_y = position.at(1);
        double target_z = position.at(2);
            
        double tmin = - ( px*(x-target_x) + py*(y-target_y) + pz*(z-target_z) );
        double dx_min = target_x - x - tmin*px;
        double dy_min = target_y - y - tmin*py;
        double dz_min = target_z - z - tmin*pz;
           
        double dxy = sqrt(dx_min*dx_min + dy_min*dy_min);
        double d = sqrt( dxy*dxy + dz_min*dz_min);
        rmsDistance += (d - averageDistance)*(d - averageDistance);
     }
     rmsDistance /= chargedConsts;
     rmsDistance = sqrt( rmsDistance );
     averageDistance /= chargedConsts;

     // and fill allthe jet variables
     std::vector<double> jpts;

     jecUnc->setJetEta(j.eta() );
     jecUnc->setJetPt(j.pt()); // here you must use the CORRECTED jet pt
     double unc = jecUnc->getUncertainty(true);
     
     jpts.push_back( j.pt() );
     jpts.push_back( j.pt()*(1.+unc));
     jpts.push_back( j.pt()*(1.-unc));
     
     jet_pt->push_back( jpts );
     jet_eta->push_back( j.eta() );
     jet_phi->push_back( j.phi() );
     for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) { 
       jet_btagInfo->at(iBtagAlgo)->push_back( j.bDiscriminator( btagAlgorithms.at(iBtagAlgo) ) );
     }
     jet_vertex_x->push_back( position.at(0) );
     jet_vertex_y->push_back( position.at(1) );
     jet_vertex_z->push_back( position.at(2) );
     jet_vertex_score->push_back( theScore );
     jet_nCons->push_back( nCons );
     jet_averageDistance->push_back( averageDistance );
     jet_rmsDistance->push_back( rmsDistance );
    /* 
     jet_constVertex_x->push_back( constVert_x ); 
     jet_constVertex_y->push_back( constVert_y ); 
     jet_constVertex_z->push_back( constVert_z ); 
     jet_const_closestVertex_dxy->push_back(constVert_closestVertex_dxy);
     jet_const_closestVertex_dz->push_back(constVert_closestVertex_dz);
     jet_const_closestVertex_d->push_back(constVert_closestVertex_d);
     jet_const_pt->push_back( const_pt );
     jet_const_eta->push_back( const_eta );
     jet_const_phi->push_back( const_phi );
     jet_const_charge->push_back( const_charge );
     jet_const_px->push_back( const_px );
     jet_const_py->push_back( const_py );
     jet_const_pz->push_back( const_pz );
     jet_const_pca0_x->push_back( const_pca0_x );
     jet_const_pca0_y->push_back( const_pca0_y );
     jet_const_pca0_z->push_back( const_pca0_z );
     */ 
  }
  }

  
   // now fill the secondary vertex information
   for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
      secVertex_x -> push_back( v.vx() );
      secVertex_y -> push_back( v.vy() );
      secVertex_z -> push_back( v.vz() );
      secVertex_ndof -> push_back( v.vertexNdof() );
      secVertex_chi2 -> push_back( v.vertexChi2() );
      secVertex_pt -> push_back( v.pt() );
      secVertex_dx -> push_back( v.vertexCovariance(0,0) );
      secVertex_dy -> push_back( v.vertexCovariance(1,1) );
      secVertex_dz -> push_back( v.vertexCovariance(2,2) );
   }


   // and fill the met
   // now also store the shifted met values considering the uncertainties
   const pat::MET &themet = mets->front();
   met->push_back( themet.corPt(pat::MET::METCorrectionLevel::Type01) ); 
   met_x->push_back( themet.corPx(pat::MET::METCorrectionLevel::Type01) );
   met_y->push_back( themet.corPy(pat::MET::METCorrectionLevel::Type01) );
   for( pat::MET::METUncertainty iMETUNC = pat::MET::METUncertainty::JetResUp; iMETUNC < pat::MET::METUncertainty::METUncertaintySize; iMETUNC = pat::MET::METUncertainty( iMETUNC + 1) ) {
     met->push_back( themet.shiftedSumEt( iMETUNC, pat::MET::METCorrectionLevel::Type01 ) );
     met_x->push_back( themet.shiftedPx( iMETUNC, pat::MET::METCorrectionLevel::Type01 ) );
     met_y->push_back( themet.shiftedPy( iMETUNC, pat::MET::METCorrectionLevel::Type01 ) );
   }

   // if is MC, get the gen level HT
   if( !isData ) {
     GenLevel_HT = 0.;
     for( size_t i = 0; i < pruned->size(); ++i ) {
       if( (*pruned)[i].status() != 23 ) continue;
       if( (*pruned)[i].pdgId() != 21 && abs((*pruned)[i].pdgId()) > 6 ) continue;
       GenLevel_HT += (*pruned)[i].pt();
     }
   }
   // and the gen info if required:
   if( storeGenData && !isData ) {
     for(size_t i = 0; i<packed->size(); i++){
        std::vector<int> mothersId;
        std::vector<int> mothersStatus;
        std::vector<double> mothersPx;
        std::vector<double> mothersPy;
        std::vector<double> mothersPz;
        std::vector<double> mothersE;
        const reco::Candidate *mother = (*packed)[i].mother(0) ; 

        while( mother ) {
          mothersId.push_back( mother->pdgId() );
          mothersStatus.push_back( mother->status() );
          mothersPx.push_back( mother->px() );
          mothersPy.push_back( mother->py() );
          mothersPz.push_back( mother->pz() );
          mothersE.push_back( mother->energy() );
          mother = (reco::GenParticle*)mother->mother(0);
        }
        mct_px -> push_back( (*packed)[i].px() );
        mct_py -> push_back( (*packed)[i].py() );
        mct_pz -> push_back( (*packed)[i].pz() );
        mct_E -> push_back( (*packed)[i].energy() );
        mct_id -> push_back( (*packed)[i].pdgId() );
        mct_status -> push_back( (*packed)[i].status() );
        mct_parentId -> push_back( mothersId );
        mct_parentStatus -> push_back( mothersStatus );
        mct_parentPx -> push_back( mothersPx );
        mct_parentPy -> push_back( mothersPy );
        mct_parentPz -> push_back( mothersPz );
        mct_parentE -> push_back( mothersE );
     }
   }


   delete jecUnc;
   // finally, write it to the tree
   nEventsAccepted += 1;
   tOutput->Fill(); 

}

std::vector<double> CalculateVertex( std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> weight, std::vector<int> charge, std::vector<double> distance, int &nConsidered, double &weightednConsidered, std::vector<double> &error, double &maxScore ) {

   nConsidered = 0;
   std::vector<double> diff_x;
   std::vector<double> diff_y;
   std::vector<double> diff_z;
   std::vector<double> score;


   for( unsigned int i = 0; i < x.size(); ++i ) {
      if( charge.at(i) == 0 ) continue;
      nConsidered += 1;
      bool knownPoint = false;
      int iKnown = -1;
      for( unsigned int i2 = 0; i2 < diff_x.size(); ++i2 ) {
        if( fabs( diff_x.at(i2) - x.at(i) ) < 1.e-10 && fabs( diff_y.at(i2) - y.at(i) ) < 1.e-10 && fabs( diff_z.at(i2) - z.at(i) ) < 1.e-10 ) {
            knownPoint = true;
            iKnown = i2;
        }
      }

      if( knownPoint ) {
        if( distance.at(i) == 0. ) score.at(iKnown) += 1.e12;
        else                       score.at(iKnown) += weight.at(i)/distance.at(i);
      }
      else {
        diff_x.push_back( x.at(i) );
        diff_y.push_back( y.at(i) );
        diff_z.push_back( z.at(i) );
        if( distance.at(i) == 0. ) score.push_back( 1.e12 );
        else                       score.push_back( weight.at(i)/distance.at(i) );
      }
   }

   double scoreMax = 0.;
   std::vector<double> mean(3, -10000.);
   for( unsigned int i = 0; i < diff_x.size(); ++i ) {
        if ( score.at(i) > scoreMax ) {
            scoreMax = score.at(i);
            maxScore = scoreMax;
            mean.at(0) = diff_x.at(i);
            mean.at(1) = diff_y.at(i);
            mean.at(2) = diff_z.at(i);
        }
    }
    return mean;
}


// ------------ method called once each job just before starting event loop  ------------
void 
LLGDVMiniAODAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LLGDVMiniAODAnalysis::endJob() 
{

  // save the metadata tree
  tMetaData->Fill();

  // save the output tree and file
  gDirectory = fOutput;
  tOutput->Write();
  tMetaData->Write();
  fOutput->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
LLGDVMiniAODAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
LLGDVMiniAODAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
LLGDVMiniAODAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
LLGDVMiniAODAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLGDVMiniAODAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LLGDVMiniAODAnalysis);

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
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TFile.h"
#include "TTree.h"
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
      edm::EDGetTokenT<reco::PFJetCollection> jetToken2_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secVtxToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<edm::TriggerResults> METFilterBits_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetToken electronToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

      // the b-tagging algorithms for which we want the numbers
      std::vector<std::string> btagAlgorithms;


      // the output file and tree
      TFile *fOutput = new TFile("RecoOutput.root", "RECREATE");
      TTree *tOutput = new TTree("RecoData", "RecoData");
      TTree *tMetaData = new TTree("MetaData", "MetaData" );

      //setup the variables for the metadata tree first:
      int nEventsProcessed = 0;
      int nEventsAccepted = 0;

      
      // the variables used for output
      // missing transverse energy
      double met = 0.;
      double met_x = 0.;
      double met_y = 0.;

      // the jet variables
      std::vector<double> *jet_eta = new std::vector<double>;
      std::vector<double> *jet_phi = new std::vector<double>;
      std::vector<double> *jet_pt = new std::vector<double>;
      std::vector<std::vector<double>* > *jet_btagInfo = new std::vector<std::vector<double>* >;

      // the muon variables
      std::vector<double> *muon_px = new std::vector<double>;
      std::vector<double> *muon_py = new std::vector<double>;
      std::vector<double> *muon_pz = new std::vector<double>;
      std::vector<double> *muon_eta = new std::vector<double>;
      std::vector<double> *muon_phi = new std::vector<double>;
      std::vector<double> *muon_iso = new std::vector<double>;
      std::vector<double> *muon_charge = new std::vector<double>;
      std::vector<bool> *muon_isTightMuon = new std::vector<bool>;
      std::vector<bool> *muon_isLooseMuon = new std::vector<bool>;

      // the electron variables
      std::vector<double> *electron_px = new std::vector<double>;
      std::vector<double> *electron_py = new std::vector<double>;
      std::vector<double> *electron_pz = new std::vector<double>;
      std::vector<double> *electron_eta = new std::vector<double>;
      std::vector<double> *electron_phi = new std::vector<double>;
      std::vector<double> *electron_iso = new std::vector<double>;
      std::vector<double> *electron_charge = new std::vector<double>;
      std::vector<bool> *electron_isVeto = new std::vector<bool>;
      std::vector<bool> *electron_isLoose = new std::vector<bool>;
      std::vector<bool> *electron_isMedium = new std::vector<bool>;
      std::vector<bool> *electron_isTight = new std::vector<bool>;
      std::vector<bool> *electron_isHEEP = new std::vector<bool>;

      // the trigger bits and names
      std::vector<int> *triggerBits = new std::vector<int>;
      std::vector<std::string> *triggerNames = new std::vector<std::string>;
      std::vector<int> *METFilterBits = new std::vector<int>;
      std::vector<std::string> *METFilterNames = new std::vector<std::string>;

      // the jet constituents 
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
      std::vector<double> *secVertex_dx = new std::vector<double>;
      std::vector<double> *secVertex_dy = new std::vector<double>;
      std::vector<double> *secVertex_dz = new std::vector<double>;


      // event metadata
      int RunNumber = 0;
      int EventNumber = 0;
      int LuminosityBlock = 0;

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
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  secVtxToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  METFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METFilters"))),
  genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
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
   triggerNames->push_back("HLT_PFMET170_NoiseCleaned_v1" );
   if( RunLeptonTriggers ) {
     triggerNames->push_back("HLT_Ele27_eta2p1_WP85_Gsf_v1");
     triggerNames->push_back("HLT_Ele32_eta2p1_WP85_Gsf_v1");
     triggerNames->push_back("HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1");
     triggerNames->push_back("HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1");
     triggerNames->push_back("HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v1");
     triggerNames->push_back("HLT_Ele27_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v1");
     triggerNames->push_back("HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet40_v1");
     triggerNames->push_back("HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v1");
     triggerNames->push_back("HLT_Ele32_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v1");
     triggerNames->push_back("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1");
     triggerNames->push_back("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1");
     triggerNames->push_back("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1");
     triggerNames->push_back("HLT_Ele95_CaloIdVT_GsfTrkIdT_v1");
     triggerNames->push_back("HLT_Ele20WP60_Ele8_Mass55_v1");
     triggerNames->push_back("HLT_Ele25WP60_SC4_Mass55_v1");
     triggerNames->push_back("HLT_Mu40_v1");
     triggerNames->push_back("HLT_Mu17_Mu8_v1");
     triggerNames->push_back("HLT_Mu17_TkMu8_v1");
     triggerNames->push_back("HLT_Mu30_TkMu11_v1");
     triggerNames->push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1");
     triggerNames->push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1");
     triggerNames->push_back("HLT_Mu25_TkMu0_dEta18_Onia_v1");
     triggerNames->push_back("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_v1");
     triggerNames->push_back("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_v1");
     triggerNames->push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v1");
     triggerNames->push_back("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1");
     triggerNames->push_back("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1");
   }
   // the btag algorithms for which we want infos
   btagAlgorithms.push_back("pfJetBProbabilityBJetTags");
   btagAlgorithms.push_back("pfJetProbabilityBJetTags");
   btagAlgorithms.push_back("pfTrackCountingHighPurBJetTags");
   btagAlgorithms.push_back("pfTrackCountingHighEffBJetTags");
   btagAlgorithms.push_back("pfCombinedInclusiveSecondaryVertexV2BJetTags");
   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
      jet_btagInfo->push_back( new std::vector<double> );
   }

   // set the output branches for the tree
   tOutput -> Branch("RecoJet_eta", &jet_eta );
   tOutput -> Branch("RecoJet_phi", &jet_phi );
   tOutput -> Branch("RecoJet_pt", &jet_pt );
   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
     std::string branchname = "RecoJet_btag_" + btagAlgorithms.at(iBtagAlgo);
     tOutput -> Branch( branchname.c_str(), &(jet_btagInfo->at(iBtagAlgo)) ); 
   }
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
   tOutput -> Branch("RecoMuon_px", &muon_px );
   tOutput -> Branch("RecoMuon_py", &muon_pz );
   tOutput -> Branch("RecoMuon_pz", &muon_py );
   tOutput -> Branch("RecoMuon_eta", &muon_eta );
   tOutput -> Branch("RecoMuon_phi", &muon_phi );
   tOutput -> Branch("RecoMuon_iso", &muon_iso);
   tOutput -> Branch("RecoMuon_charge", &muon_charge);
   tOutput -> Branch("RecoMuon_isTightMuon", &muon_isTightMuon );
   tOutput -> Branch("RecoMuon_isLooseMuon", &muon_isLooseMuon );
   tOutput -> Branch("RecoElectron_px", &electron_px );
   tOutput -> Branch("RecoElectron_py", &electron_pz );
   tOutput -> Branch("RecoElectron_pz", &electron_py );
   tOutput -> Branch("RecoElectron_eta", &electron_eta );
   tOutput -> Branch("RecoElectron_phi", &electron_phi );
   tOutput -> Branch("RecoElectron_iso", &electron_iso );
   tOutput -> Branch("RecoElectron_charge", &electron_charge );
   tOutput -> Branch("RecoElectron_isVeto", &electron_isVeto );
   tOutput -> Branch("RecoElectron_isLoose", &electron_isLoose );
   tOutput -> Branch("RecoElectron_isMedium", &electron_isMedium );
   tOutput -> Branch("RecoElectron_isTight", &electron_isTight );
   tOutput -> Branch("RecoElectron_isHEEP", &electron_isHEEP );
   tOutput -> Branch("TriggerBits", &triggerBits );
   tOutput -> Branch("TriggerNames", &triggerNames );
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
   tOutput -> Branch("RecoSecVertex_pt", &secVertex_pt );
   tOutput -> Branch("RecoSecVertex_xError", &secVertex_dx );
   tOutput -> Branch("RecoSecVertex_yError", &secVertex_dy );
   tOutput -> Branch("RecoSecVertex_zError", &secVertex_dz );
   tOutput -> Branch("MET", &met );
   tOutput -> Branch("MET_x", &met_x );
   tOutput -> Branch("MET_y", &met_y );
   tOutput -> Branch("EventNumber", &EventNumber );
   tOutput -> Branch("RunNumber", &RunNumber );
   tOutput -> Branch("LuminosityBlock", &LuminosityBlock );
  
   // the metadata tree
   tMetaData -> Branch("nEventsProcessed", &nEventsProcessed );
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
   met = 0.;
   met_x = 0.;
   met_y = 0.;
   muon_px->clear();
   muon_py->clear();
   muon_pz->clear();
   muon_eta->clear();
   muon_phi->clear();
   muon_iso->clear();
   muon_charge->clear();
   muon_isTightMuon->clear();
   muon_isLooseMuon->clear();
   electron_px->clear();
   electron_py->clear();
   electron_pz->clear();
   electron_eta->clear();
   electron_phi->clear();
   electron_iso->clear();
   electron_charge->clear();
   electron_isVeto->clear();
   electron_isLoose->clear(); 
   electron_isMedium->clear();
   electron_isTight->clear();
   electron_isHEEP->clear();
   triggerBits->clear();
   METFilterBits->clear();
   jet_eta->clear();
   jet_phi->clear();
   jet_pt->clear();
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
   secVertex_pt -> clear();
   secVertex_dx -> clear();
   secVertex_dy -> clear();
   secVertex_dz -> clear();
   
   for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) {
     jet_btagInfo->at(iBtagAlgo)->clear();
   }

   // get all the physics objects from the miniAOD
   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   Handle<edm::View<reco::GenParticle> > pruned;
   iEvent.getByToken(prunedGenToken_, pruned);

   Handle<edm::View<pat::PackedGenParticle> > packed;
   iEvent.getByToken(packedGenToken_,packed);
   
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
   
   edm::Handle<edm::TriggerResults> evMETFilterBits;
   iEvent.getByToken( METFilterBits_, evMETFilterBits );
  
   Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken( genJetToken_, genJets );
  
   Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);


   edm::EventAuxiliary aux = iEvent.eventAuxiliary();
   edm::EventID id = aux.id();
  
   EventNumber = id.event();
   RunNumber = id.run();
   LuminosityBlock = id.luminosityBlock();


   
   const edm::TriggerNames &filterNames = iEvent.triggerNames(*evMETFilterBits);
   bool passETMissFilter = true;
   for(unsigned int i = 0; i < evMETFilterBits->size(); ++i ) {
    if(    filterNames.triggerName(i) == "Flag_HBHENoiseFilter" 
        || filterNames.triggerName(i) == "Flag_CSCTightHaloFilter"
        //|| filterNames.triggerName(i) == "Flag_goodVertices"
        || filterNames.triggerName(i) == "Flag_eeBadScFilter" ) {
          if( !evMETFilterBits->accept(i) ) passETMissFilter = false;    
        }
   } 
   if( !passETMissFilter ) return;


   const edm::TriggerNames &names = iEvent.triggerNames(*evTriggerBits);
   bool passTrigger = false;
   //for(unsigned int i = 0; i < evTriggerBits->size(); ++i ) {
   // std::cout << "got trigger " << names.triggerName(i) << std::endl;
   //}
   for( unsigned int j = 0; j < triggerNames->size(); ++j ) {
    for(unsigned int i = 0; i < evTriggerBits->size(); ++i ) {
        if( names.triggerName(i) == triggerNames->at(j) ) {
          triggerBits->push_back( evTriggerBits->accept(i) ? 1 : 0 );
          if( evTriggerBits->accept(i) ) passTrigger = true;
          continue;
        }
      }
   }
   // store only events in the ntuple which pass the trigger
   if( !passTrigger ) return;
 
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

   if( firstGoodVertexIdx == -1 ) return;


   
   // muons
   // currently using the muon id taken from here:
   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
   // using loose PF combined relative isolation cut (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyObjectExperts)
   // adding a pT cut of 15 GeV
   // not using a cut that constrains the muon to the pv
   for( const pat::Muon &m : *muons ) {
      if( !m.isPFMuon() ) continue;
      if( !( m.isGlobalMuon() || m.isTrackerMuon()) ) continue;
      if( !m.isLooseMuon() ) continue;
      if( m.pt() < 15. ) continue;
      if( fabs(m.eta()) > 2.5 ) continue;
      double pfRelIso = ( m.pfIsolationR04().sumChargedHadronPt 
                        + std::max(0., m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt ) ) 
                        / m.pt();
      //if( pfRelIso > 0.2 ) continue;
      
      muon_px->push_back( m.px() );
      muon_py->push_back( m.pz() );
      muon_pz->push_back( m.py() );
      muon_phi->push_back( m.phi() );
      muon_eta->push_back( m.eta() );
      muon_iso->push_back( pfRelIso );
      muon_charge->push_back( m.charge() );
      muon_isLooseMuon->push_back( m.isLooseMuon() );
      muon_isTightMuon->push_back( m.isTightMuon(vertices->at(firstGoodVertexIdx)) );
   }
  

   // store only events without muons
   // if( muon_px->size() > 0 ) return;


   // electrons
   // implemented the electron id taken from here:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/Eg2012AnalysesSupportingMaterial#Electrons
   // and also
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyObjectExperts
   // not using a cut that constraints the electron to the PV
   // calculate the isolation according to: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
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

      electron_px->push_back( e->px() );
      electron_py->push_back( e->py() );
      electron_pz->push_back( e->pz() );
      electron_phi->push_back( e->superCluster()->phi() );
      electron_eta->push_back( e->superCluster()->eta() );
      electron_charge->push_back( e->charge() );
      electron_isVeto->push_back( (*veto_id_decisions)[e]);
      electron_isLoose->push_back( (*loose_id_decisions)[e]);
      electron_isMedium->push_back( (*medium_id_decisions)[e]);
      electron_isTight->push_back( (*tight_id_decisions)[e]);
      electron_isHEEP->push_back( (*heep_id_decisions)[e]);
      //electron_iso->push_back( iso );
   }

   // store only events without electrons
   // if( electron_px->size() > 0 ) return;


   // jets
   // using the tight selection from:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
   int ctrJet = -1;
   for( const pat::Jet &j : *jets ) {
     
     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.90 ) continue;
     if( j.numberOfDaughters() <= 1 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) continue;
     
     if( fabs(j.eta()) < 2.4 ) {
        if( j.chargedEmEnergyFraction() >= 0.9 ) continue;
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
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

     // and fill allthe jet variables
     jet_pt->push_back( j.pt() );
     jet_eta->push_back( j.eta() );
     jet_phi->push_back( j.phi() );
     for( unsigned int iBtagAlgo = 0; iBtagAlgo < btagAlgorithms.size(); ++iBtagAlgo ) { 
       jet_btagInfo->at(iBtagAlgo)->push_back( j.bDiscriminator( btagAlgorithms.at(iBtagAlgo) ) );
     } 
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
   }
 

  
   // now fill the secondary vertex information
   for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
      secVertex_x -> push_back( v.vx() );
      secVertex_y -> push_back( v.vy() );
      secVertex_z -> push_back( v.vz() );
      secVertex_ndof -> push_back( v.vertexNdof() );
      secVertex_pt -> push_back( v.pt() );
      secVertex_dx -> push_back( v.vertexCovariance(0,0) );
      secVertex_dy -> push_back( v.vertexCovariance(1,1) );
      secVertex_dz -> push_back( v.vertexCovariance(2,2) );
   }


   // and fill the met
   const pat::MET &themet = mets->front();
   met = themet.pt(); 
   met_x = themet.px();
   met_y = themet.py();

   // finally, write it to the tree
   nEventsAccepted += 1;
   tOutput->Fill(); 

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

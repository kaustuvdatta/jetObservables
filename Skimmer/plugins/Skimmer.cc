// -*- C++ -*-
//
// Package:    Analysis/QJetMass
// Class:      Skimmer
// 
/**\class Skimmer Skimmer.cc Analysis/QJetMass/plugins/Skimmer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Alejandro Gomez Espinosa
//         Created:  Wed, 06 Mar 2019 13:50:12 GMT
// Edited by: Kaustuv Datta
//        on: Mon, 11 Mar 2019
//  Based on https://github.com/cms-jet/JMEValidator/blob/28a125c1e2cc94354c09e4c5a3288d3f07f6e644/python/FrameworkConfiguration.py
//  and https://github.com/cms-jet/JMEValidator/blob/28a125c1e2cc94354c09e4c5a3288d3f07f6e644/src/JMEJetAnalyzer.cc
//

// system include files
#include <memory>
#include <vector>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

//
// class declaration
//
using namespace std;
using namespace edm;

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Skimmer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit Skimmer(const edm::ParameterSet&);
  ~Skimmer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  TTree *SimpleGentree;
      
  edm::EDGetTokenT<edm::View<reco::GenJet>> srcGenJet_;
  std::string srcGenjetNSub0p5;
  std::string srcGenjetNSub1;
  std::string srcGenjetNSub2;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau1;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau2;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau3;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau4;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau5;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau6;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau7;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau8;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau9;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub0p5tau10;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau1;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau2;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau3;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau4;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau5;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau6;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau7;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau8;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau9;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub1tau10;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau1;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau2;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau3;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau4;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau5;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau6;
  edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau7;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau8;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau9;
  // edm::EDGetTokenT< edm::ValueMap<float> > Nsub2tau10;

  vector<float> *genJetPt = new vector<float>();
  vector<float> *genJetMass = new vector<float>();
  vector<float> *genJetEta = new vector<float>();
  vector<float> *genJetPhi = new vector<float>();
  vector<float> *genJetEnergy = new vector<float>();
  vector<float> *genJetNsub0p5tau1 = new vector<float>();
  vector<float> *genJetNsub0p5tau2 = new vector<float>();
  vector<float> *genJetNsub0p5tau3 = new vector<float>();
  vector<float> *genJetNsub0p5tau4 = new vector<float>();
  vector<float> *genJetNsub0p5tau5 = new vector<float>();
  vector<float> *genJetNsub0p5tau6 = new vector<float>();
  vector<float> *genJetNsub0p5tau7 = new vector<float>();
  // vector<float> *genJetNsub0p5tau8 = new vector<float>();
  // vector<float> *genJetNsub0p5tau9 = new vector<float>();
  // vector<float> *genJetNsub0p5tau10 = new vector<float>();
  vector<float> *genJetNsub1tau1 = new vector<float>();
  vector<float> *genJetNsub1tau2 = new vector<float>();
  vector<float> *genJetNsub1tau3 = new vector<float>();
  vector<float> *genJetNsub1tau4 = new vector<float>();
  vector<float> *genJetNsub1tau5 = new vector<float>();
  vector<float> *genJetNsub1tau6 = new vector<float>();
  vector<float> *genJetNsub1tau7 = new vector<float>();
  // vector<float> *genJetNsub1tau8 = new vector<float>();
  // vector<float> *genJetNsub1tau9 = new vector<float>();
  // vector<float> *genJetNsub1tau10 = new vector<float>();
  vector<float> *genJetNsub2tau1 = new vector<float>();
  vector<float> *genJetNsub2tau2 = new vector<float>();
  vector<float> *genJetNsub2tau3 = new vector<float>();
  vector<float> *genJetNsub2tau4 = new vector<float>();
  vector<float> *genJetNsub2tau5 = new vector<float>();
  vector<float> *genJetNsub2tau6 = new vector<float>();
  vector<float> *genJetNsub2tau7 = new vector<float>();
  // vector<float> *genJetNsub2tau8 = new vector<float>();
  // vector<float> *genJetNsub2tau9 = new vector<float>();
  // vector<float> *genJetNsub2tau10 = new vector<float>();

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
Skimmer::Skimmer(const edm::ParameterSet& iConfig):
  //srcGenJets_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjets"))),
  srcGenJet_ (consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("srcGenJet"))),
  srcGenjetNSub0p5(iConfig.getParameter<std::string>("GenjetNsub0p5")),
  srcGenjetNSub1(iConfig.getParameter<std::string>("GenjetNsub1")),
  srcGenjetNSub2(iConfig.getParameter<std::string>("GenjetNsub2")),
  Nsub0p5tau1(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau1",""))),
  Nsub0p5tau2(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau2",""))),
  Nsub0p5tau3(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau3",""))),
  Nsub0p5tau4(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau4",""))),
  Nsub0p5tau5(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau5",""))),
  Nsub0p5tau6(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau6",""))),
  Nsub0p5tau7(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau7",""))),
  // Nsub0p5tau8(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau8",""))),
  // Nsub0p5tau9(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau9",""))),
  // Nsub0p5tau10(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub0p5,"tau10",""))),
  Nsub1tau1(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau1",""))),
  Nsub1tau2(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau2",""))),
  Nsub1tau3(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau3",""))),
  Nsub1tau4(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau4",""))),
  Nsub1tau5(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau5",""))),
  Nsub1tau6(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau6",""))),
  Nsub1tau7(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau7",""))),
  // Nsub1tau8(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau8",""))),
  // Nsub1tau9(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau9",""))),
  // Nsub1tau10(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub1,"tau10",""))),
  Nsub2tau1(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau1",""))),
  Nsub2tau2(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau2",""))),
  Nsub2tau3(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau3",""))),
  Nsub2tau4(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau4",""))),
  Nsub2tau5(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau5",""))),
  Nsub2tau6(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau6",""))),
  Nsub2tau7(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau7","")))
  // Nsub2tau8(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau8",""))),
  // Nsub2tau9(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau9",""))),
  // Nsub2tau10(consumes<edm::ValueMap<float>>(edm::InputTag(srcGenjetNSub2,"tau10","")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");

}


Skimmer::~Skimmer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void Skimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<edm::View<reco::GenJet>> genJets;
  iEvent.getByToken(srcGenJet_, genJets);
  //iEvent.getByToken(srcGenJets_, genJets);

  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau1 ;
  iEvent.getByToken( Nsub0p5tau1 , genjet_Nsub0p5tau1 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau2 ;
  iEvent.getByToken( Nsub0p5tau2 , genjet_Nsub0p5tau2 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau3 ;
  iEvent.getByToken( Nsub0p5tau3 , genjet_Nsub0p5tau3 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau4 ;
  iEvent.getByToken( Nsub0p5tau4 , genjet_Nsub0p5tau4 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau5 ;
  iEvent.getByToken( Nsub0p5tau5 , genjet_Nsub0p5tau5 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau6 ;
  iEvent.getByToken( Nsub0p5tau6 , genjet_Nsub0p5tau6 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau7 ;
  iEvent.getByToken( Nsub0p5tau7 , genjet_Nsub0p5tau7 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau8 ;
  // iEvent.getByToken( Nsub0p5tau8 , genjet_Nsub0p5tau8 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau9 ;
  // iEvent.getByToken( Nsub0p5tau9 , genjet_Nsub0p5tau9 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub0p5tau10 ;
  // iEvent.getByToken( Nsub0p5tau10 , genjet_Nsub0p5tau10 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau1 ;
  iEvent.getByToken( Nsub1tau1 , genjet_Nsub1tau1 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau2 ;
  iEvent.getByToken( Nsub1tau2 , genjet_Nsub1tau2 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau3 ;
  iEvent.getByToken( Nsub1tau3 , genjet_Nsub1tau3 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau4 ;
  iEvent.getByToken( Nsub1tau4 , genjet_Nsub1tau4 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau5 ;
  iEvent.getByToken( Nsub1tau5 , genjet_Nsub1tau5 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau6 ;
  iEvent.getByToken( Nsub1tau6 , genjet_Nsub1tau6 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau7 ;
  iEvent.getByToken( Nsub1tau7 , genjet_Nsub1tau7 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau8 ;
  // iEvent.getByToken( Nsub1tau8 , genjet_Nsub1tau8 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau9 ;
  // iEvent.getByToken( Nsub1tau9 , genjet_Nsub1tau9 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub1tau10 ;
  // iEvent.getByToken( Nsub1tau10 , genjet_Nsub1tau10 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau1 ;
  iEvent.getByToken( Nsub2tau1 , genjet_Nsub2tau1 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau2 ;
  iEvent.getByToken( Nsub2tau2 , genjet_Nsub2tau2 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau3 ;
  iEvent.getByToken( Nsub2tau3 , genjet_Nsub2tau3 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau4 ;
  iEvent.getByToken( Nsub2tau4 , genjet_Nsub2tau4 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau5 ;
  iEvent.getByToken( Nsub2tau5 , genjet_Nsub2tau5 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau6 ;
  iEvent.getByToken( Nsub2tau6 , genjet_Nsub2tau6 );
  edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau7 ;
  iEvent.getByToken( Nsub2tau7 , genjet_Nsub2tau7 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau8 ;
  // iEvent.getByToken( Nsub2tau8 , genjet_Nsub2tau8 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau9 ;
  // iEvent.getByToken( Nsub2tau9 , genjet_Nsub2tau9 );
  // edm::Handle< edm::ValueMap<float> > genjet_Nsub2tau10 ;
  // iEvent.getByToken( Nsub2tau10 , genjet_Nsub2tau10 );
  
  genJetPt->clear();
  genJetMass->clear();
  genJetEta->clear();
  genJetPhi->clear();
  genJetEnergy->clear();
  genJetNsub0p5tau1->clear();
  genJetNsub0p5tau2->clear();
  genJetNsub0p5tau3->clear();
  genJetNsub0p5tau4->clear();
  genJetNsub0p5tau5->clear();
  genJetNsub0p5tau6->clear();
  genJetNsub0p5tau7->clear();
  genJetNsub1tau1->clear();
  genJetNsub1tau2->clear();
  genJetNsub1tau3->clear();
  genJetNsub1tau4->clear();
  genJetNsub1tau5->clear();
  genJetNsub1tau6->clear();
  genJetNsub1tau7->clear();
  genJetNsub2tau1->clear();
  genJetNsub2tau2->clear();
  genJetNsub2tau3->clear();
  genJetNsub2tau4->clear();
  genJetNsub2tau5->clear();
  genJetNsub2tau6->clear();
  genJetNsub2tau7->clear();


  for (size_t iGenJet = 0; iGenJet < genJets -> size(); iGenJet++) {
	
    if (iGenJet>0)
      break;
    const reco::GenJet & genjet = genJets->at(iGenJet) ;

    const edm::Ptr<reco::GenJet> genjet_ptr = genJets -> ptrAt( iGenJet ) ;
    LogInfo("test") << genjet.pt() 
		    << " " << (*genjet_Nsub0p5tau1)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau2)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau3)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau4)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau5)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau6)[ genjet_ptr ]
		    << " " << (*genjet_Nsub0p5tau7)[ genjet_ptr ]

      ;


    genJetPt->push_back( genjet.pt() );
    genJetMass->push_back( genjet.mass() );
    genJetEta->push_back( genjet.eta() );
    genJetPhi->push_back( genjet.phi() );
    genJetEnergy->push_back( genjet.energy() );
    genJetNsub0p5tau1->push_back( (*genjet_Nsub0p5tau1)[ genjet_ptr ] );
    genJetNsub0p5tau2->push_back( (*genjet_Nsub0p5tau2)[ genjet_ptr ] );
    genJetNsub0p5tau3->push_back( (*genjet_Nsub0p5tau3)[ genjet_ptr ] );
    genJetNsub0p5tau4->push_back( (*genjet_Nsub0p5tau4)[ genjet_ptr ] );
    genJetNsub0p5tau5->push_back( (*genjet_Nsub0p5tau5)[ genjet_ptr ] );
    genJetNsub0p5tau6->push_back( (*genjet_Nsub0p5tau6)[ genjet_ptr ] );
    genJetNsub0p5tau7->push_back( (*genjet_Nsub0p5tau7)[ genjet_ptr ] );
    // genJetNsub0p5tau8->push_back( (*genjet_Nsub0p5tau8)[ genjet_ptr ] );
    // genJetNsub0p5tau9->push_back( (*genjet_Nsub0p5tau9)[ genjet_ptr ] );
    // genJetNsub0p5tau10->push_back( (*genjet_Nsub0p5tau10)[ genjet_ptr ] );
    genJetNsub1tau1->push_back( (*genjet_Nsub1tau1)[ genjet_ptr ] );
    genJetNsub1tau2->push_back( (*genjet_Nsub1tau2)[ genjet_ptr ] );
    genJetNsub1tau3->push_back( (*genjet_Nsub1tau3)[ genjet_ptr ] );
    genJetNsub1tau4->push_back( (*genjet_Nsub1tau4)[ genjet_ptr ] );
    genJetNsub1tau5->push_back( (*genjet_Nsub1tau5)[ genjet_ptr ] );
    genJetNsub1tau6->push_back( (*genjet_Nsub1tau6)[ genjet_ptr ] );
    genJetNsub1tau7->push_back( (*genjet_Nsub1tau7)[ genjet_ptr ] );
    // genJetNsub1tau8->push_back( (*genjet_Nsub1tau8)[ genjet_ptr ] );
    // genJetNsub1tau9->push_back( (*genjet_Nsub1tau9)[ genjet_ptr ] );
    // genJetNsub1tau10->push_back( (*genjet_Nsub1tau10)[ genjet_ptr ] );
    genJetNsub2tau1->push_back( (*genjet_Nsub2tau1)[ genjet_ptr ] );
    genJetNsub2tau2->push_back( (*genjet_Nsub2tau2)[ genjet_ptr ] );
    genJetNsub2tau3->push_back( (*genjet_Nsub2tau3)[ genjet_ptr ] );
    genJetNsub2tau4->push_back( (*genjet_Nsub2tau4)[ genjet_ptr ] );
    genJetNsub2tau5->push_back( (*genjet_Nsub2tau5)[ genjet_ptr ] );
    genJetNsub2tau6->push_back( (*genjet_Nsub2tau6)[ genjet_ptr ] );
    genJetNsub2tau7->push_back( (*genjet_Nsub2tau7)[ genjet_ptr ] );
    // genJetNsub2tau8->push_back( (*genjet_Nsub2tau8)[ genjet_ptr ] );
    // genJetNsub2tau9->push_back( (*genjet_Nsub2tau9)[ genjet_ptr ] );
    // genJetNsub2tau10->push_back( (*genjet_Nsub2tau10)[ genjet_ptr ] );
  }

  SimpleGentree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void Skimmer::beginJob() {

  edm::Service< TFileService > fileService;

  SimpleGentree = fileService->make< TTree >("SimpleGentree", "SimpleGentree");
  SimpleGentree->Branch( "genJetPt", "vector<float>", &genJetPt);
  SimpleGentree->Branch( "genJetMass", "vector<float>", &genJetMass);
  SimpleGentree->Branch( "genJetEta", "vector<float>", &genJetEta);
  SimpleGentree->Branch( "genJetPhi", "vector<float>", &genJetPhi);
  SimpleGentree->Branch( "genJetEnergy", "vector<float>", &genJetEnergy);
  SimpleGentree->Branch( "genJetNsub0p5tau1", "vector<float>", &genJetNsub0p5tau1);
  SimpleGentree->Branch( "genJetNsub0p5tau2", "vector<float>", &genJetNsub0p5tau2);
  SimpleGentree->Branch( "genJetNsub0p5tau3", "vector<float>", &genJetNsub0p5tau3);
  SimpleGentree->Branch( "genJetNsub0p5tau4", "vector<float>", &genJetNsub0p5tau4);
  SimpleGentree->Branch( "genJetNsub0p5tau5", "vector<float>", &genJetNsub0p5tau5);
  SimpleGentree->Branch( "genJetNsub0p5tau6", "vector<float>", &genJetNsub0p5tau6);
  SimpleGentree->Branch( "genJetNsub0p5tau7", "vector<float>", &genJetNsub0p5tau7);
  // SimpleGentree->Branch( "genJetNsub0p5tau8", "vector<float>", &genJetNsub0p5tau8);
  // SimpleGentree->Branch( "genJetNsub0p5tau9", "vector<float>", &genJetNsub0p5tau9);
  // SimpleGentree->Branch( "genJetNsub0p5tau10", "vector<float>", &genJetNsub0p5tau10);
  SimpleGentree->Branch( "genJetNsub1tau1", "vector<float>", &genJetNsub1tau1);
  SimpleGentree->Branch( "genJetNsub1tau2", "vector<float>", &genJetNsub1tau2);
  SimpleGentree->Branch( "genJetNsub1tau3", "vector<float>", &genJetNsub1tau3);
  SimpleGentree->Branch( "genJetNsub1tau4", "vector<float>", &genJetNsub1tau4);
  SimpleGentree->Branch( "genJetNsub1tau5", "vector<float>", &genJetNsub1tau5);
  SimpleGentree->Branch( "genJetNsub1tau6", "vector<float>", &genJetNsub1tau6);
  SimpleGentree->Branch( "genJetNsub1tau7", "vector<float>", &genJetNsub1tau7);
  // SimpleGentree->Branch( "genJetNsub1tau8", "vector<float>", &genJetNsub1tau8);
  // SimpleGentree->Branch( "genJetNsub1tau9", "vector<float>", &genJetNsub1tau9);
  // SimpleGentree->Branch( "genJetNsub1tau10", "vector<float>", &genJetNsub1tau10);
  SimpleGentree->Branch( "genJetNsub2tau1", "vector<float>", &genJetNsub2tau1);
  SimpleGentree->Branch( "genJetNsub2tau2", "vector<float>", &genJetNsub2tau2);
  SimpleGentree->Branch( "genJetNsub2tau3", "vector<float>", &genJetNsub2tau3);
  SimpleGentree->Branch( "genJetNsub2tau4", "vector<float>", &genJetNsub2tau4);
  SimpleGentree->Branch( "genJetNsub2tau5", "vector<float>", &genJetNsub2tau5);
  SimpleGentree->Branch( "genJetNsub2tau6", "vector<float>", &genJetNsub2tau6);
  SimpleGentree->Branch( "genJetNsub2tau7", "vector<float>", &genJetNsub2tau7);
  // SimpleGentree->Branch( "genJetNsub2tau8", "vector<float>", &genJetNsub2tau8);
  // SimpleGentree->Branch( "genJetNsub2tau9", "vector<float>", &genJetNsub2tau9);
  // SimpleGentree->Branch( "genJetNsub2tau10", "vector<float>", &genJetNsub2tau10);

}

// ------------ method called once each job just after ending the event loop  ------------
void Skimmer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Skimmer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Skimmer);

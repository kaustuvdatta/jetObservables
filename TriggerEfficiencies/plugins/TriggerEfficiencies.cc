// -*- C++ -*-
//
// Package:    RUNA/RUNTriggerEfficiency
// Class:      TriggerEfficiencies
// Original Author:  alejandro gomez
//         Created:  Tue, 14 Oct 2014 23:13:13 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "jetObservables/TriggerEfficiencies/plugins/CommonVariablesStructure.h"

using namespace edm;
using namespace std;

//
// constants, enums and typedefs
//

//
// class declaration
//
class TriggerEfficiencies : public EDAnalyzer {
	public:
		explicit TriggerEfficiencies(const ParameterSet&);
		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
		~TriggerEfficiencies();

	private:
		virtual void beginJob() override;
		virtual void analyze(const Event&, const EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(Run const&, EventSetup const&) override;
		//virtual void endRun(Run const&, EventSetup const&) override;
		//virtual void beginLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;
		//virtual void endLuminosityBlock(LuminosityBlock const&, EventSetup const&) override;

		// ----------member data ---------------------------
		Service<TFileService> fs_;
		TTree *PrescaleTree;
		map< string, TH1D* > histos1D_;
		map< string, TH2D* > histos2D_;
		vector< string > cutLabels;

		double cutAK8jetPt;
		double cutAK8jet1Pt;
		double cutAK8jet2Pt;
		double cutAK8jet1Mass;
		TString baseTrigger;
		vector<string> listOfTriggers;
		vector<string> triggerThresholds;

        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
        edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
		edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
		edm::EDGetTokenT<trigger::TriggerEvent> triggerEvent_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;

        unsigned int lumi = 0, run=0;
        ULong64_t event = 0;
        //int preAK8PFJet80 = 0, preAK8PFJet140 = 0, preAK8PFJet200 = 0, preAK8PFJet260 = 0, preAK8PFJet320 = 0, preAK8PFJet400 = 0, preAK8PFJet450 = 0, preAK8PFJet500 = 0, preAK8PFJet550 = 0; 
        vector<bool> triggerDesicion; 
        vector<int> triggerPrescale;
};

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerEfficiencies::TriggerEfficiencies(const ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
	triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
	triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("hltTrigger"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("recoJets")))
{
	cutAK8jetPt = iConfig.getParameter<double>("cutAK8jetPt");
	cutAK8jet1Pt = iConfig.getParameter<double>("cutAK8jet1Pt");
	cutAK8jet2Pt = iConfig.getParameter<double>("cutAK8jet2Pt");
	cutAK8jet1Mass = iConfig.getParameter<double>("cutAK8jet1Mass");
	baseTrigger = iConfig.getParameter<string>("baseTrigger");
	listOfTriggers = iConfig.getParameter<vector<string>>("listOfTriggers");
	triggerThresholds = iConfig.getParameter<vector<string>>("triggerThresholds");
}


TriggerEfficiencies::~TriggerEfficiencies()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TriggerEfficiencies::analyze(const Event& iEvent, const EventSetup& iSetup) {

    triggerDesicion.clear();
    triggerPrescale.clear();

    edm::Handle<reco::VertexCollection> vertices;
	edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
	edm::Handle<trigger::TriggerEvent> trigEvent; 
	edm::Handle<pat::JetCollection> jets;

    iEvent.getByToken(vtxToken_, vertices);
	iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
	iEvent.getByToken(triggerPrescales_, triggerPrescales);
	iEvent.getByToken(triggerEvent_,trigEvent);
	iEvent.getByToken(jetToken_, jets);

    if (vertices->empty()) return; // skip the event if no PV found
    //const reco::Vertex &PV = vertices->front();

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    //bool ORTriggers = checkORListOfTriggerBitsMiniAOD( names, triggerBits, triggerPrescales, listOfTriggers, false );
    //bool basedTriggerFired = checkTriggerBitsMiniAOD( names, triggerBits, triggerPrescales, baseTrigger, true );

    for (size_t t = 0; t < listOfTriggers.size(); t++) {
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if (TString(names.triggerName(i)).Contains(listOfTriggers[t])) { 
                if (triggerBits->accept(i)) {
                    triggerDesicion.push_back( 1 );
                    triggerPrescale.push_back( triggerPrescales->getPrescaleForIndex(i) );
                    //LogWarning("triggerbit") << names.triggerName(i) << " " <<  hltAK8PFJet80_bool << " " << hltAK8PFJet80_prescale;
                } else {
                    triggerDesicion.push_back( 0 );
                    triggerPrescale.push_back( 0 );
                }
            }
        }
    }
    event = iEvent.id().event();
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    PrescaleTree->Fill();

	/// Applying kinematic, trigger and jet ID
	pat::JetCollection JETS;
    //vector<pat::Jet> JETS;

	for (const pat::Jet &jet : *jets) {

		if ( TMath::Abs( jet.eta() ) > 2.4 ) continue;
		//string typeOfJetID = "looseJetID";	// check trigger with looser jet id
		//bool idL = jetID( jet.eta(), jet.energy(), jet.jecFactor(0), jet.neutralHadronEnergyFraction(), jet.neutralEmEnergyFraction(), jet.chargedHadronEnergyFraction(), jet.muonEnergy(), jet.chargedEmEnergyFraction(), jet.chargedMultiplicity(), jet.neutralMultiplicity(), typeOfJetID ); 

		if( ( jet.pt() > cutAK8jetPt ) ) { 
			JETS.push_back( jet );
		}
	}


    if ( JETS.size()>1 ) {

        TLorentzVector hltJet, recoJet, recoJet2;
        recoJet.SetPtEtaPhiM( JETS[0].pt(), JETS[0].eta(), JETS[0].phi(), JETS[0].mass() );
        recoJet2.SetPtEtaPhiM( JETS[1].pt(), JETS[1].eta(), JETS[1].phi(), JETS[1].mass() );
        double deltaPhi = recoJet.DeltaPhi(recoJet2);
        double ptAsym = TMath::Abs( recoJet.Pt() - recoJet2.Pt() )/( recoJet.Pt() + recoJet2.Pt() );

        if ( ( deltaPhi > 2 ) and ( ptAsym < 0.3 ) ){

            if ( ( listOfTriggers.size() == triggerDesicion.size() ) and ( triggerDesicion.size() == triggerPrescale.size() ) ){
                for (size_t t = 0; t < listOfTriggers.size(); t++) {
                    if (triggerDesicion[t]==1) {
                        auto tmp = "jet1Pt_" + listOfTriggers[t] ;
                        histos1D_[ tmp ]->Fill( JETS[0].pt() );
                        histos1D_[ tmp+"_scaled" ]->Fill( JETS[0].pt(), triggerPrescale[t] );
                        auto tmp2 = "NPV_" + listOfTriggers[t] ;
                        histos1D_[ tmp2 ]->Fill( vertices->size() );
                    }
                }
            }


            double hltPt = 0;
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
                obj.unpackPathNames(names);
                obj.unpackFilterLabels(iEvent, *triggerBits );
                for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
                    TString filterLabel = obj.filterLabels()[h];
                    for (size_t t = 0; t < listOfTriggers.size()-1; t++) {
                        //cout << filterLabel << endl;
                        TString hltObjectPt_ = "hltSinglePFJet"+ triggerThresholds[t]  +"AK8";
                        if ( filterLabel.Contains( hltObjectPt_ ) ) {
                            hltPt = obj.pt();
                            hltJet.SetPtEtaPhiM( obj.pt(), obj.eta(), obj.phi(), obj.mass()  );
                                /*cout << "\tTrigger object Pt:  pt " << obj.pt()
                                                    << ", eta " << obj.eta()
                                                    << ", phi " << obj.phi()
                                                    << ", mass " << obj.mass() << endl;*/
                            if ( ( triggerDesicion[t]==1 ) and ( hltJet.DeltaR( recoJet )<0.8 ) ){
                                auto tmpName = "jet1Pt_AK8PFJet";
                                if ( hltPt>stod(triggerThresholds[t]) ) histos1D_[ tmpName+triggerThresholds[t] +"_only" ]->Fill( JETS[0].pt() );
                                if ( hltPt>stod(triggerThresholds[t+1]) ) histos1D_[ tmpName+triggerThresholds[t+1] +"_simulated" ]->Fill( JETS[0].pt() );
                            }
                        }
                    }
                }
            }
                
        }
    }

}


// ------------ method called once each job just before starting event loop  ------------
void TriggerEfficiencies::beginJob() {

    histos1D_[ "NPV_HLT_AK8PFJet60" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet60", "NPV_HLT_AK8PFJet60", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet80" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet80", "NPV_HLT_AK8PFJet80", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet140" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet140", "NPV_HLT_AK8PFJet140", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet200" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet200", "NPV_HLT_AK8PFJet200", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet260" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet260", "NPV_HLT_AK8PFJet260", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet320" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet320", "NPV_HLT_AK8PFJet320", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet400" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet400", "NPV_HLT_AK8PFJet400", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet450" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet450", "NPV_HLT_AK8PFJet450", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet500" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet500", "NPV_HLT_AK8PFJet500", 100, 0., 100. );
    histos1D_[ "NPV_HLT_AK8PFJet550" ] = fs_->make< TH1D >( "NPV_HLT_AK8PFJet550", "NPV_HLT_AK8PFJet550", 100, 0., 100. );

    histos1D_[ "jet1Pt_HLT_AK8PFJet60" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet60", "jet1Pt_HLT_AK8PFJet60", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet80", "jet1Pt_HLT_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet140" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet140", "jet1Pt_HLT_AK8PFJet140", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet200" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet200", "jet1Pt_HLT_AK8PFJet200", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet260" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet260", "jet1Pt_HLT_AK8PFJet260", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet320" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet320", "jet1Pt_HLT_AK8PFJet320", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet400" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet400", "jet1Pt_HLT_AK8PFJet400", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet450" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet450", "jet1Pt_HLT_AK8PFJet450", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet500" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet500", "jet1Pt_HLT_AK8PFJet500", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet550" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet550", "jet1Pt_HLT_AK8PFJet550", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet60_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet60_scaled", "jet1Pt_HLT_AK8PFJet60_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet140_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet140_scaled", "jet1Pt_HLT_AK8PFJet140_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet200_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet200_scaled", "jet1Pt_HLT_AK8PFJet200_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet260_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet260_scaled", "jet1Pt_HLT_AK8PFJet260_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet320_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet320_scaled", "jet1Pt_HLT_AK8PFJet320_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet400_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet400_scaled", "jet1Pt_HLT_AK8PFJet400_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet450_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet450_scaled", "jet1Pt_HLT_AK8PFJet450_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet500_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet500_scaled", "jet1Pt_HLT_AK8PFJet500_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet550_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet550_scaled", "jet1Pt_HLT_AK8PFJet550_scaled", 150, 0., 1500. );


    histos1D_[ "jet1Pt_AK8PFJet60_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet60_only", "jet1Pt_AK8PFJet60_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet80_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet80_only", "jet1Pt_AK8PFJet80_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet140_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet140_only", "jet1Pt_AK8PFJet140_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet200_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet200_only", "jet1Pt_AK8PFJet200_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet260_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet260_only", "jet1Pt_AK8PFJet260_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet320_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet320_only", "jet1Pt_AK8PFJet320_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet400_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet400_only", "jet1Pt_AK8PFJet400_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet450_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet450_only", "jet1Pt_AK8PFJet450_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet500_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet500_only", "jet1Pt_AK8PFJet500_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet550_only" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet550_only", "jet1Pt_AK8PFJet550_only", 150, 0., 1500. );

    histos1D_[ "jet1Pt_AK8PFJet60_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet60_simulated", "jet1Pt_AK8PFJet60_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet80_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet80_simulated", "jet1Pt_AK8PFJet80_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet140_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet140_simulated", "jet1Pt_AK8PFJet140_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet200_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet200_simulated", "jet1Pt_AK8PFJet200_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet260_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet260_simulated", "jet1Pt_AK8PFJet260_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet320_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet320_simulated", "jet1Pt_AK8PFJet320_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet400_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet400_simulated", "jet1Pt_AK8PFJet400_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet450_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet450_simulated", "jet1Pt_AK8PFJet450_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet500_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet500_simulated", "jet1Pt_AK8PFJet500_simulated", 150, 0., 1500. );
    histos1D_[ "jet1Pt_AK8PFJet550_simulated" ] = fs_->make< TH1D >( "jet1Pt_AK8PFJet550_simulated", "jet1Pt_AK8PFJet550_simulated", 150, 0., 1500. );



	
	///// Sumw2 all the histos
	for( auto const& histo : histos1D_ ) histos1D_[ histo.first ]->Sumw2();
	for( auto const& histo : histos2D_ ) histos2D_[ histo.first ]->Sumw2();

    PrescaleTree = fs_->make< TTree >("TriggerObjects", "TriggerObjects");
    PrescaleTree->Branch( "run", &run, "run/I" );
    PrescaleTree->Branch( "lumi", &lumi, "lumi/I" );
    PrescaleTree->Branch( "event", &event, "event/I" );
    PrescaleTree->Branch( "triggerPrescales", "vector<int>", &triggerPrescale);


}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerEfficiencies::endJob() {

}

void TriggerEfficiencies::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

	edm::ParameterSetDescription desc;

	desc.add<double>("cutAK8jetPt", 150.);
	desc.add<double>("cutAK8jet1Pt", 500.);
	desc.add<double>("cutAK8jet2Pt", 450.);
	desc.add<double>("cutAK8jet1Mass", 60.);
	desc.add<string>("baseTrigger", "HLT_AK8PFJet200");
	desc.add<InputTag>("vertices", 	InputTag("offlineSlimmedPrimaryVertices"));
	desc.add<InputTag>("bits", 	InputTag("TriggerResults", "", "HLT"));
	desc.add<InputTag>("prescales", 	InputTag("patTrigger"));
	desc.add<InputTag>("objects", 	InputTag("slimmedPatTrigger"));
	desc.add<InputTag>("hltTrigger", 	InputTag("hltTriggerSummaryAOD","","HLT"));
	desc.add<InputTag>("recoJets", 	InputTag("slimmedJetsAK8"));
	vector<string> HLTPass;
	HLTPass.push_back("HLT_AK8PFJet400");
	desc.add<vector<string>>("listOfTriggers",	HLTPass);
	vector<string> thresholds;
	thresholds.push_back("80");
	desc.add<vector<string>>("triggerThresholds",	thresholds);

	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEfficiencies);

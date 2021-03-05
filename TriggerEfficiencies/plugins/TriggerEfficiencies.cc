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
		vector<double> triggerThresholds;

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
	triggerThresholds = iConfig.getParameter<vector<double>>("triggerThresholds");
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

	edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
	edm::Handle<trigger::TriggerEvent> trigEvent; 
	edm::Handle<pat::JetCollection> jets;

	iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
	iEvent.getByToken(triggerPrescales_, triggerPrescales);
	iEvent.getByToken(triggerEvent_,trigEvent);
	iEvent.getByToken(jetToken_, jets);

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


    if ( JETS.size()>0 ) {

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


        if ( ( listOfTriggers.size() == triggerDesicion.size() ) and ( triggerDesicion.size() == triggerPrescale.size() ) ){
            for (size_t t = 0; t < listOfTriggers.size(); t++) {
                if (triggerDesicion[t]==1) {
                    //LogWarning("pass") << listOfTriggers[t] << " " << triggerDesicion[t] << " " << triggerPrescale[t];
                    auto tmp = "jet1Pt_" + listOfTriggers[t] ;
                    histos1D_[ tmp+ "_only" ]->Fill( JETS[0].pt() );
                    if ( (triggerDesicion[0]==1) ){
                        histos1D_[ tmp+ "_AK8PFJet80" ]->Fill( JETS[0].pt() );
                        histos1D_[ tmp+ "_AK8PFJet80_scaled" ]->Fill( JETS[0].pt(), triggerPrescale[t] );
                    }
                }
            }
        }

       event = iEvent.id().event();
       run = iEvent.id().run();
       lumi = iEvent.luminosityBlock();
       PrescaleTree->Fill();


       TLorentzVector hltJet, recoJet;
       recoJet.SetPtEtaPhiM( JETS[0].pt(), JETS[0].eta(), JETS[0].phi(), JETS[0].mass() );
        if (triggerDesicion[0]==1){
            double hltPt = 0;
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
                obj.unpackPathNames(names);
                obj.unpackFilterLabels(iEvent, *triggerBits );
                for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
                    TString filterLabel = obj.filterLabels()[h];
                    //cout << filterLabel << endl;
                    TString hltObjectPt_ = "hltSinglePFJet80AK8";
                    if ( filterLabel.Contains( hltObjectPt_ ) ) {
                        hltPt = obj.pt();
                        hltJet.SetPtEtaPhiM( obj.pt(), obj.eta(), obj.phi(), obj.mass()  );
                        /*cout << "\tTrigger object Pt:  pt " << obj.pt()
                                            << ", eta " << obj.eta()
                                            << ", phi " << obj.phi()
                                            << ", mass " << obj.mass() << endl;*/
                    }
                }
            }
            if ( hltJet.DeltaR( recoJet )<0.4 ){ 
            if ( hltPt>80 ) histos1D_[ "jet1Pt_AK8PFJet80_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>140 ) histos1D_[ "jet1Pt_AK8PFJet140_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>200 ) histos1D_[ "jet1Pt_AK8PFJet200_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>260 ) histos1D_[ "jet1Pt_AK8PFJet260_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>320 ) histos1D_[ "jet1Pt_AK8PFJet320_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>400 ) histos1D_[ "jet1Pt_AK8PFJet400_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>450 ) histos1D_[ "jet1Pt_AK8PFJet450_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>500 ) histos1D_[ "jet1Pt_AK8PFJet500_simulated" ]->Fill( JETS[0].pt() );
            if ( hltPt>550 ) histos1D_[ "jet1Pt_AK8PFJet550_simulated" ]->Fill( JETS[0].pt() );

            }

        }

    }


            //for (size_t t = 0; t < triggerDesicion.size(); t++) {
            //    cout << t << " " << triggerDesicion[t] << " " << triggerPrescale[t] << endl;
           // }


}


// ------------ method called once each job just before starting event loop  ------------
void TriggerEfficiencies::beginJob() {


    histos1D_[ "jet1Pt_HLT_AK8PFJet80_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet80_only", "jet1Pt_HLT_AK8PFJet80_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet140_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet140_only", "jet1Pt_HLT_AK8PFJet140_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet200_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet200_only", "jet1Pt_HLT_AK8PFJet200_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet260_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet260_only", "jet1Pt_HLT_AK8PFJet260_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet320_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet320_only", "jet1Pt_HLT_AK8PFJet320_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet400_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet400_only", "jet1Pt_HLT_AK8PFJet400_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet450_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet450_only", "jet1Pt_HLT_AK8PFJet450_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet500_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet500_only", "jet1Pt_HLT_AK8PFJet500_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet550_only" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet550_only", "jet1Pt_HLT_AK8PFJet550_only", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet80_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet80_AK8PFJet80", "jet1Pt_HLT_AK8PFJet80_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet140_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet140_AK8PFJet80", "jet1Pt_HLT_AK8PFJet140_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet200_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet200_AK8PFJet80", "jet1Pt_HLT_AK8PFJet200_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet260_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet260_AK8PFJet80", "jet1Pt_HLT_AK8PFJet260_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet320_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet320_AK8PFJet80", "jet1Pt_HLT_AK8PFJet320_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet400_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet400_AK8PFJet80", "jet1Pt_HLT_AK8PFJet400_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet450_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet450_AK8PFJet80", "jet1Pt_HLT_AK8PFJet450_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet500_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet500_AK8PFJet80", "jet1Pt_HLT_AK8PFJet500_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet550_AK8PFJet80" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet550_AK8PFJet80", "jet1Pt_HLT_AK8PFJet550_AK8PFJet80", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet80_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet80_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet80_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet140_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet140_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet140_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet200_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet200_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet200_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet260_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet260_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet260_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet320_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet320_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet320_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet400_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet400_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet400_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet450_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet450_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet450_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet500_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet500_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet500_AK8PFJet80_scaled", 150, 0., 1500. );
    histos1D_[ "jet1Pt_HLT_AK8PFJet550_AK8PFJet80_scaled" ] = fs_->make< TH1D >( "jet1Pt_HLT_AK8PFJet550_AK8PFJet80_scaled", "jet1Pt_HLT_AK8PFJet550_AK8PFJet80_scaled", 150, 0., 1500. );


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
	desc.add<InputTag>("bits", 	InputTag("TriggerResults", "", "HLT"));
	desc.add<InputTag>("prescales", 	InputTag("patTrigger"));
	desc.add<InputTag>("objects", 	InputTag("slimmedPatTrigger"));
	desc.add<InputTag>("hltTrigger", 	InputTag("hltTriggerSummaryAOD","","HLT"));
	desc.add<InputTag>("recoJets", 	InputTag("slimmedJetsAK8"));
	vector<string> HLTPass;
	HLTPass.push_back("HLT_AK8PFJet400");
	desc.add<vector<string>>("listOfTriggers",	HLTPass);
	vector<double> thresholds;
	thresholds.push_back(80.);
	desc.add<vector<double>>("triggerThresholds",	thresholds);

	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEfficiencies);

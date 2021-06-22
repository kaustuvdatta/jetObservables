// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

namespace Rivet {



    // Nsubjetiness measurement
    class CMS_2021_PAS_SMP_21_XXX: public Analysis {
        public:

            /// @name Constructors etc.
            //@{

            /// Constructor
            CMS_2021_PAS_SMP_21_XXX()
            : Analysis("CMS_2021_PAS_SMP_21_XXX"),
                _softdrop(fjcontrib::SoftDrop(0, 0.1, 0.8) ), // parameters are beta, zcut, R0
                _nsub_0p5_1(fjcontrib::Nsubjettiness(1,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(0.5, 0.8))),
                _nsub_0p5_2(fjcontrib::Nsubjettiness(2,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(0.5, 0.8))),
                _nsub_0p5_3(fjcontrib::Nsubjettiness(3,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(0.5, 0.8))),
                _nsub_0p5_4(fjcontrib::Nsubjettiness(4,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(0.5, 0.8))),
                _nsub_1_1(fjcontrib::Nsubjettiness(1,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8))),
                _nsub_1_2(fjcontrib::Nsubjettiness(2,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8))),
                _nsub_1_3(fjcontrib::Nsubjettiness(3,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8))),
                _nsub_1_4(fjcontrib::Nsubjettiness(4,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8))),
                _nsub_2_1(fjcontrib::Nsubjettiness(1,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(2., 0.8))),
                _nsub_2_2(fjcontrib::Nsubjettiness(2,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(2., 0.8))),
                _nsub_2_3(fjcontrib::Nsubjettiness(3,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(2., 0.8))),
                _nsub_2_4(fjcontrib::Nsubjettiness(4,   fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(2., 0.8))),
                _nsub21(fjcontrib::NsubjettinessRatio(2,1,   fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8))),
                _nsub32(fjcontrib::NsubjettinessRatio(3,2,   fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1., 0.8)))
            {    }

            //@}


            /// @name Analysis methods
            //@{

            /// Book histograms and initialise projections before the run
            void init() {
                // define a projection that keeps all the particles up to |eta|=5
                const FinalState fs(Cuts::abseta < 5.);

                // use FastJet, anti-kt(R=0.8) to do the clustering
                declare(FastJets(fs, FastJets::ANTIKT, 0.8), "JetsAK8");

                // Histograms
                book( _h_ungroomedJet0_nSub[0][1], "unfoldHistoJet2_tau_0p5_1" );
                book( _h_ungroomedJet1_nSub[0][1], "unfoldHistoJet1_tau_0p5_1" );
                book( _h_ungroomedJet0_nSub[0][2], "unfoldHistoJet2_tau_0p5_2" );
                book( _h_ungroomedJet1_nSub[0][2], "unfoldHistoJet1_tau_0p5_2" );
                book( _h_ungroomedJet0_nSub[0][3], "unfoldHistoJet2_tau_0p5_3" );
                book( _h_ungroomedJet1_nSub[0][3], "unfoldHistoJet1_tau_0p5_3" );
                book( _h_ungroomedJet0_nSub[0][4], "unfoldHistoJet2_tau_0p5_4" );
                book( _h_ungroomedJet1_nSub[0][4], "unfoldHistoJet1_tau_0p5_4" );
                book( _h_ungroomedJet0_nSub[1][1], "unfoldHistoJet2_tau_1_1" );
                book( _h_ungroomedJet1_nSub[1][1], "unfoldHistoJet1_tau_1_1" );
                book( _h_ungroomedJet0_nSub[1][2], "unfoldHistoJet2_tau_1_2" );
                book( _h_ungroomedJet1_nSub[1][2], "unfoldHistoJet1_tau_1_2" );
                book( _h_ungroomedJet0_nSub[1][3], "unfoldHistoJet2_tau_1_3" );
                book( _h_ungroomedJet1_nSub[1][3], "unfoldHistoJet1_tau_1_3" );
                book( _h_ungroomedJet0_nSub[1][4], "unfoldHistoJet2_tau_1_4" );
                book( _h_ungroomedJet1_nSub[1][4], "unfoldHistoJet1_tau_1_4" );
                book( _h_ungroomedJet0_nSub[2][1], "unfoldHistoJet2_tau_2_1" );
                book( _h_ungroomedJet1_nSub[2][1], "unfoldHistoJet1_tau_2_1" );
                book( _h_ungroomedJet0_nSub[2][2], "unfoldHistoJet2_tau_2_2" );
                book( _h_ungroomedJet1_nSub[2][2], "unfoldHistoJet1_tau_2_2" );
                book( _h_ungroomedJet0_nSub[2][3], "unfoldHistoJet2_tau_2_3" );
                book( _h_ungroomedJet1_nSub[2][3], "unfoldHistoJet1_tau_2_3" );
                book( _h_ungroomedJet0_nSub[2][4], "unfoldHistoJet2_tau_2_4" );
                book( _h_ungroomedJet1_nSub[2][4], "unfoldHistoJet1_tau_2_4" );
                book( _h_ungroomedJet0_nSub21, "unfoldHistoJet2_tau21" );
                book( _h_ungroomedJet1_nSub21, "unfoldHistoJet1_tau21" );
                book( _h_ungroomedJet0_nSub32, "unfoldHistoJet2_tau32" );
                book( _h_ungroomedJet1_nSub32, "unfoldHistoJet1_tau32" );
            }

            /// Perform the per-event analysis
            void analyze(const Event& event) {

                // Look at events with >= 2 jets
                auto jetsAK8 = applyProjection<FastJets>(event, "JetsAK8").jetsByPt(Cuts::pT > 200*GeV and Cuts::abseta < 2.5);
                if (jetsAK8.size() < 2) vetoEvent;

                // Get the leading two jets
                const fastjet::PseudoJet& j0 = jetsAK8[0].pseudojet();
                const fastjet::PseudoJet& j1 = jetsAK8[1].pseudojet();

                // Calculate delta phi and the pt asymmetry
                double deltaPhi = Rivet::deltaPhi( j0.phi(), j1.phi() );
                double ptasym = (j0.pt() - j1.pt()) / (j0.pt() + j1.pt());
                if (deltaPhi < 2.0 ) vetoEvent;
                if (ptasym > 0.3) vetoEvent;

                vector<PseudoJet> dijets = {j0, j1};
                std::sort(dijets.begin(), dijets.end(),
                        [] (const PseudoJet & A, const PseudoJet & B)
                        { return fabs(A.rapidity()) < fabs(B.rapidity()); }
                );

                // Now run the substructure algs...
                //fastjet::PseudoJet sd0 = _softdrop(j0);
                //fastjet::PseudoJet sd1 = _softdrop(j1);

                _h_ungroomedJet0_nSub[0][1]->fill( _nsub_0p5_1(dijets[0]) );
                _h_ungroomedJet0_nSub[0][2]->fill( _nsub_0p5_2(dijets[0]) );
                _h_ungroomedJet0_nSub[0][3]->fill( _nsub_0p5_3(dijets[0]) );
                _h_ungroomedJet0_nSub[0][4]->fill( _nsub_0p5_4(dijets[0]) );
                _h_ungroomedJet0_nSub[1][1]->fill( _nsub_1_1(dijets[0]) );
                _h_ungroomedJet0_nSub[1][2]->fill( _nsub_1_2(dijets[0]) );
                _h_ungroomedJet0_nSub[1][3]->fill( _nsub_1_3(dijets[0]) );
                _h_ungroomedJet0_nSub[1][4]->fill( _nsub_1_4(dijets[0]) );
                _h_ungroomedJet0_nSub[2][1]->fill( _nsub_2_1(dijets[0]) );
                _h_ungroomedJet0_nSub[2][2]->fill( _nsub_2_2(dijets[0]) );
                _h_ungroomedJet0_nSub[2][3]->fill( _nsub_2_3(dijets[0]) );
                _h_ungroomedJet0_nSub[2][4]->fill( _nsub_2_4(dijets[0]) );
                _h_ungroomedJet0_nSub21->fill( _nsub21(dijets[0]) );
                _h_ungroomedJet0_nSub32->fill( _nsub32(dijets[0]) );

                _h_ungroomedJet1_nSub[0][1]->fill( _nsub_0p5_1(dijets[1]) );
                _h_ungroomedJet1_nSub[0][2]->fill( _nsub_0p5_2(dijets[1]) );
                _h_ungroomedJet1_nSub[0][3]->fill( _nsub_0p5_3(dijets[1]) );
                _h_ungroomedJet1_nSub[0][4]->fill( _nsub_0p5_4(dijets[1]) );
                _h_ungroomedJet1_nSub[1][1]->fill( _nsub_1_1(dijets[1]) );
                _h_ungroomedJet1_nSub[1][2]->fill( _nsub_1_2(dijets[1]) );
                _h_ungroomedJet1_nSub[1][3]->fill( _nsub_1_3(dijets[1]) );
                _h_ungroomedJet1_nSub[1][4]->fill( _nsub_1_4(dijets[1]) );
                _h_ungroomedJet1_nSub[2][1]->fill( _nsub_2_1(dijets[1]) );
                _h_ungroomedJet1_nSub[2][2]->fill( _nsub_2_2(dijets[1]) );
                _h_ungroomedJet1_nSub[2][3]->fill( _nsub_2_3(dijets[1]) );
                _h_ungroomedJet1_nSub[2][4]->fill( _nsub_2_4(dijets[1]) );
                _h_ungroomedJet1_nSub21->fill( _nsub21(dijets[1]) );
                _h_ungroomedJet1_nSub32->fill( _nsub32(dijets[1]) );
            }


            /// Normalise histograms etc., after the run
            void finalize() {
                /*/ Normalize the normalized cross section histograms to unity,
                for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
                    normalize(_h_ungroomedJetMass_dj[i][1]);
                    normalize(_h_sdJetMass_dj[i][1]);
                }
                // Normalize the absolute cross section histograms to xs * lumi.
                for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
                    scale(_h_ungroomedJetMass_dj[i][0],   crossSection()/picobarn / sumOfWeights() / (ptBins_dj[i+1]-ptBins_dj[i]) );
                    scale(_h_sdJetMass_dj[i][0],          crossSection()/picobarn / sumOfWeights() / (ptBins_dj[i+1]-ptBins_dj[i]) );
                }*/
            }
        //@}


        private:

            /// @name FastJet grooming tools (configured in constructor init list)
            //@{
            const fjcontrib::SoftDrop _softdrop;
            const fjcontrib::Nsubjettiness _nsub_0p5_1;
            const fjcontrib::Nsubjettiness _nsub_0p5_2;
            const fjcontrib::Nsubjettiness _nsub_0p5_3;
            const fjcontrib::Nsubjettiness _nsub_0p5_4;
            const fjcontrib::Nsubjettiness _nsub_1_1;
            const fjcontrib::Nsubjettiness _nsub_1_2;
            const fjcontrib::Nsubjettiness _nsub_1_3;
            const fjcontrib::Nsubjettiness _nsub_1_4;
            const fjcontrib::Nsubjettiness _nsub_2_1;
            const fjcontrib::Nsubjettiness _nsub_2_2;
            const fjcontrib::Nsubjettiness _nsub_2_3;
            const fjcontrib::Nsubjettiness _nsub_2_4;
            const fjcontrib::NsubjettinessRatio _nsub21;
            const fjcontrib::NsubjettinessRatio _nsub32;
            //@}


            /// @name Histograms
            //@{
            Histo1DPtr _h_ungroomedJet0pt, _h_ungroomedJet1pt;
            Histo1DPtr _h_sdJet0pt, _h_sdJet1pt;
            // Here, store both the absolute (index 0) and normalized (index 1) cross sections.
            static const int betaNsub=3;
            static const int indexNsub=4;
            Histo1DPtr _h_ungroomedJet0_nSub[betaNsub][indexNsub];
            Histo1DPtr _h_ungroomedJet1_nSub[betaNsub][indexNsub];
            Histo1DPtr _h_ungroomedJet0_nSub21, _h_ungroomedJet1_nSub21,_h_ungroomedJet0_nSub32, _h_ungroomedJet1_nSub32;
            //@}


    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(CMS_2021_PAS_SMP_21_XXX);


}

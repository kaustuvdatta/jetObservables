# Selection that requires >=1 AK8 jet (w/ pT>200 GeV,y<|1.7| and m_SD>50 GeV for skimming here), 
# ==1 high pT muon (pT>55) in the event, 
# 0 further high pT muons in system passing the above's requirements at gen/reco-level based criteria
# Veto on any further non-prompt or loose leptons (muons/electrons) w/ pT>15 GeV  in event (loose/cutBasedId criteria at reco level MC and data, pdgId at gen level),
# MET pT> 50 GeV, a leptonic hemisphere W (mu+nu 4-vecs) system with pT > 150 GeV, 
# exactly one AK4 in leptonic/hadronic hemisphere w/ pT>30 GeV, 
# AK4 in leptonic hemisphere is required to be the leading (in pT) AK4 close to the, i.e., within deltaR(mu,AK4)<1.6 
# and it must also be b-tagged/originate from a b (medium WP DeepJet a la DeepFlavB/ghost hadron flvaour matched for reco/gen)
# additional pT requirements on MET (MET_T1_pt_nom/genMET>30 GeV) and pT of leptonically decaying W (4-momentum addition of mu + MET)
# otherwise event is discarded since the recoiling leptonic top isn't reconstructed well.
# require geometric separations between measurement candidate jet (AK8) and the leptonic hemisphere objects in additoin to above considerations



import ROOT
import math, os, sys
import numpy as np
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METXYCorrector import METXYCorr_Met_MetPhi

ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

class nSubProd(Module):

    def __init__(self, sysSource=[], leptonSF={}, year='2017', isMC=True, 
                 onlyUnc='', evtSelection='_WtopSel', isSigMC=False, 
                 effMapMode=False, controlHistos=False, applyMETxyCorrections=True):

        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.constCheckPlots = False

        self.applyMETxyCorrections=applyMETxyCorrections

        if '_constituentJES' in self.onlyUnc and (self.onlyUnc.endswith(('charged', 'neutral', 'photon'))):

            if not('photon' in self.onlyUnc):
                self.constJESVariation = 0.05 if 'neutral' in self.onlyUnc else 0.01
            else:
                self.constJESVariation = 0.03

        self.controlHistos=controlHistos
        
        self.isSigMC = isSigMC
        self.alsoDeltaRMatchedB = False # track events where the hadronic hemisphere b jets in reco and gen  are also deltaR matched just as the AK8's are

        self.nrecoEvents = 0
        self.recoWeightSum = 0
        self.ngenEvents = 0
        self.genWeightSum = 0
        self.evtCounter=0
        self.nTruerecoEvents = 0
        self.nAccepgenEvents = 0

        self.evtSelection = evtSelection
        self.recoLevel=0
        self.fakes=0
        self.miss=0
        self.genLevel=0
        self.response=0
        self.ufo=0
        self.ufoResponse=0
        self.ufoFake=0
        self.ufoMiss=0
        self.recoLevelW=0
        self.fakesW=0
        self.missW=0
        self.genLevelW=0
        self.responseW=0
        self.ufoW=0
        self.recoLeveltop=0
        self.fakestop=0
        self.misstop=0
        self.genLeveltop=0
        self.responsetop=0
        self.ufotop=0
        self.effMapOnly=effMapMode

        ### Kinematics Cuts AK8Jets ###
        self.minAK8JetPt = 170.  ### this is the basic minimum, not the final
        self.maxJetAK8Rap = 1.7

        ### Cut on minLeadAK8JetPtW and minLeadJetMass (using SD mass now and not inv. mass) to store both W/top-like events (and for the CR in the gap in mSD window b/w both) in same skims
        self.minLeadJetMass = 50.      #30. reducing to check if inv. mass cut is reasonable in test skims
        self.minLeadAK8JetPtW = 200.
        self.minSDMassW = 65. 
        self.maxSDMassW = 120.
        self.minLeadAK8JetPtTop= 400.
        self.minSDMassTop = 140.
        self.maxSDMassTop = 300.
        
        self.METCutWtop = 30. #50; loosening to check effs of selection
        self.minLeptonicWPt = 100. #150; loosening to check effs of selection


        ### Kinematics Cuts AK4 Jets ###
        self.minJetPt = 25.
        self.maxJetEta = 2.4

        if self.year=='2017':
            self.minBDisc = 0.3040   ### L: 0.0532, M: 0.3040, T: 0.7476, for DeepJet (ie, DeepFlavB);https://btv-wiki.docs.cern.ch/ScaleFactors/UL2017/
            #self.jetIDFlag = 2 #tight in 17/18, 3 is tight in 2016 as per https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#nanoAOD%20Flags~!!!!!! this is outdated, all are 2 now
        elif self.year=='2018':
            self.minBDisc = 0.2783   ### L: 0.0490, M: 0.2783, T: 0.7476; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
            #self.jetIDFlag = 2 #tight in 17/18, 3 is tight in 2016
        elif self.year=='2016_preVFP':
            self.minBDisc = 0.2598   ### L: 0.0508, M: 0.2598, T: 0.6502; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016preVFP/
            #self.jetIDFlag = 2 #tight in 17/18, 3 is tight in 2016   as per          
        elif self.year=='2016':
            self.minBDisc = 0.2489   ### L: 0.0480, M: 0.2489, T: 0.6377; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016postVFP/
            #self.jetIDFlag = 2 #tight in 17/18, 3 is tight in 2016
        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 55. 
        self.minLooseMuonPt = 15. 
        self.maxMuonEta = 2.4
        #self.minMuonMETPt = 50.

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 15.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]
        self.maxElectronEta = 2.4

        #const.-level energy scsale variation test plots    
        self.ptBins = [0, 1, 2, 3, 4, 5, 10, 15]
        self.pT_bins_names = ['0_1', '1_2', '2_3', '3_4', '4_5', '5_10', '10_15', 'gt15']

        #overall event weights, updated in functions below as required per nominal/systematics runs of the skimmer
        self.totalRecoWeight = 1.
        
        # Nominal values of systematic weights handled by the following class-level variables for events passing selections
        # Only storing for some systematics, since storing the nominals for the others don't make sense 
        # (ie, eg, nominal for isr/fsrWeight is 1, and PSWeights store correlated (w_var/w_nom) weights a la ISRup/FSRnom,ISRnom/FSRnom,ISRdown/FSRnom, ISRnom/FSRdown) 

        self.puWeight = 1. # storing nom. puWeight for reco
        self.evtGenWeight = 1. # storing event.genWeight for selected reco+gen events
        self.pdfWeight = 1. 
        self.pdfWeightAll = np.ones((103,), dtype=np.float32)*1.
        self.l1PreFireWeight = 1.
        self.leptonWeight = 1.

        self.pdfWeightUp = 1.
        self.puWeightUp = 1.
        self.l1PreFireWeightUp = 1.
        self.isrWeightUp = 1.
        self.fsrWeightUp = 1.
        self.leptonWeightAllUp = 1.
        self.leptonWeightISOUp = 1.
        self.leptonWeightIDUp = 1.
        self.leptonWeightTrigUp = 1.
        self.leptonWeightRecoEffUp = 1.

        self.pdfWeightDown = 1.
        self.puWeightDown = 1.
        self.l1PreFireWeightDown = 1.
        self.isrWeightDown = 1.
        self.fsrWeightDown = 1.
        self.leptonWeightAllDown = 1.
        self.leptonWeightISODown = 1.
        self.leptonWeightIDDown = 1.
        self.leptonWeightTrigDown = 1.
        self.leptonWeightRecoEffDown = 1.

        self.dummy = 0

        ### Defining nsubjetiness basis
        self.maxTau = 6
        self.nSub_labels = {
                        "_tau_0p25_1": [0., 1., 1000  ],
                        "_tau_0p25_2": [0., 1., 1000  ],
                        "_tau_0p25_3": [0., 1., 1000  ],
                        "_tau_0p25_4": [0., 1., 1000  ],
                        "_tau_0p25_5": [0., 1., 1000  ],
                        "_tau_0p25_6": [0., 1., 1000  ],
                        
                        "_tau_0p5_1": [0., 0.9, 900  ],
                        "_tau_0p5_2": [0., 0.9, 900  ],
                        "_tau_0p5_3": [0., 0.8, 800  ],
                        "_tau_0p5_4": [0., 0.8, 800  ],
                        "_tau_0p5_5": [0., 0.8, 800  ],
                        "_tau_0p5_6": [0., 0.8, 800  ],

                        "_tau_1_1": [ 0., 0.9, 900  ],
                        "_tau_1_2": [ 0., 0.7, 700  ],
                        "_tau_1_3": [ 0., 0.5, 500  ],
                        "_tau_1_4": [ 0., 0.5, 500  ],
                        "_tau_1_5": [ 0., 0.5, 500  ],
                        "_tau_1_6": [ 0., 0.5, 500  ],

                        "_tau_1p5_1": [ 0., 0.9, 900  ],
                        "_tau_1p5_2": [ 0., 0.7, 700  ],
                        "_tau_1p5_3": [ 0., 0.5, 500  ],
                        "_tau_1p5_4": [ 0., 0.5, 500  ],
                        "_tau_1p5_5": [ 0., 0.5, 500  ],
                        "_tau_1p5_6": [ 0., 0.5, 500  ],

                        "_tau_2_1": [ 0., 0.9, 900  ],
                        "_tau_2_2": [ 0., 0.5, 500  ],
                        "_tau_2_3": [ 0., 0.5, 500  ],
                        "_tau_2_4": [ 0., 0.5, 500  ],
                        "_tau_2_5": [ 0., 0.5, 500  ],
                        "_tau_2_6": [ 0., 0.5, 500  ]
                }

        self.tauPrefixes = ['_tau_0p25_','_tau_0p5_','_tau_1_','_tau_1p5_','_tau_2_']
        self.nSub0p25 = ROOT.NsubjettinessWrapper( 0.25, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) 
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub1p5 = ROOT.NsubjettinessWrapper( 1.5, 0.8, 0, 0 ) 
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )

        self.nSub1_OP_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 ) ##### needed for genjet tau21 or tau32, WTA_kT=3, OP_kT=6
        self.nSub1_WTA_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 3 )

        ### Softdrop quantities
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, 0.)#self.minAK8JetPt

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")

        ### Helpers
        #self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
        self.nJet = [ 'Jet']#, 'sdJet' ]
        #if self.runSDVariables: self.nJet = [ 'sdJet']

        ### Uncertainties
        self.sysSource = ['_nom'] +[ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') and not self.effMapOnly]
      
        # pu/leptonWeights used only to change reco event weight (up/down) without application to gen weight, 
        # others applied to modify self.genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)] ##### this is now done in the histogramming #####
        # weights+up/down variations are saved for accepted events, to be applied in histogramming in off-crab postprocessing
        self.sysWeightList = ( '_pu', '_btag', '_lepton', '_l1', '_pdf', '_isr', '_fsr') if not self.effMapOnly else ()
        self.sysrecoWeightList = ( '_pu', '_btag', '_lepton', '_l1') if not self.effMapOnly else ()
        self.sysgenWeightList = ( '_pdf', '_isr', '_fsr' ) if not self.effMapOnly else ()
        
        #############################################################################

    def beginJob(self, histFile, histDirName):

        Module.beginJob(self, histFile, histDirName)

        ### Booking histograms
        
        if 'Wtop' in self.evtSelection: 
            selList = ['_WtopSel']#[ '_'+x+i for i in ['_WtopSel'] for x in selList] #['_WSel','_topSel']
        elif self.evtSelection.endswith('WSel') and not 'top' in self.evtSelection: 
            selList = ['_WSel']
        elif self.evtSelection.endswith('topSel'): 
            selList = ['_topSel']
        else: 
            print ("Selection %s cannot be handled with this skimmer; what funny business are you even up to?"%self.evtSelection)
            print ("############################################ SELECTION ERROR ############################################")

        print("Working in the following selections:", selList)

        

        if (self.onlyUnc==''):
        
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   25, 0, 25) )
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            allSel = selList + ['_WSel','_topSel'] if self.evtSelection.startswith('_Wtop') else selList
            #### general selection
            if self.controlHistos:
                for isel in [ '_noSelnoWeight', '_noSel' ] + allSel:
                    
                    self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                    
                
                    self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   10, 0, 10) )
                    self.addP4Hists( 'muons', isel )
                    self.addP4Hists( 'eles', isel )
                    
                    self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   10, 0, 10) )
                    self.addP4Hists( 'AK8jets', isel )
                    
                    self.addObject( ROOT.TH1F('nAK4jets'+isel,   ';number of AK4 jets',   35, 0, 35) )
                    self.addObject( ROOT.TH1F('nAK4bjets'+isel,   ';number of b-tagged AK4s',   6, 0, 6) )
                    self.addObject( ROOT.TH1F('nAK4leptjets'+isel,   ';number of lept. hem. AK4s',   6, 0, 6) )
                    self.addObject( ROOT.TH1F('nAK4bleptjets'+isel,   ';number of b-tagged AK4s in lept. hem.',   6, 0, 6) )
                    self.addObject( ROOT.TH1F('dR_Mu_AK4blept'+isel, '; #Delta R(#mu, b-tagged AK4 in lept. hem.)', 40, 0., 2.))
                    self.addObject( ROOT.TH1F('dPhi_AK8_AK4blept'+isel, '; #Delta #phi(AK8, b-tagged AK4 in lept. hem.)', 40, 0., 2.))

                    self.addP4Hists( 'AK4jets', isel )
                    self.addP4Hists( 'AK4btaggedjets', isel )
                    self.addP4Hists( 'AK4leptjets', isel )
                    self.addP4Hists( 'AK4leptbjet', isel )
                    self.addP4Hists( 'leptonicW', isel )
                    self.addP4Hists( 'leptonicTop', isel )
                                    
                    self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )
                    self.addObject( ROOT.TH1F('Mtt'+isel,   '; m_{t#bar{t}}(GeV)',   200, 0, 2000) )
                    self.addObject( ROOT.TH1F('HT'+isel,   ';HT (GeV)',   200, 0, 2000) )

            
                
            if self.isMC:
                for isel in [ '_noSel' ] + allSel:
                    if self.controlHistos:
                        self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   10, 0, 10) )
                        
                        self.addP4Hists( 'genmuons', isel )
                        self.addP4Hists( 'geneles', isel )
                        
                        self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   10, 0, 10) )
                        self.addP4Hists( 'AK8genjets', isel )
                        
                        self.addObject( ROOT.TH1F('ngenAK4jets'+isel,   ';number of AK4 genjets',   35, 0, 35) ) #ngenAK4jets_noSel,ngenAK4bjets_noSel,ngenAK4bleptjets_noSel
                        self.addObject( ROOT.TH1F('ngenAK4bjets'+isel,   ';number of b-flav. matched AK4 genjets',   6, 0, 6) )
                        self.addObject( ROOT.TH1F('ngenAK4leptjets'+isel,   ';number of lept. hem. AK4 genjets',   6, 0, 6) )
                        self.addObject( ROOT.TH1F('ngenAK4bleptjets'+isel,   ';number of b-flav. matched AK4 genjets in lept. hem.',   6, 0, 6) )
                        self.addObject( ROOT.TH1F('dR_genMu_genAK4blept'+isel, '; #Delta R(gen #mu, b-matched gen AK4 in lept. hem.)', 40, 0., 2.))
                        self.addObject( ROOT.TH1F('dPhi_genAK8_genAK4blept'+isel, '; #Delta #phi(gen AK8, b-matched gen AK4 in lept. hem.)', 40, 0., 2.))
        
                        self.addP4Hists( 'AK4genjets', isel )
                        self.addP4Hists( 'AK4bmatchedgenjets', isel )
                        self.addP4Hists( 'AK4genleptjets', isel )
                        self.addP4Hists( 'AK4genleptbjet', isel )
                        self.addP4Hists( 'genleptonicW', isel )
                        self.addP4Hists( 'genleptonicTop', isel )
                        
                        self.addObject( ROOT.TH1F('genMETPt'+isel,   ';gen MET (GeV)',   200, 0, 2000) )
                        self.addObject( ROOT.TH1F('genMtt'+isel,   '; gen m_{t#bar{t}}(GeV)',   200, 0, 2000) )
                        self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )                

                

                self.addObject( ROOT.TH1F('gen_dR_charged',   '; #Delta R(const_{ch}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('gen_dRap_dPhi_charged',   '; #Delta y(const_{ch}, genAK8); #Delta #phi(const_{ch}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('gen_dR_chargedHadrons',   '; #Delta R(h_{ch}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('gen_dRap_dPhi_chargedHadrons',   '; #Delta y(h_{ch}, genAK8); #Delta #phi(h_{ch}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('gen_dR_neutral',   '; #Delta R(h_{neut.}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('gen_dRap_dPhi_neutral',   '; #Delta y(h_{neut.}, genAK8); #Delta #phi(h_{neut.}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('gen_dR_photons',   '; #Delta R(#gamma, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('gen_dRap_dPhi_photons',   '; #Delta y(#gamma, genAK8); #Delta #phi(#gamma, genAK8)', 100, -1., 1., 100, -1., 1.) )
                

                self.addObject( ROOT.TH1F('reco_dR_charged',   '; #Delta R(const_{ch}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('reco_dRap_dPhi_charged',   '; #Delta y(const_{ch}, recoAK8); #Delta #phi(const_{ch}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('reco_dR_chargedHadrons',   '; #Delta R(h_{ch}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('reco_dRap_dPhi_chargedHadrons',   '; #Delta y(h_{ch}, recoAK8); #Delta #phi(h_{ch}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('reco_dR_neutral',   '; #Delta R(h_{neut.}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('reco_dRap_dPhi_neutral',   '; #Delta y(h_{neut.}, recoAK8); #Delta #phi(h_{neut.}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('reco_dR_photons',   '; #Delta R(#gamma, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('reco_dRap_dPhi_photons',   '; #Delta y(#gamma, recoAK8); #Delta #phi(#gamma, recoAK8)', 100, -1., 1., 100, -1., 1.) )


                self.addObject( ROOT.TH1F('gen_ncharged',   '; n const_{ch}in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('gen_nchargedHadrons',   '; n h_{ch} in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('gen_nneutral',   '; n h_{neut.} in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('gen_nphotons',   '; n #gamma in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('reco_ncharged',   '; n const_{ch} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('reco_nchargedHadrons',   '; n h_{ch} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('reco_nneutral',   '; n h_{neut.} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('reco_nphotons',   '; n #gamma in recoAK8', 100, 0, 100) )


                self.addObject( ROOT.TH1F('gen_pT_charged',   '; p_{T,ch}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('gen_pT_chargedHadrons',   '; p_{T,ch had}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('gen_pT_neutral',   '; p_{T,neut.}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('gen_pT_photons',   '; p_{T,#gamma}', 200, 0, 100.) )

                self.addObject( ROOT.TH1F('reco_pT_charged',   '; p_{T,ch}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('reco_pT_chargedHadrons',   '; p_{T,ch had}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('reco_pT_neutral',   '; p_{T,neut.}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('reco_pT_photons',   '; p_{T,#gamma}', 200, 0, 100.) )


                self.addObject( ROOT.TH1F('reco_dR_SJ1',   '; #Delta R(AK8, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR_SJ2',   '; #Delta R(AK8, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2_SJ1',   '; #Delta R2(AK8, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2_SJ2',   '; #Delta R2(AK8, AK8 SJ2)', 100, 0, 1.) )


                self.addObject( ROOT.TH1F('reco_dRSJ1_charged',   '; #Delta R(const_{ch}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ1_chargedHadrons',   '; #Delta R(h_{ch}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ1_neutral',   '; #Delta R(h_{neut.}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ1_photons',   '; #Delta R(#gamma, AK8 SJ1)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('reco_dRSJ2_charged',   '; #Delta R(const_{ch}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ2_chargedHadrons',   '; #Delta R(h_{ch}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ2_neutral',   '; #Delta R(h_{neut.}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dRSJ2_photons',   '; #Delta R(#gamma, AK8 SJ2)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('reco_dR2SJ1_charged',   '; #Delta R2(const_{ch}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ1_chargedHadrons',   '; #Delta R2(h_{ch}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ1_neutral',   '; #Delta R2(h_{neut.}, AK8 SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ1_photons',   '; #Delta R2(#gamma, AK8 SJ1)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('reco_dR2SJ2_charged',   '; #Delta R2(const_{ch}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ2_chargedHadrons',   '; #Delta R2(h_{ch}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ2_neutral',   '; #Delta R2(h_{neut.}, AK8 SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('reco_dR2SJ2_photons',   '; #Delta R2(#gamma, AK8 SJ2)', 100, 0, 1.) )


                self.addObject( ROOT.TH2F('genjetPt_dR_charged',   '; #Delta R(const_{ch}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('genjetPt_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('genjetPt_dR_neutral',   '; #Delta R(h_{neut.}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('genjetPt_dR_photons',   '; #Delta R(#gamma, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )

                self.addObject( ROOT.TH2F('recojetPt_dR_charged',   '; #Delta R(const_{ch}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('recojetPt_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('recojetPt_dR_neutral',   '; #Delta R(h_{neut.}, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('recojetPt_dR_photons',   '; #Delta R(#gamma, AK8); p_{T} AK8', 100, 0, 1., 120, 0, 1200.) )

                self.addObject( ROOT.TH2F('genjetmSD_dR_charged',   '; #Delta R(const_{ch}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('genjetmSD_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('genjetmSD_dR_neutral',   '; #Delta R(h_{neut.}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('genjetmSD_dR_photons',   '; #Delta R(#gamma, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )

                self.addObject( ROOT.TH2F('recojetmSD_dR_charged',   '; #Delta R(const_{ch}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('recojetmSD_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('recojetmSD_dR_neutral',   '; #Delta R(h_{neut.}, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('recojetmSD_dR_photons',   '; #Delta R(#gamma, AK8); m_{SD} AK8', 100, 0, 1., 30, 50, 300.) )


                self.addObject( ROOT.TH2F('gen_dR_Pt_charged',   '; p_{T} const_{ch};  #Delta R(const_{ch}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('gen_dR_Pt_chargedHadrons',   '; p_{T} h_{ch};  #Delta R(h_{ch}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('gen_dR_Pt_neutral',   '; p_{T} h_{neut.};  #Delta R(h_{neut.}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('gen_dR_Pt_photons',   '; p_{T} #gamma;  #Delta R(#gamma, AK8); ',  200, 0, 100., 100, 0, 1.,) )

                self.addObject( ROOT.TH2F('reco_dR_Pt_charged',   '; p_{T} const_{ch};  #Delta R(const_{ch}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('reco_dR_Pt_chargedHadrons',   '; p_{T} h_{ch};  #Delta R(h_{ch}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('reco_dR_Pt_neutral',   '; p_{T} h_{neut.};  #Delta R(h_{neut.}, AK8); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('reco_dR_Pt_photons',   '; p_{T} #gamma;  #Delta R(#gamma, AK8); ',  200, 0, 100., 100, 0, 1.,) )


                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('gen_dR_chargedHadrons_{}'.format(bin_name), '; #Delta R(h_{{ch}}, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('gen_dRap_dPhi_chargedHadrons_{}'.format(bin_name), '; #Delta y(h_{{ch}}, genAK8, pT {}); #Delta #phi(h_{{ch}}, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('reco_dR_chargedHadrons_{}'.format(bin_name), '; #Delta R(h_{{ch}}, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('reco_dRap_dPhi_chargedHadrons_{}'.format(bin_name), '; #Delta y(h_{{ch}}, recoAK8, pT {}); #Delta #phi(h_{{ch}}, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))

                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('gen_dR_neutral_{}'.format(bin_name), '; #Delta R(h_{{neut.}}, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('gen_dRap_dPhi_neutral_{}'.format(bin_name), '; #Delta y(h_{{neut.}}, genAK8, pT {}); #Delta #phi(h_{{neut.}}, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('reco_dR_neutral_{}'.format(bin_name), '; #Delta R(h_{{neut.}}, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('reco_dRap_dPhi_neutral_{}'.format(bin_name), '; #Delta y(h_{{neut.}}, recoAK8, pT {}); #Delta #phi(h_{{neut.}}, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))

                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('gen_dR_photons_{}'.format(bin_name), '; #Delta R(#gamma, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('gen_dRap_dPhi_photons_{}'.format(bin_name), '; #Delta y(#gamma, genAK8, pT {}); #Delta #phi(#gamma, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('reco_dR_photons_{}'.format(bin_name), '; #Delta R(#gamma, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('reco_dRap_dPhi_photons_{}'.format(bin_name), '; #Delta y(#gamma, recoAK8, pT {}); #Delta #phi(#gamma, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))


    #############################################################################
    def addP4Hists(self, s, t ):


        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        
        if not 'met' in s.lower():
            self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -2.5, 2.5 ) )
            self.addObject( ROOT.TH1F(s+'_y'+t, s+';y', 100, -2.5, 2.5 ) )
        
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14259) )
        
        if not 'mu' in s.lower() or 'met' in s.lower(): 
            self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 500) )
        
            if 'ak8' in s.lower():
                self.addObject( ROOT.TH1F(s+'_msoftdrop'+t,s+';mass (GeV)', 100, 0, 500) )
            

    #############################################################################
    def DrRapPhi(self,va,vb ):
        dy = va.Rapidity()-vb.Rapidity()
        dphi = va.DeltaPhi(vb)
        return ROOT.TMath.Sqrt(dy*dy+dphi*dphi)

    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        print (wrappedOutputTree, outputFile, inputFile, inputTree, self.sysSource)
        
        self.out = wrappedOutputTree
        
        tmplistAK8=[]
        tmplistAK4=[]
        tmplistLeptAK8=[]
        tmplistMuon=[]
        tmplistLeptW=[]
        tmplistLeptTop=[]

        tmplistMETPt=[]

        # need to keep branches for selected AK8s, selected muon, selected(reconstructed) leptonicW,
        # MET pT, selected/b-tagged leptonic hemisphere AK4

        if not self.isMC: 
            
            self.out.branch('totalRecoWeight_nom', "F" )
            self.out.branch('passRecoSel_nom', "I" )
            self.out.branch('nRecoBtags_nom',  "I") 
            self.out.branch('nRecoAK4s_nom',  "I") 
            #self.out.branch('nRecoCrackElectrons',  "I") 
            self.out.branch('nRecoLeptBtags_nom',  "I") 
            #self.out.branch('nRecoAK8s_nom',  "I") 
            self.out.branch('recoSelectedEventNumber_nom', "L" )
            self.out.branch('good_nPVs_nom', "F" )
            
            self.out.branch('selRecoLeptHemDeltaR_nom', "F" ) # dR(muon, accepted AK4 b-jet in lept hem)
            self.out.branch('selRecoLeptHemDeltaPhi_nom', "F" ) # dphi(muon, accepted AK4 b-jet in lept hem)
            self.out.branch('selRecoLeptHemDeltaRap_nom', "F" ) # dy(muon, accepted AK4 b-jet in lept hem)
            self.out.branch('selRecoLeptHemAK8DeltaPhi_nom', "F" ) # dphi(AK8, accepted AK4 b-jet in lept hem)
            self.out.branch('selRecoLeptHemAK8DeltaR_nom', "F" ) # dphi(AK8, accepted AK4 b-jet in lept hem)
            self.out.branch('passHLT_Mu50_nom', "I" ) 
            self.out.branch('passHLT_TkMu50_nom', "I" ) 
            self.out.branch('passHLT_OldMu100_nom', "I" ) 
            self.out.branch('passHLT_TkMu100_nom', "I" ) 
            

            

            tmplistAK8.append('selRecoJets_nom')
            tmplistMuon.append('selRecoMu_nom')
            tmplistLeptW.append('selRecoLeptW_nom')
            tmplistLeptTop.append('selRecoLeptTop_nom')
            tmplistMETPt.append('selRecoMET_nom')
            if self.applyMETxyCorrections:
                tmplistLeptTop.append('selRecoLeptTop_uncorrected_nom')
                tmplistLeptW.append('selRecoLeptW_uncorrected_nom')
                tmplistMETPt.append('selRecoMET_uncorrected_nom')

            # We require exactly one leptonic hemisphere AK4 (b-tagged),
            # so keeping that objects info in a branch, 
            tmplistAK4.append('selRecoAK4bjetLeptHem_nom') # to store leptonic hem., lead b-tagged AK4's
            tmplistLeptAK8.append('selRecoAK8jetLeptHem_nom')
        elif self.isMC: 
            for sys in self.sysSource:

                if not('const' in sys):

                    self.out.branch('totalRecoWeight'+sys, "F" )
                    self.out.branch('passRecoSel'+sys, "I" )
                    self.out.branch('nRecoBtags'+sys,  "I") 
                    self.out.branch('nRecoAK4s'+sys,  "I") 
                    #self.out.branch('nRecoCrackElectrons'+sys,  "I") 
                    #self.out.branch('nRecoAK8s'+sys,  "I") 
                    self.out.branch('nRecoLeptBtags'+sys,  "I") 
                    self.out.branch('recoSelectedEventNumber'+sys, "L" )
                    self.out.branch('good_nPVs'+sys, "F" )
                                    
                    self.out.branch('leptonWeightNom'+sys,  "F")
                    self.out.branch('puWeightNom'+sys,  "F")
                    self.out.branch('l1prefiringWeightNom'+sys,  "F")
                    self.out.branch('pdfWeightNom'+sys, "F" ) 

                    self.out.branch('FlagDeltaRMatchedBjets'+sys, "I")
                    
                    self.out.branch( 'accepgenSelectedEventNumber'+sys, "L" ) 
                    self.out.branch( 'truerecoSelectedEventNumber'+sys, "L" )
                    
                    # For tagging leptonic hemisphere top, holds info on geometric separations between b-jet and muon
                    self.out.branch('selRecoLeptHemAK8DeltaPhi'+sys, "F" ) 
                    self.out.branch('selRecoLeptHemAK8DeltaR'+sys, "F" ) 

        
                    self.out.branch('evtGenWeight'+sys, "F" )  # nominal generator weight for evt.
                    self.out.branch('passGenSel'+sys, "I" )  # nominal generator weight for evt.
                    self.out.branch('genSelectedEventNumber'+sys, "L" )


                if sys.startswith('_nom'): 
                    
                    #nGenBtags=1 for evt. cat. 0; 2 for evt. cat. 1; >2 for evt. cat. 1 as well, 0 for evt. cat. -1
                    self.out.branch('nGenBtags'+sys,  "I") #for hadron flavour ghost-matched b's
                    self.out.branch('nGenAK4s'+sys,  "I") #for all AK4's close to muon
                    #self.out.branch('nGenCrackElectrons'+sys,  "I") 
                    #self.out.branch('nGenAK8s'+sys,  "I") 
                    self.out.branch('nGenLeptBtags'+sys,  "I") #for hadron flavour ghost-matched b's

                    #systematics on jes/jer studies on measurement AK8 are not changing what goes into these following branches
                    
                    
                    self.out.branch('selGenLeptHemDeltaR'+sys, "F")
                    self.out.branch('selGenLeptHemDeltaPhi'+sys, "F")
                    self.out.branch('selGenLeptHemDeltaRap'+sys, "F")
                
                    self.out.branch('selGenLeptHemAK8DeltaPhi'+sys, "F")
                    self.out.branch('selGenLeptHemAK8DeltaR'+sys, "F")
                    self.out.branch('passHLT_Mu50_nom', "I" ) 
                    self.out.branch('passHLT_TkMu50_nom', "I" ) 
                    self.out.branch('passHLT_OldMu100_nom', "I" ) 
                    self.out.branch('passHLT_TkMu100_nom', "I" ) 

                    
                self.out.branch('selRecoLeptHemDeltaR'+sys,"F")
                self.out.branch('selRecoLeptHemDeltaPhi'+sys,"F")
                self.out.branch('selRecoLeptHemDeltaRap'+sys,"F")
                
                if self.isSigMC and (self.onlyUnc=='') and not (self.effMapOnly): 
                    # only storing in nominal case for sigMC
                    # leaving functionality for non nominal cases in case (ie, switch off onlyUnc req.)
                    self.out.branch( 'pdfWeightAll'+sys, 'F', 103 )
                    self.out.branch( 'pdfWeightUp'+sys, "F" )
                    self.out.branch( 'pdfWeightDown'+sys, "F" )

                    self.out.branch( 'isrWeightUp'+sys, "F" )
                    self.out.branch( 'isrWeightDown'+sys, "F" )

                    self.out.branch( 'fsrWeightUp'+sys, "F" )
                    self.out.branch( 'fsrWeightDown'+sys, "F" )

                    self.out.branch( 'puWeightUp'+sys, "F" )
                    self.out.branch( 'puWeightDown'+sys, "F" )

                    self.out.branch( 'l1prefiringWeightUp'+sys,  "F")
                    self.out.branch( 'l1prefiringWeightDown'+sys,  "F")


                    self.out.branch( 'leptonWeightAllUp'+sys, "F" )
                    self.out.branch( 'leptonWeightAllDown'+sys, "F" )

                    self.out.branch( 'leptonWeightISOUp'+sys, "F" )
                    self.out.branch( 'leptonWeightISODown'+sys, "F" )

                    self.out.branch( 'leptonWeightIDUp'+sys, "F" )
                    self.out.branch( 'leptonWeightIDDown'+sys, "F" )

                    self.out.branch( 'leptonWeightTrigUp'+sys, "F" )
                    self.out.branch( 'leptonWeightTrigDown'+sys, "F" )

                    self.out.branch( 'leptonWeightRecoEffUp'+sys, "F" )
                    self.out.branch( 'leptonWeightRecoEffDown'+sys, "F" )
                   
            # since gen selection will always be unchanged we only store selGen objects for the nominal selection case
            # recojet kinematics may vary, due to systematics like jes/jer, so accepgen/true reco quantities may also vary, thus we store them separately for each sys
            # masks from weight branches required help to keep this consistent after skimming

            tmplistAK8_reco = [ 'selRecoJets'+sys for sys in self.sysSource ]
            tmplistAK8_truereco = ['trueRecoJets'+sys for sys in self.sysSource]
            tmplistAK8_gen = ['selGenJets'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            tmplistAK8_accepgen = ['accepGenJets'+sys for sys in self.sysSource ]
            for x in tmplistAK8_reco+tmplistAK8_truereco+tmplistAK8_gen+tmplistAK8_accepgen: tmplistAK8.append(x)

            #if sys.startswith('_nom'):
            tmplistmu_reco =  [ 'selRecoMu'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistmu_gen  =  [ 'selGenMu'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            for x in tmplistmu_reco+tmplistmu_gen : tmplistMuon.append(x)

            tmplistLeptW_reco =  [ 'selRecoLeptW'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistLeptW_gen  =  [ 'selGenLeptW'+sys for sys in self.sysSource  if sys.startswith('_nom')]
            for x in tmplistLeptW_reco+tmplistLeptW_gen: tmplistLeptW.append(x)

            tmplistLeptTop_reco =  [ 'selRecoLeptTop'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistLeptTop_gen  =  [ 'selGenLeptTop'+sys for sys in self.sysSource  if sys.startswith('_nom')]
            for x in tmplistLeptTop_reco+tmplistLeptTop_gen: tmplistLeptTop.append(x)

            tmplistMETPt_reco =  [ 'selRecoMET'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistMETPt_gen  =  [ 'selGenMET'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            for x in tmplistMETPt_reco+tmplistMETPt_gen: tmplistMETPt.append(x)

            if self.applyMETxyCorrections:
                tmplistLeptW_reco =  [ 'selRecoLeptW_uncorrected'+sys for sys in self.sysSource if not('const' in sys)]
                tmplistLeptTop_reco =  [ 'selRecoLeptTop_uncorrected'+sys for sys in self.sysSource if not('const' in sys)]
                tmplistMETPt_reco =  [ 'selRecoMET_uncorrected'+sys for sys in self.sysSource if not('const' in sys)]
                for x in tmplistLeptW_reco: tmplistLeptW.append(x)
                for x in tmplistLeptTop_reco: tmplistLeptTop.append(x)
                for x in tmplistMETPt_reco: tmplistMETPt.append(x)

            
            # Store leading reco/gen AK4 (b-tagged/hadron flavour matched) jet if any in the hadronic hemisphere
            #requring >=1 btag, but with exactly one b-tag in the hadronic hemisphere
            tmplistAK4bjetLeptHem_reco = ['selRecoAK4bjetLeptHem'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistAK4bjetLeptHem_gen = [ 'selGenAK4bjetLeptHem'+sys for sys in self.sysSource if sys.startswith('_nom') ] 
            for x in tmplistAK4bjetLeptHem_reco+tmplistAK4bjetLeptHem_gen: tmplistAK4.append(x)
            
            tmplistLeptAK8_reco = ['selRecoAK8jetLeptHem'+sys for sys in self.sysSource if not('const' in sys)]
            tmplistLeptAK8_gen = [ 'selGenAK8jetLeptHem'+sys for sys in self.sysSource if sys.startswith('_nom') ] 
            for x in tmplistLeptAK8_reco+tmplistLeptAK8_gen: tmplistLeptAK8.append(x)
            
        print ("Stored AK8 jet branches:", tmplistAK8)
        for iJ in tmplistAK8:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F') 
            self.out.branch(iJ+'_eta',  'F') 
            self.out.branch(iJ+'_y',  'F') 
            self.out.branch(iJ+'_phi',  'F') 
            self.out.branch(iJ+'_mass',  'F') 
            self.out.branch(iJ+'_msoftdrop',  'F') 
            #self.out.branch(iJ+'_mSD',  'F') 
            if not('gen' in iJ.lower()): 
                #self.out.branch(iJ+'_pt_nom',  'F')# a la nanoAOD-tools recorrection with JECs in postproc. 
                #self.out.branch(iJ+'_msoftdrop_nom',  'F')  # a la nanoAOD-tools, stored for posterity as its own branch in custom flat, unjagged selRecoJets collection
                #self.out.branch(iJ+'_msoftdrop_new',  'F')  # a la nanoAOD-tools, stored for posterity as its own branch in custom flat, unjagged selRecoJets collection
                #self.out.branch(iJ+'_msoftdrop_corr_JMR',  'F')                 
                self.out.branch(iJ+'_pt_raw',  'F') 
                self.out.branch(iJ+'_mass_raw',  'F') 
                self.out.branch(iJ+'_mSD',  'F') #AK8 JERC corrected (no subjet corrections beyond fixing defn of raw subjets in nanoAOD-tools as per v1 of changes a la PR281 from Fikri)
                self.out.branch(iJ+'_msoftdrop_nom_PUPPICorred',  'F')  # a la nanoAOD-tools, stored for posterity as its own branch in custom flat, unjagged selRecoJets collection
                self.out.branch(iJ+'_msoftdrop_corr_PUPPI',  'F') 
                self.out.branch(iJ+'_msoftdrop_raw',  'F') 
                self.out.branch(iJ+'_msoftdrop_corr_JMS',  'F') 
                self.out.branch(iJ+'_msoftdrop_corr_subjetJEC',  'F') 

                self.out.branch(iJ+'_JECfactor',  'F') 
                self.out.branch(iJ+'_JERfactor',  'F') 
                self.out.branch(iJ+'_JMSfactor',  'F') 

            self.out.branch(iJ+'_tau21',  'F') 
            self.out.branch(iJ+'_tau32',  'F') 
            self.out.branch(iJ+'_tau21_WTA',  'F') 
            self.out.branch(iJ+'_tau32_WTA',  'F') 
            self.out.branch(iJ+'_tau21_exkT',  'F') 
            self.out.branch(iJ+'_tau32_exkT',  'F') 
            #self.out.branch(iJ+'_genMatched',  'O')# 
            
            for x in self.nSub_labels:
                self.out.branch(iJ+x, 'F')
        
        print ("Stored AK4 jet branches:", tmplistAK4)
        for iJ in tmplistAK4:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F') 
            self.out.branch(iJ+'_eta',  'F') 
            self.out.branch(iJ+'_y',  'F') 
            self.out.branch(iJ+'_phi',  'F') 
            self.out.branch(iJ+'_mass',  'F') 
            if not('gen' in iJ.lower()):
                self.out.branch(iJ+'_jetId', 'I') 
                self.out.branch(iJ+'_btagDeepFlavB', 'F') 

                if self.isMC: 
                    self.out.branch(iJ+'_hadronFlavour', 'I') 
                
            

        print ("Stored mu branches:", tmplistMuon)
        for iMu in tmplistMuon:
            self.out.branch('n'+iMu,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iMu+'_pt',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_eta',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_y',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_phi',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_mass',  'F')#, lenVar='n'+iMu)
            if not('gen' in iMu.lower()):
                self.out.branch(iMu+'_ptRel', 'F')# 
                self.out.branch(iMu+'_tkRelIso', 'F')# #if not(self.isMC):
                self.out.branch(iMu+'_p', 'F')# 
            
        print ("Stored MET branches:", tmplistMETPt)
        for i in tmplistMETPt:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            
        print ("Stored leptonic W branches:", tmplistLeptW)
        for i in tmplistLeptW:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_eta',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_y',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mass',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mt1',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mt',  'F')#, lenVar='n'+i)

        print ("Stored leptonic top branches:", tmplistLeptTop)
        for i in tmplistLeptTop:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_eta',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_y',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mass',  'F')#, lenVar='n'+i)
            #self.out.branch(i+'_mt1',  'F')#, lenVar='n'+i)
            #self.out.branch(i+'_mt',  'F')#, lenVar='n'+i)

        print ("Stored leptonic hem. AK8 branches:", tmplistLeptAK8)
        for i in tmplistLeptAK8:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_eta',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_y',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mass',  'F')#, lenVar='n'+i)
            #self.out.branch(i+'_mt',  'F')#, lenVar='n'+i)


        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #if not self.onlyTrees:
        if self.isMC and self.onlyUnc=='':

            self.genLevel = self.response+self.miss

            print("evtCounter","nrecoEvents","ngenEvents","n_AccepgenEvents","n_TruerecoEvents")
            print(self.evtCounter,self.nrecoEvents,self.ngenEvents,self.nAccepgenEvents,self.nTruerecoEvents)
            print("recoWeightSum=",self.recoWeightSum)
            print("genWeightSum=",self.genWeightSum)

            #self.genLevelW = self.responseW+self.missW
            #self.genLeveltop = self.responsetop+self.misstop
            '''
            getattr( self, 'cutflow_test' ).SetBinContent( 1, self.recoLevel )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 1, "allReco" )
            getattr( self, 'cutflow_test' ).SetBinContent( 2, self.response )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 2, "allTrueRecoAccepGen" )
            getattr( self, 'cutflow_test' ).SetBinContent( 3, self.fakes )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 3, "allFakes" )
            getattr( self, 'cutflow_test' ).SetBinContent( 4, self.miss )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 4, "allMiss" )
            getattr( self, 'cutflow_test' ).SetBinContent( 5, self.genLevel )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 5, "allGen" )

            getattr( self, 'cutflow_test' ).SetBinContent( 7, self.ufo )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 7, "allUfo" )
            getattr( self, 'cutflow_test' ).SetBinContent( 8, self.ufoResponse )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 8, "ufoResponse" )
            getattr( self, 'cutflow_test' ).SetBinContent( 9, self.ufoFake )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 9, "ufoFake" )
            getattr( self, 'cutflow_test' ).SetBinContent( 10, self.ufoMiss )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 10, "ufoMiss" )
            '''
    #############################################################################
    def getleptonSF(self, lepton, leptonP4 ): #figure out how to uproot this approach with a cup of Coffe(a)

        leptonP4eta = abs(leptonP4.eta)
        leptonP = ROOT.TMath.Sqrt(leptonP4.p4().Px()**2 + leptonP4.p4().Py()**2 + leptonP4.p4().Pz()**2)

        SFFileTrigger = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['Trigger'][0] )
        histoSFTrigger = SFFileTrigger.Get( self.leptonSFhelper[lepton]['Trigger'][1] )
        SFTrigger = histoSFTrigger.GetBinContent( histoSFTrigger.GetXaxis().FindBin( leptonP4eta ), histoSFTrigger.GetYaxis().FindBin( leptonP4.pt ) )
        SFTriggerUp = histoSFTrigger.GetBinErrorUp( histoSFTrigger.GetXaxis().FindBin( leptonP4eta ), histoSFTrigger.GetYaxis().FindBin( leptonP4.pt ) )+SFTrigger
        SFTriggerDown = -histoSFTrigger.GetBinErrorLow( histoSFTrigger.GetXaxis().FindBin( leptonP4eta ), histoSFTrigger.GetYaxis().FindBin( leptonP4.pt ) )+SFTrigger

        SFFileID = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ID'][0] )
        histoSFID = SFFileID.Get( self.leptonSFhelper[lepton]['ID'][1] )
        histoSFID_X = histoSFID.GetXaxis().FindBin( leptonP4eta)
        histoSFID_Y = histoSFID.GetYaxis().FindBin( leptonP4.pt )
        SFID = histoSFID.GetBinContent( histoSFID_X, histoSFID_Y )
        SFID = SFID if SFID>0 else 1
        SFIDUp = histoSFID.GetBinErrorUp( histoSFID_X, histoSFID_Y )+SFID
        SFIDDown = -histoSFID.GetBinErrorLow( histoSFID_X, histoSFID_Y )+SFID

        SFFileISO = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ISO'][0] )
        histoSFISO = SFFileISO.Get( self.leptonSFhelper[lepton]['ISO'][1] )
        histoSFISO_X = histoSFISO.GetXaxis().FindBin( leptonP4eta )
        histoSFISO_Y = histoSFISO.GetYaxis().FindBin( leptonP4.pt )
        SFISO = histoSFISO.GetBinContent( histoSFISO_X, histoSFISO_Y )
        SFISO = SFISO if SFISO>0 else 1
        SFISOUp = histoSFISO.GetBinErrorUp( histoSFISO_X, histoSFISO_Y )+SFISO
        SFISODown = -histoSFISO.GetBinErrorLow( histoSFISO_X, histoSFISO_Y )+SFISO
        
        SFFileRecoEff = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['RecoEff'][0] )
        histoSFRecoEff = SFFileRecoEff.Get( self.leptonSFhelper[lepton]['RecoEff'][1] )
        histoSFRecoEff_X = histoSFRecoEff.GetXaxis().FindBin( leptonP4eta )
        histoSFRecoEff_Y = histoSFRecoEff.GetYaxis().FindBin( leptonP )
        SFRecoEff = histoSFRecoEff.GetBinContent( histoSFRecoEff_X, histoSFRecoEff_Y )
        SFRecoEff = SFRecoEff if SFRecoEff>0 else 1
        SFRecoEffUp = histoSFRecoEff.GetBinErrorUp( histoSFRecoEff_X, histoSFRecoEff_Y )+SFRecoEff
        SFRecoEffDown = -histoSFRecoEff.GetBinErrorLow( histoSFRecoEff_X, histoSFRecoEff_Y )+SFRecoEff

        #print (SFTrigger * SFID * SFISO), SFTrigger , SFID , SFISO, leptonP4.pt, leptonP4.eta
        return [[SFTrigger,SFTriggerUp,SFTriggerDown], [SFID,SFIDUp,SFIDDown] , [SFISO,SFISOUp,SFISODown], [SFRecoEff,SFRecoEffUp,SFRecoEffDown]]

    
    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''
        self.isMC = event.run == 1

        if self.isMC:
            # effectively a preselection
            muons = list(Collection(event, 'Muon'))
            genLeptons = list(Collection( event, 'GenDressedLepton' ))
            
            recoMuons = [ x for x in muons if x.pt > 50 and abs(x.eta) < self.maxMuonEta ]
            recoMuons.sort(key=lambda x: x.pt, reverse=True) 
            
            genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > 50 and abs(x.eta) < self.maxMuonEta ]
            genMuons.sort(key=lambda x: x.pt, reverse=True)
            
            if len(recoMuons)>=1 or len(genMuons)>=1:
                pass
            else:
                return False


        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGendR_Mu_leptAK4b, selGendPhi_Mu_leptAK4b, selGendRap_Mu_leptAK4b, selGendPhi_AK8_leptAK4b, selGendR_AK8_leptAK4b, selGenMuons, selGenAK4bjets, selGenAK4bleptjets, selGenAK4leptjets, selGenLeptW, selGenLeptTop, selGenJets, selGenMET, nselGenAK4 = self.genSelection(event)  #, nselGenCrackEl
        passRecoSel, iRecoSel, selRecodR_Mu_leptAK4b, selRecodPhi_Mu_leptAK4b, selRecodRap_Mu_leptAK4b, selRecodPhi_AK8_leptAK4b, selRecodR_AK8_leptAK4b, selRecoMuons, selRecoAK4bjets, selRecoAK4bleptjets, selRecoAK4leptjets, selRecoLeptW, selRecoLeptW_uncorr, selRecoLeptTop, selRecoLeptTop_uncorr, selRecoJets, selRecoMET, selRecoMET_uncorr, nselRecoAK4 = self.recoSelection(event)  #, nselRecoCrackEl
        #print (selRecodR_Mu_leptAK4b)


        if self.isMC:#in lieu of using this in postproc presel cuts, since we don't want to lose gen events...
            PV = (getattr(event, "PV_npvsGood") > 0)

            if not ('2016' in self.year):
                METFilters = (
                                (getattr(event, "Flag_goodVertices") == 1) and
                                (getattr(event, "Flag_globalSuperTightHalo2016Filter") == 1) and
                                (getattr(event, "Flag_HBHENoiseFilter") == 1) and
                                (getattr(event, "Flag_HBHENoiseIsoFilter") == 1) and
                                (getattr(event, "Flag_EcalDeadCellTriggerPrimitiveFilter") == 1) and
                                (getattr(event, "Flag_BadPFMuonFilter") == 1) and
                                (getattr(event, "Flag_BadPFMuonDzFilter") == 1) and
                                (getattr(event, "Flag_eeBadScFilter") == 1) and
                                (getattr(event, "Flag_ecalBadCalibFilter") == 1)
                            )
                
            else:
                METFilters = (
                                (getattr(event, "Flag_goodVertices") == 1) and
                                (getattr(event, "Flag_globalSuperTightHalo2016Filter") == 1) and
                                (getattr(event, "Flag_HBHENoiseFilter") == 1) and
                                (getattr(event, "Flag_HBHENoiseIsoFilter") == 1) and
                                (getattr(event, "Flag_EcalDeadCellTriggerPrimitiveFilter") == 1) and
                                (getattr(event, "Flag_BadPFMuonFilter") == 1) and
                                (getattr(event, "Flag_BadPFMuonDzFilter") == 1) and
                                (getattr(event, "Flag_eeBadScFilter") == 1)
                            )
                
            cuts = PV and METFilters# and Triggers

            if cuts and passGenSel:
                pass
            elif not(cuts):
                for sys in self.sysSource:
                    passRecoSel[sys]=False
                
                if not(passGenSel):
                    return False
                else:
                    pass


        if not (self.isMC) and not(passRecoSel['_nom']): 
            #self.totalRecoWeight=0.
            return False

        #elif self.isMC and (not(self.onlyUnc) or (self.onlyUnc and 'const' in self.onlyUnc)):
        if (self.isMC or self.isSigMC) and ((self.onlyUnc=='') or (self.onlyUnc!='' and 'const' in self.onlyUnc)):

            if not(passGenSel) and not(passRecoSel['_nom']): 
                #self.totalRecoWeight=0.
                #self.evtGenWeight=0.
                return False
                

        elif (self.isMC or self.isSigMC) and ((self.onlyUnc!='') and not('const' in self.onlyUnc)):
            sysPassFlag=False
            for sys in self.sysSource:
                if passRecoSel[sys] or passGenSel: 
                    sysPassFlag=True
                    break
            
            if not(sysPassFlag):#==False: 
                #self.totalRecoWeight=0.
                #self.evtGenWeight=0.
                    
                return False
                
        self.evtCounter+=1

        for sys in self.sysSource:

            if 'const' in sys: 
                s='_nom'
            else:
                s=sys

            
            if self.isMC: 
                self.puWeight = event.puWeight
                self.evtGenWeight = event.genWeight
                self.l1PreFireWeight = event.L1PreFiringWeight_Nom

                
                
                if len(selRecoMuons[s])>0 and passRecoSel[s]: 

                    if len(selRecoMuons[s])>1: print ("!!!!!!!!!!!!!!!Warning, extra muons leaking, check selection!!!!!!!!!!!!!!!!!!!")
                    leptonSFs = [[1.,1.,1.], [1.,1.,1.], [1.,1.,1.], [1.,1.,1.]] #dummy, done now in hist producer #self.getleptonSF( "muon", selRecoMuons[s][0] )
                    self.leptonWeight = 1.
                else: 

                    leptonSFs = [[0,0,0], [0,0,0], [0,0,0], [0,0,0]]
                    self.leptonWeight = 0.
                
                self.totalRecoWeight = self.evtGenWeight*self.puWeight*self.leptonWeight*self.l1PreFireWeight#*self.btaggingWeight

                if not sys.startswith(('_jes', '_jer')):#self.sysWeightList):
    
                    selRecoJets[sys] = selRecoJets['_nom']
                    selRecoAK4bleptjets[sys] = selRecoAK4bleptjets['_nom']
                    selRecodR_Mu_leptAK4b[sys] = selRecodR_Mu_leptAK4b['_nom']
                    selRecodPhi_AK8_leptAK4b[sys] = selRecodPhi_AK8_leptAK4b['_nom']
                    selRecodR_AK8_leptAK4b[sys] = selRecodR_AK8_leptAK4b['_nom']
                    iRecoSel[sys] = iRecoSel['_nom']
                    passRecoSel[sys] = passRecoSel['_nom']
                    selRecoMuons[sys] = selRecoMuons['_nom']
                    selRecoLeptW[sys] = selRecoLeptW['_nom']
                    selRecoLeptTop[sys] = selRecoLeptTop['_nom']
                    selRecoMET[sys] = selRecoMET['_nom']
                    selRecoLeptW_uncorr[sys] = selRecoLeptW_uncorr['_nom']
                    selRecoLeptTop_uncorr[sys] = selRecoLeptTop_uncorr['_nom']
                    selRecoMET_uncorr[sys] = selRecoMET_uncorr['_nom']

                    # PDF sets for RunIISummer20UL samples seem to be NNPDF31_nnlo_as_0.118_mc_hessian_pdfas (pdfid=325300),
                    # structure of the pdf set's array of 103 members a la: https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
                    # general structure of pdf weights = w_i,var/w_nom, where i runeth over 0-102, where [0]=>nominal from avg. over replicas, 
                    # [1-100]=> PDF eigenvectors of the covariance matrix in the parameter space, 
                    # [101,102]=> central value for (forced positive definite) a_s=[0.116,0.120]
                    if self.isSigMC:
                        

                        ############## PDF and gen weights for nominal+wtUnc runs
                        pdfWeights =  getattr( event, 'LHEPdfWeight' ) #convert to a simple numpy array to extract replica weights and thereby the variations
                        self.pdfWeightAll = np.array([pdfWeights[i] for i in range(pdfWeights.GetSize())], dtype=np.float32)
                        self.pdfWeight = pdfWeights[0]
                        
                        # 'Naively' defining pdf uncertainty symmetrically by taking the envelope defined as per Eq. 15 of https://link.springer.com/article/10.1140/epjc/s10052-015-3318-8#Sec29 
                        # a la sigma_up/dn=sigma=np.sqrt(np.sum([(arr[i]-arr[0])**2. for i in np.arange(1,103)])) from symmhessian diff calc; 
                        # however i'm not sure if this is also applicable when eig vect contributions are stored as weights, still in tests seems to give less conservative sigma estimate than 
                        # doing max(arr[i]-arr[0]) and adding that in as the up/down variation), 
                        # so saving all pdf weights anyway in case I need to unfold 1+102 times with reweighted response matrices...

                        pdf_std = np.sqrt(np.sum([(self.pdfWeightAll[i]-self.pdfWeightAll[0])**2. for i in range(1,103)]))
                        self.pdfWeightUp = self.pdfWeight + pdf_std #central + 1 sigma
                        self.pdfWeightDown = self.pdfWeight - pdf_std #central - 1 sigma
                        self.isrWeightDown = getattr( event, 'PSWeight' )[2]
                        self.fsrWeightDown = getattr( event, 'PSWeight' )[3]
                        self.isrWeightUp = getattr( event, 'PSWeight' )[0]
                        self.fsrWeightUp = getattr( event, 'PSWeight' )[1]
                        self.puWeightUp = event.puWeightUp
                        self.puWeightDown = event.puWeightDown          
                        self.l1PreFireWeightUp = event.L1PreFiringWeight_Up
                        self.l1PreFireWeightDown = event.L1PreFiringWeight_Dn


                        ############# Naive all sources up/down variations of lepton SF  for lepton weight calc.
                        self.leptonWeightAllUp = np.prod([ i[1] for i in leptonSFs]) if not self.leptonWeight==0 else 0
                        self.leptonWeightAllDown = np.prod([ i[2] for i in leptonSFs]) if not self.leptonWeight==0 else 0

                        ############# Separated, all sources up/down variations of lepton SF for lepton weight variation calc.

                        if self.leptonWeight==0:

                            self.leptonWeightTrigUp = 0.
                            self.leptonWeightTrigDown = 0.

                            self.leptonWeightIDUp = 0.
                            self.leptonWeightIDDown = 0.

                            self.leptonWeightISOUp = 0.
                            self.leptonWeightISODown = 0.

                            self.leptonWeightRecoEffUp = 0.
                            self.leptonWeightRecoEffDown = 0.
                        
                        else:

                            tempWtTrig=np.prod([i[0] for i in leptonSFs])/leptonSFs[0][0] #divide out nominal trig SF
                            tempWtID=np.prod([i[0] for i in leptonSFs])/leptonSFs[1][0] #divide out nominal ID SF
                            tempWtISO=np.prod([i[0] for i in leptonSFs])/leptonSFs[2][0] #divide out nominal ISO SF
                            tempWtRecoEff=np.prod([i[0] for i in leptonSFs])/leptonSFs[3][0] #divide out nominal reco eff SF

                            #multiply in appropriate variations
                            self.leptonWeightTrigUp = tempWtTrig*leptonSFs[0][1] 
                            self.leptonWeightTrigDown = tempWtTrig*leptonSFs[0][2]                            

                            self.leptonWeightIDUp = tempWtID*leptonSFs[1][1]
                            self.leptonWeightIDDown = tempWtID*leptonSFs[1][2]

                            self.leptonWeightISOUp = tempWtISO*leptonSFs[2][1]
                            self.leptonWeightISODown = tempWtISO*leptonSFs[2][2]

                            self.leptonWeightRecoEffUp = tempWtRecoEff*leptonSFs[3][1]
                            self.leptonWeightRecoEffDown = tempWtRecoEff*leptonSFs[3][2]



                else: 
                    self.totalRecoWeight = self.evtGenWeight*self.puWeight*self.leptonWeight*self.l1PreFireWeight#*self.btaggingWeight
            else:
                self.totalRecoWeight=1.
                self.leptonWeight = 1. #np.prod(leptonSFs)
                self.puWeight = 1.
                self.l1PreFireWeight = 1. 

            genJet = OrderedDict() if (self.isMC) else None
            recoJet = OrderedDict()
            
            if self.isMC and not( passRecoSel[s]): 
                #fill in 'missing' gen jets, ie those going into RM underflow since event didn't pass the nominal reco seleciton

                if passGenSel:
                    #### Gen Misses
                    if self.onlyUnc=='' or ((self.onlyUnc!='') and ('nom' in sys)): 
                        self.ngenEvents+=1


                    self.miss=self.miss+1
                
                    # inp. to createNsubBasis: ( AK8 jet p4, evt. pointer, particle collection, isGen)
                    genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True  )

                    # inp. to fillBranches: ( evt. pointer, obj name/type, gen/reco jet, isDummy, sys name)
                    # if not a dummy fill set an isDummy flag (arg name is actually dummy)=False
                    if sys.startswith('_nom'): 
                
                        # In this case this is a selected gen jet, so fill dummy=False
                        self.fillAK8Branches( event, 'selGenJets'+sys, genJet, len(selGenJets), False, sys )
                        self.fillOtherBranches( event, 'selGenAK4bjetLeptHem'+sys, selGenAK4bleptjets[0], 1)#, False, sys )
                        if len(selGenJets)>1: self.fillOtherBranches( event, 'selGenAK8jetLeptHem'+sys, selGenJets[1], 1)#, False, sys )
                        self.fillOtherBranches( event, 'selGenMu'+sys, selGenMuons[0], 1 )#, False, sys )
                        self.fillOtherBranches( event, 'selGenLeptW'+sys, selGenLeptW[0], 1)#, False, sys )
                        self.fillOtherBranches( event, 'selGenLeptTop'+sys, selGenLeptTop[0], 1)#, False, sys )
                        self.fillOtherBranches( event, 'selGenMET'+sys, selGenMET, 1)#, False, sys )
                        
                        self.out.fillBranch( 'selGenLeptW'+sys+'_mt', ROOT.TMath.Sqrt(2*selGenMuons[0].pt*selGenMET.Pt()*(1.-ROOT.TMath.Cos(selGenMuons[0].p4().DeltaPhi(selGenMET)))) )#, False, sys )
                
                    # filling for all jer/jes systematics separately from nominal case
                    # since in those cases the selected gen jet and smeared selected reco jet might not pass the deltaR match
                    # In this case this is not an accepted (for filling the centre of the RM) gen jet, so fill dummy=True to track misses
                    self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, len(selGenJets), True, s ) #Will use dummies in accepgen arrays as masks to selGen to extract fakes' info

                else: 
                    self.evtGenWeight=0.
                self.totalRecoWeight=0.
            
            elif passRecoSel[s]:
                self.recoLevel = self.recoLevel+1
                if self.onlyUnc=='' or ((self.onlyUnc!='') and ('nom' in sys)): 
                       
                    self.nrecoEvents+=1

                #print("Filling reco:",sys,s,True if ('const' in self.onlyUnc and 'const' in sys ) else False,'' if not('const' in self.onlyUnc  and 'const' in sys) else sys.split(self.onlyUnc)[1])
                
                recoJet['Jet'] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands', constJES=True if ('const' in self.onlyUnc and 'const' in sys ) else False, varUpDown = '' if not('const' in self.onlyUnc  and 'const' in sys) else sys.split(self.onlyUnc)[1] )#tmpRecoJets[sys][0] 
                
                self.fillAK8Branches( event, 'selRecoJets'+sys, recoJet, len(selRecoJets[s]), False, s ) #not a dummy fill so dummy=False
                self.fillOtherBranches( event, 'selRecoAK4bjetLeptHem'+s, selRecoAK4bleptjets[s][0], 1, s)#, False, sys )
                
                if len(selRecoJets[s])>1: self.fillOtherBranches( event, 'selRecoAK8jetLeptHem'+s, selRecoJets[s][1], 1, s)#, False, sys )

                self.fillOtherBranches( event, 'selRecoMu'+s, selRecoMuons[s][0], 1)#, False, sys )
                
                self.fillOtherBranches( event, 'selRecoLeptW_uncorrected'+s, selRecoLeptW_uncorr[s][0], 1)#, False, sys ) 
                self.fillOtherBranches( event, 'selRecoLeptW'+s, selRecoLeptW[s][0], 1)#, False, sys ) 

                
                self.fillOtherBranches( event, 'selRecoLeptTop_uncorrected'+s, selRecoLeptTop_uncorr[s][0], 1)#, False, sys ) 
                self.fillOtherBranches( event, 'selRecoLeptTop'+s, selRecoLeptTop[s][0], 1)#, False, sys ) 

                
                self.fillOtherBranches( event, 'selRecoMET_uncorrected'+s, selRecoMET_uncorr[s], 1)#, False, sys ) 
                self.fillOtherBranches( event, 'selRecoMET'+s, selRecoMET[s], 1)#, False, sys ) 

                
                self.out.fillBranch( 'selRecoLeptW_uncorrected'+s+'_mt', ROOT.TMath.Sqrt(2*selRecoMuons[s][0].pt*selRecoMET_uncorr[s].Pt()*(1.-ROOT.TMath.Cos(selRecoMuons[s][0].p4().DeltaPhi(selRecoMET_uncorr[s])))))  
                self.out.fillBranch( 'selRecoLeptW'+s+'_mt', ROOT.TMath.Sqrt(2*selRecoMuons[s][0].pt*selRecoMET[s].Pt()*(1.-ROOT.TMath.Cos(selRecoMuons[s][0].p4().DeltaPhi(selRecoMET[s])))))  

                

                if self.isMC: 
                    
                    deltaRmatch = False

                    if not passGenSel:  ##### fake reco
                        
                        self.fakes = self.fakes+1
                        
                        #if self.onlyUnc=='' or ((self.onlyUnc!='') and ('nom' in sys)): 

                        #    #self.nCombFakerecoEvents+=1

                        self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, len(selRecoJets[s]), True, s ) #True=dummy fill, to track fake reco, similar strategy as with gen
                        #self.evtGenWeight=0.

                    else:

                        if self.onlyUnc=='' or ((self.onlyUnc!='') and ('nom' in sys)): 

                            self.ngenEvents+=1
                        
                        genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                        
                        # consider event as missing gen if not deltaR matched 
                        if sys.startswith('_nom'): 
                            self.fillAK8Branches( event, 'selGenJets'+sys, genJet, len(selGenJets), False, sys )
                            self.fillOtherBranches( event, 'selGenAK4bjetLeptHem'+sys, selGenAK4bleptjets[0], 1)#, False, sys )
                            if len(selGenJets)>1: self.fillOtherBranches( event, 'selGenAK8jetLeptHem'+sys, selGenJets[1], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenMu'+sys, selGenMuons[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenLeptW'+sys, selGenLeptW[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenLeptTop'+sys, selGenLeptTop[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenMET'+sys, selGenMET, 1)#, False, sys )
                            self.out.fillBranch( 'selGenLeptW'+sys+'_mt', ROOT.TMath.Sqrt(2*selGenMuons[0].pt*selGenMET.Pt()*(1.-ROOT.TMath.Cos(selGenMuons[0].p4().DeltaPhi(selGenMET)))))

                        #need to delta R match for true reco and accepted gen
                        if self.DrRapPhi( recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() ) < 0.4 :
                            deltaRmatch = True
                            self.response= self.response+1
                        
                        self.alsoDeltaRMatchedB=1 if self.DrRapPhi( selRecoAK4bleptjets[s][0].p4(),selGenAK4bleptjets[0].p4())<0.2 else 0

                        
                        if deltaRmatch:
                            # fill only if deltaR matched, for eventual response matrix filling
                            self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, len(selGenJets), False, s )    
                            self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, len(selRecoJets[s]), False, s )
                            self.nAccepgenEvents+=1
                            self.nTruerecoEvents+=1  
                        else:
                            self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, len(selGenJets), True, s )    #fill dummy=True to indicate fakes and misses 
                            self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, len(selRecoJets[s]), True, s )
                            self.fakes = self.fakes+1
                            self.miss = self.miss+1

            #recoleptHemBjets=[]                    
            recoleptHemBjets = [x for x in selRecoAK4bjets[s] if not (x in selRecoAK4bleptjets[s])]
            #self.recoEventCategory=len(selRecoAK4bjets) if ((len(selRecoAK4bjets)>0 and len(selRecoAK4bleptjets[sys])>0)) else -1
            
            if sys.startswith('_nom') and self.isMC:
                genleptHemBjets=[]
                genleptHemBjets = [x for x in selGenAK4bjets if not( x in selGenAK4bleptjets)] 
                #self.genEventCategory=len(selGenAK4bjets) if ((len(selGenAK4bjets)>0 and len(selGenAK4bleptjets)>0)) else -1
            
            if not self.isMC:
                #self.out.fillBranch( 'recoEventCategory'+sys, self.recoEventCategory if passRecoSel[sys] else -1  )
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood'))# if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoBtags'+sys, len(selRecoAK4bjets[s]) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoAK4s'+sys, nselRecoAK4[s] if passRecoSel[sys] else 0)
                #self.out.fillBranch( 'nRecoCrackElectrons'+sys, nselRecoCrackEl if passRecoSel[sys] else 0)
                #self.out.fillBranch( 'nRecoAK8s'+sys, len(selRecoJets[sys]) if passRecoSel[sys] else 0)

                self.out.fillBranch( 'nRecoLeptBtags'+sys, len(recoleptHemBjets) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'selRecoLeptHemDeltaR'+sys, selRecodR_Mu_leptAK4b[s] if passRecoSel[sys] else 929.)
                self.out.fillBranch( 'selRecoLeptHemDeltaPhi'+sys, selRecodPhi_Mu_leptAK4b[s] if passRecoSel[sys] else 929.)
                self.out.fillBranch( 'selRecoLeptHemDeltaRap'+sys, selRecodRap_Mu_leptAK4b[s] if passRecoSel[sys] else 929.)

                self.out.fillBranch( 'selRecoLeptHemAK8DeltaPhi'+sys, selRecodPhi_AK8_leptAK4b[s] if passRecoSel[sys] else 929.)
                self.out.fillBranch( 'selRecoLeptHemAK8DeltaR'+sys, selRecodR_AK8_leptAK4b[s] if passRecoSel[sys] else 929.)


                self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                
                self.out.fillBranch( 'passRecoSel'+sys, 1 if passRecoSel[sys] else 0)
                self.out.fillBranch( 'passHLT_Mu50'+sys, 1 if getattr(event, 'HLT_Mu50')==1 else 0)# and passRecoSel[sys] else 0 ) 
                
                try: 
                    self.out.fillBranch( 'passHLT_TkMu50'+sys, 1 if getattr(event, 'HLT_TkMu50')==1 else 0)# and passRecoSel[sys] else 0 ) 
                except RuntimeError:
                    self.out.fillBranch( 'passHLT_TkMu50'+sys, 0 ) #set unavailable triggers to 0, even enters histo if pass reco sel and mu50 only at that point

                try: 
                    self.out.fillBranch( 'passHLT_OldMu100'+sys, 1 if getattr(event, 'HLT_OldMu100')==1 else 0)# and passRecoSel[sys] else 0 ) 
                except RuntimeError:
                    self.out.fillBranch( 'passHLT_OldMu100'+sys, 0 ) 

                try: 
                    self.out.fillBranch( 'passHLT_TkMu100'+sys, 1 if getattr(event, 'HLT_TkMu100')==1 else 0)# and passRecoSel[sys] else 0 ) 
                except RuntimeError:
                    self.out.fillBranch( 'passHLT_TkMu100'+sys, 0 )                     


            else:
                #self.out.fillBranch( 'recoEventCategory'+sys, self.recoEventCategory if passRecoSel[sys] else -1)

                if not('const' in sys):
                    self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') if passRecoSel[sys] else 0)
                    self.out.fillBranch( 'nRecoBtags'+sys, len(selRecoAK4bjets[s]) if passRecoSel[sys] else 0)
                    self.out.fillBranch( 'nRecoAK4s'+sys, nselRecoAK4[s] if passRecoSel[sys] else 0)
                    #self.out.fillBranch( 'nRecoCrackElectrons'+sys, nselRecoCrackEl if passRecoSel[sys] else 0)
                    #self.out.fillBranch( 'nRecoAK8s'+sys, len(selRecoJets[sys]) if passRecoSel[sys] else 0)

                    self.out.fillBranch( 'nRecoLeptBtags'+sys, len(recoleptHemBjets) if passRecoSel[sys] else 0)
                    self.out.fillBranch( 'selRecoLeptHemAK8DeltaPhi'+sys, selRecodPhi_AK8_leptAK4b[sys] if passRecoSel[sys] else 929.)                       
                    self.out.fillBranch( 'selRecoLeptHemAK8DeltaR'+sys, selRecodR_AK8_leptAK4b[sys] if passRecoSel[sys] else 929.)                       
                    #self.out.fillBranch( 'FlagRecoLeptHemBjet'+sys, 1 if (passRecoSel[sys] and len(recoleptHemBjets)>0) else 0)
                    self.out.fillBranch( 'FlagDeltaRMatchedBjets'+sys, self.alsoDeltaRMatchedB if (passRecoSel[sys] and passGenSel) else 0)
                    
                    self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'passRecoSel'+sys, 1 if passRecoSel[sys] else 0)
                    
                    self.out.fillBranch( 'puWeightNom'+sys, self.puWeight if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'l1prefiringWeightNom'+sys, self.l1PreFireWeight if passRecoSel[sys] else 0.)
                    #self.out.fillBranch( 'btagWeightNom'+sys, self.btaggingWeight if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightNom'+sys, self.leptonWeight if passRecoSel[sys] else 0.)

                    #self.out.fillBranch( 'genEventCategory'+sys, self.genEventCategory if passGenSel or passRecoSel[sys] else -1)
                    self.out.fillBranch( 'evtGenWeight'+sys, self.evtGenWeight if passGenSel or passRecoSel[sys] else 0.) 
                    self.out.fillBranch( 'passGenSel'+sys, 1 if passGenSel else 0)#1 if passGenSel or passRecoSel[sys] else 0
                    
                    self.out.fillBranch( 'selRecoLeptHemDeltaR'+sys, selRecodR_Mu_leptAK4b[sys] if passRecoSel[sys] else 929.)            
                    self.out.fillBranch( 'selRecoLeptHemDeltaPhi'+sys, selRecodPhi_Mu_leptAK4b[sys] if passRecoSel[sys] else 929.)            
                    self.out.fillBranch( 'selRecoLeptHemDeltaRap'+sys, selRecodRap_Mu_leptAK4b[sys] if passRecoSel[sys] else 929.) 

                if sys.startswith('_nom'): 
                    self.out.fillBranch( 'passHLT_Mu50'+sys, 1 if getattr(event, 'HLT_Mu50')==1 else 0)# and passRecoSel[sys] else 0 ) 
                                    
                    try: 
                        self.out.fillBranch( 'passHLT_TkMu50'+sys, 1 if getattr(event, 'HLT_TkMu50')==1 else 0)# and passRecoSel[sys] else 0 ) 
                    except RuntimeError:
                        self.out.fillBranch( 'passHLT_TkMu50'+sys, 0 ) 

                    try: 
                        self.out.fillBranch( 'passHLT_OldMu100'+sys, 1 if getattr(event, 'HLT_OldMu100')==1 else 0)# and passRecoSel[sys] else 0 ) 
                    except RuntimeError:
                        self.out.fillBranch( 'passHLT_OldMu100'+sys, 0 ) 

                    try: 
                        self.out.fillBranch( 'passHLT_TkMu100'+sys, 1 if getattr(event, 'HLT_TkMu100')==1 else 0)# and passRecoSel[sys] else 0 ) 
                    except RuntimeError:
                        self.out.fillBranch( 'passHLT_TkMu100'+sys, 0 )   

                    self.out.fillBranch( 'nGenBtags'+sys, len(selGenAK4bjets) if passGenSel else 0)
                    self.out.fillBranch( 'nGenAK4s'+sys, nselGenAK4 if passGenSel else 0)
                    #self.out.fillBranch( 'nGenCrackElectrons'+sys, nselGenCrackEl if passGenSel else 0)
                    #self.out.fillBranch( 'nGenAK8s'+sys, len(selGenJets) if passGenSel else 0)
                    self.out.fillBranch( 'nGenLeptBtags'+sys, len(genleptHemBjets) if passGenSel else 0)
                    #self.out.fillBranch( 'FlagGenLeptHemBjet'+sys, 1 if (passGenSel and len(genleptHemBjets)>0) else 0)
                    
                    
                    self.out.fillBranch( 'selGenLeptHemDeltaR'+sys, selGendR_Mu_leptAK4b if passGenSel or passRecoSel[sys] else 929.)
                    self.out.fillBranch( 'selGenLeptHemDeltaPhi'+sys, selGendPhi_Mu_leptAK4b if passGenSel or passRecoSel[sys] else 929.)
                    self.out.fillBranch( 'selGenLeptHemDeltaRap'+sys, selGendRap_Mu_leptAK4b if passGenSel or passRecoSel[sys] else 929.)
                    self.out.fillBranch( 'selGenLeptHemAK8DeltaPhi'+sys, selGendPhi_AK8_leptAK4b if passGenSel or passRecoSel[sys] else 929.)
                    self.out.fillBranch( 'selGenLeptHemAK8DeltaR'+sys, selGendR_AK8_leptAK4b if passGenSel or passRecoSel[sys] else 929.)
                
                    

                if self.isSigMC and self.onlyUnc=='' and not (self.effMapOnly):
                    self.out.fillBranch( 'pdfWeightNom'+sys, self.pdfWeight if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'pdfWeightAll'+sys, self.pdfWeightAll if passGenSel or passRecoSel[sys] else np.zeros((103,),dtype=np.float32))
                    self.out.fillBranch( 'pdfWeightUp'+sys, self.pdfWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'pdfWeightDown'+sys, self.pdfWeightDown if passGenSel or passRecoSel[sys] else 0.)
                    
                    self.out.fillBranch( 'isrWeightUp'+sys, self.isrWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'isrWeightDown'+sys, self.isrWeightDown if passGenSel or passRecoSel[sys] else 0.)
                    
                    self.out.fillBranch( 'fsrWeightUp'+sys, self.fsrWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'fsrWeightDown'+sys, self.fsrWeightDown if passGenSel or passRecoSel[sys] else 0.)
                    
                    self.out.fillBranch( 'puWeightUp'+sys, self.puWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'puWeightDown'+sys, self.puWeightDown if passRecoSel[sys] else 0.)

                    self.out.fillBranch( 'l1prefiringWeightUp'+sys, self.l1PreFireWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'l1prefiringWeightDown'+sys, self.l1PreFireWeightDown if passRecoSel[sys] else 0.)

                    #self.out.fillBranch( 'btagWeightUp'+sys, self.puWeightUp if passRecoSel[sys] else 0.)
                    #self.out.fillBranch( 'btagWeightDown'+sys, self.puWeightDown if passRecoSel[sys] else 0.)

                    self.out.fillBranch( 'leptonWeightAllUp'+sys, self.leptonWeightAllUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightAllDown'+sys, self.leptonWeightAllDown if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightISOUp'+sys, self.leptonWeightISOUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightISODown'+sys, self.leptonWeightISODown if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightIDUp'+sys, self.leptonWeightIDUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightIDDown'+sys, self.leptonWeightIDDown if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightTrigUp'+sys, self.leptonWeightTrigUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightTrigDown'+sys, self.leptonWeightTrigDown if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightRecoEffUp'+sys, self.leptonWeightRecoEffUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'leptonWeightRecoEffDown'+sys, self.leptonWeightRecoEffDown if passRecoSel[sys] else 0.)

        self.recoWeightSum += (self.totalRecoWeight if passRecoSel['_nom'] else 0.)
        self.genWeightSum += (self.evtGenWeight if passGenSel or passRecoSel['_nom'] else 0.)# or passRecoSel['_nom'] else 0.)

        
        if self.evtCounter%500==0 or self.evtCounter==1:# and self.onlyUnc=='':
            print("event.event","evtCounter","nrecoEvents","ngenEvents", "passGenSel", "passRecoSel['_nom']") #"n_CAccepgenEvents","n_CTruerecoEvents")
            print(event.event,self.evtCounter,self.nrecoEvents,self.ngenEvents, passGenSel, passRecoSel['_nom']) #self.nAccepgenEvents,self.nTruerecoEvents)
            print("recoWeightSum=",self.recoWeightSum)
            print("genWeightSum=",self.genWeightSum)
            #if passRecoSel['_nom']:
            #    print('reco mSDs', getattr(recoJet['Jet']['jet'],'mSD_nom'), 
            #          getattr(recoJet['Jet']['jet'],'msoftdrop_nom'),
            #          getattr(recoJet['Jet']['jet'],'msoftdrop'),
            #          getattr(recoJet['Jet']['jet'],'msoftdrop_raw'), 
            #          getattr(recoJet['Jet']['jet'],'msoftdrop_nom_PUPPICorred'),#p4().Pt(),recoJet['Jet']['jet'].p4().Eta(),recoJet['Jet']['jet'].p4().Phi(),recoJet['Jet']['jet'].p4().M())
            #          )#getattr(recoJet['Jet']['jet'],'mSD_nom')) #p4().Pt(),recoJet['Jet']['jet'].p4().Eta(),recoJet['Jet']['jet'].p4().Phi(),recoJet['Jet']['jet'].p4().M())
            #    print('reco masses', getattr(recoJet['Jet']['jet'],'mass_nom'), 
            #          getattr(recoJet['Jet']['jet'],'mass'))

            #if passGenSel:
            #    print('gen mSDs', getattr(genJet['Jet']['jet'],'msoftdrop'))#,
            #    print('gen masses', getattr(genJet['Jet']['jet'],'mass'))


        return True


    ################ Selection of objects for reco/gen done below ################ 
    # track n-btags, hence track gen/recoEvtCategories (ie, only had hem b tag or both lep+had hem)
    # change class-level flags below to figure out what to fill and when
    # ie, flags based on btags decides which evt sel of the above two the event falls in
    # also, fill the had hem b-tagged/flavourMatched AK4 branches for later W/top selecting
     

    #############################################################################
    def recoSelection( self, event, sysUnc=[]  ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' )) #anti-kT R=0.8 PF-PUPPI jets
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet')) #anti-kT R=0.4 PF-CHS jets
        
        ########### Lepton selection ###############
        # only loose selection for electrons, to veto on (cutbased id == 1 for veto, 2 for loose and ...)
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 1 and abs(x.eta)<self.maxElectronEta ] 

        # applying effectively the isHighPt selection on muons here (since that flag is not explicitly available in nano)
        # dRapils of highPt selectors used implicitly provided in: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#HighPt_Muon
        # (the isGlobal,isTracker,dxy and dz cuts explicitly used previously are redundant in principle since isHighPt and, therefore I assume highPtId==2, requires this) 
        # The highPtId is used since the UL recommendations (e.g., 2017): https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017#High_pT_above_120_GeV 
        # strongly recommend it, I just use that and the lower eff. (at high pT) tightID req. from before is removed, along with the redundant selectors
        # We are also requiring a tracker-based relative isolation <0.3, where nano documentation says dR=0.3 for highPt, (instead of a PF ISO requirement)
        # and UL/high pT recommendations (2nd link above) state that PF based quantities (like isolation) are not to be used unless really necessary
        # Further, relevant to the high Pt ID, 'TrackerISO' is imposed on top of the high pT ID for the centrally provided medium/hight pT muon ISO efficiencies 
        # so a tracker based isolation is used (which also peaks at <=0.1, roughly even as low as 0.05 for muons passing the highPtID selection and kinematic reqs. on pT&eta). 
        # Further for ID/ISO effs. @high pT medium effs are meant to be extended, and in the latter a high Pt WP is avialable, so highPtID is the only way to go!

        #update to tight tkRelIso cut,(<0.05 is tight as per: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection#HighPt_Muon), updating muon SF files as necessary to deal with this
        #The tracker based relative isolation is computed using tracks from the leading PV in the event, which are in a cone of size DeltaR < 0.3. * sum(pT(tracker tracks from PV))/pT(\mu)
        #highPtId==2 => global highPtId (which is then also including the tracker)
        recoMuons = [x for x in muons if x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta and x.highPtId==2 and x.tkRelIso<0.05]
        #  and abs(x.dxy)<0.2 and abs(x.dz)<0.5  and x.isGlobal] #redundant

        #loose selection for muons, to veto on
        recoVetoMuons = [x for x in muons if x.pt > self.minLooseMuonPt and abs(x.eta) < self.maxMuonEta and x.looseId==1 and not(x in recoMuons)]

        
        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoVetoMuons.sort(key=lambda x:x.pt,reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)
        
        nleptons = len(recoMuons)+len(recoElectrons)+len(recoVetoMuons)

        #if nleptons>=2: 
        leptonVetoFlag = nleptons>=2 
        #if leptonVetoFlag: print(event.event,'reco',leptonVetoFlag,nleptons,len(recoMuons),len(recoElectrons),len(recoLooseMuons))
        ############################################

        
        
        if len(AK8jets)!=0:
            for ijets in AK8jets: 
                #jecFactor_nom = getattr( ijets, 'pt_nom' ) / getattr( ijets, 'pt_raw' )

                ijets.rapidity = ijets.p4().Rapidity()  # add rapidity attr. 
                #ijets.mSD_nom = ijets.msoftdrop_raw*jecFactor_nom # add mSD_nom attr., corrected explicitly 


        recoAK8jets = {}
        passSel = {}
        iSel = {}
        recodR_Mu_leptAK4b = {}
        recodRap_Mu_leptAK4b = {}
        recodPhi_Mu_leptAK4b = {}
        recodPhi_AK8_leptAK4b = {}
        recodR_AK8_leptAK4b = {}
        recoAK4jets = {}
        recoAK4bjets = {}
        recoAK4bLeptJets = {}
        recoAK4LeptJets = {}

        recoMuons_dict = {} 
        recoleptW_dict = {} 
        recoleptW_uncorr_dict = {} 
        recoleptTop_dict = {} 
        recoleptTop_uncorr_dict = {} 
        MET = {} 
        MET_uncorr = {} 
        MET_dict = {}
        MET_uncorr_dict = {}
        nAK4_dict = {} 


        for sys in self.sysSource:
            met = Object(event, 'MET_T1')    

             
            MET_Pt_corrected, MET_Phi_corrected = METXYCorr_Met_MetPhi(getattr( met, 'pt'+ (sys if not(sys.startswith(('_nom', '_const'))) else '') ),
                                                                       getattr( met, 'phi'+ (sys if not(sys.startswith(('_nom', '_const'))) else '') ),
                                                                       event.run, 
                                                                       self.year, 
                                                                       self.isMC, 
                                                                       getattr(event, 'PV_npvs'), 
                                                                       True, False) #isUL, ispuppi
            

            if sys.startswith(self.sysWeightList) or 'const' in sys: sys = '_nom'

            ################### Basic AK8 jet selection ############################
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs(x.rapidity) < self.maxJetAK8Rap and (x.jetId >=2)]
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)

            AK8HT = sum( [ getattr( x, 'pt'+sys ) for x in recoAK8jets[sys] ] )

            if len(recoAK8jets[sys])!=0:
                #not used in final selection, prev. attempt at constructing softdrop mass based on JERC corrected AK8 kinematics, now we use JEC (re)corrected subjets's p4().M() and AK8 JER smearing factor to
                #obtain the msoftdrop+'sys' (sys=  nom, jerUp/Down...) via the fatjetuncertainties script in nanoAOD-tools with UL patches I've finally introduced, since they weren't centrally made available
                for jet in recoAK8jets[sys]:
                    jercFactor = getattr( jet, 'pt'+sys ) / getattr( jet, 'pt_raw' ) 
                    setattr(jet, 'mSD'+sys,  ( getattr( jet, 'msoftdrop_raw' ) * jercFactor )  )

            ########################################################################

            ################### MET  #######################
            MET[sys] = ROOT.TLorentzVector()
            MET[sys].SetPtEtaPhiM(MET_Pt_corrected, 0., MET_Phi_corrected, 0.)

            MET_uncorr[sys] = ROOT.TLorentzVector()
            MET_uncorr[sys].SetPtEtaPhiM(getattr( met, 'pt'+ (sys if not(sys.startswith(('_nom', '_const'))) else '') ), 0., getattr( met, 'phi'+ (sys if not(sys.startswith(('_nom', '_const'))) else '') ), 0.)
            ################################################


            ############ Basic AK4 selection ###############
            recoAK4jets[sys] = [ x for x in jets if getattr( x, 'pt'+sys ) > self.minJetPt and abs(x.eta) < self.maxJetEta and (x.jetId>=2)]# and x.btagDeepFlavB > self.minBDisc]  
            #puID not included for now (but, ref. in case needed: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL)
            recoAK4jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)
            nAK4 = len(recoAK4jets[sys])

            recoAK4bjets[sys] = [ x for x in recoAK4jets[sys] if x.btagDeepFlavB > self.minBDisc]
            recoAK4bjets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)

            ################################################
            
            ################################################################ Applying selection ########################################################################
            
            # Select the event as passing the top selection (ie, in a boosted semileptonic ttbar topology):
            # >=1 boosted heavy AK8 (so as not to reject events where the leptonic hem. is reconstructed as an AK8, 
            # which does slightly increases chances of letting in other ttbar channels...),
            # ==1 high pT mu, ==1 btag (w/ exactly one jet in the lept hem., requiring it to also be b-tagged
            # apart from being the closest AK4 to the mu), no further leptons passing loose/highPt criteria, 
            # large MET, large lept. W pT. The lept. W and b-cand define the hemisphere and the 'tag' in the TnP

            if not(leptonVetoFlag):
                if not self.effMapOnly:
                    passSel[sys], iSel[sys], recodR_Mu_leptAK4b[sys], recodPhi_Mu_leptAK4b[sys], recodRap_Mu_leptAK4b[sys], recodPhi_AK8_leptAK4b[sys], recodR_AK8_leptAK4b[sys], recoAK4bLeptJets[sys], recoAK4LeptJets[sys]  = self.WtopSelection( False, event, 
                                                                                                                                                                                                                           recoMuons, recoElectrons, 
                                                                                                                                                                                                                           recoAK4jets[sys], recoAK4bjets[sys], 
                                                                                                                                                                                                                           recoAK8jets[sys], MET[sys], sys) 

                    
                else:
                    passSel[sys], iSel[sys], recodR_Mu_leptAK4b[sys], recodPhi_Mu_leptAK4b[sys], recodRap_Mu_leptAK4b[sys], recodPhi_AK8_leptAK4b[sys], recodR_AK8_leptAK4b[sys], recoAK4bLeptJets[sys], recoAK4LeptJets[sys]  = self.WtopSelection_effMap( False, event, 
                                                                                                                                                                                                                              recoMuons, recoElectrons, 
                                                                                                                                                                                                                              recoAK4jets[sys], recoAK4bjets[sys], 
                                                                                                                                                                                                                              recoAK8jets[sys], MET[sys], sys)                             
                
            else:
                passSel[sys], iSel[sys], recodR_Mu_leptAK4b[sys], recodPhi_Mu_leptAK4b[sys], recodRap_Mu_leptAK4b[sys], recodPhi_AK8_leptAK4b[sys], recodR_AK8_leptAK4b[sys], recoAK4bLeptJets[sys], recoAK4LeptJets[sys]  = False, None, 929, 929, 929, 929, 929, [], []



            ############################################################################################################################################################
            
            
            
            if passSel[sys]:

                if len(recoMuons)>0:

                    recoleptWd = [recoMuons[0].p4()+MET[sys] ]
                    recoleptWd_uncorr = [recoMuons[0].p4()+MET_uncorr[sys] ]
                    recoleptTopd = [recoMuons[0].p4()+MET[sys]+recoAK4bLeptJets[sys][0].p4() ]
                    recoleptTopd_uncorr = [recoMuons[0].p4()+MET_uncorr[sys]+recoAK4bLeptJets[sys][0].p4() ]

                    recoMuons_dict[sys] = recoMuons
                    recoleptW_dict[sys] = recoleptWd
                    recoleptW_uncorr_dict[sys] = recoleptWd_uncorr
                    recoleptTop_dict[sys] = recoleptTopd
                    recoleptTop_uncorr_dict[sys] = recoleptTopd_uncorr

                MET_dict[sys] = MET[sys]
                MET_uncorr_dict[sys] = MET_uncorr[sys]
                nAK4_dict[sys] = nAK4

            else:

                recoMuons_dict[sys] = []
                recoleptW_dict[sys] = []
                recoleptW_uncorr_dict[sys] = []
                recoleptTop_dict[sys] = []
                recoleptTop_uncorr_dict[sys] = []
                MET_dict[sys] = []
                MET_uncorr_dict[sys] = []



        ######################### Weights ############################

        # all these weights are to plot nominal, control histos below
        # weights for array masks are done in analyzer
        # deprecated so commenting out
        weight=1.
        
        if self.isMC:  
            

            #### Lepton Weights (also reconstruct leptonic W and top objects if possible)####
            #if self.isMC:
            if len(recoMuons)>0:
                
                self.leptonWeight = 1
            else: 
                

                leptonSFs = [[0, 0, 0,],[0, 0, 0,],[0, 0, 0,],[0, 0, 0,]]
                self.leptonWeight=0.
             
            ##########################################################


            self.puWeight = event.puWeight
            self.l1PreFireWeight = event.L1PreFiringWeight_Nom
            self.evtGenWeight = event.genWeight 

            #### Applying all reco weights for control histos ####, not being plotted anymore in skimmer by default though
            weight =  self.evtGenWeight  * self.puWeight *  self.l1PreFireWeight * self.leptonWeight #btagweights['_nom'] # * self.topweight * self.btaggingWeight            

        else:
            leptonSFs = [1, 1, 1, 1.]
            
            weight = 1.
            ##############################################################

        if self.controlHistos:

            # filling test histos
            if not self.onlyUnc:

                #### Checking not selected without weights
                getattr( self, 'nPVs_noSelnoWeight' ).Fill( getattr( event, 'PV_npvsGood') )
                getattr( self, 'nleps_noSelnoWeight' ).Fill( nleptons )
                
                for imuon in recoMuons:
                    getattr( self, 'muons_pt_noSelnoWeight' ).Fill( imuon.pt )
                    getattr( self, 'muons_eta_noSelnoWeight' ).Fill( imuon.eta )
                    getattr( self, 'muons_y_noSelnoWeight' ).Fill( imuon.p4().Rapidity() )
                    getattr( self, 'muons_phi_noSelnoWeight' ).Fill( imuon.phi )

                for iele in recoElectrons:
                    getattr( self, 'eles_pt_noSelnoWeight' ).Fill( iele.pt )
                    getattr( self, 'eles_eta_noSelnoWeight' ).Fill( iele.eta )
                    getattr( self, 'eles_y_noSelnoWeight' ).Fill( iele.p4().Rapidity() )
                    getattr( self, 'eles_phi_noSelnoWeight' ).Fill( iele.phi )

                getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recoAK8jets['_nom']) )
                getattr( self, 'HT_noSelnoWeight' ).Fill( AK8HT )

                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt_nom )
                    getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                    getattr( self, 'AK8jets_y_noSelnoWeight' ).Fill( ijet.rapidity )
                    getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                    getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.mass_nom)#if self.evtSelection.startswith('W') else ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI )
                    getattr( self, 'AK8jets_msoftdrop_noSelnoWeight' ).Fill( ijet.msoftdrop)
                
                
                getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4jets) )
                getattr( self, 'nAK4bjets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
                getattr( self, 'nAK4leptjets_noSelnoWeight' ).Fill( len(recoAK4LeptJets['_nom']) )
                getattr( self, 'nAK4bleptjets_noSelnoWeight' ).Fill( len(recoAK4bLeptJets['_nom']) )
                getattr( self, 'dR_Mu_AK4blept_noSelnoWeight' ).Fill( recodR_Mu_leptAK4b['_nom'] )
                getattr( self, 'dPhi_AK8_AK4blept_noSelnoWeight' ).Fill( recodPhi_AK8_leptAK4b['_nom'] )

                for ijet in recoAK4jets:
                    getattr( self, 'AK4jets_pt_noSelnoWeight' ).Fill( ijet.pt )
                    getattr( self, 'AK4jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                    getattr( self, 'AK4jets_y_noSelnoWeight' ).Fill( ijet.p4().Rapidity() )
                    getattr( self, 'AK4jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                    getattr( self, 'AK4jets_mass_noSelnoWeight' ).Fill( ijet.mass )

                for ijet in recoAK4bjets:
                    getattr( self, 'AK4btaggedjets_pt_noSelnoWeight' ).Fill( ijet.pt )
                    getattr( self, 'AK4btaggedjets_eta_noSelnoWeight' ).Fill( ijet.eta )
                    getattr( self, 'AK4btaggedjets_y_noSelnoWeight' ).Fill( ijet.p4().Rapidity())
                    getattr( self, 'AK4btaggedjets_phi_noSelnoWeight' ).Fill( ijet.phi )
                    getattr( self, 'AK4btaggedjets_mass_noSelnoWeight' ).Fill( ijet.mass )

                for ijet in recoAK4LeptJets['_nom']:
                    getattr( self, 'AK4leptjets_pt_noSelnoWeight' ).Fill( ijet.pt )
                    getattr( self, 'AK4leptjets_eta_noSelnoWeight' ).Fill( ijet.eta )
                    getattr( self, 'AK4leptjets_y_noSelnoWeight' ).Fill( ijet.p4().Rapidity() )
                    getattr( self, 'AK4leptjets_phi_noSelnoWeight' ).Fill( ijet.phi )
                    getattr( self, 'AK4leptjets_mass_noSelnoWeight' ).Fill( ijet.mass )

                if len( recoAK4bLeptJets['_nom'])>=1:
                    getattr( self, 'AK4leptbjet_pt_noSelnoWeight' ).Fill( recoAK4bLeptJets['_nom'][0].pt )
                    getattr( self, 'AK4leptbjet_eta_noSelnoWeight' ).Fill( recoAK4bLeptJets['_nom'][0].eta )
                    getattr( self, 'AK4leptbjet_y_noSelnoWeight' ).Fill( recoAK4bLeptJets['_nom'][0].p4().Rapidity())
                    getattr( self, 'AK4leptbjet_phi_noSelnoWeight' ).Fill( recoAK4bLeptJets['_nom'][0].phi )
                    getattr( self, 'AK4leptbjet_mass_noSelnoWeight' ).Fill( recoAK4bLeptJets['_nom'][0].mass )



                getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )
                for lW in recoleptW:
                    getattr( self, 'leptonicW_pt_noSelnoWeight' ).Fill( lW.Pt() )
                    getattr( self, 'leptonicW_eta_noSelnoWeight' ).Fill( lW.Eta() )
                    getattr( self, 'leptonicW_y_noSelnoWeight' ).Fill( lW.Rapidity() )
                    getattr( self, 'leptonicW_phi_noSelnoWeight' ).Fill( lW.Phi() )
                    getattr( self, 'leptonicW_mass_noSelnoWeight' ).Fill( lW.M() )
                for lTop in recoleptTop:
                    getattr( self, 'leptonicTop_pt_noSelnoWeight' ).Fill( lTop.Pt() )
                    getattr( self, 'leptonicTop_eta_noSelnoWeight' ).Fill( lTop.Eta() )
                    getattr( self, 'leptonicTop_y_noSelnoWeight' ).Fill( lTop.Rapidity() )
                    getattr( self, 'leptonicTop_phi_noSelnoWeight' ).Fill( lTop.Phi() )
                    getattr( self, 'leptonicTop_mass_noSelnoWeight' ).Fill( lTop.M() )
                    

                #### Checking not selected with weights
                getattr( self, 'nPVs_noSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
                getattr( self, 'nleps_noSel' ).Fill( nleptons , weight)
                
                for imuon in recoMuons:
                    getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt , weight)
                    getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta , weight)
                    getattr( self, 'muons_y_noSel' ).Fill( imuon.p4().Rapidity() , weight)
                    getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi , weight)

                for iele in recoElectrons:
                    getattr( self, 'eles_pt_noSel' ).Fill( iele.pt , weight)
                    getattr( self, 'eles_eta_noSel' ).Fill( iele.eta , weight)
                    getattr( self, 'eles_y_noSel' ).Fill( iele.p4().Rapidity() , weight)
                    getattr( self, 'eles_phi_noSel' ).Fill( iele.phi , weight)

                getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets['_nom']) , weight)
                getattr( self, 'HT_noSel' ).Fill( AK8HT , weight)

                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt_nom , weight)
                    getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta , weight)
                    getattr( self, 'AK8jets_y_noSel' ).Fill( ijet.rapidity , weight)
                    getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi , weight)
                    getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.mass_nom, weight)# if self.evtSelection.startswith('W') else ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI , weight)
                    getattr( self, 'AK8jets_msoftdrop_noSel' ).Fill( ijet.msoftdrop, weight)

                getattr( self, 'nAK4jets_noSel' ).Fill( len(recoAK4jets), weight )
                getattr( self, 'nAK4bjets_noSel' ).Fill( len(recoAK4bjets) , weight)
                getattr( self, 'nAK4leptjets_noSel' ).Fill( len(recoAK4LeptJets['_nom']), weight )
                getattr( self, 'nAK4bleptjets_noSel' ).Fill( len(recoAK4bLeptJets['_nom']), weight )
                getattr( self, 'dR_Mu_AK4blept_noSel' ).Fill( recodR_Mu_leptAK4b['_nom'], weight )
                getattr( self, 'dPhi_AK8_AK4blept_noSel' ).Fill( recodPhi_AK8_leptAK4b['_nom'], weight )
                
                for ijet in recoAK4jets:
                    getattr( self, 'AK4jets_pt_noSel' ).Fill( ijet.pt , weight)
                    getattr( self, 'AK4jets_eta_noSel' ).Fill( ijet.eta , weight)
                    getattr( self, 'AK4jets_y_noSel' ).Fill( ijet.p4().Rapidity() , weight)
                    getattr( self, 'AK4jets_phi_noSel' ).Fill( ijet.phi , weight)
                    getattr( self, 'AK4jets_mass_noSel' ).Fill( ijet.mass , weight)
                    
                for ijet in recoAK4bjets:
                    getattr( self, 'AK4btaggedjets_pt_noSel' ).Fill( ijet.pt , weight)
                    getattr( self, 'AK4btaggedjets_eta_noSel' ).Fill( ijet.eta , weight)
                    getattr( self, 'AK4btaggedjets_y_noSel' ).Fill( ijet.p4().Rapidity() , weight)
                    getattr( self, 'AK4btaggedjets_phi_noSel' ).Fill( ijet.phi , weight)
                    getattr( self, 'AK4btaggedjets_mass_noSel' ).Fill( ijet.mass , weight)

                for ijet in recoAK4LeptJets['_nom']:
                    getattr( self, 'AK4leptjets_pt_noSel' ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4leptjets_eta_noSel' ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4leptjets_y_noSel' ).Fill( ijet.p4().Rapidity(), weight )
                    getattr( self, 'AK4leptjets_phi_noSel' ).Fill( ijet.phi, weight )
                    getattr( self, 'AK4leptjets_mass_noSel' ).Fill( ijet.mass, weight )

                if len( recoAK4bLeptJets['_nom'])>=1:
                    getattr( self, 'AK4leptbjet_pt_noSel' ).Fill( recoAK4bLeptJets['_nom'][0].pt , weight)
                    getattr( self, 'AK4leptbjet_eta_noSel' ).Fill( recoAK4bLeptJets['_nom'][0].eta , weight)
                    getattr( self, 'AK4leptbjet_y_noSel' ).Fill( recoAK4bLeptJets['_nom'][0].p4().Rapidity(), weight)
                    getattr( self, 'AK4leptbjet_phi_noSel' ).Fill( recoAK4bLeptJets['_nom'][0].phi , weight)
                    getattr( self, 'AK4leptbjet_mass_noSel' ).Fill( recoAK4bLeptJets['_nom'][0].mass , weight)


                getattr( self, 'METPt_noSel' ).Fill( MET.Pt() , weight)
                
                for lW in recoleptW:
                    getattr( self, 'leptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                    getattr( self, 'leptonicW_eta_noSel' ).Fill( lW.Eta() , weight)
                    getattr( self, 'leptonicW_y_noSel' ).Fill( lW.Rapidity() , weight)
                    getattr( self, 'leptonicW_phi_noSel' ).Fill( lW.Phi() , weight)
                    getattr( self, 'leptonicW_mass_noSel' ).Fill( lW.M() , weight)
                    #getattr( self, 'leptonicWMT_noSel').Fill( lW.Mt() , weight)

                for lTop in recoleptTop:
                    getattr( self, 'leptonicTop_pt_noSel' ).Fill( lTop.Pt(), weight )
                    getattr( self, 'leptonicTop_eta_noSel' ).Fill( lTop.Eta() , weight)
                    getattr( self, 'leptonicTop_y_noSel' ).Fill( lTop.Rapidity() , weight)
                    getattr( self, 'leptonicTop_phi_noSel' ).Fill( lTop.Phi() , weight)
                    getattr( self, 'leptonicTop_mass_noSel' ).Fill( lTop.M() , weight)
                    #getattr( self, 'leptonicWMT_noSel').Fill( lW.Mt() , weight)
                


                #### Checking nominal Wtop selection with weights
                if ('_nom' in iSel.keys()) and passSel['_nom'] and iSel['_nom']:
                    
                    # basic reco histos
                    getattr( self, 'nPVs'+iSel['_nom'] ).Fill( getattr( event, 'PV_npvsGood'), weight )
                    getattr( self, 'nleps'+iSel['_nom'] ).Fill( nleptons , weight)
                    
                    for imuon in recoMuons:
                        getattr( self, 'muons_pt'+iSel['_nom'] ).Fill( imuon.pt , weight)
                        getattr( self, 'muons_eta'+iSel['_nom'] ).Fill( imuon.eta , weight)
                        getattr( self, 'muons_y'+iSel['_nom'] ).Fill( imuon.p4().Rapidity() , weight)
                        getattr( self, 'muons_phi'+iSel['_nom'] ).Fill( imuon.phi , weight)

                    for iele in recoElectrons:
                        getattr( self, 'eles_pt'+iSel['_nom'] ).Fill( iele.pt , weight)
                        getattr( self, 'eles_eta'+iSel['_nom'] ).Fill( iele.eta , weight)
                        getattr( self, 'eles_y'+iSel['_nom'] ).Fill( iele.p4().Rapidity() , weight)
                        getattr( self, 'eles_phi'+iSel['_nom'] ).Fill( iele.phi , weight)

                    getattr( self, 'nAK8jets'+iSel['_nom'] ).Fill( len(recoAK8jets['_nom']) , weight)
                    getattr( self, 'HT'+iSel['_nom'] ).Fill( AK8HT , weight)

                    for ijet in recoAK8jets['_nom']:
                        getattr( self, 'AK8jets_pt'+iSel['_nom'] ).Fill( ijet.pt_nom , weight)
                        getattr( self, 'AK8jets_eta'+iSel['_nom'] ).Fill( ijet.eta , weight)
                        getattr( self, 'AK8jets_y'+iSel['_nom'] ).Fill( ijet.rapidity , weight)
                        getattr( self, 'AK8jets_phi'+iSel['_nom'] ).Fill( ijet.phi , weight)
                        getattr( self, 'AK8jets_mass'+iSel['_nom'] ).Fill( ijet.mass_nom, weight)
                        getattr( self, 'AK8jets_msoftdrop'+iSel['_nom'] ).Fill( ijet.msoftdrop, weight)
                    
                    getattr( self, 'nAK4jets'+iSel['_nom'] ).Fill( len(recoAK4jets), weight )
                    getattr( self, 'nAK4bjets'+iSel['_nom'] ).Fill( len(recoAK4bjets) , weight)
                    getattr( self, 'nAK4leptjets'+iSel['_nom'] ).Fill( len(recoAK4LeptJets['_nom']), weight )
                    getattr( self, 'nAK4bleptjets'+iSel['_nom'] ).Fill( len(recoAK4bLeptJets['_nom']), weight )
                    getattr( self, 'dR_Mu_AK4blept'+iSel['_nom']  ).Fill( recodR_Mu_leptAK4b['_nom'] , weight)
                    getattr( self, 'dPhi_AK8_AK4blept'+iSel['_nom'] ).Fill( recodPhi_AK8_leptAK4b['_nom'], weight )

                    for ijet in recoAK4jets:
                        getattr( self, 'AK4jets_pt'+iSel['_nom'] ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4jets_eta'+iSel['_nom'] ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4jets_y'+iSel['_nom'] ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4jets_phi'+iSel['_nom'] ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4jets_mass'+iSel['_nom'] ).Fill( ijet.mass , weight)
                        
                    for ijet in recoAK4bjets:
                        getattr( self, 'AK4btaggedjets_pt'+iSel['_nom'] ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4btaggedjets_eta'+iSel['_nom'] ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4btaggedjets_y'+iSel['_nom'] ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4btaggedjets_phi'+iSel['_nom'] ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4btaggedjets_mass'+iSel['_nom'] ).Fill( ijet.mass , weight)

                    for ijet in recoAK4LeptJets['_nom']:
                        getattr( self, 'AK4leptjets_pt'+iSel['_nom'] ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4leptjets_eta'+iSel['_nom'] ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4leptjets_y'+iSel['_nom'] ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4leptjets_phi'+iSel['_nom'] ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4leptjets_mass'+iSel['_nom'] ).Fill( ijet.mass, weight )

                    if len( recoAK4bLeptJets['_nom'])>=1:
                        getattr( self, 'AK4leptbjet_pt'+iSel['_nom'] ).Fill( recoAK4bLeptJets['_nom'][0].pt , weight)
                        getattr( self, 'AK4leptbjet_eta'+iSel['_nom'] ).Fill( recoAK4bLeptJets['_nom'][0].eta , weight)
                        getattr( self, 'AK4leptbjet_y'+iSel['_nom'] ).Fill( recoAK4bLeptJets['_nom'][0].p4().Rapidity(), weight)
                        getattr( self, 'AK4leptbjet_phi'+iSel['_nom'] ).Fill( recoAK4bLeptJets['_nom'][0].phi , weight)
                        getattr( self, 'AK4leptbjet_mass'+iSel['_nom'] ).Fill( recoAK4bLeptJets['_nom'][0].mass , weight)

                    getattr( self, 'METPt'+iSel['_nom'] ).Fill( MET.Pt() , weight)
                    
                    for lW in recoleptW:
                        getattr( self, 'leptonicW_pt'+iSel['_nom'] ).Fill( lW.Pt(), weight )
                        getattr( self, 'leptonicW_eta'+iSel['_nom'] ).Fill( lW.Eta() , weight)
                        getattr( self, 'leptonicW_y'+iSel['_nom'] ).Fill( lW.Rapidity() , weight)
                        getattr( self, 'leptonicW_phi'+iSel['_nom'] ).Fill( lW.Phi() , weight)
                        getattr( self, 'leptonicW_mass'+iSel['_nom'] ).Fill( lW.M() , weight)
                        #getattr( self, 'leptonicWMT'+iSel['_nom']).Fill( lW.Mt() , weight)

                    for lTop in recoleptTop:
                        getattr( self, 'leptonicTop_pt'+iSel['_nom'] ).Fill( lTop.Pt(), weight )
                        getattr( self, 'leptonicTop_eta'+iSel['_nom'] ).Fill( lTop.Eta() , weight)
                        getattr( self, 'leptonicTop_y'+iSel['_nom'] ).Fill( lTop.Rapidity() , weight)
                        getattr( self, 'leptonicTop_phi'+iSel['_nom'] ).Fill( lTop.Phi() , weight)
                        getattr( self, 'leptonicTop_mass'+iSel['_nom'] ).Fill( lTop.M() , weight)
                        #getattr( self, 'leptonicWMT'+iSel['_nom']).Fill( lW.Mt() , weight)

            

        return passSel, iSel, recodR_Mu_leptAK4b, recodPhi_Mu_leptAK4b, recodRap_Mu_leptAK4b, recodPhi_AK8_leptAK4b, recodR_AK8_leptAK4b, recoMuons_dict, recoAK4bjets, recoAK4bLeptJets, recoAK4LeptJets, recoleptW_dict, recoleptW_uncorr_dict, recoleptTop_dict, recoleptTop_uncorr_dict, recoAK8jets, MET_dict, MET_uncorr_dict, nAK4_dict#, nCrackEl

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing gen information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genJets = list(Collection( event, 'GenJet'))
        genParticles = Collection(event, 'GenPart')
        genmet = Object( event, 'GenMET')
        
        tops =  [x for x in genParticles if x.pdgId==6 and x.statusFlags==14]
        tops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)
        antiTops =  [x for x in genParticles if x.pdgId==-6 and x.statusFlags==14]
        antiTops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)

        #### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and abs(x.eta)<self.maxElectronEta ]# ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        #allGenElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and (abs(x.eta)<2.5)]# or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #and abs(x.eta)<self.maxMuonEta ]
        #allGenElectrons.sort(key=lambda x:x.pt, reverse=True)
        #allCrackElectrons = [x for x in allGenElectrons if not(x in genElectrons)]
        #nCrackEl = len(allCrackElectrons)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        
        genVetoMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minLooseMuonPt and abs(x.eta) < self.maxMuonEta and not(x in genMuons)]
        genVetoMuons.sort(key=lambda x:x.pt, reverse=True)
        
        #leptonVetoFlag=False
        ngenLeptons = len(genMuons)+len(genElectrons)+len(genVetoMuons)
        leptonVetoFlag = ngenLeptons>=2 
        #if leptonVetoFlag: print(event.event,'gen',leptonVetoFlag,ngenLeptons,len(genMuons),len(genElectrons),len(genLooseMuons))
        ##################################################

        #### MET (not sure if needed)
        genMET = ROOT.TLorentzVector()
        genMET.SetPtEtaPhiM(genmet.pt, 0., genmet.phi, 0)
        ##################################################

        #### Basic AK4 jet selection
        genAK4jets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta]# and abs(x.hadronFlavour)==5 ] 
        genAK4jets.sort(key=lambda x:x.pt,reverse=True)
        nAK4 = len(genAK4jets)
        genAK4bjets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta and abs(x.hadronFlavour)==5 ]  # hadron flavour ghost matched (to B's) jets as a proxy for b-jets 
        genAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        #### Basic AK8 jet selection 
        
        if len(genJetsAK8)!=0: 
            for ijets in genJetsAK8: 
                ijets.rapidity = ijets.p4().Rapidity()
                ijets.msoftdrop = self.getGenJetAK8softdropmass( AK8jet=ijets, event=event, PFCollection='GenCands', isGen=True)
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.rapidity) < self.maxJetAK8Rap]
        genAK8HT = sum( [ x.pt for x in genAK8jets ] )
        ##################################################    
        
        
        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        if not(leptonVetoFlag):
            if not self.effMapOnly:
                passSel, iSel, gendR_Mu_leptAK4b, gendPhi_Mu_leptAK4b, gendRap_Mu_leptAK4b, gendPhi_AK8_leptAK4b, gendR_AK8_leptAK4b, genAK4bleptjets, genAK4leptjets = self.WtopSelection( True, event, 
                                                                                                                                                                        genMuons, genElectrons, 
                                                                                                                                                                        genAK4jets, genAK4bjets, 
                                                                                                                                                                        genAK8jets, genMET, '' )
            else:
                passSel, iSel, gendR_Mu_leptAK4b, gendPhi_Mu_leptAK4b, gendRap_Mu_leptAK4b, gendPhi_AK8_leptAK4b, gendR_AK8_leptAK4b, genAK4bleptjets, genAK4leptjets = self.WtopSelection_effMap( True, event, 
                                                                                                                                                                           genMuons, genElectrons, 
                                                                                                                                                                           genAK4jets, genAK4bjets, 
                                                                                                                                                                           genAK8jets, genMET, '' )
        else:

            passSel, iSel, gendR_Mu_leptAK4b, gendPhi_Mu_leptAK4b, gendRap_Mu_leptAK4b, gendPhi_AK8_leptAK4b, gendR_AK8_leptAK4b, genAK4bleptjets, genAK4leptjets  = False, None, 929, 929, 929, 929, 929, [], []
        

        #### reconstruct leptonic W and top if possible
        if len(genMuons)>0:
            genleptW = [genMuons[0].p4()+genMET]
            genleptTop = []
            if len(genAK4bleptjets)>0:
                genleptTop = [genMuons[0].p4()+genMET+genAK4bleptjets[0].p4()]
        else:
            genleptW = [] #genMuons[0].p4()+genMET]
            genleptTop = []

        #### Weight
        weight = event.genWeight
        #self.totalGenWeight = weight

        if self.controlHistos:

            ##### Filling histograms
            if not self.onlyUnc: 
                #### Checking no selection
                getattr( self, 'ngenleps_noSel' ).Fill( len(genMuons)+len(genElectrons), weight )
                
                for imuon in genMuons:
                    getattr( self, 'genmuons_pt_noSel' ).Fill( imuon.pt, weight )
                    getattr( self, 'genmuons_eta_noSel' ).Fill( imuon.eta, weight )
                    getattr( self, 'genmuons_y_noSel' ).Fill( imuon.p4().Rapidity(), weight )
                    getattr( self, 'genmuons_phi_noSel' ).Fill( imuon.phi, weight )
                
                for iele in genElectrons:
                    getattr( self, 'geneles_pt_noSel' ).Fill( iele.pt, weight )
                    getattr( self, 'geneles_eta_noSel' ).Fill( iele.eta, weight )
                    getattr( self, 'geneles_y_noSel' ).Fill( iele.p4().Rapidity(), weight )
                    getattr( self, 'geneles_phi_noSel' ).Fill( iele.phi, weight )
                
                getattr( self, 'ngenAK8jets_noSel' ).Fill( len(genAK8jets), weight )
                getattr( self, 'genHT_noSel' ).Fill( genAK8HT, weight )
                
                for ijet in genAK8jets:
                    getattr( self, 'AK8genjets_pt_noSel' ).Fill( ijet.pt, weight )
                    getattr( self, 'AK8genjets_eta_noSel' ).Fill( ijet.eta, weight )
                    getattr( self, 'AK8genjets_y_noSel' ).Fill( ijet.rapidity, weight )
                    getattr( self, 'AK8genjets_phi_noSel' ).Fill( ijet.phi, weight )
                    getattr( self, 'AK8genjets_mass_noSel' ).Fill( ijet.mass, weight )
                    getattr( self, 'AK8genjets_msoftdrop_noSel' ).Fill( ijet.msoftdrop, weight )
                
                '''
                AK4genjets
                AK4bmatchedgenjets
                AK4genleptjets
                AK4genleptbjet
                '''

                getattr( self, 'ngenAK4jets_noSel' ).Fill( len(genAK4jets), weight )
                getattr( self, 'ngenAK4bjets_noSel' ).Fill( len(genAK4bjets), weight )
                getattr( self, 'ngenAK4leptjets_noSel' ).Fill( len(genAK4leptjets), weight )
                getattr( self, 'ngenAK4bleptjets_noSel' ).Fill( len(genAK4bleptjets), weight )
                getattr( self, 'dR_genMu_genAK4blept_noSel' ).Fill( gendR_Mu_leptAK4b, weight )
                getattr( self, 'dPhi_genAK8_genAK4blept_noSel' ).Fill( gendPhi_AK8_leptAK4b, weight )
                
                for ijet in genAK4jets:
                    getattr( self, 'AK4genjets_pt_noSel' ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4genjets_eta_noSel' ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4genjets_y_noSel' ).Fill( ijet.p4().Rapidity(), weight )
                    getattr( self, 'AK4genjets_phi_noSel' ).Fill( ijet.phi, weight )
                    getattr( self, 'AK4genjets_mass_noSel' ).Fill( ijet.mass, weight )

                for ijet in genAK4bjets:
                    getattr( self, 'AK4bmatchedgenjets_pt_noSel' ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4bmatchedgenjets_eta_noSel' ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4bmatchedgenjets_y_noSel' ).Fill( ijet.p4().Rapidity(), weight )
                    getattr( self, 'AK4bmatchedgenjets_phi_noSel' ).Fill( ijet.phi, weight )
                    getattr( self, 'AK4bmatchedgenjets_mass_noSel' ).Fill( ijet.mass, weight )

                for ijet in genAK4leptjets:
                    getattr( self, 'AK4genleptjets_pt_noSel' ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4genleptjets_eta_noSel' ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4genleptjets_y_noSel' ).Fill( ijet.p4().Rapidity(), weight )
                    getattr( self, 'AK4genleptjets_phi_noSel' ).Fill( ijet.phi, weight )
                    getattr( self, 'AK4genleptjets_mass_noSel' ).Fill( ijet.mass, weight )

                if len( genAK4bleptjets)>=1:
                    getattr( self, 'AK4genleptbjet_pt_noSel' ).Fill( genAK4bleptjets[0].pt, weight )
                    getattr( self, 'AK4genleptbjet_eta_noSel' ).Fill( genAK4bleptjets[0].eta, weight )
                    getattr( self, 'AK4genleptbjet_y_noSel' ).Fill( genAK4bleptjets[0].p4().Rapidity(), weight )
                    getattr( self, 'AK4genleptbjet_phi_noSel' ).Fill( genAK4bleptjets[0].phi, weight )
                    getattr( self, 'AK4genleptbjet_mass_noSel' ).Fill( genAK4bleptjets[0].mass, weight )
                
                getattr( self, 'genMETPt_noSel' ).Fill( genMET.Pt(), weight )

                for lW in genleptW:
                    getattr( self, 'genleptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                    getattr( self, 'genleptonicW_eta_noSel' ).Fill( lW.Eta(), weight )
                    getattr( self, 'genleptonicW_y_noSel' ).Fill( lW.Rapidity(), weight )
                    getattr( self, 'genleptonicW_phi_noSel' ).Fill( lW.Phi(), weight )
                    getattr( self, 'genleptonicW_mass_noSel' ).Fill( lW.M(), weight )
                    #getattr( self, 'genleptonicWMT_noSel' ).Fill( lW.Mt(), weight )
                
                for lTop in genleptTop:
                    getattr( self, 'genleptonicTop_pt_noSel' ).Fill( lTop.Pt(), weight )
                    getattr( self, 'genleptonicTop_eta_noSel' ).Fill( lTop.Eta(), weight )
                    getattr( self, 'genleptonicTop_y_noSel' ).Fill( lTop.Rapidity(), weight )
                    getattr( self, 'genleptonicTop_phi_noSel' ).Fill( lTop.Phi(), weight )
                    getattr( self, 'genleptonicTop_mass_noSel' ).Fill( lTop.M(), weight )
                    #getattr( self, 'genleptonicTopMT_noSel' ).Fill( lW.Mt(), weight )
                
                
                if passSel and iSel:
                    #print ("Filling gen histos")
                    #### Checking nominal Wtop selection with weights
                    getattr( self, 'ngenleps'+iSel ).Fill( len(genMuons)+len(genElectrons), weight )
                    for imuon in genMuons:
                        getattr( self, 'genmuons_pt'+iSel ).Fill( imuon.pt, weight )
                        getattr( self, 'genmuons_eta'+iSel ).Fill( imuon.eta, weight )
                        getattr( self, 'genmuons_y'+iSel ).Fill( imuon.p4().Rapidity(), weight )
                        getattr( self, 'genmuons_phi'+iSel ).Fill( imuon.phi, weight )

                    for iele in genElectrons:
                        getattr( self, 'geneles_pt'+iSel ).Fill( iele.pt, weight )
                        getattr( self, 'geneles_eta'+iSel ).Fill( iele.eta, weight )
                        getattr( self, 'geneles_y'+iSel ).Fill( iele.p4().Rapidity(), weight )
                        getattr( self, 'geneles_phi'+iSel ).Fill( iele.phi, weight )

                    getattr( self, 'ngenAK8jets'+iSel ).Fill( len(genAK8jets), weight )
                    getattr( self, 'genHT'+iSel ).Fill( genAK8HT, weight )

                    for ijet in genAK8jets:
                        getattr( self, 'AK8genjets_pt'+iSel ).Fill( ijet.pt, weight )
                        getattr( self, 'AK8genjets_eta'+iSel ).Fill( ijet.eta, weight )
                        getattr( self, 'AK8genjets_y'+iSel ).Fill( ijet.rapidity, weight )
                        getattr( self, 'AK8genjets_phi'+iSel ).Fill( ijet.phi, weight )
                        getattr( self, 'AK8genjets_mass'+iSel ).Fill( ijet.mass, weight )
                        getattr( self, 'AK8genjets_msoftdrop'+iSel ).Fill( ijet.msoftdrop, weight )

                    getattr( self, 'ngenAK4jets'+iSel ).Fill( len(genAK4jets), weight )
                    getattr( self, 'ngenAK4bjets'+iSel ).Fill( len(genAK4bjets), weight )
                    getattr( self, 'ngenAK4leptjets'+iSel ).Fill( len(genAK4leptjets), weight )
                    getattr( self, 'ngenAK4bleptjets'+iSel ).Fill( len(genAK4bleptjets), weight )
                    getattr( self, 'dR_genMu_genAK4blept'+iSel).Fill( gendR_Mu_leptAK4b, weight )
                    getattr( self, 'dPhi_genAK8_genAK4blept'+iSel ).Fill( gendPhi_AK8_leptAK4b, weight )

                    for ijet in genAK4jets:
                        getattr( self, 'AK4genjets_pt'+iSel  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4genjets_eta'+iSel  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4genjets_y'+iSel  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4genjets_phi'+iSel  ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4genjets_mass'+iSel  ).Fill( ijet.mass, weight )

                    for ijet in genAK4bjets:
                        getattr( self, 'AK4bmatchedgenjets_pt'+iSel  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4bmatchedgenjets_eta'+iSel  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4bmatchedgenjets_y'+iSel  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4bmatchedgenjets_phi'+iSel ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4bmatchedgenjets_mass'+iSel ).Fill( ijet.mass, weight )

                    for ijet in genAK4leptjets:
                        getattr( self, 'AK4genleptjets_pt'+iSel ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4genleptjets_eta'+iSel ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4genleptjets_y'+iSel ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4genleptjets_phi'+iSel ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4genleptjets_mass'+iSel ).Fill( ijet.mass, weight )

                    if len(genAK4bleptjets)>=1:
                        getattr( self, 'AK4genleptbjet_pt'+iSel ).Fill( genAK4bleptjets[0].pt, weight )
                        getattr( self, 'AK4genleptbjet_eta'+iSel ).Fill( genAK4bleptjets[0].eta, weight )
                        getattr( self, 'AK4genleptbjet_y'+iSel ).Fill( genAK4bleptjets[0].p4().Rapidity(), weight )
                        getattr( self, 'AK4genleptbjet_phi'+iSel ).Fill( genAK4bleptjets[0].phi, weight )
                        getattr( self, 'AK4genleptbjet_mass'+iSel ).Fill( genAK4bleptjets[0].mass, weight )

                    getattr( self, 'genMETPt'+iSel ).Fill( genMET.Pt(), weight )
                    
                    for lW in genleptW:
                        getattr( self, 'genleptonicW_pt'+iSel  ).Fill( lW.Pt(), weight )
                        getattr( self, 'genleptonicW_eta'+iSel  ).Fill( lW.Eta(), weight )
                        getattr( self, 'genleptonicW_y'+iSel  ).Fill( lW.Rapidity(), weight )
                        getattr( self, 'genleptonicW_phi'+iSel  ).Fill( lW.Phi(), weight )
                        getattr( self, 'genleptonicW_mass'+iSel  ).Fill( lW.M(), weight )
                        #getattr( self, 'genleptonicWMT'+iSel ).Fill( lW.Mt(), weight )

                    for lTop in genleptTop:
                        getattr( self, 'genleptonicTop_pt'+iSel  ).Fill( lTop.Pt(), weight )
                        getattr( self, 'genleptonicTop_eta'+iSel  ).Fill( lTop.Eta(), weight )
                        getattr( self, 'genleptonicTop_y'+iSel  ).Fill( lTop.Rapidity(), weight )
                        getattr( self, 'genleptonicTop_phi'+iSel  ).Fill( lTop.Phi(), weight )
                        getattr( self, 'genleptonicTop_mass'+iSel  ).Fill( lTop.M(), weight )
                        #getattr( self, 'genleptonicTopMT_noSel' ).Fill( lW.Mt(), weight )
            

        return passSel, iSel, gendR_Mu_leptAK4b, gendPhi_Mu_leptAK4b, gendRap_Mu_leptAK4b, gendPhi_AK8_leptAK4b, gendR_AK8_leptAK4b, genMuons, genAK4bjets, genAK4bleptjets, genAK4leptjets, genleptW, genleptTop, genAK8jets, genMET, nAK4#, nCrackEl

    #############################################################################
    def WtopSelection( self, isGen, event, muons, electrons, AK4jets, AK4bjets, AK8jets, MET, ptLabel):
    
        if (len(muons)==1) and (len(electrons)==0) and (len(AK8jets)>=1) and MET.Pt()>self.METCutWtop and (len(AK4jets)>=1):# and (len(AK4bjets)>=1)              
            
            # reconstruct leptonic hemisphere objects for selection/'tag' (a la TnP) here 
            leptW=muons[0].p4()+MET
            
            AK4lepjets = [x for x in AK4jets if self.DrRapPhi(muons[0].p4(),x.p4())<1.6]# and (self.DrRapPhi(x.p4(), muons[0].p4() )>0.4) ]
            AK4lepjets.sort(key=lambda x: getattr(x,'pt'+ptLabel),reverse=True) 

            #leptTop=leptW+AK4leptbjets[0].p4() if len(AK4leptbjets)>=1 else None
            #and (leptTop and leptTop.M()<AK8jets[0].mass)
            if not(len(AK4lepjets)>=1): 
                return False, None, 929, 929, 929, 929, 929, [], []

            isLeadJetBCandFlag = False
            if not isGen: 
                #AK4leptbjets = [x for x in AK4lepjets if x.btagDeepFlavB > self.minBDisc] #stupid [0:1] to ensure it's obvious that the b-tag is required on the leading  AK4
                #AK4leptbjets.sort(key=lambda x:x.pt,reverse=True) 
                if AK4lepjets[0].btagDeepFlavB>self.minBDisc: isLeadJetBCandFlag=True
            else: 
                #AK4leptbjets = [x for x in AK4lepjets if abs(x.hadronFlavour)==5] #stupid [0:1] to ensure it's obvious that the b-tag is required on the leading  AK4
                #AK4leptbjets.sort(key=lambda x:x.pt,reverse=True) 
                if abs(AK4lepjets[0].hadronFlavour)==5: isLeadJetBCandFlag=True

            isBCandHadronicallyIsolated = not any( self.DrRapPhi(AK4lepjets[0].p4(), x.p4()) < 0.4 and getattr(AK4lepjets[0],'pt'+ptLabel) != getattr(x,'pt'+ptLabel) for x in AK4jets )

            if isLeadJetBCandFlag and isBCandHadronicallyIsolated and leptW.Pt()>self.minLeptonicWPt and AK8jets[0].p4().DeltaPhi(muons[0].p4())>2. and self.DrRapPhi(AK8jets[0].p4(),AK4lepjets[0].p4())>1.2: #2~=2pi/3 #len(AK4leptbjets)==1

                # to reconstruct the hadronic hemisphere: ensure separation of lead AK8 jet in the event 
                # from the leptonic hemisphere mu and the leading leptonic hemisphere AK4, which also must be b-tagged,  
                # additionally, ensure hadronic hemisphere's topology is boosted, with the leptonic W pT cut 
                # (ie, the leptonic top's decay products )
                # This is effectively also requiring that the measurement AK8 is the leading AK8 in the event,
                # and is reconstructed in the in had. hem., given any other AK8 (say from the leptW+lept.b)
                # would be expected to have less energy/pT ()

                leadJetpT = getattr( AK8jets[0], 'pt'+ptLabel )
                leadJetmass = getattr( AK8jets[0], 'mass'+ptLabel ) #msoftdrop
                AK4leptbjets = AK4lepjets[0:1]


                # found leptonic hemisphere AK4s within dR<1.6 of muon, check if the leading such AK4 (ak4leptjets[0]) is b-tagged (or hadr. flavour matched) else discard event (so if >1 b-tags in lept hem event still accepted)
                # event kept if there is exactly one AK4 in the leptonic hemisphere, close to the mu, and it is b-tagged
                
                dR_mu_leptAK4b = self.DrRapPhi(muons[0].p4(),AK4leptbjets[0].p4()) 

                dPhi_mu_leptAK4b = AK4leptbjets[0].p4().DeltaPhi(muons[0].p4()) 
                dPhi_AK8_leptAK4b = AK4leptbjets[0].p4().DeltaPhi(AK8jets[0].p4()) 
                dR_AK8_leptAK4b = self.DrRapPhi(AK8jets[0].p4(),AK4leptbjets[0].p4()) 

                dRap_mu_leptAK4b = AK4leptbjets[0].p4().Rapidity() - muons[0].p4().Rapidity()

                if self.evtSelection.startswith('_Wtop'):
                    if leadJetmass>self.minLeadJetMass and leadJetpT>self.minLeadAK8JetPtW:# and (len(AK4hadbjets)==1):
                        # apply a lowish inv. mass and pT cut to get rid of extra events, 
                        # selectively keeping some info on the lept. had hem b-candidate to return, to do control plots with them
                        
                        return True, '_WtopSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets #store dPhi between AK4 and AK8 to check if they're far enough apart in the azimuth
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets

                elif 'WSel' in self.evtSelection: 
                    if (((leadJetmass<self.maxSDMassW) and (leadJetmass>self.minSDMassW)) and (leadJetpT>self.minLeadAK8JetPtW)):
                        
                        return True, '_WSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets

                elif '_topSel' in self.evtSelection:
                    if ((leadJetmass > self.minSDMassTop) and (leadJetpT> self.minLeadAK8JetPtTop)):# and (leadJetmass/(1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop)
                        return True, '_topSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets   
            else: 
                #if len(AK4leptbjets)==0:
                return False, None, 929, 929, 929, 929, 929, [], AK4lepjets if len(AK4lepjets)!=0 else []
                #else:
                #    return False, None, 929, 929, 929, 929, AK4leptbjets, AK4lepjets if len(AK4lepjets)!=0 else []
        else: return False, None, 929, 929, 929, 929, 929, [], []

    #############################################################################
    def WtopSelection_effMap( self, isGen, event, muons, electrons, AK4jets, AK4bjets, AK8jets, MET, ptLabel):
        #removing b-tag requirements, and saving in the btagged collections only the leading leptonic hem. AK4 in dR<1.6 of the muon

        if (len(muons)==1) and (len(electrons)==0) and (len(AK8jets)>=1) and MET.Pt()>self.METCutWtop and (len(AK4jets)>=1):# and (len(AK4bjets)>=1)              
            
            # reconstruct leptonic hemisphere objects for selection
            leptW=muons[0].p4()+MET
            
            AK4lepjets = [x for x in AK4jets if self.DrRapPhi(muons[0].p4(),x.p4())<1.6 ]#and (self.DrRapPhi(x.p4(), muons[0].p4() )>0.4) ]
            AK4lepjets.sort(key=lambda x:x.pt,reverse=True) 

            #leptTop=leptW+AK4leptbjets[0].p4() if len(AK4leptbjets)>=1 else None
            #and (leptTop and leptTop.M()<AK8jets[0].mass)
            if not(len(AK4lepjets)>=1): 
                return False, None, 929, 929, 929, 929, 929, [], []

            isLeadJetBCandFlag = False 
            if not isGen: 

                #Set to true by default to allow events where the lead AK4 in the lept hem is isolated from other AK4s and is accepted as the b-cand w/ and w/o a b-tag
                if AK4lepjets[0]: isLeadJetBCandFlag=True #.btagDeepFlavB>self.minBDisc
            else: 
                if abs(AK4lepjets[0].hadronFlavour)==5: isLeadJetBCandFlag=True

            isBCandHadronicallyIsolated = not any( self.DrRapPhi(AK4lepjets[0].p4(), x.p4()) < 0.4 and AK4lepjets[0].pt != x.pt  for x in AK4jets )

            if isLeadJetBCandFlag and isBCandHadronicallyIsolated and leptW.Pt()>self.minLeptonicWPt and AK8jets[0].p4().DeltaPhi(muons[0].p4())>2. and self.DrRapPhi(AK8jets[0].p4(),AK4lepjets[0].p4())>1.2: #2~=2pi/3 #len(AK4leptbjets)==1

               
                leadJetpT = getattr( AK8jets[0], 'pt'+ptLabel )
                leadJetmass = getattr( AK8jets[0], 'mass'+ptLabel )#msoftdrop 
                AK4leptbjets = AK4lepjets[0:1]


                # found leptonic hemisphere AK4s within dR<1.6 of muon, check if the leading such AK4 (ak4leptjets[0]) is b-tagged (or hadr. flavour matched) else discard event (so if >1 b-tags in lept hem event still accepted)
                # event kept if there is exactly one AK4 in the leptonic hemisphere, close to the mu, and it is b-tagged
                
                dR_mu_leptAK4b = self.DrRapPhi(muons[0].p4(),AK4leptbjets[0].p4()) 

                dPhi_mu_leptAK4b = AK4leptbjets[0].p4().DeltaPhi(muons[0].p4()) 
                dPhi_AK8_leptAK4b = AK4leptbjets[0].p4().DeltaPhi(AK8jets[0].p4()) 
                dR_AK8_leptAK4b = self.DrRapPhi(AK8jets[0].p4(),AK4leptbjets[0].p4()) 

                dRap_mu_leptAK4b = AK4leptbjets[0].p4().Rapidity() - muons[0].p4().Rapidity()

                if self.evtSelection.startswith('_Wtop'):
                    if leadJetmass>self.minLeadJetMass and leadJetpT>self.minLeadAK8JetPtW:# and (len(AK4hadbjets)==1):
                        # apply a lowish inv. mass and pT cut to get rid of some events
                        
                        return True, '_WtopSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets

                elif 'WSel' in self.evtSelection: 
                    if (((leadJetmass<self.maxSDMassW) and (leadJetmass>self.minSDMassW)) and (leadJetpT>self.minLeadAK8JetPtW)):
                        
                        return True, '_WSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets

                elif '_topSel' in self.evtSelection:
                    if ((leadJetmass > self.minSDMassTop) and (leadJetpT> self.minLeadAK8JetPtTop)):# and (leadJetmass/(1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop)
                        return True, '_topSel', dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets
                    else: 
                        return False, None, dR_mu_leptAK4b, dPhi_mu_leptAK4b, dRap_mu_leptAK4b, dPhi_AK8_leptAK4b, dR_AK8_leptAK4b, AK4leptbjets, AK4lepjets   
            else: 
                #if len(AK4leptbjets)==0:
                return False, None, 929, 929, 929, 929, 929, [], AK4lepjets if len(AK4lepjets)!=0 else []
                #else:
                #    return False, None, 929, 929, 929, 929, AK4leptbjets, AK4lepjets if len(AK4lepjets)!=0 else []
        else: return False, None, 929, 929, 929, 929, 929, [], []
        
        

       

    #############################################################################
    def getGenJetAK8softdropmass(self, AK8jet, event, PFCollection, isGen=True ): 
        '''for gen, in this use case, but generically useful for reco if needed; uses same inputs as createNsubBasis function below'''
        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        #### Applying PUPPI weights to the PF candidates for reco and not for gen and pushing back constituents
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()

        for p in pfCands :
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            tp = tp * p.puppiWeight if not isGen else tp
            #except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        #print ("pushing back candidates")
        for x in CandsPUPPIweightedVec:
            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
        #print ("pushed back candidates")

        if isGen: #to add in softdrop mass as a variable to the selected (accep)gen jet branches
            sd_AK8jets = self.sd.result( constituents)
            if len(sd_AK8jets)>0: #stupidly, in some rare cases, this will not be true (with ptmin>=170) leading to errors, hence, switching to ptmin=0 in function calls (but, since not sure this is error, free I use this if-else block)
                ak8jet['msoftdrop'] = sd_AK8jets[0].m()
            else: ak8jet['msoftdrop'] = -1. 
        return ak8jet['msoftdrop']

    #############################################################################
    

    def createNsubBasis(self, AK8jet, event, PFCollection, isGen=False, constJES=False, varUpDown='',noPlot=True ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use
        type = 'reco' if not(isGen) else 'gen'

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        if not(noPlot) and not(isGen) and not(self.onlyUnc.startswith('_je')):
            subjets = list(Collection(event, 'SubJet' )) if not(isGen) else []#list(Collection(event, 'SubGenJetAK8' )) 
            subjet1 = subjets[AK8jet.subJetIdx1] if not(isGen) else []
            subjet2 = subjets[AK8jet.subJetIdx2] if not(isGen) else [] 

            dRSJ1_AK8 = self.DrRapPhi(AK8jet.p4(),subjet1.p4()) if not (isGen) else None
            dRSJ2_AK8 = self.DrRapPhi(AK8jet.p4(),subjet2.p4()) if not (isGen) else None

            getattr(self, '%s_dR_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
            getattr(self, '%s_dR_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)        
            if dRSJ1_AK8<0.6: getattr(self, '%s_dR2_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
            if dRSJ2_AK8<0.6: getattr(self, '%s_dR2_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)

        #### Run calculations of NSub bases and store for ungroomed AK8jets (default in CMS)

        #### Applying PUPPI weights to the PF candidates for reco and not for gen and pushing back constituents
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()
        #if constJES:
        if not(noPlot):

            charges = []
            pIDs = []
            dRs = []
            dRaps = []
            dPhis = []
            pTs = []
            for p in pfCands :
                #if p.p4().M()<0.: #to check on -ve mass electrons in constituents
                #    p.p4().M()=0.
                tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                tp = tp * p.puppiWeight if not isGen else tp
                #except RuntimeError: tp = tp    ### for genjets
                CandsPUPPIweightedVec.push_back(tp)
                #if constJES:
                charges.append(p.charge)
                pIDs.append(p.pdgId)
                dRs.append(self.DrRapPhi( AK8jet.p4(), p.p4() ))
                dRaps.append(AK8jet.p4().Rapidity() - p.p4().Rapidity() )
                dPhis.append(AK8jet.p4().DeltaPhi(p.p4()) )
                pTs.append(p.p4().Pt())

            #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
            #print ("pushing back candidates")
            if not(constJES):
                #print(AK8jet,isGen,constJES,central,noPlot)

                nch = 0
                npho = 0
                nchHad = 0
                nneutHad = 0

                #for ix,x in enumerate(CandsPUPPIweightedVec):
                for c,pID,dR,dRap,dPhi,pT, x in zip(charges, pIDs, dRs, dRaps, dPhis, pTs, CandsPUPPIweightedVec):

                    puppi_dR = self.DrRapPhi( AK8jet.p4(), x )
                    if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                    #print ("pushed back candidates")
                    if not (self.onlyUnc.startswith('_je')):
                        if isGen and dR>0.8: continue

                        if not(isGen) and puppi_dR>0.8: continue

                        if pT > 15:
                            bin_name = 'gt15'
                        elif pT > 10:
                            bin_name = '10_15'
                        elif pT > 5:
                            bin_name = '5_10'
                        elif pT > 4:
                            bin_name = '4_5'
                        elif pT > 3:
                            bin_name = '3_4'
                        elif pT > 2:
                            bin_name = '2_3'
                        elif pT > 1:
                            bin_name = '1_2'
                        else:
                            bin_name = '0_1'


                        if pID==22:
                            if not(noPlot): 
                                getattr(self, '%s_dR_photons'%(type)).Fill(dR if isGen else puppi_dR)
                                getattr(self, '%s_dRap_dPhi_photons'%(type)).Fill(dRap, dPhi)

                                getattr(self, '%s_dR_photons_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                getattr(self, '%s_dRap_dPhi_photons_%s' % (type, bin_name)).Fill(dRap, dPhi)


                                if not(isGen):
                                    getattr(self, '%s_dRSJ1_photons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                    if dRSJ1_AK8<0.6: getattr(self, '%s_dR2SJ1_photons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                    getattr(self, '%s_dRSJ2_photons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))
                                    if dRSJ2_AK8<0.6: getattr(self, '%s_dR2SJ2_photons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))

                                getattr(self, '%s_pT_photons'%(type)).Fill(pT)
                                getattr(self, '%sjetPt_dR_photons'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())
                                getattr(self, '%sjetmSD_dR_photons'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) 
                                getattr(self, '%s_dR_Pt_photons'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())
                            npho+=1

                        elif c==0 and pID!=22:
                            if not(noPlot): 
                                getattr(self, '%s_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR)
                                getattr(self, '%s_dRap_dPhi_neutral'%(type)).Fill(dRap, dPhi)
                                getattr(self, '%s_dR_neutral_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                getattr(self, '%s_dRap_dPhi_neutral_%s' % (type, bin_name)).Fill(dRap, dPhi)

                                if not(isGen):
                                    getattr(self, '%s_dRSJ1_neutral'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                    if dRSJ1_AK8<0.6: getattr(self, '%s_dR2SJ1_neutral'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                    getattr(self, '%s_dRSJ2_neutral'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))
                                    if dRSJ2_AK8<0.6: getattr(self, '%s_dR2SJ2_neutral'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))

                                getattr(self, '%s_pT_neutral'%(type)).Fill(pT)
                                getattr(self, '%sjetPt_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())
                                getattr(self, '%sjetmSD_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) 
                                getattr(self, '%s_dR_Pt_neutral'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())
                            nneutHad+=1

                        elif abs(c)>0:# and abs(pID)!=11:# and abs(pID)!=13:
                            if not(noPlot): 
                                getattr(self, '%s_dR_charged'%(type)).Fill(dR if isGen else puppi_dR)
                                getattr(self, '%s_dRap_dPhi_charged'%(type)).Fill(dRap, dPhi)
                                #getattr(self, '%s_dR_charged_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                #getattr(self, '%s_dRap_dPhi_charged_%s' % (type, bin_name)).Fill(dRap, dPhi)

                                if not(isGen):
                                    getattr(self, '%s_dRSJ1_charged'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                    if dRSJ1_AK8<0.6: getattr(self, '%s_dR2SJ1_charged'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                    getattr(self, '%s_dRSJ2_charged'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies
                                    if dRSJ2_AK8<0.6: getattr(self, '%s_dR2SJ2_charged'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies

                                getattr(self, '%s_pT_charged'%(type)).Fill(pT)
                                getattr(self, '%sjetPt_dR_charged'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())#x is dR, y is pT
                                getattr(self, '%sjetmSD_dR_charged'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) #x is dR, y is pT
                                getattr(self, '%s_dR_Pt_charged'%(type)).Fill(pT, dR if isGen else puppi_dR)  
                            nch+=1
                            if abs(pID)!=11 and abs(pID)!=13:
                                if not(noPlot): 
                                    getattr(self, '%s_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, '%s_dRap_dPhi_chargedHadrons'%(type)).Fill(dRap, dPhi)
                                    getattr(self, '%s_dR_chargedHadrons_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, '%s_dRap_dPhi_chargedHadrons_%s' % (type, bin_name)).Fill(dRap, dPhi)

                                    if not(isGen):
                                        
                                        getattr(self, '%s_dRSJ1_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        if dRSJ1_AK8<0.6: getattr(self, '%s_dR2SJ1_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        
                                        getattr(self, '%s_dRSJ2_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        if dRSJ2_AK8<0.6: getattr(self, '%s_dR2SJ2_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies

                                    getattr(self, '%s_pT_chargedHadrons'%(type)).Fill(pT)
                                    getattr(self, '%sjetPt_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())#x is dR, y is pT
                                    getattr(self, '%sjetmSD_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) #x is dR, y is pT
                                    getattr(self, '%s_dR_Pt_chargedHadrons'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())#x is dR, y is pT
                                    
                                nchHad+=1
                if not (self.onlyUnc.startswith('_je')):
                    if not(noPlot):# and not(self.onlyUnc):
                        if central: 
                            getattr(self, '%s_nphotons'%(type)).Fill(npho)
                            getattr(self, '%s_nneutral'%(type)).Fill(nneutHad)
                            getattr(self, '%s_ncharged'%(type)).Fill(nch)
                            getattr(self, '%s_nchargedHadrons'%(type)).Fill(nchHad)

        else:
            if constJES:
                charges = []
                pIDs = []
            for p in pfCands :
                #if p.p4().M()<0.: #to check on -ve mass electrons in constituents
                #    p.p4().M()=0.
                tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                tp = tp * p.puppiWeight if not isGen else tp
                #except RuntimeError: tp = tp    ### for genjets
                CandsPUPPIweightedVec.push_back(tp)
                if constJES:
                    charges.append(p.charge)
                    pIDs.append(p.pdgId)

            #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
            #print ("pushing back candidates")
            if not(constJES):
                for x in CandsPUPPIweightedVec:
                    if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                    #print ("pushed back candidates")
            else:

                modifier= 1. if 'Up' in varUpDown else -1.
                for c,pID,x in zip(charges, pIDs, CandsPUPPIweightedVec):
                    if 'photon' in self.onlyUnc:

                        if pID==22:
                            x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                            x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))

                            if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                        else: 
                            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                    
                    elif ('neutral' in self.onlyUnc):

                         if c==0 and pID!=22:
                            x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                            x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))
                            if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                         else:
                            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                
                    elif ('charged' in self.onlyUnc):
                        
                         if c<0 or c>0:
                            x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                            x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))
                            if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                         else:
                            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)         

        #### Computing n-subjetiness basis from PF PUPPI constituents
        nsub0p25 = self.nSub0p25.getTau( self.maxTau, constituents )
        
        nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
        
        nsub1 = self.nSub1.getTau( self.maxTau, constituents )
        
        nsub1p5 = self.nSub1p5.getTau( self.maxTau, constituents )
        
        nsub2 = self.nSub2.getTau( self.maxTau, constituents )
        
        nsub1_OP_kT = self.nSub1_OP_kT.getTau( 3, constituents )
        nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, constituents )

        ### default in CMS OP_kT https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        try: ak8jet['tau21'] = nsub1_OP_kT[1]/nsub1_OP_kT[0]
        except ZeroDivisionError: ak8jet['tau21'] = -1
        try: ak8jet['tau32'] = nsub1_OP_kT[2]/nsub1_OP_kT[1]
        except ZeroDivisionError: ak8jet['tau32'] = -1

        try: ak8jet['tau21_WTA'] = nsub1_WTA_kT[1]/nsub1_WTA_kT[0]
        except ZeroDivisionError: ak8jet['tau21_WTA'] = -1
        try: ak8jet['tau32_WTA'] = nsub1_WTA_kT[2]/nsub1_WTA_kT[1]
        except ZeroDivisionError: ak8jet['tau32_WTA'] = -1

        #### filling histos and branches with nsub basis
        for tauN in range(self.maxTau):
            ak8jet['0p25'+str(tauN+1)] = nsub0p25[tauN]
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['1p5'+str(tauN+1)] = nsub1p5[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]

        try: ak8jet['tau21_exkT'] = nsub1[1]/nsub1[0]
        except ZeroDivisionError: ak8jet['tau21_exkT'] = -1
        try: ak8jet['tau32_exkT'] = nsub1[2]/nsub1[1]
        except ZeroDivisionError: ak8jet['tau32_exkT'] = -1

        #if isGen: #to add in softdrop mass as a variable to the selected (accep)gen jet branches
        #    sd_AK8jets = self.sd.result( constituents)
        #    if len(sd_AK8jets)>0: #stupidly, in some rare cases, this will not be true (with ptmin>=170) leading to errors, hence, switching to ptmin=0 in function calls (but, since not sure this is error, free I use this if-else block)
        #        ak8jet['msoftdrop'] = sd_AK8jets[0].m()
        #    else: ak8jet['msoftdrop'] = -1. 
        '''
        #### Computing Softdrop jets
        if self.runSDVariables:
            sdAK8jets = self.sd.result( constituents ) #CandsPUPPIweightedVec )

            ### Storing good jet as list for later use
            #if len(sdAK8jets)>0:

            ak8jet['sdjet'] = sdAK8jets[0]

            # Cluster only the particles near the appropriate jet to save time
            sd_constituents =  ROOT.vector("TLorentzVector")()

            for x in ak8jet['sdjet'].constituents():
                sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
                if abs(ak8jet['sdjet'].delta_R( x )) < 0.8:
                    sd_constituents.push_back(sd_constits)
            sd_nsub0p5 = self.nSub0p5.getTau( self.maxTau, sd_constituents )
            sd_nsub1 = self.nSub1.getTau( self.maxTau, sd_constituents )
            sd_nsub2 = self.nSub2.getTau( self.maxTau, sd_constituents )
            sd_nsub1_OP_kT = self.nSub1_OP_kT.getTau( 3, sd_constituents )

            try: ak8jet['sdtau21'] = sd_nsub1_OP_kT[1]/sd_nsub1_OP_kT[0]
            except ZeroDivisionError: ak8jet['sdtau21'] = -1
            try: ak8jet['sdtau32'] = sd_nsub1_OP_kT[2]/sd_nsub1_OP_kT[1]
            except ZeroDivisionError: ak8jet['sdtau32'] = -1

            for tauN in range(self.maxTau):
                ak8jet['sd0p5'+str(tauN+1)] = sd_nsub0p5[tauN]
                ak8jet['sd1'+str(tauN+1)] = sd_nsub1[tauN]
                ak8jet['sd2'+str(tauN+1)] = sd_nsub2[tauN]
        else:
            ak8jet['sdjet'] = ROOT.TLorentzVector( )
            ak8jet['sdtau21'] = -1
            ak8jet['sdtau32'] = -1
            for tauN in range(self.maxTau):
                ak8jet['sd0p5'+str(tauN+1)] = -1
                ak8jet['sd1'+str(tauN+1)] = -1
                ak8jet['sd2'+str(tauN+1)] = -1
        '''
        return ak8jet

    
    #############################################################################
    def fillAK8Branches( self, event, jetLabel, jetInfo, length, dummy=False, sys='_nom' ): 
        #use dummy to fill dummy values into branches for pass/nonpass (gen/reco)sel to maintain correspondence between events in different branches
        dummyFill=-929
        #### Filling branch with passAK8jet info after selection
        #print (jetLabel+"_pt", [ getattr( iJ['jet'], 'pt'+jetLabel.split('Jets')[1] ) )
        
        if not dummy:
            #if event.event%2==True: print ([ iJ['jet'].eta  )
            
            if 'selgenjets' in jetLabel.lower(): self.out.fillBranch('genSelectedEventNumber'+sys, event.event)
            elif 'accepgenjets' in jetLabel.lower(): self.out.fillBranch('accepgenSelectedEventNumber'+sys, event.event)
            elif 'selrecojets' in jetLabel.lower(): self.out.fillBranch('recoSelectedEventNumber'+sys, event.event)
            elif 'truerecojets' in jetLabel.lower(): self.out.fillBranch('truerecoSelectedEventNumber'+sys, event.event)
            
            self.out.fillBranch( 'n'+jetLabel, length )
            c=0
            for i,iJ in jetInfo.items():
                if c==0:
                    if 'reco' in jetLabel.lower():

                        self.out.fillBranch(jetLabel+"_pt_raw",  getattr(iJ['jet'], 'pt_raw')  )
                        self.out.fillBranch(jetLabel+"_mass_raw",  getattr(iJ['jet'], 'mass_raw')  )
                        
                        #JER/C recorrected branches

                        self.out.fillBranch(jetLabel+"_pt",  getattr(iJ['jet'], 'pt'+sys)  )
                        self.out.fillBranch(jetLabel+"_mass",  getattr(iJ['jet'], 'mass'+sys)  )

                        self.out.fillBranch(jetLabel+"_msoftdrop",  getattr(iJ['jet'], 'msoftdrop'+sys) ) # default, new fully corrected softdrop from nanoAOD-tools (with new raw, SJ JEC recorrection) and energy resolution smearing of AK8 applied 

                        self.out.fillBranch(jetLabel+"_mSD",  getattr(iJ['jet'], 'mSD'+sys) )

                        #self.out.fillBranch(jetLabel+"_msoftdrop",  getattr(iJ['jet'], 'msoftdrop') ) # nano default, w/ outdated AK4PUPPI JECs at GT for PAT step, not used anymore for cuts

                        self.out.fillBranch(jetLabel+"_msoftdrop_nom_PUPPICorred",  getattr(iJ['jet'], 'msoftdrop_nom_PUPPICorred') ) # prev. incorrect calc of msoftdrop_nom in nanoAOD-tools, now fixed, as a result of fixing raw, but puppi corrs not for UL, keeping in for posterity/commparisons

                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI",  getattr(iJ['jet'], 'msoftdrop_corr_PUPPI') ) # PUPPI corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_raw",  getattr(iJ['jet'], 'msoftdrop_raw') )
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMS",  getattr(iJ['jet'], 'msoftdrop_corr_JMS') ) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_subjetJEC",  getattr(iJ['jet'], 'msoftdrop_corr_subjetJEC') ) # default subjet JEC corrected softdrop's corr factor vs. raw

                        self.out.fillBranch(jetLabel+"_JERfactor",  getattr(iJ['jet'], 'corr_JER') if self.isMC else 1. ) # JER corr factor
                        self.out.fillBranch(jetLabel+"_JMSfactor",  getattr(iJ['jet'], 'corr_JMS') if self.isMC else 1. ) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_JECfactor",  getattr(iJ['jet'], 'corr_JEC') ) #  if self.isMC else 1. getattr(iJ['jet'], 'pt'+sys)/getattr(iJ['jet'], 'pt_raw') ) # JEC corr factor


                        
                    elif 'gen' in jetLabel.lower():
                        self.out.fillBranch(jetLabel+"_pt",  iJ['jet'].pt  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass",  iJ['jet'].mass  )
                        self.out.fillBranch(jetLabel+"_msoftdrop",  iJ['jet'].msoftdrop  )
                        

                    self.out.fillBranch(jetLabel+"_eta",  iJ['jet'].eta  )
                    self.out.fillBranch(jetLabel+"_y",  iJ['jet'].rapidity  )
                    self.out.fillBranch(jetLabel+"_phi",  iJ['jet'].phi  )

                    self.out.fillBranch(jetLabel+"_tau21",  iJ['tau21']  )
                    self.out.fillBranch(jetLabel+"_tau32",  iJ['tau32']  )
                    self.out.fillBranch(jetLabel+"_tau21_WTA",  iJ['tau21_WTA']  )
                    self.out.fillBranch(jetLabel+"_tau32_WTA",  iJ['tau32_WTA']  )
                    self.out.fillBranch(jetLabel+"_tau21_exkT",  iJ['tau21_exkT']  )
                    self.out.fillBranch(jetLabel+"_tau32_exkT",  iJ['tau32_exkT']  )
            
                    for tauN in range(1, self.maxTau+1):
                        for pref in self.tauPrefixes:
                            beta=pref.split('_')[2]
                            #if 'reco' in jetLabel.lower(): print (event.event, jetLabel+pref+str(tauN), iJ[beta+str(tauN)] )
                            self.out.fillBranch(jetLabel+pref+str(tauN), iJ[beta+str(tauN)]  )
                c+=1
        
        else:
            #fill dummies
            
            if 'accepgenjets' in jetLabel.lower(): self.out.fillBranch('accepgenSelectedEventNumber'+sys, dummyFill)
            elif 'truerecojets' in jetLabel.lower(): self.out.fillBranch('truerecoSelectedEventNumber'+sys, dummyFill)
        
            self.out.fillBranch( 'n'+jetLabel, dummyFill)

            c=0
            for i,iJ in jetInfo.items():
                if c==0:
                    if 'reco' in jetLabel.lower():


                        self.out.fillBranch(jetLabel+"_pt_raw",  dummyFill)
                        self.out.fillBranch(jetLabel+"_mass_raw", dummyFill )

                        #JER/C recorrected branches

                        self.out.fillBranch(jetLabel+"_pt",  dummyFill )
                        self.out.fillBranch(jetLabel+"_mass",  dummyFill )

                        self.out.fillBranch(jetLabel+"_msoftdrop",  dummyFill ) # default, new fully corrected softdrop from nanoAOD-tools (with new raw, SJ JEC recorrection) and energy resolution smearing of AK8 applied

                        self.out.fillBranch(jetLabel+"_mSD",  dummyFill )

                        #self.out.fillBranch(jetLabel+"_msoftdrop",  dummyFill ) # nano default, w/ outdated AK4PUPPI JECs at GT for PAT step, not used anymore for cuts

                        self.out.fillBranch(jetLabel+"_msoftdrop_nom_PUPPICorred",  dummyFill ) # prev. incorrect calc of msoftdrop_nom in nanoAOD-tools, now fixed, as a result of fixing raw, but puppi corrs not for UL, keeping in for posterity/commparisons

                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI",  dummyFill ) # PUPPI corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_raw",  dummyFill )
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMS",  dummyFill) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_subjetJEC",  dummyFill ) # default subjet JEC corrected softdrop's corr factor vs. raw

                        self.out.fillBranch(jetLabel+"_JERfactor",  dummyFill ) # JER corr factor
                        self.out.fillBranch(jetLabel+"_JMSfactor",  dummyFill ) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_JECfactor",  dummyFill ) #  if self.isMC else 1. getattr(iJ['jet'], 'pt'+sys)/getattr(iJ['jet'], 'pt_raw') ) # JEC corr factor


                    elif 'gen' in jetLabel.lower():
                        self.out.fillBranch(jetLabel+"_pt", dummyFill  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass", dummyFill  )
                        self.out.fillBranch(jetLabel+"_msoftdrop", dummyFill  )
                        
                    self.out.fillBranch(jetLabel+"_eta", dummyFill  )
                    self.out.fillBranch(jetLabel+"_y", dummyFill  )
                    self.out.fillBranch(jetLabel+"_phi", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau21", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau32", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau21_WTA", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau32_WTA", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau21_exkT", dummyFill  )
                    self.out.fillBranch(jetLabel+"_tau32_exkT", dummyFill  )

                    for tauN in range(1, self.maxTau+1):
                        for pref in self.tauPrefixes:
                            self.out.fillBranch(jetLabel+pref+str(tauN),  dummyFill  )
                c+=1


    #############################################################################
    def fillOtherBranches( self, event, objectLabel, objectInfo, lenVar, sys='' ): 
                    
        self.out.fillBranch( 'n'+objectLabel, lenVar )
    
        if 'ak4' in objectLabel.lower() or 'mu' in objectLabel.lower() or 'ak8' in objectLabel.lower():
            self.out.fillBranch(objectLabel+"_pt",  ( getattr(objectInfo, 'pt' + sys) )  if ( ('ak' in objectLabel.lower()) or ('jet' in objectLabel.lower()) )  else objectInfo.pt )
            self.out.fillBranch(objectLabel+"_eta",  objectInfo.eta  )
            self.out.fillBranch(objectLabel+"_y",  objectInfo.p4().Rapidity()  )
            self.out.fillBranch(objectLabel+"_phi",  objectInfo.phi  )
            self.out.fillBranch(objectLabel+"_mass",  ( getattr(objectInfo, 'mass' + sys) )  if ( ('ak' in objectLabel.lower()) or ('jet' in objectLabel.lower()) )  else objectInfo.mass )
            if not('gen' in objectLabel.lower()):
                if ('mu' in objectLabel.lower()):
                    #leptonP = ROOT.TMath.Sqrt(objectInfo.p4().Px()**2 + objectInfo.p4().Py()**2 + objectInfo.p4().Pz()**2)
                    self.out.fillBranch(objectLabel+"_p",  ROOT.TMath.Sqrt(objectInfo.p4().Px()**2 + objectInfo.p4().Py()**2 + objectInfo.p4().Pz()**2))
                    self.out.fillBranch(objectLabel+"_ptRel",  objectInfo.jetPtRelv2  )
                    self.out.fillBranch(objectLabel+"_tkRelIso",  objectInfo.tkRelIso  )

                elif 'ak4' in objectLabel.lower():
                    self.out.fillBranch(objectLabel+"_jetId",  objectInfo.jetId)
                    self.out.fillBranch(objectLabel+"_btagDeepFlavB",  objectInfo.btagDeepFlavB  )
                    if self.isMC: 
                        self.out.fillBranch(objectLabel+"_hadronFlavour",  objectInfo.hadronFlavour  )


                    
        else:
            #for MET and leptonic W/top
            self.out.fillBranch(objectLabel+"_pt",  objectInfo.Pt()  )
            self.out.fillBranch(objectLabel+"_phi",  objectInfo.Phi()  )

            if not 'met' in objectLabel.lower(): #for leptonic W
                self.out.fillBranch(objectLabel+"_eta",  objectInfo.Eta()  )
                self.out.fillBranch(objectLabel+"_y",  objectInfo.Rapidity()  )
                self.out.fillBranch(objectLabel+"_mass",  objectInfo.M()  )
                if not('lepttop' in objectLabel.lower()): self.out.fillBranch(objectLabel+"_mt1",  objectInfo.Mt()  )
                #self.out.fillBranch(objectLabel+"_mt",  objectInfo.Mt()  )
            
            
    """
    def createNsubBasis(self, AK8jet, event, PFCollection, isGen=False, constJES=False, varUpDown='' ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        #### Run calculations of NSub bases and store for ungroomed AK8jets (default in CMS)

        #### Applying PUPPI weights to the PF candidates for reco and not for gen and pushing back constituents
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()
        if constJES:
            charges = []
            pIDs = []
        for p in pfCands :
            #if p.p4().M()<0.: #to check on -ve mass electrons in constituents
            #    p.p4().M()=0.
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            tp = tp * p.puppiWeight if not isGen else tp
            #except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)
            if constJES:
                charges.append(p.charge)
                pIDs.append(p.pdgId)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        #print ("pushing back candidates")
        if not(constJES):
            for x in CandsPUPPIweightedVec:
                if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                #print ("pushed back candidates")
        else:

            modifier= 1. if 'Up' in varUpDown else -1.
            for c,pID,x in zip(charges, pIDs, CandsPUPPIweightedVec):
                if 'photon' in self.onlyUnc:

                    if pID==22:
                        x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                        x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))

                        if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                    else: 
                        if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
                
                elif ('neutral' in self.onlyUnc):

                     if c==0 and pID!=22:
                        x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                        x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))
                        if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                     else:
                        if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
            
                elif ('charged' in self.onlyUnc):
                    
                     if c<0 or c>0:
                        x_new = ROOT.TLorentzVector()#x.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                        x_new.SetPtEtaPhiM(x.Pt()*(1+modifier*self.constJESVariation), x.Eta(), x.Phi(), x.M()*(1+modifier*self.constJESVariation))
                        if self.DrRapPhi( AK8jet.p4(), x_new ) < 0.8: constituents.push_back(x_new)
                     else:
                        if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)            


        #### Computing n-subjetiness basis from PF PUPPI constituents
        nsub0p25 = self.nSub0p25.getTau( self.maxTau, constituents )
        
        nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
        
        nsub1 = self.nSub1.getTau( self.maxTau, constituents )
        
        nsub1p5 = self.nSub1p5.getTau( self.maxTau, constituents )
        
        nsub2 = self.nSub2.getTau( self.maxTau, constituents )
        
        nsub1_OP_kT = self.nSub1_OP_kT.getTau( 3, constituents )
        nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, constituents )

        ### default in CMS OP_kT https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        try: ak8jet['tau21'] = nsub1_OP_kT[1]/nsub1_OP_kT[0]
        except ZeroDivisionError: ak8jet['tau21'] = -1
        try: ak8jet['tau32'] = nsub1_OP_kT[2]/nsub1_OP_kT[1]
        except ZeroDivisionError: ak8jet['tau32'] = -1

        try: ak8jet['tau21_WTA'] = nsub1_WTA_kT[1]/nsub1_WTA_kT[0]
        except ZeroDivisionError: ak8jet['tau21_WTA'] = -1
        try: ak8jet['tau32_WTA'] = nsub1_WTA_kT[2]/nsub1_WTA_kT[1]
        except ZeroDivisionError: ak8jet['tau32_WTA'] = -1

        #### filling histos and branches with nsub basis
        for tauN in range(self.maxTau):
            ak8jet['0p25'+str(tauN+1)] = nsub0p25[tauN]
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['1p5'+str(tauN+1)] = nsub1p5[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]

        try: ak8jet['tau21_exkT'] = nsub1[1]/nsub1[0]
        except ZeroDivisionError: ak8jet['tau21_exkT'] = -1
        try: ak8jet['tau32_exkT'] = nsub1[2]/nsub1[1]
        except ZeroDivisionError: ak8jet['tau32_exkT'] = -1

        #if isGen: #to add in softdrop mass as a variable to the selected (accep)gen jet branches
        #    sd_AK8jets = self.sd.result( constituents)
        #    if len(sd_AK8jets)>0: #stupidly, in some rare cases, this will not be true (with ptmin>=170) leading to errors, hence, switching to ptmin=0 in function calls (but, since not sure this is error, free I use this if-else block)
        #        ak8jet['msoftdrop'] = sd_AK8jets[0].m()
        #    else: ak8jet['msoftdrop'] = -1. 
        '''
        #### Computing Softdrop jets
        if self.runSDVariables:
            sdAK8jets = self.sd.result( constituents ) #CandsPUPPIweightedVec )

            ### Storing good jet as list for later use
            #if len(sdAK8jets)>0:

            ak8jet['sdjet'] = sdAK8jets[0]

            # Cluster only the particles near the appropriate jet to save time
            sd_constituents =  ROOT.vector("TLorentzVector")()

            for x in ak8jet['sdjet'].constituents():
                sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
                if abs(ak8jet['sdjet'].delta_R( x )) < 0.8:
                    sd_constituents.push_back(sd_constits)
            sd_nsub0p5 = self.nSub0p5.getTau( self.maxTau, sd_constituents )
            sd_nsub1 = self.nSub1.getTau( self.maxTau, sd_constituents )
            sd_nsub2 = self.nSub2.getTau( self.maxTau, sd_constituents )
            sd_nsub1_OP_kT = self.nSub1_OP_kT.getTau( 3, sd_constituents )

            try: ak8jet['sdtau21'] = sd_nsub1_OP_kT[1]/sd_nsub1_OP_kT[0]
            except ZeroDivisionError: ak8jet['sdtau21'] = -1
            try: ak8jet['sdtau32'] = sd_nsub1_OP_kT[2]/sd_nsub1_OP_kT[1]
            except ZeroDivisionError: ak8jet['sdtau32'] = -1

            for tauN in range(self.maxTau):
                ak8jet['sd0p5'+str(tauN+1)] = sd_nsub0p5[tauN]
                ak8jet['sd1'+str(tauN+1)] = sd_nsub1[tauN]
                ak8jet['sd2'+str(tauN+1)] = sd_nsub2[tauN]
        else:
            ak8jet['sdjet'] = ROOT.TLorentzVector( )
            ak8jet['sdtau21'] = -1
            ak8jet['sdtau32'] = -1
            for tauN in range(self.maxTau):
                ak8jet['sd0p5'+str(tauN+1)] = -1
                ak8jet['sd1'+str(tauN+1)] = -1
                ak8jet['sd2'+str(tauN+1)] = -1
        '''
        return ak8jet
    """
    #############################################################################
    def matchRecoGenParticle( self, event, recoJet ):

        genParticles = Collection(event, "GenPart")

        quarksFromW = GenQuarkFromW( genParticles )
        bquarksFromTop = GenBquarkFromTop( genParticles )
        allQuarksFromWtop = quarksFromW + bquarksFromTop

        listMatched = []
        for q in allQuarksFromWtop:
            #print event.event, q.pdgId, genParticles[q.genPartIdxMother].pdgId, q.p4().Pt()
            if recoJet.p4().DeltaR( q.p4() )<0.3: listMatched.append( True )
            else: listMatched.append( False )

        if (len(listMatched)==4) and all(listMatched[:2]) and not all(listMatched[2:]): boosted = 2  ## only boostedW
        elif (len(listMatched)==4) and all(listMatched[:2]) and any(listMatched[2:]): boosted = 4      ## boosted Top
        else: boosted = 0
        #print listMatched, boosted

        return boosted

    
    #############################################################################
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)

    #############################################################################






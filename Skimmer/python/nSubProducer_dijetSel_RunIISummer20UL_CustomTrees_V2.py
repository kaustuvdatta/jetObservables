# Changelog Oct '22: 
# 1a) updated to correct PS/pdfWeight variation implementation, 
# 1b) added in 2016 trigger turn ons calculated by Alejandro
# 2) added better functionality for storing the nominal selection trees (to be updated for systematics, needs a bit more refactoring)
# 3) udpated DeltaR's to use the y-phi metric when calculating separations from the massive AK8
# 4) Setting ROOT histos to Sumw2 mode by default
# 5) Next step, refactor the code and then only work with ROOT trees, and histos will only be made live during analysis and not on the crab post-processing side
# Changelog Jan/Feb '23
# Refactoring to produce custom trees complete, run with crab using the multicrab_nSubProducer_dijetSel.py in ../test/ which calls the jetObservablesnSubProducer_dijetSel.py which calls this skimmer module
# run skimmer with onlyUnc='', onlyTrees=False and isMC=True/False (isMC=isSigMC=True)  to produce control histograms and nanoskims with nominal (+wt. variations for signal MC)
# use onlyUnc='jer'/'jesCorr'/'jesUncorr', isSigMC=isMC=True, to run with systematics for jer, correlated/uncorrelated jes variations on signal MC samples


import ROOT
import math, os, sys
import numpy as np
import pandas as pd
import root_numpy
from root_numpy import *
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
#import numpy.testing as npt

class nSubProd(Module):

    def __init__(self, sysSource=[], year='2017', isMC=True, onlyUnc='', onlyTrees=False, isSigMC=False, controlPlotsOnly=False, leptonVeto=False, perTriggerPassSel=False):
        self.writeHistFile=True
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.isSigMC = isSigMC
        self.onlyTrees = onlyTrees
        self.constCheckPlots = False

        if '_constituentJES' in self.onlyUnc and (self.onlyUnc.endswith(('charged', 'neutral', 'photon'))):

            if not('photon' in self.onlyUnc):
                self.constJESVariation = 0.05 if 'neutral' in self.onlyUnc else 0.01
            else:
                self.constJESVariation = 0.03

        self.runSDVariables = False
        self.perTriggerPassSel = perTriggerPassSel

        self.nrecoEvents = 0
        self.ngenEvents = 0
        self.nCombTruerecoEvents = 0
        self.nCombAccepgenEvents = 0

        self.nTruerecoEvents = 0
        self.nAccepgenEvents = 0


        self.FlagBadPFCands = False
        self.FlagBadPFCandsCounter=0
        self.evtCounter=0

        self.recoLevel=0
        self.fakes=0
        self.miss=0
        self.genLevel=0
        self.response=0
        self.ufo=0
        self.ufoResponse=0
        self.ufoFake=0
        self.ufoMiss=0
        self.leptonVeto = leptonVeto
        ### Kinematics Cuts AK8Jets ###
        self.minLeadAK8JetPtDijet = 200.
        self.minAK8JetPt = 170.  ### this is the basic minimum, not the final
        #self.maxJetAK8Eta = 2.4
        self.maxAK8JetRap = 1.7
        self.diffPt = [ 200, 10000 ] # 350, 500, 750, 1000,

        ### Kinenatic Cuts Muons ###
        self.minLooseMuonPt = 30.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 20.
        self.maxElectronEta = 2.5
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]
        self.ptBins = [0, 1, 2, 3, 4, 5, 10, 15]
        self.pT_bins_names = ['0_1', '1_2', '2_3', '3_4', '4_5', '5_10', '10_15', 'gt15']

        #self.range1ElectronEta = [0,1.442]
        #self.range2ElectronEta = [1.56,2.4]

        #overall event weights, updated in functions below as required per nominal/systematics runs of the skimmer
        self.totalRecoWeight = 1.
        #self.totalGenWeight = 1.
        
        #self.totalWeight = 1.
        #self.triggerWeight = 1.
        
        # Nominal values of systematic weights handled by the following class-level variables for events passing selections
        # Only storing for some systematics, since storing the nominals for the others don't make sense 
        # (ie, eg, nominal for isr/fsrWeight is 1, and PSWeights store correlated (w_var/w_nom) weights a la ISRup/FSRnom,ISRnom/FSRnom,ISRdown/FSRnom, ISRnom/FSRdown) 
        self.puWeight = 1. # storing nom. puWeight for reco
        self.evtGenWeight = 1. # storing event.genWeight for selected events
        self.pdfWeight = 1. # storing nominal pdfWeight, to see if it deviates from 1 
        #self.psWeight = 1.
        self.pdfWeightAll = np.ones((103,), dtype=np.float32)*1.
        self.pdfWeightUp = 1.
        self.puWeightUp = 1.
        self.isrWeightUp = 1.
        self.fsrWeightUp = 1.
        self.pdfWeightDown = 1.
        self.puWeightDown = 1.
        self.isrWeightDown = 1.
        self.fsrWeightDown = 1.

        self.l1PreFireWeight = 1.
        self.l1PreFireWeightUp = 1.
        self.l1PreFireWeightDown = 1.


        self.eventCategory = -1
        self.dummy = 0

       
        self.triggerTable = OrderedDict()
        #not used if skimmer isn't run with events selected exclusively per trigger, except to have the trigger names available for filling certain branches  (a la passHLT... )

        self.triggerTable['AK8PFJet80'] = {
                                            '2016': [117.649,191.368],
                                            '2016_preVFP': [121.938,192.909],
                                            '2017': [146.47,220.122],
                                            '2018': [141.586,224.871],
                                           }
        self.triggerTable['AK8PFJet140'] = {
                                            '2016': [191.368,263.426],
                                            '2016_preVFP': [192.909,266.569],
                                            '2017': [220.122,294.317],
                                            '2018': [224.871,299.834],
                                           }
        self.triggerTable['AK8PFJet200'] = {
                                            '2016': [263.426,335.765],
                                            '2016_preVFP': [266.569,338.288],
                                            '2017': [294.317,363.808],
                                            '2018': [299.834,372.496],
                                           }
        self.triggerTable['AK8PFJet260'] = {
                                            '2016': [335.765,409.866],
                                            '2016_preVFP': [338.288,411.402],
                                            '2017': [363.808,435.629],
                                            '2018': [372.496,443.459],
                                           }
        self.triggerTable['AK8PFJet320'] = {
                                            '2016': [409.866,509.961],
                                            '2016_preVFP': [411.402,510.262],
                                            '2017': [435.629,530.812],
                                            '2018': [443.459,540.081],
                                           }
        self.triggerTable['AK8PFJet400'] = {
                                            '2016': [509.961,560.765],
                                            '2016_preVFP': [510.262,568.293],
                                            '2017': [530.812,580.872],
                                            '2018': [540.081,585.585],
                                           }
        self.triggerTable['AK8PFJet450'] = {
                                            '2016': [560.765, 6500.],#634.024],
                                            '2016_preVFP': [568.293, 6500.],#631.693],
                                            '2017': [580.872,640.092],
                                            '2018': [585.585,645.01],
                                           }
        if '2017' in self.year or '2018' in self.year:
        
            self.triggerTable['AK8PFJet500'] = {
                                                #'2016': [634.024,6500.0],
                                                #'2016_preVFP': [631.693,6500.0],
                                                '2017': [640.092,6500.0],#,700.13],
                                                '2018': [645.01,6500.0],#,701.393],
                                               }
        #if '2017' in self.year or '2018' in self.year:
        #    self.triggerTable['AK8PFJet550'] = {
        #                                        '2017': [700.13,6500.0],
        #                                        '2018': [701.393,6500.0],
        #                                       }


        ### Defining nsubjetiness basis
        self.maxTau = 5

        self.nSub_labels = {
                        "_tau_0p25_1": [0., 1., 1000  ],
                        "_tau_0p25_2": [0., 1., 1000  ],
                        "_tau_0p25_3": [0., 1., 1000  ],
                        "_tau_0p25_4": [0., 1., 1000  ],
                        "_tau_0p25_5": [0., 1., 1000  ],
                        
                        "_tau_0p5_1": [0., 0.9, 900  ],
                        "_tau_0p5_2": [0., 0.9, 900  ],
                        "_tau_0p5_3": [0., 0.8, 800  ],
                        "_tau_0p5_4": [0., 0.8, 800  ],
                        "_tau_0p5_5": [0., 0.8, 800  ],

                        "_tau_1_1": [ 0., 0.9, 900  ],
                        "_tau_1_2": [ 0., 0.7, 700  ],
                        "_tau_1_3": [ 0., 0.5, 500  ],
                        "_tau_1_4": [ 0., 0.5, 500  ],
                        "_tau_1_5": [ 0., 0.5, 500  ],

                        "_tau_1p5_1": [ 0., 0.9, 900  ],
                        "_tau_1p5_2": [ 0., 0.7, 700  ],
                        "_tau_1p5_3": [ 0., 0.5, 500  ],
                        "_tau_1p5_4": [ 0., 0.5, 500  ],
                        "_tau_1p5_5": [ 0., 0.5, 500  ],

                        "_tau_2_1": [ 0., 0.9, 900  ],
                        "_tau_2_2": [ 0., 0.5, 500  ],
                        "_tau_2_3": [ 0., 0.5, 500  ],
                        "_tau_2_4": [ 0., 0.5, 500  ],
                        "_tau_2_5": [ 0., 0.5, 500  ]
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
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, 0.)#self.minAK8JetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")
        self.controlPlotsOnly = controlPlotsOnly #hackey way to get around long processing times for new control plotss

        ### Uncertainties
        self.sysSource = ['_nom'] + ([ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ] if not(self.controlPlotsOnly) else [])
        #if onlyUnc: self.sysSource = [ onlyUnc+i for i in [ 'Up', 'Down' ] ]
        
        # puWeights used only to change reco event weight (up/down) without application to gen weight, 
        # others applied to modify self.genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)] 
        # weights+up/down variations are saved for accepted events, to be applied in histogramming
          
        self.sysWeightList = ( '_pu', '_pdf', '_isr', '_fsr', '_l1' ) ## '_ps'
        self.sysrecoWeightList = ( '_pu', '_l1' ) 
        self.sysgenWeightList = ( '_pdf', '_isr', '_fsr' ) #, '_ps'



    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        tauBins = 100

        ### Booking histograms
        selList = [ '_'+x+'_dijetSel' for x in self.triggerTable  ] if not self.isMC else [ '_dijetSel' ]
        #selList = ([ '_'+x+'_dijetSel' for x in self.triggerTable  ] + [ '_'+x+'_weight_dijetSel' for x in self.triggerTable  ]) if not self.isMC else [ '_dijetSel' ]
        
        
        #if not self.isMC:
        #    for isel in [ '_only'+x+'_dijetSel' for x in self.triggerTable  ]:
        #        self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
        #        self.addP4Hists( 'LeadPtJetAK8', isel )
        #        self.addP4Hists( 'SubleadPtJetAK8', isel )
        #        self.addP4Hists( 'CentralJetAK8', isel )
        #        self.addP4Hists( 'ForwardJetAK8', isel )

        if not (self.onlyUnc.startswith('_je')):# or self.controlPlotsOnly ):
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   12, 0, 12) ) 
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            if not self.isMC:
                for isel in [ '_only'+x+'_dijetSel' for x in self.triggerTable  ]:
                    self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                    self.addP4Hists( 'LeadPtJetAK8', isel )
                    self.addP4Hists( 'SubleadPtJetAK8', isel )
                    self.addP4Hists( 'CentralJetAK8', isel )
                    self.addP4Hists( 'ForwardJetAK8', isel )

            if self.isMC and self.constCheckPlots and not(self.onlyUnc.startswith('_je')):

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


                self.addObject( ROOT.TH1F('Fgen_dR_charged',   '; #Delta R(const_{ch}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Fgen_dRap_dPhi_charged',   '; #Delta y(const_{ch}, genAK8); #Delta #phi(const_{ch}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Fgen_dR_chargedHadrons',   '; #Delta R(h_{ch}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Fgen_dRap_dPhi_chargedHadrons',   '; #Delta y(h_{ch}, genAK8); #Delta #phi(h_{ch}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Fgen_dR_neutral',   '; #Delta R(h_{neut.}, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Fgen_dRap_dPhi_neutral',   '; #Delta y(h_{neut.}, genAK8); #Delta #phi(h_{neut.}, genAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Fgen_dR_photons',   '; #Delta R(#gamma, genAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Fgen_dRap_dPhi_photons',   '; #Delta y(#gamma, genAK8); #Delta #phi(#gamma, genAK8)', 100, -1., 1., 100, -1., 1.) )
                

                self.addObject( ROOT.TH1F('Freco_dR_charged',   '; #Delta R(const_{ch}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Freco_dRap_dPhi_charged',   '; #Delta y(const_{ch}, recoAK8); #Delta #phi(const_{ch}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Freco_dR_chargedHadrons',   '; #Delta R(h_{ch}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Freco_dRap_dPhi_chargedHadrons',   '; #Delta y(h_{ch}, recoAK8); #Delta #phi(h_{ch}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Freco_dR_neutral',   '; #Delta R(h_{neut.}, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Freco_dRap_dPhi_neutral',   '; #Delta y(h_{neut.}, recoAK8); #Delta #phi(h_{neut.}, recoAK8)', 100, -1., 1., 100, -1., 1.) )

                self.addObject( ROOT.TH1F('Freco_dR_photons',   '; #Delta R(#gamma, recoAK8)', 100, 0, 1.) )
                self.addObject( ROOT.TH2F('Freco_dRap_dPhi_photons',   '; #Delta y(#gamma, recoAK8); #Delta #phi(#gamma, recoAK8)', 100, -1., 1., 100, -1., 1.) )


                self.addObject( ROOT.TH1F('Fgen_ncharged',   '; n const_{ch}in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Fgen_nchargedHadrons',   '; n h_{ch} in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Fgen_nneutral',   '; n h_{neut.} in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Fgen_nphotons',   '; n #gamma in genAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Freco_ncharged',   '; n const_{ch} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Freco_nchargedHadrons',   '; n h_{ch} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Freco_nneutral',   '; n h_{neut.} in recoAK8', 100, 0, 100) )
                self.addObject( ROOT.TH1F('Freco_nphotons',   '; n #gamma in recoAK8', 100, 0, 100) )


                self.addObject( ROOT.TH1F('Fgen_pT_charged',   '; p_{T,ch}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Fgen_pT_chargedHadrons',   '; p_{T,ch had}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Fgen_pT_neutral',   '; p_{T,neut.}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Fgen_pT_photons',   '; p_{T,#gamma}', 200, 0, 100.) )

                self.addObject( ROOT.TH1F('Freco_pT_charged',   '; p_{T,ch}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Freco_pT_chargedHadrons',   '; p_{T,ch had}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Freco_pT_neutral',   '; p_{T,neut.}', 200, 0, 100.) )
                self.addObject( ROOT.TH1F('Freco_pT_photons',   '; p_{T,#gamma}', 200, 0, 100.) )


                self.addObject( ROOT.TH1F('Freco_dR_SJ1',   '; #Delta R(AK8_F, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR_SJ2',   '; #Delta R(AK8_F, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2_SJ1',   '; #Delta R2(AK8_F, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2_SJ2',   '; #Delta R2(AK8_F, AK8_F SJ2)', 100, 0, 1.) )


                self.addObject( ROOT.TH1F('Freco_dRSJ1_charged',   '; #Delta R(const_{ch}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ1_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ1_neutral',   '; #Delta R(h_{neut.}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ1_photons',   '; #Delta R(#gamma, AK8_F SJ1)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('Freco_dRSJ2_charged',   '; #Delta R(const_{ch}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ2_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ2_neutral',   '; #Delta R(h_{neut.}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dRSJ2_photons',   '; #Delta R(#gamma, AK8_F SJ2)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('Freco_dR2SJ1_charged',   '; #Delta R2(const_{ch}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ1_chargedHadrons',   '; #Delta R2(h_{ch}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ1_neutral',   '; #Delta R2(h_{neut.}, AK8_F SJ1)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ1_photons',   '; #Delta R2(#gamma, AK8_F SJ1)', 100, 0, 1.) )

                self.addObject( ROOT.TH1F('Freco_dR2SJ2_charged',   '; #Delta R2(const_{ch}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ2_chargedHadrons',   '; #Delta R2(h_{ch}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ2_neutral',   '; #Delta R2(h_{neut.}, AK8_F SJ2)', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('Freco_dR2SJ2_photons',   '; #Delta R2(#gamma, AK8_F SJ2)', 100, 0, 1.) )


                self.addObject( ROOT.TH2F('FgenjetPt_dR_charged',   '; #Delta R(const_{ch}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FgenjetPt_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FgenjetPt_dR_neutral',   '; #Delta R(h_{neut.}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FgenjetPt_dR_photons',   '; #Delta R(#gamma, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )

                self.addObject( ROOT.TH2F('FrecojetPt_dR_charged',   '; #Delta R(const_{ch}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FrecojetPt_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FrecojetPt_dR_neutral',   '; #Delta R(h_{neut.}, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )
                self.addObject( ROOT.TH2F('FrecojetPt_dR_photons',   '; #Delta R(#gamma, AK8_F); p_{T} AK8_F', 100, 0, 1., 120, 0, 1200.) )


                self.addObject( ROOT.TH2F('FgenjetmSD_dR_charged',   '; #Delta R(const_{ch}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FgenjetmSD_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FgenjetmSD_dR_neutral',   '; #Delta R(h_{neut.}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FgenjetmSD_dR_photons',   '; #Delta R(#gamma, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )

                self.addObject( ROOT.TH2F('FrecojetmSD_dR_charged',   '; #Delta R(const_{ch}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FrecojetmSD_dR_chargedHadrons',   '; #Delta R(h_{ch}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FrecojetmSD_dR_neutral',   '; #Delta R(h_{neut.}, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )
                self.addObject( ROOT.TH2F('FrecojetmSD_dR_photons',   '; #Delta R(#gamma, AK8_F); m_{SD} AK8_F', 100, 0, 1., 30, 50, 300.) )


                self.addObject( ROOT.TH2F('Fgen_dR_Pt_charged',   '; p_{T} const_{ch};  #Delta R(const_{ch}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Fgen_dR_Pt_chargedHadrons',   '; p_{T} h_{ch};  #Delta R(h_{ch}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Fgen_dR_Pt_neutral',   '; p_{T} h_{neut.};  #Delta R(h_{neut.}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Fgen_dR_Pt_photons',   '; p_{T} #gamma;  #Delta R(#gamma, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )

                self.addObject( ROOT.TH2F('Freco_dR_Pt_charged',   '; p_{T} const_{ch};  #Delta R(const_{ch}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Freco_dR_Pt_chargedHadrons',   '; p_{T} h_{ch};  #Delta R(h_{ch}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Freco_dR_Pt_neutral',   '; p_{T} h_{neut.};  #Delta R(h_{neut.}, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )
                self.addObject( ROOT.TH2F('Freco_dR_Pt_photons',   '; p_{T} #gamma;  #Delta R(#gamma, AK8_F); ',  200, 0, 100., 100, 0, 1.,) )


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

                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('Fgen_dR_chargedHadrons_{}'.format(bin_name), '; #Delta R(h_{{ch}}, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Fgen_dRap_dPhi_chargedHadrons_{}'.format(bin_name), '; #Delta y(h_{{ch}}, genAK8, pT {}); #Delta #phi(h_{{ch}}, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('Freco_dR_chargedHadrons_{}'.format(bin_name), '; #Delta R(h_{{ch}}, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Freco_dRap_dPhi_chargedHadrons_{}'.format(bin_name), '; #Delta y(h_{{ch}}, recoAK8, pT {}); #Delta #phi(h_{{ch}}, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))

                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('Fgen_dR_neutral_{}'.format(bin_name), '; #Delta R(h_{{neut.}}, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Fgen_dRap_dPhi_neutral_{}'.format(bin_name), '; #Delta y(h_{{neut.}}, genAK8, pT {}); #Delta #phi(h_{{neut.}}, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('Freco_dR_neutral_{}'.format(bin_name), '; #Delta R(h_{{neut.}}, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Freco_dRap_dPhi_neutral_{}'.format(bin_name), '; #Delta y(h_{{neut.}}, recoAK8, pT {}); #Delta #phi(h_{{neut.}}, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))

                for bin_name in self.pT_bins_names:
                    self.addObject(ROOT.TH1F('Fgen_dR_photons_{}'.format(bin_name), '; #Delta R(#gamma, genAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Fgen_dRap_dPhi_photons_{}'.format(bin_name), '; #Delta y(#gamma, genAK8, pT {}); #Delta #phi(#gamma, genAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))
                    self.addObject(ROOT.TH1F('Freco_dR_photons_{}'.format(bin_name), '; #Delta R(#gamma, recoAK8, pT {})'.format(bin_name), 100, 0, 1.))
                    self.addObject(ROOT.TH2F('Freco_dRap_dPhi_photons_{}'.format(bin_name), '; #Delta y(#gamma, recoAK8, pT {}); #Delta #phi(#gamma, recoAK8)'.format(bin_name), 100, -1., 1., 100, -1., 1.))

    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  ';p_{T} (GeV)',   200, 0, 3000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, ';#eta', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_y'+t, ';y', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, ';#phi', 100, -3.5, 3.5) )
        if not ('mu' in s.lower() or 'electron' in s.lower() or 'met' in s.lower()): 
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
        tmplist=[]
        tmplist2=[]
        if not self.isMC:

            self.out.branch( 'eventCategory_nom',  "I")
            self.out.branch( 'totalRecoWeight_nom', "F" )
            self.out.branch( 'passRecoSel_nom', "I" )
            self.out.branch( 'recoSelectedEventNumber_nom', "L" )
            self.out.branch( 'good_nPVs_nom', "F" )
            
            self.out.branch( 'nRecoLeptons_nom', "I" )        
            
            self.out.branch( 'pt_asymm_nom', "F" )
            self.out.branch( 'delta_phi_nom', "F" )
            self.out.branch( 'delta_R_nom', "F" )
            

            tmplist.append( 'selRecoJets_nom')
            tmplist.append( 'selRecoJetsF_nom')
            tmplist.append( 'selRecoLeadingJets_nom')
            tmplist.append( 'selRecoSubleadingJets_nom'  )
            
            for x in self.triggerTable.keys():
                self.out.branch( 'passHLT_'+x, "I" )
            if not('2016' in self.year): self.out.branch( 'passHLT_AK8PFJet550',"I" )# 1 if getattr(event, 'HLT_AK8PFJet550')==1 and passRecoSel[sys] else 0)  

        elif self.isMC:
            for sys in self.sysSource: 
                self.out.branch( 'eventCategory'+sys,  "I")

                self.out.branch( 'totalRecoWeight'+sys, "F" )
                self.out.branch( 'passRecoSel'+sys, "I" )
                self.out.branch( 'puWeightNom'+sys, "F" )
                self.out.branch( 'l1prefiringWeightNom'+sys, "F")

                self.out.branch( 'pt_asymm'+sys, "F" )
                self.out.branch( 'delta_phi'+sys, "F" )
                self.out.branch( 'delta_R'+sys, "F" )

                if sys.endswith('nom'):
                    self.out.branch( 'gen_pt_asymm'+sys, "F" )
                    self.out.branch( 'gen_delta_phi'+sys, "F" )
                    self.out.branch( 'gen_delta_R'+sys, "F" )
                    self.out.branch( 'good_nPVs'+sys, "F" )
                    self.out.branch( 'nRecoLeptons'+sys, "I" )
                    self.out.branch( 'nGenLeptons'+sys, "I" )


                self.out.branch( 'evtGenWeight'+sys, "F" )
                self.out.branch( 'passGenSel'+sys, "I" )
                self.out.branch( 'recoSelectedEventNumber'+sys, "L" )
                self.out.branch( 'genSelectedEventNumber'+sys, "L" )
                self.out.branch( 'accepgenSelectedEventNumber'+sys, "L" )
                self.out.branch( 'truerecoSelectedEventNumber'+sys, "L" )

                if self.isSigMC and not self.onlyUnc:
                    self.out.branch( 'pdfWeightNom'+sys, "F" )
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


            tmplist_reco = [ 'selRecoJets'+sys for sys in self.sysSource ]
            tmplist_truereco = ['trueRecoJets'+sys for sys in self.sysSource]
            tmplist_gen = ['selGenJets'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            tmplist_accepgen = ['accepGenJets'+sys for sys in self.sysSource]

            for x in tmplist_reco+tmplist_truereco +tmplist_gen+tmplist_accepgen : tmplist.append(x)


            tmplist_reco = [ 'selRecoJetsF'+sys for sys in self.sysSource ]
            tmplist_gen = ['selGenJetsF'+sys for sys in self.sysSource if sys.startswith('_nom')] 

            for x in tmplist_reco+tmplist_gen: tmplist.append(x)


            #if self.controlPlotsOnly: #hack to run only creating some new branches
            tmplist_reco2 = [ 'selRecoLeadingJets'+sys for sys in self.sysSource ]
            tmplist_gen2 = ['selGenLeadingJets'+sys for sys in self.sysSource if sys.startswith('_nom')]
            tmplist_reco3 = [ 'selRecoSubleadingJets'+sys for sys in self.sysSource ]
            tmplist_gen3 = ['selGenSubleadingJets'+sys for sys in self.sysSource if sys.startswith('_nom')]
            for x in tmplist_reco2+tmplist_reco3+tmplist_gen2+tmplist_gen3: tmplist.append(x)

            
            tmplist_truereco = ['trueRecoJetsF'+sys for sys in self.sysSource]
            tmplist_accepgen = ['accepGenJetsF'+sys for sys in self.sysSource]
            #    
            for x in tmplist_truereco+tmplist_accepgen: tmplist.append(x)

            for x in self.triggerTable.keys():
                self.out.branch( 'passHLT_'+x, "I" )
            if not('2016' in self.year): self.out.branch( 'passHLT_AK8PFJet550',"I" )

        print ("Stored jet branches:", tmplist)
        for iJ in tmplist:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_y',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_msoftdrop',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_mSD',  'F')# hack, new branch to store recorrected msoftdrop and circumvent name clashes or overwriting nanoAOD-tools branches

            #if not(self.controlPlotsOnly):#avoid going through the bulk of the info being saved and processed otherwise
            if not('gen' in iJ.lower()): 
                self.out.branch(iJ+'_pt_raw',  'F') 
                #self.out.branch(iJ+'_pt_nom',  'F')# a la nanoAOD-tools recorrection with JECs in postproc. 
                self.out.branch(iJ+'_msoftdrop_nom',  'F')  # a la nanoAOD-tools, stored for posterity as its own branch in custom flat, unjagged selRecoJets collection
                self.out.branch(iJ+'_msoftdrop_corr_PUPPI',  'F') 
                self.out.branch(iJ+'_msoftdrop_raw',  'F') 
                self.out.branch(iJ+'_msoftdrop_corr_JMS',  'F') 
                self.out.branch(iJ+'_msoftdrop_corr_JMR',  'F') 
                self.out.branch(iJ+'_msoftdrop_JMAR',  'F') 
                self.out.branch(iJ+'_JERfactor',  'F') 
                self.out.branch(iJ+'_JECfactor',  'F') 

            self.out.branch(iJ+'_tau21',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_exkT',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_exkT',  'F')#, lenVar='n'+iJ)
            
            for x in list(self.nSub_labels.keys()):
                self.out.branch(iJ+x, 'F')#, lenVar='n'+iJ )
        
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #if not self.onlyTrees:
        if self.isMC and not self.onlyUnc:
            self.genLevel = self.response+self.miss
        '''
        if self.isMC and not self.onlyUnc:
            self.genLevel = self.response+self.miss

            #getattr( self, 'cutflow_test' ).SetBinContent( 1, self.recoLevel )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 1, "allReco" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 2, self.response )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 2, "allTrueRecoAccepGen" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 3, self.fakes )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 3, "allFakes" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 4, self.miss )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 4, "allMiss" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 5, self.genLevel )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 5, "allGen" )

            #getattr( self, 'cutflow_test' ).SetBinContent( 7, self.ufo )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 7, "allUfo" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 8, self.ufoResponse )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 8, "ufoResponse" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 9, self.ufoFake )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 9, "ufoFake" )
            #getattr( self, 'cutflow_test' ).SetBinContent( 10, self.ufoMiss )
            #getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 10, "ufoMiss" )
        '''


    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''
        self.isMC = event.run == 1
        #print(event.event)
        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenJets, gen_deltaR, gen_ptAsymm, gen_deltaPhi, gen_nleptons = self.genSelection( event)#gen_ncrackElectrons,selGenMuons, selGenElectrons,
        passRecoSel, iRecoSel, selRecoJets, reco_deltaR, reco_ptAsymm, reco_deltaPhi, reco_nleptons = self.recoSelection( event )#reco_ncrackElectrons,selRecoMuons, selRecoElectrons,
        #print(event.event)


        if not( self.isMC):
            if not passRecoSel['_nom']: 
                self.totalRecoWeight=0.
                self.out.fillBranch( 'recoSelectedEventNumber_nom', -1  )
                self.out.fillBranch( 'eventCategory_nom', -1  )
                self.out.fillBranch( 'good_nPVs_nom', getattr( event, 'PV_npvsGood') )
                self.out.fillBranch( 'pt_asymm_nom', -929  )
                self.out.fillBranch( 'delta_phi_nom', 929  )
                self.out.fillBranch( 'delta_R_nom', -929  )
                self.out.fillBranch( 'totalRecoWeight_nom', 0.)
                self.out.fillBranch( 'passRecoSel_nom', 0)
                self.out.fillBranch( 'nRecoLeptons_nom', 0)

                return False

        elif self.isMC and (not(self.onlyUnc) or (self.onlyUnc and 'const' in self.onlyUnc)):
            if not( passGenSel) and not( passRecoSel['_nom']): 
                self.totalRecoWeight=0.
                self.evtGenWeight=0.
                self.out.fillBranch( 'recoSelectedEventNumber_nom', -1  )
                self.out.fillBranch( 'genSelectedEventNumber_nom', -1  )
                self.out.fillBranch( 'truerecoSelectedEventNumber_nom', -1  )
                self.out.fillBranch( 'accepgenSelectedEventNumber_nom', -1  )
                self.out.fillBranch( 'eventCategory_nom', -1)
                self.out.fillBranch( 'good_nPVs_nom', getattr( event, 'PV_npvsGood') )

                self.out.fillBranch( 'pt_asymm_nom', -929  )
                self.out.fillBranch( 'delta_phi_nom', 929  )
                self.out.fillBranch( 'delta_R_nom', -929  )
                
                self.out.fillBranch( 'nRecoLeptons_nom', 0)

                self.out.fillBranch( 'nGenLeptons_nom', 0)

                self.out.fillBranch( 'gen_pt_asymm_nom', -929  )
                self.out.fillBranch( 'gen_delta_phi_nom', 929  )
                self.out.fillBranch( 'gen_delta_R_nom', -929  )

                self.out.fillBranch( 'totalRecoWeight_nom', 0.)
                self.out.fillBranch( 'passRecoSel_nom', 0)
                self.out.fillBranch( 'puWeightNom_nom', 0.)
                self.out.fillBranch( 'l1prefiringWeightNom_nom', 0.)
                self.out.fillBranch( 'evtGenWeight_nom', 0.) 
                self.out.fillBranch( 'passGenSel_nom', 0)

                if self.isSigMC and not self.onlyUnc:
                    self.out.fillBranch( 'pdfWeightNom_nom', 0.)
                    self.out.fillBranch( 'pdfWeightAll_nom', np.zeros((103,),dtype=np.float32))
                    self.out.fillBranch( 'pdfWeightUp_nom', 0.)
                    self.out.fillBranch( 'pdfWeightDown_nom', 0.)
                
                    self.out.fillBranch( 'isrWeightUp_nom', 0.)
                    self.out.fillBranch( 'isrWeightDown_nom', 0.)
                    
                    self.out.fillBranch( 'fsrWeightUp_nom', 0.)
                    self.out.fillBranch( 'fsrWeightDown_nom', 0.)
                    
                    self.out.fillBranch( 'puWeightUp_nom', 0.)
                    self.out.fillBranch( 'puWeightDown_nom', 0.)

                    self.out.fillBranch( 'l1prefiringWeightUp_nom', 0.)
                    self.out.fillBranch( 'l1prefiringWeightDown_nom', 0.)

                return False

        elif self.isSigMC and self.onlyUnc:

            sysPassFlag=False
            for sys in self.sysSource:
                if passRecoSel[sys] or passGenSel: 
                    sysPassFlag=True
                    break
            if sysPassFlag==False: 
                for sys in self.sysSource:        
                    #if sys=='_nom': continue
                    self.totalRecoWeight=0.
                    self.evtGenWeight=0.
                    self.out.fillBranch( 'eventCategory'+sys, -1)
                    self.out.fillBranch( 'recoSelectedEventNumber'+sys, -1  )
                    self.out.fillBranch( 'genSelectedEventNumber'+sys, -1  )
                    self.out.fillBranch( 'truerecoSelectedEventNumber'+sys, -1  )
                    self.out.fillBranch( 'accepgenSelectedEventNumber'+sys, -1  )
                    if sys.endswith('nom'):
                        self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood'))
                        self.out.fillBranch( 'gen_pt_asymm'+sys, -929  )
                        self.out.fillBranch( 'gen_delta_phi'+sys, 929  )
                        self.out.fillBranch( 'gen_delta_R'+sys, -929  )
                        self.out.fillBranch( 'nRecoLeptons'+sys, 0)
                        self.out.fillBranch( 'nGenLeptons'+sys, 0)
                

                    self.out.fillBranch( 'pt_asymm'+sys, -929  )
                    self.out.fillBranch( 'delta_phi'+sys, 929  )
                    self.out.fillBranch( 'delta_R'+sys, -929  )
                

                    self.out.fillBranch( 'totalRecoWeight'+sys, 0.)
                    self.out.fillBranch( 'passRecoSel'+sys, 0)
                    
                    self.out.fillBranch( 'puWeightNom'+sys, 0.)
                    self.out.fillBranch( 'l1prefiringWeightNom'+sys, 0.)
                    
                    self.out.fillBranch( 'evtGenWeight'+sys, 0.) 
                    self.out.fillBranch( 'passGenSel'+sys, 0) 
                    
                return False
        
        self.evtCounter+=1.

        for sys in self.sysSource:

            #if not self.isMC: 
            self.totalRecoWeight=1.
            if self.isMC:
                
                self.puWeight = event.puWeight
                self.evtGenWeight=event.genWeight
                try: 
                    self.l1PreFireWeight = event.L1PreFiringWeight_Nom                
                    self.totalRecoWeight = event.genWeight*event.puWeight*event.L1PreFiringWeight_Nom 
                except RuntimeError:  # to deal with 2018 QCD_Pt samples not having the l1 weight branches for some reason
                    self.l1PreFireWeight = 1.
                    self.totalRecoWeight = event.genWeight*event.puWeight*1.#event.L1PreFiringWeight_Nom 
                 
                if not sys.startswith(('_jes', '_jer')):
                   
                    selRecoJets[sys] = selRecoJets['_nom']
                    iRecoSel[sys] = iRecoSel['_nom']
                    passRecoSel[sys] = passRecoSel['_nom']

                    # PDF sets for most/all(to be checked if true) RunIISummer20UL samples seem to be NNPDF31_nnlo_as_0.118_mc_hessian_pdfas (pdfid=325300), definitely true for QCD_HT signal MC for dijets
                    # structure of the pdf set's array of 103 members a la: https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
                    # general structure of pdf weights = w_i,var/w_nom, where i runeth over 0-102, where [0]=>nominal from avg. over replicas, 
                    # [1-100]=> PDF eigenvectors of the covariance matrix in the parameter space, 
                    # [101,102]=> central value for (forced positive definite) a_s=[0.116,0.120]
                    if self.isSigMC:
                        pdfWeights =  getattr( event, 'LHEPdfWeight' ) #convert to a simple numpy array to extract replica weights and thereby the variations
                        self.pdfWeightAll = np.array([pdfWeights[i] for i in range(pdfWeights.GetSize())], dtype=np.float32)
                        self.pdfWeight = pdfWeights[0]
                        
                        # Naively defining pdf uncertainty symmetrically by taking the envelope defined as per Eq. 15 of https://link.springer.com/article/10.1140/epjc/s10052-015-3318-8#Sec29 
                        # a la sigma_up/dn=sigma=np.sqrt(np.sum([(arr[i]-arr[0])**2. for i in np.arange(1,103)])) from symmhessian diff calc; 
                        # however i'm not sure if this is also applicable when eig vect contributions are stored as weights, still in tests seems to give less conservative sigma estimate than doing max(arr[i]-arr[0]) 
                        # and adding that in as the up/down variation), so saving all pdf weights anyway in case I need to unfold 1+102 times with reweighted response matrices

                        pdf_std = np.sqrt(np.sum([(self.pdfWeightAll[i]-self.pdfWeightAll[0])**2. for i in range(1,103)]))
                        self.pdfWeightUp = self.pdfWeight + pdf_std #central+sigma
                        self.pdfWeightDown = self.pdfWeight - pdf_std #central-sigma
                        self.isrWeightDown = getattr( event, 'PSWeight' )[2]
                        self.fsrWeightDown = getattr( event, 'PSWeight' )[3]
                        self.isrWeightUp = getattr( event, 'PSWeight' )[0]
                        self.fsrWeightUp = getattr( event, 'PSWeight' )[1]
                        
                        self.puWeightUp = event.puWeightUp
                        self.puWeightDown = event.puWeightDown
                        try:
                            self.l1PreFireWeightUp = event.L1PreFiringWeight_Up
                            self.l1PreFireWeightDown = event.L1PreFiringWeight_Dn
                        except RuntimeError:
                            self.l1PreFireWeightUp = 1.
                            self.l1PreFireWeightDown = 1.
                    
                else: 
                    #could remove this else block, doesn't do anything (will remove when running final check)
                    try: self.totalRecoWeight = event.genWeight*event.puWeight*event.L1PreFiringWeight_Nom #self.totalRecoWeight*self.puWeight
                    except RuntimeError: self.totalRecoWeight = event.genWeight*event.puWeight*1.
            else: 
                self.totalRecoWeight=1.        
            
            genJet = OrderedDict() if (self.isMC) else None
            recoJet = OrderedDict()    
            

            #if self.controlPlotsOnly:
            genJetLeading = OrderedDict() if (self.isMC) else None
            recoJetLeading = OrderedDict()
            genJetSubleading = OrderedDict() if (self.isMC) else None
            recoJetSubleading = OrderedDict()
        
            genJetF = OrderedDict() if (self.isMC) else None
            recoJetF = OrderedDict()

            tmpRecoJets = OrderedDict()
            tmpGenJets = OrderedDict() if self.isMC else None
            tmpRecoJets0 = OrderedDict()
            tmpGenJets0 = OrderedDict() if self.isMC else None


            if self.isMC and not( passRecoSel[sys]):# and not self.onlyUnc.startswith(self.sysWeightList):

                if passGenSel:# and not self.onlyTrees:

                    #### Misses
                    self.miss = self.miss+1
                    self.eventCategory = 3
                    tmpGenJets={}
                    tmpGenJets0={}
                    #if not(self.controlPlotsOnly):
                    tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True ) # add PUPPI weights, and, if not(self.controlPlotsOnly) nsub bases as well
                    tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )

                        
                    #if self.controlPlotsOnly:
                    genJetLeading['Jet'] = tmpGenJets[0]
                    genJetSubleading['Jet'] = tmpGenJets[1]

                    if abs(tmpGenJets[0]['jet'].rapidity) > abs(tmpGenJets[1]['jet'].rapidity):

                        tmpGenJets0[0] = self.createNsubBasis( AK8jet = selGenJets[0], event = event, PFCollection = 'GenCands', isGen = True, constJES = False, varUpDown = '', central = False, noPlot = False) if self.constCheckPlots else tmpGenJets[0]# add PUPPI weights, and, if not(self.controlPlotsOnly) nsub bases as well
                        tmpGenJets0[1] = self.createNsubBasis( AK8jet = selGenJets[1], event = event, PFCollection = 'GenCands', isGen = True, constJES = False, varUpDown = '', central = True, noPlot = False) if self.constCheckPlots else tmpGenJets[1]

                        genJet['Jet'] = tmpGenJets0[1]#self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                        genJetF['Jet'] = tmpGenJets0[0]
                    else: 
                        tmpGenJets0[0] = self.createNsubBasis( AK8jet = selGenJets[0], event = event, PFCollection = 'GenCands', isGen = True, constJES = False, varUpDown = '', central = True, noPlot = False) if self.constCheckPlots else tmpGenJets[0]# add PUPPI weights, and, if not(self.controlPlotsOnly) nsub bases as well
                        tmpGenJets0[1] = self.createNsubBasis( AK8jet = selGenJets[1], event = event, PFCollection = 'GenCands', isGen = True, constJES = False, varUpDown = '', central = False, noPlot = False) if self.constCheckPlots else tmpGenJets[1]

                        genJet['Jet'] = tmpGenJets0[0]#self.createNsubBasis( selGenJets[0], event, 'GenCands', True  )
                        genJetF['Jet'] = tmpGenJets0[1]

                    #self.ngenEvents+=1
                    if sys.startswith('_nom'): 

                        self.fillBranches( event, 'selGenJets'+sys, genJet, False, sys ) 
                        #if self.controlPlotsOnly:
                        self.fillBranches( event, 'selGenLeadingJets'+sys, genJetLeading, False, sys ) 
                        self.fillBranches( event, 'selGenSubleadingJets'+sys, genJetSubleading, False, sys ) 
                        #else:

                        self.fillBranches( event, 'selGenJetsF'+sys, genJetF, False, sys ) 
                    self.fillBranches( event, 'accepGenJets'+sys, genJet, True, sys )

                    self.fillBranches( event, 'accepGenJetsF'+sys, genJetF, True, sys )

                else: self.evtGenWeight=0.
                self.totalRecoWeight=0.

            if passRecoSel[sys]:  #### Detector level dist.
                
                self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel
                self.eventCategory = 1
                #self.nrecoEvents+=1

                tmpRecoJets[sys]={}
                tmpRecoJets0[sys]={}

                tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                tmpRecoJets[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )

                #if self.controlPlotsOnly:
                recoJetLeading['Jet'] = tmpRecoJets[sys][0]
                recoJetSubleading['Jet'] = tmpRecoJets[sys][1]
                    

                

                if passGenSel:
                    tmpGenJets={}
                    tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                    tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                    #if self.controlPlotsOnly:
                    genJetLeading['Jet'] = tmpGenJets[0]
                    genJetSubleading['Jet'] = tmpGenJets[1]
    
                   



                if abs(tmpRecoJets[sys][0]['jet'].rapidity) > abs(tmpRecoJets[sys][1]['jet'].rapidity):
                    tmpRecoJets0[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands', False, False, '', False, False) if self.constCheckPlots else tmpRecoJets[sys][0]
                    tmpRecoJets0[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands', False, False, '', True, False) if self.constCheckPlots else tmpRecoJets[sys][1]

                    recoJet['Jet'] = tmpRecoJets0[sys][1]
                    recoJetF['Jet'] = tmpRecoJets0[sys][0]#self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )
                else: 
                    tmpRecoJets0[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands', False, False, '', True, False) if self.constCheckPlots else tmpRecoJets[sys][0]
                    tmpRecoJets0[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands', False, False, '', False, False) if self.constCheckPlots else tmpRecoJets[sys][1]

                    recoJet['Jet'] = tmpRecoJets0[sys][0]#self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )                    
                    recoJetF['Jet'] = tmpRecoJets0[sys][1]


                self.fillBranches( event, 'selRecoJets'+sys, recoJet, False, sys ) #not a dummy fill so dummy=False
                self.fillBranches( event, 'selRecoJetsF'+sys, recoJetF, False, sys ) #not a dummy fill so dummy=False
                #if (self.controlPlotsOnly):
                self.fillBranches( event, 'selRecoLeadingJets'+sys, recoJetLeading, False, sys ) #not a dummy fill so dummy=False
                self.fillBranches( event, 'selRecoSubleadingJets'+sys, recoJetSubleading, False, sys ) #not a dummy fill so dummy=False

            
                if self.isMC:
                    deltaRmatch=False
                    if not passGenSel:  ##### fakes
                        self.fakes = self.fakes+1
                        self.eventCategory = 2
                        #if not(self.controlPlotsOnly):
                        self.fillBranches( event, 'trueRecoJets'+sys, recoJet, True, sys ) #True=dummyFill
                        self.fillBranches( event, 'trueRecoJetsF'+sys, recoJetF, True, sys ) #True=dummyFill
                        self.evtGenWeight=0.

                    else:        ##### go to matrix
                        
                        tmpGenJets={}
                        tmpGenJets0={}
                        tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                        tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )

                        #if self.controlPlotsOnly:
                        genJetLeading['Jet'] = tmpGenJets[0]
                        genJetSubleading['Jet'] = tmpGenJets[1]
                        #self.ngenEvents+=1

                        if abs(tmpGenJets[0]['jet'].rapidity) > abs(tmpGenJets[1]['jet'].rapidity):

                            tmpGenJets0[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True, False, '', False, False)  if self.constCheckPlots else tmpGenJets[0] # add PUPPI weights, and, if not(self.controlPlotsOnly) nsub bases as well
                            tmpGenJets0[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True, False, '', True, False)  if self.constCheckPlots else tmpGenJets[1] 

                            genJet['Jet'] = tmpGenJets0[1]#self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                            genJetF['Jet'] = tmpGenJets0[0]
                        else: 
                            tmpGenJets0[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True, False, '', True, False)  if self.constCheckPlots else tmpGenJets[0] # add PUPPI weights, and, if not(self.controlPlotsOnly) nsub bases as well
                            tmpGenJets0[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True, False, '', False, False)  if self.constCheckPlots else tmpGenJets[1] 

                            genJet['Jet'] = tmpGenJets0[0]#self.createNsubBasis( selGenJets[0], event, 'GenCands', True  )
                            genJetF['Jet'] = tmpGenJets0[1]


                        if sys.startswith('_nom'): 
                            self.fillBranches( event, 'selGenJets'+sys, genJet, False, sys ) 
                            self.fillBranches( event, 'selGenJetsF'+sys, genJetF, False, sys ) 
                            #if self.controlPlotsOnly:
                            self.fillBranches( event, 'selGenLeadingJets'+sys, genJetLeading, False, sys ) 
                            self.fillBranches( event, 'selGenSubleadingJets'+sys, genJetSubleading, False, sys ) 
                            

                        # small redundancy here since reco/genJet['Jet'] is already one of the two tmpReco/GenJets:
                        # basically, checking that the measurement jets as stored in the recoJet and genJet objects are in whack with one another, 
                        # while also checking that both of the pT ordered tmp jets in gen/reco are deltaR matched, since if not it's not safe to continue
                        # this block is to ensure that the measurement jet is pristinely selected in an event where the entire dijet system is consistent b/w gen- and reco-level 
                        # doing all this so that it doesn't need to be rechecked later via storing the subleading jet too (too much computational overhead for QCD samples)
                        if ( self.DrRapPhi( recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() ) < 0.4 ) and ( self.DrRapPhi( recoJetF['Jet']['jet'].p4(), genJetF['Jet']['jet'].p4()) < 0.4 ):#self.DrRapPhi( tmpRecoJets[sys][0]['jet'].p4(), tmpGenJets[0]['jet'].p4() ) < 0.4 and self.DrRapPhi( tmpRecoJets[sys][1]['jet'].p4(), tmpGenJets[1]['jet'].p4() ) < 0.4):
                            deltaRmatch = True
                            self.response= self.response+1
                            self.eventCategory = 4
                            #if not(self.controlPlotsOnly):

                        #if ( self.DrRapPhi( recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() ) < 0.4 ):
                        #    #self.nAccepgenEvents+=1
                        #    #self.nTruerecoEvents+=1    

                        if deltaRmatch:
                            self.fillBranches( event, 'accepGenJets'+sys, genJet, False, sys )    
                            self.fillBranches( event, 'accepGenJetsF'+sys, genJetF, False, sys )    
                            self.fillBranches( event, 'trueRecoJets'+sys, recoJet, False, sys )
                            self.fillBranches( event, 'trueRecoJetsF'+sys, recoJetF, False, sys )
                            #self.nCombAccepgenEvents+=1
                            #self.nCombTruerecoEvents+=1

                        
                        
                        #if not deltaR matched, fill branches with dummies to posteriorly obtain misreconstructed gen and fake reco jets
                        else:
                            self.fillBranches( event, 'accepGenJets'+sys, genJet, True, sys )    #fill dummy=True to indicate fakes and misses 
                            self.fillBranches( event, 'accepGenJetsF'+sys, genJetF, True, sys )    #fill dummy=True to indicate fakes and misses 
                            self.fillBranches( event, 'trueRecoJets'+sys, recoJet, True, sys )
                            self.fillBranches( event, 'trueRecoJetsF'+sys, recoJetF, True, sys )
                            
                            self.fakes = self.fakes+1
                            self.eventCategory = 2
                            self.miss = self.miss+1
                            self.eventCategory = 3



            if not self.isMC:
                self.out.fillBranch( 'eventCategory'+sys, self.eventCategory )
                self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'passRecoSel'+sys, 1 if passRecoSel[sys] else 0)
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood'))# if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'pt_asymm'+sys, reco_ptAsymm[sys] if passRecoSel[sys] else -929.)
                self.out.fillBranch( 'delta_R'+sys, reco_deltaR[sys] if passRecoSel[sys] else -929.)
                self.out.fillBranch( 'delta_phi'+sys, reco_deltaPhi[sys] if passRecoSel[sys] else 929.)
                self.out.fillBranch( 'nRecoLeptons'+sys, reco_nleptons)
                

                for x in self.triggerTable.keys():
                    self.out.fillBranch( 'passHLT_'+x, 1 if getattr(event, 'HLT_'+x)==1 and passRecoSel[sys] else 0)

                if not('2016' in self.year): self.out.fillBranch( 'passHLT_AK8PFJet550', 1 if getattr(event, 'HLT_AK8PFJet550')==1 and passRecoSel[sys] else 0)
            
            #if self.isMC: #and not sys.startswith(('_jes', '_jer')):
            else:
                self.out.fillBranch( 'eventCategory'+sys, self.eventCategory )
                self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'passRecoSel'+sys, 1 if passRecoSel[sys] else 0)
                self.out.fillBranch( 'puWeightNom'+sys, self.puWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'l1prefiringWeightNom'+sys, self.l1PreFireWeight if passRecoSel[sys] else 0.)

                self.out.fillBranch( 'pt_asymm'+sys, reco_ptAsymm[sys] if passRecoSel[sys] else -929.)
                self.out.fillBranch( 'delta_R'+sys, reco_deltaR[sys] if passRecoSel[sys] else -929.)
                self.out.fillBranch( 'delta_phi'+sys, reco_deltaPhi[sys] if passRecoSel[sys] else 929.)

                if sys.endswith('nom'):
                    self.out.fillBranch( 'gen_pt_asymm'+sys, gen_ptAsymm if passGenSel else -929.)
                    self.out.fillBranch( 'gen_delta_R'+sys, gen_deltaR if passGenSel else -929.)
                    self.out.fillBranch( 'gen_delta_phi'+sys, gen_deltaPhi if passGenSel else 929.)
                    self.out.fillBranch( 'nGenLeptons'+sys, gen_nleptons)
                    
                    self.out.fillBranch( 'nRecoLeptons'+sys, reco_nleptons)
                    
                    self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') )#if passRecoSel[sys] else 0.)

                #if sys.startswith('_nom'): 
                self.out.fillBranch( 'evtGenWeight'+sys, self.evtGenWeight if passGenSel or passRecoSel[sys] else 0.) 
                self.out.fillBranch( 'passGenSel'+sys, 1 if passGenSel else 0)#1 if passGenSel or passRecoSel[sys] else 0 
                    
                if self.isSigMC and not self.onlyUnc:
                    self.out.fillBranch( 'pdfWeightNom'+sys, self.pdfWeight if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'pdfWeightAll'+sys, self.pdfWeightAll if passGenSel or passRecoSel[sys] else np.zeros((103,),dtype=np.float32))
                    self.out.fillBranch( 'pdfWeightUp'+sys, self.pdfWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'pdfWeightDown'+sys, self.pdfWeightDown if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'isrWeightUp'+sys, self.isrWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'isrWeightDown'+sys, self.isrWeightDown if passGenSel or passRecoSel[sys]  else 0.)
                    self.out.fillBranch( 'fsrWeightUp'+sys, self.fsrWeightUp if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'fsrWeightDown'+sys, self.fsrWeightDown if passGenSel or passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'puWeightUp'+sys, self.puWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'puWeightDown'+sys, self.puWeightDown if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'l1prefiringWeightUp'+sys, self.l1PreFireWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'l1prefiringWeightDown'+sys, self.l1PreFireWeightDown if passRecoSel[sys] else 0.)

        #if self.evtCounter%100==0:
        #    print("evtCounter","nrecoEvents","ngenEvents","nCombAccepgenEvents","nCombTruerecoEvents","nAccepgenEvents","nTruerecoEvents")
        #    print(self.evtCounter,self.nrecoEvents,self.ngenEvents,self.nCombAccepgenEvents,self.nCombTruerecoEvents,self.nAccepgenEvents,self.nTruerecoEvents)

        return True


    #############################################################################
    def recoSelection( self, event, sysUnc=[] ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))

        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and abs(x.p4().Eta()) < self.maxElectronEta]#and  #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minLooseMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.mediumId and x.miniPFRelIso_all<0.2 ] #and x.pfIsoId>=2 and abs(x.dxy)<0.2 and abs(x.dz)<0.5

        #recoCrackElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))]#and abs(x.p4().Eta()) < self.maxElectronEta ]]
        #nCrackEl = len(recoCrackElectrons)

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons) 

        ##################################################
        if not len(AK8jets)==0:
            for ijets in AK8jets: 
                ijets.rapidity = ijets.p4().Rapidity()
        recoAK8jets = {}
        passSel = {}
        ptAsymm = {}
        deltaPhi = {}
        deltaR = {}
        iSel = {}
        for sys in self.sysSource:
            #if sys.startswith(self.sysWeightList): sys = '_nom'
            #### Basic AK8 jet selection
            #JetID (https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD): For UL16-18 samples jetId==2 means: pass tight ID, fail tightLepVeto, jetId==6 means: pass tight and tightLepVeto ID, so we keep the >2 (ie, 6) 
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs( x.rapidity ) < self.maxAK8JetRap and (x.jetId > 2)] #self.maxJetAK8Eta
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ),reverse=True)

            if len(recoAK8jets[sys])!=0:

                for jet in recoAK8jets[sys]:
                    jecFactor = getattr( jet, 'pt'+sys ) / getattr( jet, 'pt_raw' ) 
                    setattr(jet, 'mSD'+sys,  ( getattr( jet, 'msoftdrop_raw' ) * jecFactor )  )


            ##################################################

            ##### Applying selection
            #print ("in reco sel", event.event, recoMuons, recoElectrons, recoAK8jets[sys], sys)
            #print (self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets[sys], sys ))
            passSel[sys], ptAsymm[sys], deltaPhi[sys], deltaR[sys] = self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets[sys], sys )
            iSel[sys] = '_dijetSel' if passSel[sys] else None
        
        if not self.isMC:
            if passSel['_nom']:
               
                for itrigger, itrValues in self.triggerTable.items():
                    if ( getattr(event, 'HLT_'+itrigger)==1 ):

                        getattr( self, 'nPVs_only'+itrigger+'_dijetSel' ).Fill( getattr( event, 'PV_npvsGood') )
                        #if len(recoAK8jets['_nom'])>0:
                        getattr( self, 'LeadPtJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].pt_nom )   ### pt_nom here to ensure data process
                        getattr( self, 'LeadPtJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].eta )
                        getattr( self, 'LeadPtJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].rapidity )
                        getattr( self, 'LeadPtJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].phi )
                        getattr( self, 'LeadPtJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                        getattr( self, 'LeadPtJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( getattr(recoAK8jets['_nom'][0], 'mSD_nom') )
                        
                        getattr( self, 'SubleadPtJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].pt_nom )   ### pt_nom here to ensure data process
                        getattr( self, 'SubleadPtJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].eta )
                        getattr( self, 'SubleadPtJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].rapidity )
                        getattr( self, 'SubleadPtJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].phi )
                        getattr( self, 'SubleadPtJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                        getattr( self, 'SubleadPtJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( getattr(recoAK8jets['_nom'][1], 'mSD_nom') )
                        
                        
                        if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                            c=0
                            f=1
                        else:
                            c=1
                            f=0

                        getattr( self, 'CentralJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][c].pt_nom )   ### pt_nom here to ensure data process
                        getattr( self, 'CentralJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][c].eta )
                        getattr( self, 'CentralJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][c].rapidity )
                        getattr( self, 'CentralJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][c].phi )
                        getattr( self, 'CentralJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][c].mass_nom )
                        getattr( self, 'CentralJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( getattr(recoAK8jets['_nom'][c], 'mSD_nom') )

                        getattr( self, 'ForwardJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][f].pt_nom )   ### pt_nom here to ensure data process
                        getattr( self, 'ForwardJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][f].eta )
                        getattr( self, 'ForwardJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][f].rapidity )
                        getattr( self, 'ForwardJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][f].phi )
                        getattr( self, 'ForwardJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][f].mass_nom )
                        getattr( self, 'ForwardJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( getattr(recoAK8jets['_nom'][f], 'mSD_nom') )
                    
                        if self.perTriggerPassSel:
                            if ( itrValues[self.year][0] < recoAK8jets['_nom'][0].pt_nom < itrValues[self.year][1] ) :
                                iSel['_nom'] = '_'+itrigger+'_dijetSel'
                                passSel['_nom'] = True
                                break
                            else:
                                passSel['_nom'] = False
                                iSel['_nom'] = None
                                #self.triggerWeight = 0
            else:
                passSel['_nom'] = False
                iSel['_nom'] = None

        #### Weight #########
        weight=1.

        if self.isMC:
            try: 
                l1PreFireWeight = event.L1PreFiringWeight_Nom                
            except RuntimeError:  # to deal with 2018 QCD_Pt samples not having the l1 weight branches for some reason
                l1PreFireWeight = 1.

            self.puWeight = event.puWeight
            self.l1PreFireWeight = l1PreFireWeight
            self.evtGenWeight = event.genWeight 

            weight = self.puWeight * self.l1PreFireWeight * self.evtGenWeight
            
        else:
            weight = 1 #self.triggerWeight
        #self.totalRecoWeight = weight

        ################################################## Add mass/mSD and rapidity histos (then factorise ttbar skimmer and this as well)
        '''
        if not self.onlyUnc:
            #### Checking no selection without weights
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
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt_nom )
                getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK8jets_y_noSelnoWeight' ).Fill( ijet.rapidity )
                getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.mass_nom )
                getattr( self, 'AK8jets_msoftdrop_noSelnoWeight' ).Fill( ijet.msoftdrop )
            if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                getattr( self, 'LeadingPtAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                getattr( self, 'LeadingPtAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                getattr( self, 'LeadingPtAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                getattr( self, 'LeadingPtAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                getattr( self, 'LeadingPtAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop )
                if  len(recoAK8jets['_nom'])>1:
                    getattr( self, 'SubleadingPtAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                    getattr( self, 'SubleadingPtAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                    getattr( self, 'SubleadingPtAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                    getattr( self, 'SubleadingPtAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                    getattr( self, 'SubleadingPtAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                    getattr( self, 'SubleadingPtAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop )

                    if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                        getattr( self, 'CentralAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                        getattr( self, 'CentralAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                        getattr( self, 'CentralAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                        getattr( self, 'CentralAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                        getattr( self, 'CentralAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop )

                        getattr( self, 'ForwardAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                        getattr( self, 'ForwardAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                        getattr( self, 'ForwardAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                        getattr( self, 'ForwardAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                        getattr( self, 'ForwardAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop )
                    else: 
                        getattr( self, 'CentralAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                        getattr( self, 'CentralAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                        getattr( self, 'CentralAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                        getattr( self, 'CentralAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                        getattr( self, 'CentralAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop )

                        getattr( self, 'ForwardAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                        getattr( self, 'ForwardAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                        getattr( self, 'ForwardAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                        getattr( self, 'ForwardAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                        getattr( self, 'ForwardAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop )

            #### Checking no selection with weights
            getattr( self, 'nPVs_noSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
            getattr( self, 'nleps_noSel' ).Fill( nleptons, weight )
            for imuon in recoMuons:
                getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt, weight )
                getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta, weight )
                getattr( self, 'muons_y_noSel' ).Fill( imuon.p4().Rapidity(), weight )
                getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi, weight )
            for iele in recoElectrons:
                getattr( self, 'eles_pt_noSel' ).Fill( iele.pt, weight )
                getattr( self, 'eles_eta_noSel' ).Fill( iele.eta, weight )
                getattr( self, 'eles_y_noSel' ).Fill( iele.p4().Rapidity(), weight )
                getattr( self, 'eles_phi_noSel' ).Fill( iele.phi, weight )
            getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets['_nom']), weight )
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt_nom, weight )
                getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8jets_y_noSel' ).Fill( ijet.rapidity, weight )
                getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.mass_nom, weight )
                getattr( self, 'AK8jets_msoftdrop_noSel' ).Fill( ijet.msoftdrop, weight )

            if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                getattr( self, 'LeadingPtAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                getattr( self, 'LeadingPtAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                getattr( self, 'LeadingPtAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                getattr( self, 'LeadingPtAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                getattr( self, 'LeadingPtAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )

                if  len(recoAK8jets['_nom'])>1:
                    getattr( self, 'SubleadingPtAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'SubleadingPtAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'SubleadingPtAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity , weight )
                    getattr( self, 'SubleadingPtAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi , weight )
                    getattr( self, 'SubleadingPtAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'SubleadingPtAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )

                    if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                        getattr( self, 'CentralAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                        getattr( self, 'CentralAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                        getattr( self, 'CentralAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                        getattr( self, 'CentralAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                        getattr( self, 'CentralAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )

                        getattr( self, 'ForwardAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                        getattr( self, 'ForwardAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                        getattr( self, 'ForwardAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                        getattr( self, 'ForwardAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi, weight )
                        getattr( self, 'ForwardAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )
                    else: 
                        getattr( self, 'CentralAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                        getattr( self, 'CentralAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                        getattr( self, 'CentralAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                        getattr( self, 'CentralAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi, weight )
                        getattr( self, 'CentralAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )

                        getattr( self, 'ForwardAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                        getattr( self, 'ForwardAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                        getattr( self, 'ForwardAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                        getattr( self, 'ForwardAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                        getattr( self, 'ForwardAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )
                        
            getattr( self, 'recoPtAsymm_noSel' ).Fill( ptAsymm['_nom'], weight )
            getattr( self, 'recoDeltaPhi_noSel' ).Fill( deltaPhi['_nom'], weight )

            #reweight = self.totalRecoWeight
            #### Creating Nsub basis, filling histos and creating branches IF passSel
            if ('_nom' in iSel.keys()) and passSel['_nom'] and iSel['_nom']:

                #### Basic reco histos
                getattr( self, 'nPVs'+iSel['_nom'] ).Fill( getattr( event, 'PV_npvsGood'), weight )
                getattr( self, 'nleps'+iSel['_nom'] ).Fill( len(recoMuons)+len(recoElectrons), weight )
                for imuon in recoMuons:
                    getattr( self, 'muons_pt'+iSel['_nom'] ).Fill( imuon.pt, weight )
                    getattr( self, 'muons_eta'+iSel['_nom'] ).Fill( imuon.eta, weight )
                    getattr( self, 'muons_y'+iSel['_nom'] ).Fill( imuon.p4().Rapidity(), weight )
                    getattr( self, 'muons_phi'+iSel['_nom'] ).Fill( imuon.phi, weight )
                for iele in recoElectrons:
                    getattr( self, 'eles_pt'+iSel['_nom'] ).Fill( iele.pt, weight )
                    getattr( self, 'eles_eta'+iSel['_nom'] ).Fill( iele.eta, weight )
                    getattr( self, 'eles_y'+iSel['_nom'] ).Fill( iele.p4().Rapidity(), weight )
                    getattr( self, 'eles_phi'+iSel['_nom'] ).Fill( iele.phi, weight )
                getattr( self, 'nAK8jets'+iSel['_nom'] ).Fill( len(recoAK8jets['_nom']), weight )
                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt'+iSel['_nom'] ).Fill( ijet.pt_nom, weight )
                    getattr( self, 'AK8jets_eta'+iSel['_nom'] ).Fill( ijet.eta, weight )
                    getattr( self, 'AK8jets_y'+iSel['_nom'] ).Fill( ijet.rapidity, weight )
                    getattr( self, 'AK8jets_phi'+iSel['_nom'] ).Fill( ijet.phi, weight )
                    getattr( self, 'AK8jets_mass'+iSel['_nom'] ).Fill( ijet.mass_nom, weight )
                    getattr( self, 'AK8jets_msoftdrop'+iSel['_nom'] ).Fill( ijet.msoftdrop, weight )
                #if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                getattr( self, 'LeadingPtAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                getattr( self, 'LeadingPtAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                getattr( self, 'LeadingPtAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                getattr( self, 'LeadingPtAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                getattr( self, 'LeadingPtAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )

                getattr( self, 'SubleadingPtAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                getattr( self, 'SubleadingPtAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                getattr( self, 'SubleadingPtAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity , weight )
                getattr( self, 'SubleadingPtAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi , weight )
                getattr( self, 'SubleadingPtAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                getattr( self, 'SubleadingPtAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )

                if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                    getattr( self, 'CentralAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                    getattr( self, 'CentralAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                    getattr( self, 'CentralAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                    getattr( self, 'CentralAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                    getattr( self, 'CentralAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                    getattr( self, 'CentralAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )

                    getattr( self, 'ForwardAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'ForwardAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'ForwardAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                    getattr( self, 'ForwardAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi, weight )
                    getattr( self, 'ForwardAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'ForwardAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )
                else: 
                    getattr( self, 'CentralAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'CentralAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'CentralAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                    getattr( self, 'CentralAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi, weight )
                    getattr( self, 'CentralAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'CentralAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop, weight )

                    getattr( self, 'ForwardAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                    getattr( self, 'ForwardAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                    getattr( self, 'ForwardAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                    getattr( self, 'ForwardAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                    getattr( self, 'ForwardAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                    getattr( self, 'ForwardAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop, weight )

                getattr( self, 'recoPtAsymm'+iSel['_nom'] ).Fill( ptAsymm['_nom'], weight )
                getattr( self, 'recoDeltaPhi'+iSel['_nom'] ).Fill( deltaPhi['_nom'], weight )
        '''
        return passSel, iSel, recoAK8jets, deltaR, ptAsymm, deltaPhi, nleptons#, nCrackEl, recoMuons,recoElectrons,

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing reco information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genParticles = Collection(event, 'GenPart')

        ### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt > self.minLooseElectronPt and abs(x.eta)<self.maxElectronEta ]#and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))]#
        #genCrackElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt > self.minLooseElectronPt and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))]#and abs(x.p4().Eta()) < self.maxElectronEta ]]
        #nCrackEl = len(genCrackElectrons)

        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minLooseMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        #print ("in gensel", event.event, genMuons, genElectrons, genJetsAK8, '' )
        nleptons = len(genMuons)+len(genElectrons) 

        ##################################################

        #### Basic AK8 jet selection
        if not len(genJetsAK8)==0:
            for ijets in genJetsAK8: 
                ijets.rapidity = ijets.p4().Rapidity()
                ijets.mSD = self.getGenJetAK8softdropmass( AK8jet=ijets, event=event, PFCollection='GenCands', isGen=True)
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.rapidity) < self.maxAK8JetRap]#maxJetAK8Eta 
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        #print ("in gensel", event.event, genMuons, genElectrons, genAK8jets, '' )
        passSel, ptAsymm, deltaPhi, deltaR = self.dijetSelection( event, genMuons, genElectrons, genAK8jets, '' )
        iSel = '_dijetSel' if passSel else None
        #### Weight
        weight = event.genWeight
        '''
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
            for ijet in genAK8jets:
                getattr( self, 'AK8genjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK8genjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8genjets_y_noSel' ).Fill( ijet.rapidity, weight )
                getattr( self, 'AK8genjets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8genjets_mass_noSel' ).Fill( ijet.mass, weight )
                #if not self.FlagBadPFCands: 
                getattr( self, 'AK8genjets_msoftdrop_noSel' ).Fill( ijet.msoftdrop, weight )

            if not len(genAK8jets)==0:
                getattr( self, 'LeadingPtAK8genjet_pt_noSel' ).Fill( genAK8jets[0].pt, weight )
                getattr( self, 'LeadingPtAK8genjet_eta_noSel' ).Fill( genAK8jets[0].eta, weight )
                getattr( self, 'LeadingPtAK8genjet_y_noSel' ).Fill( genAK8jets[0].rapidity, weight )
                getattr( self, 'LeadingPtAK8genjet_phi_noSel' ).Fill( genAK8jets[0].phi, weight )
                getattr( self, 'LeadingPtAK8genjet_mass_noSel' ).Fill( genAK8jets[0].mass, weight )
                #if not self.FlagBadPFCands: 
                getattr( self, 'LeadingPtAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[0].msoftdrop, weight )
                
                if len(genAK8jets)>1:
                    getattr( self, 'SubleadingPtAK8genjet_pt_noSel' ).Fill( genAK8jets[1].pt, weight )
                    getattr( self, 'SubleadingPtAK8genjet_eta_noSel' ).Fill( genAK8jets[1].eta, weight )
                    getattr( self, 'SubleadingPtAK8genjet_y_noSel' ).Fill( genAK8jets[1].rapidity , weight )
                    getattr( self, 'SubleadingPtAK8genjet_phi_noSel' ).Fill( genAK8jets[1].phi , weight )
                    getattr( self, 'SubleadingPtAK8genjet_mass_noSel' ).Fill( genAK8jets[1].mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'SubleadingPtAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[1].msoftdrop, weight )

                    if abs(genAK8jets[0].rapidity)<abs(genAK8jets[1].rapidity):
                        getattr( self, 'CentralAK8genjet_pt_noSel' ).Fill( genAK8jets[0].pt, weight )
                        getattr( self, 'CentralAK8genjet_eta_noSel' ).Fill( genAK8jets[0].eta, weight )
                        getattr( self, 'CentralAK8genjet_y_noSel' ).Fill( genAK8jets[0].rapidity, weight )
                        getattr( self, 'CentralAK8genjet_phi_noSel' ).Fill( genAK8jets[0].phi, weight )
                        getattr( self, 'CentralAK8genjet_mass_noSel' ).Fill( genAK8jets[0].mass, weight )
                        #if not self.FlagBadPFCands: 
                        getattr( self, 'CentralAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[0].msoftdrop, weight )

                        getattr( self, 'ForwardAK8genjet_pt_noSel' ).Fill( genAK8jets[1].pt, weight )
                        getattr( self, 'ForwardAK8genjet_eta_noSel' ).Fill( genAK8jets[1].eta, weight )
                        getattr( self, 'ForwardAK8genjet_y_noSel' ).Fill( genAK8jets[1].rapidity, weight )
                        getattr( self, 'ForwardAK8genjet_phi_noSel' ).Fill( genAK8jets[1].phi, weight )
                        getattr( self, 'ForwardAK8genjet_mass_noSel' ).Fill( genAK8jets[1].mass, weight )
                        getattr( self, 'ForwardAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[1].msoftdrop, weight )
                    else: 
                        getattr( self, 'CentralAK8genjet_pt_noSel' ).Fill( genAK8jets[1].pt, weight )
                        getattr( self, 'CentralAK8genjet_eta_noSel' ).Fill( genAK8jets[1].eta, weight )
                        getattr( self, 'CentralAK8genjet_y_noSel' ).Fill( genAK8jets[1].rapidity, weight )
                        getattr( self, 'CentralAK8genjet_phi_noSel' ).Fill( genAK8jets[1].phi, weight )
                        getattr( self, 'CentralAK8genjet_mass_noSel' ).Fill( genAK8jets[1].mass, weight )
                        #if not self.FlagBadPFCands: 
                        getattr( self, 'CentralAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[1].msoftdrop, weight )

                        getattr( self, 'ForwardAK8genjet_pt_noSel' ).Fill( genAK8jets[0].pt, weight )
                        getattr( self, 'ForwardAK8genjet_eta_noSel' ).Fill( genAK8jets[0].eta, weight )
                        getattr( self, 'ForwardAK8genjet_y_noSel' ).Fill( genAK8jets[0].rapidity, weight )
                        getattr( self, 'ForwardAK8genjet_phi_noSel' ).Fill( genAK8jets[0].phi, weight )
                        getattr( self, 'ForwardAK8genjet_mass_noSel' ).Fill( genAK8jets[0].mass, weight )
                        #if not self.FlagBadPFCands: 
                        getattr( self, 'ForwardAK8genjet_msoftdrop_noSel' ).Fill( genAK8jets[0].msoftdrop, weight )

            getattr( self, 'genPtAsymm_noSel' ).Fill( ptAsymm, weight )
            getattr( self, 'genDeltaPhi_noSel' ).Fill( deltaPhi, weight )
            if len(genAK8jets)>0: getattr(self, 'genJetAK8_partonFlavour_noSel').Fill(genAK8jets[0].partonFlavour, weight)


            ##### Filling histograms
            if passSel and iSel:
                #print (genAK8jets[0].partonFlavour)
                #### Checking nominal selection with weights
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
                for ijet in genAK8jets:
                    getattr( self, 'AK8genjets_pt'+iSel ).Fill( ijet.pt, weight )
                    getattr( self, 'AK8genjets_eta'+iSel ).Fill( ijet.eta, weight )
                    getattr( self, 'AK8genjets_y'+iSel ).Fill( ijet.rapidity, weight )
                    getattr( self, 'AK8genjets_phi'+iSel ).Fill( ijet.phi, weight )
                    getattr( self, 'AK8genjets_mass'+iSel ).Fill( ijet.mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'AK8genjets_msoftdrop'+iSel ).Fill( ijet.msoftdrop, weight )
                
                #if not len(genAK8jets)==0:
                getattr( self, 'LeadingPtAK8genjet_pt'+iSel ).Fill( genAK8jets[0].pt, weight )
                getattr( self, 'LeadingPtAK8genjet_eta'+iSel ).Fill( genAK8jets[0].eta, weight )
                getattr( self, 'LeadingPtAK8genjet_y'+iSel ).Fill( genAK8jets[0].rapidity, weight )
                getattr( self, 'LeadingPtAK8genjet_phi'+iSel ).Fill( genAK8jets[0].phi, weight )
                getattr( self, 'LeadingPtAK8genjet_mass'+iSel ).Fill( genAK8jets[0].mass, weight )
                #if not self.FlagBadPFCands: 
                getattr( self, 'LeadingPtAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[0].msoftdrop, weight )

                getattr( self, 'SubleadingPtAK8genjet_pt'+iSel ).Fill( genAK8jets[1].pt, weight )
                getattr( self, 'SubleadingPtAK8genjet_eta'+iSel ).Fill( genAK8jets[1].eta, weight )
                getattr( self, 'SubleadingPtAK8genjet_y'+iSel ).Fill( genAK8jets[1].rapidity , weight )
                getattr( self, 'SubleadingPtAK8genjet_phi'+iSel ).Fill( genAK8jets[1].phi , weight )
                getattr( self, 'SubleadingPtAK8genjet_mass'+iSel ).Fill( genAK8jets[1].mass, weight )
                #if not self.FlagBadPFCands: 
                getattr( self, 'SubleadingPtAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[1].msoftdrop, weight )

                if abs(genAK8jets[0].rapidity)<abs(genAK8jets[1].rapidity):
                    getattr( self, 'CentralAK8genjet_pt'+iSel ).Fill( genAK8jets[0].pt, weight )
                    getattr( self, 'CentralAK8genjet_eta'+iSel ).Fill( genAK8jets[0].eta, weight )
                    getattr( self, 'CentralAK8genjet_y'+iSel ).Fill( genAK8jets[0].rapidity, weight )
                    getattr( self, 'CentralAK8genjet_phi'+iSel ).Fill( genAK8jets[0].phi, weight )
                    getattr( self, 'CentralAK8genjet_mass'+iSel ).Fill( genAK8jets[0].mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'CentralAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[0].msoftdrop, weight )

                    getattr( self, 'ForwardAK8genjet_pt'+iSel ).Fill( genAK8jets[1].pt, weight )
                    getattr( self, 'ForwardAK8genjet_eta'+iSel ).Fill( genAK8jets[1].eta, weight )
                    getattr( self, 'ForwardAK8genjet_y'+iSel ).Fill( genAK8jets[1].rapidity, weight )
                    getattr( self, 'ForwardAK8genjet_phi'+iSel ).Fill( genAK8jets[1].phi, weight )
                    getattr( self, 'ForwardAK8genjet_mass'+iSel ).Fill( genAK8jets[1].mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'ForwardAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[1].msoftdrop, weight )
                else: 
                    getattr( self, 'CentralAK8genjet_pt'+iSel ).Fill( genAK8jets[1].pt, weight )
                    getattr( self, 'CentralAK8genjet_eta'+iSel ).Fill( genAK8jets[1].eta, weight )
                    getattr( self, 'CentralAK8genjet_y'+iSel ).Fill( genAK8jets[1].rapidity, weight )
                    getattr( self, 'CentralAK8genjet_phi'+iSel ).Fill( genAK8jets[1].phi, weight )
                    getattr( self, 'CentralAK8genjet_mass'+iSel ).Fill( genAK8jets[1].mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'CentralAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[1].msoftdrop, weight )

                    getattr( self, 'ForwardAK8genjet_pt'+iSel ).Fill( genAK8jets[0].pt, weight )
                    getattr( self, 'ForwardAK8genjet_eta'+iSel ).Fill( genAK8jets[0].eta, weight )
                    getattr( self, 'ForwardAK8genjet_y'+iSel ).Fill( genAK8jets[0].rapidity, weight )
                    getattr( self, 'ForwardAK8genjet_phi'+iSel ).Fill( genAK8jets[0].phi, weight )
                    getattr( self, 'ForwardAK8genjet_mass'+iSel ).Fill( genAK8jets[0].mass, weight )
                    #if not self.FlagBadPFCands: 
                    getattr( self, 'ForwardAK8genjet_msoftdrop'+iSel ).Fill( genAK8jets[0].msoftdrop, weight )

                getattr( self, 'genPtAsymm'+iSel ).Fill( ptAsymm, weight )
                getattr( self, 'genDeltaPhi'+iSel ).Fill( deltaPhi, weight )
                getattr(self, 'genJetAK8_partonFlavour'+iSel).Fill(genAK8jets[0].partonFlavour, weight)
        '''
        return passSel, iSel, genAK8jets, deltaR, ptAsymm, deltaPhi, nleptons#, nCrackEl, genMuons, genElectrons,

    #############################################################################
    def dijetSelection( self, event, muons, electrons, AK8jets, ptLabel ):

        #print ("in main dijetSel", event.event, muons, electrons, AK8jets, ptLabel)

        if self.leptonVeto and (len(muons)+len(electrons))!=0:
            return False, 999, -999, -999

        else:
            if len(AK8jets)>1:# and abs(AK8jets[0].eta)<self.maxLeadAK8JetEta):#(len(muons)+len(electrons))==0 and
                
                jet1Pt = getattr( AK8jets[0], 'pt'+ptLabel )
                jet2Pt = getattr( AK8jets[1], 'pt'+ptLabel )
                #print ("in main dijetSel", event.event, muons, electrons, AK8jets, ptLabel, jet1Pt, jet2Pt)

                tmpJet1 = ROOT.TLorentzVector()
                tmpJet1.SetPtEtaPhiM( getattr( AK8jets[0], 'pt'+ptLabel ), AK8jets[0].eta, AK8jets[0].phi, getattr( AK8jets[0], 'mass'+ptLabel )  )
                tmpJet2 = ROOT.TLorentzVector()
                tmpJet2.SetPtEtaPhiM( getattr( AK8jets[1], 'pt'+ptLabel ), AK8jets[1].eta, AK8jets[1].phi, getattr( AK8jets[1], 'mass'+ptLabel )  )
                
                deltaR_flag = True
                # dR_val = 929.
                #print("DijetSel Loop flag",len(AK8jets))
                for i in range(1,len(AK8jets)):
                    tmpJet = ROOT.TLorentzVector()
                    tmpJet.SetPtEtaPhiM( getattr( AK8jets[i], 'pt'+ptLabel ), AK8jets[i].eta, AK8jets[i].phi, getattr( AK8jets[i], 'mass'+ptLabel )  )
                    #cleaning cut
                    if self.DrRapPhi(tmpJet1,tmpJet)<1.6:# or abs(tmpJet1.DeltaPhi(tmpJet)<2.): #2*jet radius separation + deltaphi sep. between leading AK8 jet and any other fatjets in the event
                        deltaR_flag=False
                        break

                #print("Post-dijetSel Loop flag",len(AK8jets))

                if not deltaR_flag: 
                    return False, 999, -999, -999

                else:
                    #do calcs only when necessary
                    ptAsymm = ( jet1Pt - jet2Pt ) / (jet1Pt + jet2Pt)
                    deltaPhi = tmpJet1.DeltaPhi( tmpJet2 )
                    deltaR = self.DrRapPhi(tmpJet1,tmpJet2)
                    if (jet1Pt>self.minAK8JetPt) and (abs(deltaPhi)>2.):# and (ptAsymm<0.3)  and
                        return True, ptAsymm, deltaPhi, deltaR
                    else: 
                        return False, 999, -999, -999
            
            else: 
                return False, 999, -999, -999

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
            #if p.p4().M()<0.:# and not(np.isclose(p.p4().M(),0.,rtol=10**(-9),atol=10**(-9))): 
            #    self.FlagBadPFCands=True
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
                ak8jet['mSD'] = sd_AK8jets[0].m()
            else: ak8jet['mSD'] = -1. 
        return ak8jet['mSD']


    #############################################################################

    def createNsubBasis(self, AK8jet, event, PFCollection, isGen=False, constJES=False, varUpDown='', central=False, noPlot=True ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        #print(AK8jet,isGen,constJES,central,noPlot)

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

        
        
            if central:
                getattr(self, '%s_dR_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
                getattr(self, '%s_dR_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)        
                if dRSJ1_AK8<0.6: getattr(self, '%s_dR2_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
                if dRSJ2_AK8<0.6: getattr(self, '%s_dR2_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)
            else:
                getattr(self, 'F%s_dR_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
                getattr(self, 'F%s_dR_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)        
                if dRSJ1_AK8<0.6: getattr(self, 'F%s_dR2_SJ1'%(type)).Fill(dRSJ1_AK8)#dR if isGen else puppi_dR)
                if dRSJ2_AK8<0.6: getattr(self, 'F%s_dR2_SJ2'%(type)).Fill(dRSJ2_AK8)#dR if isGen else puppi_dR)

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

                        if central==True:

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

                        else:
                            if pID==22:
                                if not(noPlot): 
                                    getattr(self, 'F%s_dR_photons'%(type)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, 'F%s_dRap_dPhi_photons'%(type)).Fill(dRap, dPhi)
                                    getattr(self, 'F%s_dR_photons_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, 'F%s_dRap_dPhi_photons_%s' % (type, bin_name)).Fill(dRap, dPhi)


                                    if not(isGen):
                                        getattr(self, 'F%s_dRSJ1_photons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                        if dRSJ1_AK8<0.6: getattr(self, 'F%s_dR2SJ1_photons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                        getattr(self, 'F%s_dRSJ2_photons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))
                                        if dRSJ2_AK8<0.6: getattr(self, 'F%s_dR2SJ2_photons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))

                                    getattr(self, 'F%s_pT_photons'%(type)).Fill(pT)
                                    getattr(self, 'F%sjetPt_dR_photons'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())
                                    getattr(self, 'F%sjetmSD_dR_photons'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) 
                                    getattr(self, 'F%s_dR_Pt_photons'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())
                                npho+=1

                            elif c==0 and pID!=22:
                                if not(noPlot): 
                                    getattr(self, 'F%s_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, 'F%s_dRap_dPhi_neutral'%(type)).Fill(dRap, dPhi)
                                    getattr(self, 'F%s_dR_neutral_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, 'F%s_dRap_dPhi_neutral_%s' % (type, bin_name)).Fill(dRap, dPhi)

                                    if not(isGen):
                                        getattr(self, 'F%s_dRSJ1_neutral'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                        if dRSJ1_AK8<0.6: getattr(self, 'F%s_dR2SJ1_neutral'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))
                                        getattr(self, 'F%s_dRSJ2_neutral'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))
                                        if dRSJ2_AK8<0.6: getattr(self, 'F%s_dR2SJ2_neutral'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))

                                    getattr(self, 'F%s_pT_neutral'%(type)).Fill(pT)
                                    getattr(self, 'F%sjetPt_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())
                                    getattr(self, 'F%sjetmSD_dR_neutral'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) 
                                    getattr(self, 'F%s_dR_Pt_neutral'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())
                                nneutHad+=1

                            elif abs(c)>0:# and abs(pID)!=11:# and abs(pID)!=13:
                                if not(noPlot): 
                                    getattr(self, 'F%s_dR_charged'%(type)).Fill(dR if isGen else puppi_dR)
                                    getattr(self, 'F%s_dRap_dPhi_charged'%(type)).Fill(dRap, dPhi)
                                    #getattr(self, 'F%s_dR_charged_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                    #getattr(self, 'F%s_dRap_dPhi_charged_%s' % (type, bin_name)).Fill(dRap, dPhi)


                                    if not(isGen):
                                        getattr(self, 'F%s_dRSJ1_charged'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        if dRSJ1_AK8<0.6: getattr(self, 'F%s_dR2SJ1_charged'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        getattr(self, 'F%s_dRSJ2_charged'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies
                                        if dRSJ2_AK8<0.6: getattr(self, 'F%s_dR2SJ2_charged'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies

                                    getattr(self, 'F%s_pT_charged'%(type)).Fill(pT)
                                    getattr(self, 'F%sjetPt_dR_charged'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())#x is dR, y is pT
                                    getattr(self, 'F%sjetmSD_dR_charged'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) #x is dR, y is pT
                                    getattr(self, 'F%s_dR_Pt_charged'%(type)).Fill(pT, dR if isGen else puppi_dR)  
                                nch+=1
                                if abs(pID)!=11 and abs(pID)!=13:
                                    if not(noPlot): 
                                        getattr(self, 'F%s_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR)
                                        getattr(self, 'F%s_dRap_dPhi_chargedHadrons'%(type)).Fill(dRap, dPhi)

                                        getattr(self, 'F%s_dR_chargedHadrons_%s' % (type, bin_name)).Fill(dR if isGen else puppi_dR)
                                        getattr(self, 'F%s_dRap_dPhi_chargedHadrons_%s' % (type, bin_name)).Fill(dRap, dPhi)

                                        if not(isGen):
                                            
                                            getattr(self, 'F%s_dRSJ1_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                            if dRSJ1_AK8<0.6: getattr(self, 'F%s_dR2SJ1_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet1.p4(),x))#0.6 for 3rd v of radial distribution studies
                                            
                                            getattr(self, 'F%s_dRSJ2_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies
                                            if dRSJ2_AK8<0.6: getattr(self, 'F%s_dR2SJ2_chargedHadrons'%(type)).Fill(self.DrRapPhi(subjet2.p4(),x))#0.6 for 3rd v of radial distribution studies

                                        getattr(self, 'F%s_pT_chargedHadrons'%(type)).Fill(pT)
                                        getattr(self, 'F%sjetPt_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR,AK8jet.p4().Pt())#x is dR, y is pT
                                        getattr(self, 'F%sjetmSD_dR_chargedHadrons'%(type)).Fill(dR if isGen else puppi_dR, getattr(AK8jet, 'mSD' + ('_nom' if not(isGen) else '' ))) #x is dR, y is pT
                                        getattr(self, 'F%s_dR_Pt_chargedHadrons'%(type)).Fill(pT, dR if isGen else puppi_dR)     #,AK8jet.p4().Pt())#x is dR, y is pT
                                        
                                    nchHad+=1 
                if not (self.onlyUnc.startswith('_je')):
                    if not(noPlot):# and not(self.onlyUnc):
                        if central: 
                            getattr(self, '%s_nphotons'%(type)).Fill(npho)
                            getattr(self, '%s_nneutral'%(type)).Fill(nneutHad)
                            getattr(self, '%s_ncharged'%(type)).Fill(nch)
                            getattr(self, '%s_nchargedHadrons'%(type)).Fill(nchHad)

                        else:
                            getattr(self, 'F%s_nphotons'%(type)).Fill(npho)
                            getattr(self, 'F%s_nneutral'%(type)).Fill(nneutHad)
                            getattr(self, 'F%s_ncharged'%(type)).Fill(nch)
                            getattr(self, 'F%s_nchargedHadrons'%(type)).Fill(nchHad)


            
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

    """
    def createNsubBasis(self, AK8jet, event, PFCollection, isGen=False ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        #### Run calculations of NSub bases and store for ungroomed AK8jets (default in CMS)

        #### Applying PUPPI weights to the PF candidates for reco and not for gen and pushing back constituents
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()

        for p in pfCands :
            #if p.p4().M()<0.:# and not(np.isclose(p.p4().M(),0.,rtol=10**(-9),atol=10**(-9))): 
            #    self.FlagBadPFCands=True
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            tp = tp * p.puppiWeight if not isGen else tp
            #except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        #print ("pushing back candidates")
        for x in CandsPUPPIweightedVec:
            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
        #print ("pushed back candidates")

        #if not(self.controlPlotsOnly):
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

        #### filing branches with nsub basis
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
        #    if len(sd_AK8jets)>0: #stupidly, in some rare cases, this will not be true (with ptmin>=170) leading to errors, hence the else (and also switching to ptmin=0.)
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
    def fillBranches( self, event, jetLabel, jetInfo, dummy=False, sys='_nom' ): 
        #use dummy to fill dummy values into branches for pass/nonpass (gen/reco)sel to maintain correspondence between events in different branches
        dummyFill=-929
        #### Filling branch with passAK8jet info after selection
        #print (jetLabel+"_pt", [ getattr( iJ['jet'], 'pt'+jetLabel.split('Jets')[1] ) )
        
        if not dummy:
            #if event.event%2==True: print ([ iJ['jet'].eta  )
            
            if 'selgenjets' in jetLabel.lower(): self.out.fillBranch('genSelectedEventNumber'+sys, event.event)
            elif 'accepgenjets' in jetLabel.lower() : self.out.fillBranch('accepgenSelectedEventNumber'+sys, event.event)
            elif 'selrecojets' in jetLabel.lower(): self.out.fillBranch('recoSelectedEventNumber'+sys, event.event)
            elif 'truerecojets' in jetLabel.lower() : self.out.fillBranch('truerecoSelectedEventNumber'+sys, event.event)
            
            if not('leading' in jetLabel.lower()): self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
            c=0
            for i,iJ in jetInfo.items():
                if c==0:
                    if 'reco' in jetLabel.lower():

                        self.out.fillBranch(jetLabel+"_mSD",  getattr(iJ['jet'], 'mSD'+sys) )

                        self.out.fillBranch(jetLabel+"_msoftdrop",  getattr(iJ['jet'], 'msoftdrop') ) #use uncorrected softdrop mass since the JMAR corrected ones are relevant only at the mass peak and existing corrections are outdated as is
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI",  getattr(iJ['jet'], 'msoftdrop_corr_PUPPI') ) # PUPPI corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_raw",  getattr(iJ['jet'], 'msoftdrop_raw') ) # prev. incorrect calc of msoftdrop_raw in nanoAOD-tools, now fixed, changes below recorrected branches for msoftdrop via nanoAOD-tools
                        self.out.fillBranch(jetLabel+"_msoftdrop_nom",  getattr(iJ['jet'], 'msoftdrop_nom') ) # prev. incorrect calc of msoftdrop_nom in nanoAOD-tools, now fixed, as a result of fixing raw
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMS",  getattr(iJ['jet'], 'msoftdrop_corr_JMS') ) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMR",  getattr(iJ['jet'], 'msoftdrop_corr_JMR') ) # JMR corr factor
                        self.out.fillBranch(jetLabel+"_JERfactor",  getattr(iJ['jet'], 'corr_JER') if self.isMC else 1. ) # JER corr factor
                        self.out.fillBranch(jetLabel+"_JECfactor",  getattr(iJ['jet'], 'pt'+sys)/getattr(iJ['jet'], 'pt_raw') ) # JEC corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_JMAR",  getattr(iJ['jet'], 'msoftdrop'+sys) ) # default fully corrected softdrop from nanoAOD-tools (with new raw fix) also stored for comparison

                        self.out.fillBranch(jetLabel+"_pt",  getattr(iJ['jet'], 'pt'+sys)  )
                        self.out.fillBranch(jetLabel+"_pt_raw",  getattr(iJ['jet'], 'pt_raw')  )
                        #self.out.fillBranch(jetLabel+"_pt_nom",  getattr(iJ['jet'], 'pt_nom')  )
                        self.out.fillBranch(jetLabel+"_mass",  getattr(iJ['jet'], 'mass'+sys)  )
                        
                    elif 'gen' in jetLabel.lower():
                        self.out.fillBranch(jetLabel+"_pt",  iJ['jet'].pt  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass",  iJ['jet'].mass  )
                        self.out.fillBranch(jetLabel+"_mSD",  iJ['jet'].mSD  )
                        

                    self.out.fillBranch(jetLabel+"_eta",  iJ['jet'].eta  )
                    self.out.fillBranch(jetLabel+"_y",  iJ['jet'].rapidity  )
                    self.out.fillBranch(jetLabel+"_phi",  iJ['jet'].phi  )

                    #if not(self.controlPlotsOnly):
                    self.out.fillBranch(jetLabel+"_tau21",  iJ['tau21']  )
                    self.out.fillBranch(jetLabel+"_tau32",  iJ['tau32']  )
                    self.out.fillBranch(jetLabel+"_tau21_WTA",  iJ['tau21_WTA']  )
                    self.out.fillBranch(jetLabel+"_tau32_WTA",  iJ['tau32_WTA']  )
                    self.out.fillBranch(jetLabel+"_tau21_exkT", iJ['tau21_exkT']  )
                    self.out.fillBranch(jetLabel+"_tau32_exkT", iJ['tau32_exkT']  )

                    """
                    self.out.branch(iJ+'_mSD',  'F')# hack, new branch to store recorrected msoftdrop and circumvent name clashes or overwriting nanoAOD-tools branches

                    #if not(self.controlPlotsOnly):#avoid going through the bulk of the info being saved and processed otherwise
                    if not('gen' in iJ.lower()): 
                        self.out.branch(iJ+'_pt_raw',  'F') 
                        #self.out.branch(iJ+'_pt_nom',  'F')# a la nanoAOD-tools recorrection with JECs in postproc. 
                        self.out.branch(iJ+'_msoftdrop_nom',  'F')  # a la nanoAOD-tools, stored for posterity as its own branch in custom flat, unjagged selRecoJets collection
                        self.out.branch(iJ+'_msoftdrop_corr_PUPPI',  'F') 
                        self.out.branch(iJ+'_msoftdrop_raw',  'F') 
                        self.out.branch(iJ+'_msoftdrop_corr_JMS',  'F') 
                        self.out.branch(iJ+'_msoftdrop_corr_JMR',  'F') 
                        self.out.branch(iJ+'_msoftdrop_JMAR',  'F') 
                        self.out.branch(iJ+'_msoftdrop_corr_JER',  'F') 
                        self.out.branch(iJ+'_JECfactor',  'F') 
                    """
        
                    for tauN in range(1, self.maxTau+1):
                        for pref in self.tauPrefixes:
                            beta=pref.split('_')[2]
                            #if 'reco' in jetLabel.lower(): print (event.event, jetLabel+pref+str(tauN), iJ[beta+str(tauN)] )
                            self.out.fillBranch(jetLabel+pref+str(tauN), iJ[beta+str(tauN)]  )
                c+=1
        
        else:
            #fill dummies
            
            #if not(self.controlPlotsOnly):
            if 'accepgenjets' in jetLabel.lower(): self.out.fillBranch('accepgenSelectedEventNumber'+sys, dummyFill)
            elif 'truerecojets' in jetLabel.lower(): self.out.fillBranch('truerecoSelectedEventNumber'+sys, dummyFill)
        
            if not('leading' in jetLabel.lower()): self.out.fillBranch( 'n'+jetLabel, dummyFill)

            c=0
            for i,iJ in jetInfo.items():
                if c==0:

                    if 'reco' in jetLabel.lower():

                        self.out.fillBranch(jetLabel+"_mSD",  dummyFill )

                        self.out.fillBranch(jetLabel+"_msoftdrop",  dummyFill ) #use uncorrected softdrop mass since the JMAR corrected ones are relevant only at the mass peak and existing corrections are outdated as is
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI",  dummyFill ) # PUPPI corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_raw",  dummyFill ) # prev. incorrect calc of msoftdrop_raw in nanoAOD-tools, now fixed, changes below recorrected branches for msoftdrop via nanoAOD-tools
                        self.out.fillBranch(jetLabel+"_msoftdrop_nom",  dummyFill ) # prev. incorrect calc of msoftdrop_nom in nanoAOD-tools, now fixed, as a result of fixing raw
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMS",  dummyFill ) # JMS corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_JMR",  dummyFill ) # JMR corr factor
                        self.out.fillBranch(jetLabel+"_JERfactor",  dummyFill ) # JER corr factor
                        self.out.fillBranch(jetLabel+"_JECfactor",  dummyFill ) # JEC corr factor
                        self.out.fillBranch(jetLabel+"_msoftdrop_JMAR",  dummyFill ) # default fully corrected softdrop from nanoAOD-tools (with new raw fix) also stored for comparison

                        self.out.fillBranch(jetLabel+"_pt",  dummyFill  )
                        self.out.fillBranch(jetLabel+"_pt_raw",  dummyFill  )
                        #self.out.fillBranch(jetLabel+"_pt_nom",  getattr(iJ['jet'], 'pt_nom')  )
                        self.out.fillBranch(jetLabel+"_mass",  dummyFill  )
                        
                    elif 'gen' in jetLabel.lower():
                        self.out.fillBranch(jetLabel+"_pt",  dummyFill  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass",  dummyFill  )
                        self.out.fillBranch(jetLabel+"_mSD",  dummyFill  )

                    self.out.fillBranch(jetLabel+"_eta", dummyFill  )
                    self.out.fillBranch(jetLabel+"_y", dummyFill  )
                    self.out.fillBranch(jetLabel+"_phi", dummyFill  )
                    #if not(self.controlPlotsOnly):
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
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)


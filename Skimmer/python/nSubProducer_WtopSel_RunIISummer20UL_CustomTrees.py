# Selection that requires >=1 AK8 jet (w/ pT>200 GeV,y<|1.7|), 1 high pT muon (pT>55) in the event, 
# 0 further leptons in system (including veto on any loose electrons w/ pT>40 GeV  in event)
# MET pT> 50 GeV, a leptonic hemisphere W (->mu nu) system with pT > 150 GeV, 
# AK4s in leptonic/hadronic hemisphere w/ pT>30 GeV, 
# AK4s in leptonic hemisphere identified as 0.4<deltaR(mu,AK4)<1.6

# ensure separation of hadronic AK8 from leptonic hemisphere via dPhi(AK8, mu)>2. 
# & dR(AK8, leading lept. AK4)>1.6,
# AK4's in hadronic hemisphere then identified by dR(AK8, AK4)<1.6 
# [ (0.8<dR(AK8, AK4)<1.6, PREV. w) /dR(AK8, AK4)<0.8, PREV top]
# require btag on leading AK4 in hadronic hemisphere (new as per theorists 11/22)
# in histogramming separate between W/top with dR(had hem AK4,AK8), softdrop mass and pT
# a la: W==> 0.8<dR(AK8, AK4)<1.6, 60.<m_SD<120., pT>200 GeV
#     top==> dR(AK8, AK4)<0.8, 140.<m_SD<220., pT>400 GeV 




import ROOT
import math, os, sys
import numpy as np
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

class nSubProd(Module):

    def __init__(self, sysSource=[], leptonSF={}, year='2017', isMC=True, onlyUnc='', onlyTrees=False, evtSelection='_WtopSel', isSigMC=False):

        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.onlyTrees = False#onlyTrees, hardcoded since option not used, need to erase from code's logic but not a problem
        self.isSigMC = isSigMC
        self.alsoDeltaRMatchedB = False # track events where the hadronic hemisphere b jets in reco and gen  are also deltaR matched just as the AK8's are

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

        ### Kinematics Cuts AK8Jets ###
        self.minAK8JetPt = 170  ### this is the basic minimum, not the final
        self.maxJetAK8Rap = 1.7

        ### Cuts for selections between W/top (previously), now used more inclusively (except in control histos) 
        ### ie, cut on minLeadAK8JetPtW and minSDMassW to store all event passing that criteria to store both W/top in same skims
        ### and, to navigate whether we store the PUPPI corrected (W/unmerged top-cand) or uncorrected (merged top-cand) softdrop mass in control histos
         
        self.minLeadAK8JetPtW = 200.
        self.minSDMassW = 60.#60.#   ### looser to pick low bins
        self.maxSDMassW = 120. #120.#  ### looser to pick higher bins
        self.minLeadAK8JetPtTop= 400.
        self.minSDMassTop = 140.
        self.maxSDMassTop = 240.
        
        self.minleadJetMass = 55. #inv. mass
        self.minleadJetpT = 200. 
        self.METCutWtop = 50.
        self.minLeptonicWPt = 150.


        ### Kinematics Cuts AK4 Jets ###
        self.minJetPt = 30.
        self.maxJetEta = 2.4
        if self.year=='2017':
            self.minBDisc = 0.3040   ### L: 0.0532, M: 0.3040, T: 0.7476, for DeepJet (ie, DeepFlavB);https://btv-wiki.docs.cern.ch/ScaleFactors/UL2017/

        elif self.year=='2018':
            self.minBDisc = 0.2783   ### L: 0.0490, M: 0.2783, T: 0.7476; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/

        elif self.year=='2016_preVFP':
            self.minBDisc = 0.2598   ### L: 0.0508, M: 0.2598, T: 0.6502; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016preVFP/

        elif self.year=='2016':
            self.minBDisc = 0.2489   ### L: 0.0480, M: 0.2489, T: 0.6377; https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016postVFP/

        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 55.
        self.maxMuonEta = 2.4
        self.minMuonMETPt = 50.

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 40.
        self.minTightElectronPt = 120.
        self.minElectronMETPt = 80.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.4]

        #overall event weights, updated in functions below as required per nominal/systematics runs of the skimmer
        self.totalRecoWeight = 1.
        #self.totalGenWeight = 1.
        
        # Nominal values of systematic weights handled by the following class-level variables for events passing selections
        # Only storing for some systematics, since storing the nominals for the others don't make sense 
        # (ie, eg, nominal for isr/fsrWeight is 1, and PSWeights store correlated (w_var/w_nom) weights a la ISRup/FSRnom,ISRnom/FSRnom,ISRdown/FSRnom, ISRnom/FSRdown) 

        #self.topreweight = topreweight 
        #self.topweight = 1.

        self.puWeight = 1. # storing nom. puWeight for reco
        self.evtGenWeight = 1. # storing event.genWeight for selected events
        self.pdfWeight = 1. # storing nominal pdfWeight, to see if it deviates from 1 
        #self.psWeight = 1.
        self.pdfWeightAll = np.ones((103,), dtype=np.float32)*1.
        self.btaggingWeight = 1.
        self.l1PreFireWeight = 1.
        #self.nGenbs = 0.
        #self.nRecobs = 0.
        self.leptonWeight = 1.

        self.pdfWeightUp = 1.
        self.puWeightUp = 1.
        self.l1PreFireWeightUp = 1.
        self.isrWeightUp = 1.
        self.fsrWeightUp = 1.
        self.btaggingWeightUp = 1.
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
        self.btaggingWeightDown = 1.
        self.leptonWeightAllDown = 1.
        self.leptonWeightISODown = 1.
        self.leptonWeightIDDown = 1.
        self.leptonWeightTrigDown = 1.
        self.leptonWeightRecoEffDown = 1.

        #evt categeories: 0 for 1 had. hem. AK4 b-tag with 1 lept hem btag in event
            #                 1 for 1 had. hem. b-tag
            #                -1 undefined/no selection
        self.eventCategory = -1
        self.recoEventCategory = -1
        self.genEventCategory = -1
        self.dummy = 0

        ### Defining nsubjetiness basis
        self.maxTau = 5
        self.nSub_labels = {
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
                        "_tau_2_1": [ 0., 0.9, 900  ],
                        "_tau_2_2": [ 0., 0.5, 500  ],
                        "_tau_2_3": [ 0., 0.5, 500  ],
                        "_tau_2_4": [ 0., 0.5, 500  ],
                        "_tau_2_5": [ 0., 0.5, 500  ]
                }
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes, WTA_kT=3, OP_kT=6
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
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
        self.sysSource = ['_nom'] +[ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
      
        # pu/btag/leptonWeights used only to change reco event weight (up/down) without application to gen weight, 
        # others applied to modify self.genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)] ##### this is now done in the histogramming #####
        # weights+up/down variations are saved for accepted events, to be applied in histogramming in off-crab postprocessing
        self.sysWeightList = ( '_pu', '_btag', '_lepton', '_pdf', '_isr', '_fsr')
        self.sysrecoWeightList = ( '_pu', '_btag', '_lepton')
        self.sysgenWeightList = ( '_pdf', '_isr', '_fsr' ) 
        
        '''
        #deprecated, used previously on top cands for top pT reweightings
        self.FLAGS = [
        "isPrompt",
        "isDecayedLeptonHadron",
        "isTauDecayProduct",
        "isPromptTauDecayProduct",
        "isDirectTauDecayProduct",
        "isDirectPromptTauDecayProduct",
        "isDirectHadronDecayProduct",
        "isHardProcess",
        "fromHardProcess",
        "isHardProcessTauDecayProduct",
        "isDirectHardProcessTauDecayProduct",
        "fromHardProcessBeforeFSR",
        "isFirstCopy",
        "isLastCopy",
        "isLastCopyBeforeFSR",
        ]
        """bit-packed statusFlags interpretations.  Use `GenParticle.hasFlags` to query""" 
        '''
        #############################################################################

    def beginJob(self, histFile, histDirName):

        Module.beginJob(self, histFile, histDirName)

        #tauBins = 1000 #====> bins decided in histogramming
        
        ### Booking histograms
        #selList = ['1hadbtag','2hem2btag']
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

        

        if not self.onlyUnc:
        
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   25, 0, 25) )
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            #self.addObject( ROOT.TH1F('Lepweight',   ';LepWeight',   20, 0, 2) )
            #self.addObject( ROOT.TH1F('Btagweight',   ';BtagWeight',   25, 0, 2) )
            #self.addObject( ROOT.TH1F('Topweight',   ';Topweight',   25, 0, 2) )
            allSel = selList + ['_WSel','_topSel'] if self.evtSelection.startswith('_Wtop') else selList
            #### general selection
            for isel in [ '_noSelnoWeight', '_noSel' ] + allSel:
                self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   10, 0, 10) )
                
                self.addP4Hists( 'muons', isel )
                self.addP4Hists( 'eles', isel )
                
                self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   10, 0, 10) )
                self.addP4Hists( 'AK8jets', isel )
                
                self.addObject( ROOT.TH1F('nAK4jets'+isel,   ';number of AK4 jets',   20, 0, 20) )
                #self.addObject( ROOT.TH1F('nAK4hadjets'+isel,   ';number of had. hem. AK4 jets',   40, 0, 20) )
                #self.addObject( ROOT.TH1F('nAK4lepjets'+isel,   ';number of lept. hem. AK4 jets',   40, 0, 20) )
                self.addObject( ROOT.TH1F('nAK4bjets'+isel,   ';number of b-tagged AK4s',   10, 0, 10) )
                self.addObject( ROOT.TH1F('nAK4bhadjets'+isel,   ';number of b-tagged had. hem. AK4s',   10, 0, 5) )
                self.addObject( ROOT.TH1F('dR_AK8_AK4bhad'+isel,   ';#Delta R( reco AK8, b-tagged AK4 in hadr. hem. )',   40, 0., 2.) )
                #self.addObject( ROOT.TH1F('nAK4blepjets'+isel,   ';number of b-tagged lept. hem. AK4s',   40, 0, 20) )
                
                self.addP4Hists( 'AK4jets', isel )
                self.addP4Hists( 'AK4btaggedjets', isel )
                self.addP4Hists( 'leptonicW', isel )
                #self.addP4Hists( 'AK4btaggedhadjets', isel )
                
                self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )
                #self.addObject( ROOT.TH1F('leptonicWMT'+isel,   ';leptonic W m_T(GeV)',   400, 0, 4000) )
                #self.addObject( ROOT.TH1F('Mtt'+isel,   '; m_{t#bar{t}}(GeV)',   400, 0, 4000) )
                self.addObject( ROOT.TH1F('HT'+isel,   ';HT (GeV)',   200, 0, 2000) )
                #self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

            if self.isMC:
                for isel in [ '_noSel' ] + allSel:
                    self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   10, 0, 10) )
                    
                    self.addP4Hists( 'genmuons', isel )
                    self.addP4Hists( 'geneles', isel )
                    
                    self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   10, 0, 10) )
                    self.addP4Hists( 'AK8genjets', isel )
                    
                    self.addObject( ROOT.TH1F('ngenAK4jets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                    #self.addObject( ROOT.TH1F('ngenAK4hadjets'+isel,   ';number of had. hem. AK4 genjets',   20, 0, 20) )
                    #self.addObject( ROOT.TH1F('ngenAK4lepjets'+isel,   ';number of lept. hem. AK4 genjets',   20, 0, 20) )
                    self.addObject( ROOT.TH1F('ngenAK4bjets'+isel,   ';number of b-flav. matched AK4 genjets',   10, 0, 10) )
                    self.addObject( ROOT.TH1F('ngenAK4bhadjets'+isel,   ';number of b-flav. matched had. hem. AK4 genjets',   10, 0, 5) )
                    self.addObject( ROOT.TH1F('dR_genAK8_genAK4bhad'+isel,   ';#Delta R( gen AK8, gen AK4 b-jet in hadr. hem. )',   40, 0, 2.) )
                    #self.addObject( ROOT.TH1F('ngenAK4blepjets'+isel,   ';number of b-flav. matched lept. hem. AK4 genjets',   20, 0, 20) )
                    
                    self.addP4Hists( 'AK4genjets', isel )
                    self.addP4Hists( 'AK4bmatchedgenjets', isel )
                    self.addP4Hists( 'genleptonicW', isel )
                    #self.addP4Hists( 'AK4bmatchedlepgenjets', isel )
                    
                    self.addObject( ROOT.TH1F('genMETPt'+isel,   ';gen MET (GeV)',   200, 0, 2000) )
                    #self.addObject( ROOT.TH1F('genleptonicWMT'+isel,   ';gen leptonic W m_T(GeV)',   400, 0, 4000) )
                    #self.addObject( ROOT.TH1F('genMtt'+isel,   '; gen m_{t#bar{t}}(GeV)',   400, 0, 4000) )
                    self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )                
            
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
        #if self.onlyTrees:
        self.out = wrappedOutputTree
        #self.out.branch('triggerWeight',  "F")
        tmplistAK8=[]
        tmplistAK4=[]
        tmplistMuon=[]
        tmplistLeptWPt=[]
        tmplistMETPt=[]

        # need to keep branches for selected AK8s, selected muon, selected(reconstructed) leptonicW,
        # MET pT, selected/b-tagged hadronic hemisphere AK4

        if not self.isMC: 
            
            self.out.branch('recoEventCategory_nom',  "I") 
            self.out.branch('totalRecoWeight_nom', "F" )
            #nRecoBtags=1 for evt. cat. 0; 2 for evt. cat. 1; >2 for evt. cat. 1 as well, 0 for evt. cat. -1
            # the nRecoBtags branch allow us to check if any other AK4 was a b-tagged jet except for the one we do the dR selections on
            self.out.branch('nRecoBtags_nom',  "I") 
            self.out.branch('nRecoHadBtags_nom',  "I") 
            self.out.branch('nRecoLepBtags_nom',  "I") 
            self.out.branch('recoSelectedEventNumber_nom', "L" )
            self.out.branch('good_nPVs_nom', "F" )
            self.out.branch('FlagRecoLeptHemBjet_nom', "I")
            self.out.branch('selRecoHadHemDeltaR_nom', "F" ) # dR(lead AK8, lead b-tagged AK4 in had. hem.)
            
            tmplistAK8.append('selRecoJets_nom')
            tmplistMuon.append('selRecoMu_nom')
            tmplistLeptWPt.append('selRecoLeptW_nom')
            tmplistMETPt.append('selRecoMET_nom')

            # We require exactly one b-tagged (leading) hadronic hemisphere AK4,
            # so keeping that objects info in a branch, 
            tmplistAK4.append('selRecoAK4bjetHadHem_nom') # to store hadronic hem., lead b-tagged AK4's

        elif self.isMC: 
            for sys in self.sysSource:

                #evt categeories: 0 for 1 had. hem. AK4 b-tag with 1 lept hem btag in event
                #                 1 for 1 had. hem. b-tag
                #                -1 for no selection 
                self.out.branch('recoEventCategory'+sys,  "I") 
                self.out.branch('totalRecoWeight'+sys, "F" )
                #nRecoBtags=1 for evt. cat. 0; 2 for evt. cat. 1; >2 for evt. cat. 1 as well, 0 for evt. cat. -1
                self.out.branch('nRecoBtags'+sys,  "I") 
                self.out.branch('nRecoHadBtags'+sys,  "I") 
                self.out.branch('nRecoLepBtags'+sys,  "I") 
                self.out.branch('recoSelectedEventNumber'+sys, "L" )
                self.out.branch('good_nPVs'+sys, "F" )
                                
                self.out.branch('leptonWeightNom'+sys,  "F")
                self.out.branch('btagWeightNom'+sys,  "F")
                self.out.branch('puWeightNom'+sys,  "F")
                self.out.branch('l1prefiringWeightNom'+sys,  "F")
                self.out.branch('pdfWeightNom'+sys, "F" ) 

                self.out.branch('FlagDeltaRMatchedBjets'+sys, "I")
                self.out.branch('FlagRecoLeptHemBjet'+sys, "I")
                
                self.out.branch( 'accepgenSelectedEventNumber'+sys, "L" ) 
                self.out.branch( 'truerecoSelectedEventNumber'+sys, "L" )
                
                # For tagging W vs. top, holds info on deltaR between b-jet and AK8jet
                self.out.branch('selRecoHadHemDeltaR'+sys,"F")

                if sys.startswith('_nom'): 
                    self.out.branch('genEventCategory'+sys,  "I") 
                    self.out.branch('evtGenWeight'+sys, "F" )  # nominal generator weight for evt.
                    self.out.branch('genSelectedEventNumber'+sys, "L" )

                    #nGenBtags=1 for evt. cat. 0; 2 for evt. cat. 1; >2 for evt. cat. 1 as well, 0 for evt. cat. -1
                    self.out.branch('nGenBtags'+sys,  "I") #for hadron flavour ghost-matched b's
                    self.out.branch('nGenHadBtags'+sys,  "I") #for hadron flavour ghost-matched b's
                    self.out.branch('nGenLepBtags'+sys,  "I") #for hadron flavour ghost-matched b's
                    
                    self.out.branch('FlagGenLeptHemBjet'+sys, "I")

                    # For tagging W vs. top, holds info on deltaR between b-jet and AK8jet
                    self.out.branch('selGenHadHemDeltaR'+sys, "F")

                
                if self.isSigMC and not self.onlyUnc: 
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

                    self.out.branch('l1prefiringWeightUp'+sys,  "F")
                    self.out.branch('l1prefiringWeightDown'+sys,  "F")

                    self.out.branch( 'btagWeightUp'+sys, "F" )
                    self.out.branch( 'btagWeightDown'+sys, "F" )

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
            tmplistAK8_accepgen = ['accepGenJets'+sys for sys in self.sysSource]
            for x in tmplistAK8_reco+tmplistAK8_truereco+tmplistAK8_gen+tmplistAK8_accepgen: tmplistAK8.append(x)

            #if sys.startswith('_nom'):
            tmplistmu_reco =  [ 'selRecoMu'+sys for sys in self.sysSource if sys.startswith('_nom')]
            tmplistmu_gen  =  [ 'selGenMu'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            for x in tmplistmu_reco+tmplistmu_gen : tmplistMuon.append(x)

            tmplistLeptWPt_reco =  [ 'selRecoLeptW'+sys for sys in self.sysSource if sys.startswith('_nom')]
            tmplistLeptWPt_gen  =  [ 'selGenLeptW'+sys for sys in self.sysSource  if sys.startswith('_nom')]
            for x in tmplistLeptWPt_reco+tmplistLeptWPt_gen: tmplistLeptWPt.append(x)

            tmplistMETPt_reco =  [ 'selRecoMET'+sys for sys in self.sysSource if sys.startswith('_nom')]
            tmplistMETPt_gen  =  [ 'selGenMET'+sys for sys in self.sysSource if sys.startswith('_nom')] 
            for x in tmplistMETPt_reco+tmplistMETPt_gen: tmplistMETPt.append(x)

            '''
            tmplistAK4hadHem_reco = [ 'selRecoAK4sHadHem'+sys for sys in self.sysSource ]
            tmplistAK4hadHem_gen = ['selGenAK4sHadHem'+sys for sys in self.sysSource if sys.startswith('_nom') ] 
            tmplistAK4lepHem_reco = [ 'selRecoAK4sLepHem'+sys for sys in self.sysSource ]
            tmplistAK4lepHem_gen = ['selGenAK4sLepHem'+sys for sys in self.sysSource if sys.startswith('_nom') ]
            for x in tmplistAK4hadHem_reco+tmplistAK4hadHem_gen+tmplistAK4lepHem_reco+tmplistAK4lepHem_gen: tmplistAK4.append(x)
            
            '''
            # Store leading reco/gen AK4 (b-tagged/hadron flavour matched) jet if any in the hadronic hemisphere
            #requring >=1 btag, but with exactly one b-tag in the hadronic hemisphere
            tmplistAK4bjethadHem_reco = ['selRecoAK4bjetHadHem'+sys for sys in self.sysSource ]
            tmplistAK4bjethadHem_gen = [ 'selGenAK4bjetHadHem'+sys for sys in self.sysSource if sys.startswith('_nom') ] 
            for x in tmplistAK4bjethadHem_reco+tmplistAK4bjethadHem_gen: tmplistAK4.append(x)
            
 
            
        print ("Stored AK8 jet branches:", tmplistAK8)
        for iJ in tmplistAK8:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_y',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_msoftdrop',  'F')#, lenVar='n'+iJ)
            if not('gen' in iJ.lower()): self.out.branch(iJ+'_msoftdrop_corr_PUPPI',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_exkT',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_exkT',  'F')#, lenVar='n'+iJ)
            #self.out.branch(iJ+'_genMatched',  'O')#, lenVar='n'+iJ)
            
            for x in self.nSub_labels:
                self.out.branch(iJ+x, 'F')#, lenVar='n'+iJ )
        
        print ("Stored AK4 jet branches:", tmplistAK4)
        for iJ in tmplistAK4:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_y',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F')#, lenVar='n'+iJ)
            #if not('gen' in ij.lower()) and ('bjet' in i.lower()):
            #    self.out.branch(iJ+'_btagVal',  'F', lenVar='n'+iJ)
            #    self.out.branch(iJ+'_btagSF',  'F', lenVar='n'+iJ)
            #    self.out.branch(iJ+'_btagSFUp',  'F', lenVar='n'+iJ)
            #    self.out.branch(iJ+'_btagSFDown',  'F', lenVar='n'+iJ)
            #elif ('gen' in ij.lower()) and ('bjet' in i.lower()):
            #    self.out.branch(iJ+'_hadronFlavourMatched',  'F', lenVar='n'+iJ)


        print ("Stored mu branches:", tmplistMuon)
        for iMu in tmplistMuon:
            self.out.branch('n'+iMu,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iMu+'_pt',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_eta',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_y',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_phi',  'F')#, lenVar='n'+iMu)
            self.out.branch(iMu+'_mass',  'F')#, lenVar='n'+iMu)
            #self.out.branch(iMu+'_mindRAK4mu')#,  'F', lenVar='n'+iMu)
            #if (not 'gen' in iMu.lower()):
            #    #self.out.branch(iMu+'_pTrelCustom')#,  'F', lenVar='n'+iMu)
            #    #self.out.branch(iMu+'_pTrelNano')#,  'F', lenVar='n'+iMu)
            #    #self.out.branch(iMu+'_miniPFRelIsoAll')#,  'F', lenVar='n'+iMu)
            #    self.out.branch(iMu+'_leptonSF','F', 4) #to store in format [SFTrigger(nom/up/down), SFID(nom/up/down), SFISO(nom/up/down), SFRecoEff(nom/up/down)]
            #    self.out.branch(iMu+'_leptonSFUp','F', 4)
            #    self.out.branch(iMu+'_leptonSFDown','F', 4)

        print ("Stored MET branches:", tmplistMETPt)
        for i in tmplistMETPt:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            
        print ("Stored leptonic W branches:", tmplistLeptWPt)
        for i in tmplistLeptWPt:
            self.out.branch('n'+i,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(i+'_pt',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_eta',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_y',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_phi',  'F')#, lenVar='n'+i)
            self.out.branch(i+'_mass',  'F')#, lenVar='n'+i)

        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        if not self.onlyTrees:
            if self.isMC and not self.onlyUnc:

                self.genLevel = self.response+self.miss
                #self.genLevelW = self.responseW+self.missW
                #self.genLeveltop = self.responsetop+self.misstop

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

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGendR_AK8_hadAK4b, selGenMuons, selGenAK4bjets, selGenAK4bhadjets, selGenLeptW, selGenJets, selGenMET = self.genSelection(event) 
        passRecoSel, iRecoSel, selRecodR_AK8_hadAK4b, selRecoMuons, selRecoAK4bjets, selRecoAK4bhadjets, selRecoLeptW, selRecoJets, selRecoMET = self.recoSelection(event) 
        #print (selRecodR_AK8_hadAK4b)
        if not self.isMC:
            if not passRecoSel['_nom']: 
                self.totalRecoWeight=0.
                return False
        else:
            if (not passGenSel) and (not passRecoSel['_nom']): 
                self.totalRecoWeight=0.
                self.evtGenWeight=0.
                return False

        for sys in self.sysSource:

            if not self.isMC: 
                self.totalRecoWeight=1.
                #leptonSFs = [1, 1, 1, 1.]
                self.leptonWeight = 1. #np.prod(leptonSFs)
                self.btaggingWeight = 1.
                self.l1PreFireWeight = 1. 

            if self.isMC: 
                self.puWeight = event.puWeight
                self.evtGenWeight = event.genWeight
                self.l1PreFireWeight = event.L1PreFiringWeight_Nom

                #fix b-tagging weight!!!!!!!!!!!!!!!!
                ######################################
                
                #### b-tagging Weights #####
                bTagSFs=[]
                w=0.
                wup=0
                wdn=0
                if len(selRecoAK4bjets)>0: 
                    w=1.
                    bTagSFs = [x.btagSF_deepjet_M for x in selRecoAK4bjets]
                    for i in bTagSFs:
                        w *= i  
                self.btaggingWeight = w 
                
                ######################################
                #self.btaggingWeight = 1. ############# Dummy
                #self.leptonWeight = 1. ############# Dummy
                
                if len(selRecoMuons)>0 and passRecoSel['_nom']: 

                    if len(selRecoMuons)>1: print ("!!!!!!!!!!!!!!!Warning, extra muons leaking, check selection!!!!!!!!!!!!!!!!!!!")
                    leptonSFs = self.getleptonSF( "muon", selRecoMuons[0] )
                    self.leptonWeight = np.prod([i[0] for i in leptonSFs])
                else: 

                    leptonSFs = [0, 0, 0, 0]
                    self.leptonWeight = np.prod(leptonSFs)
                
                self.totalRecoWeight = self.evtGenWeight*self.puWeight*self.leptonWeight*self.btaggingWeight*self.l1PreFireWeight

                if not sys.startswith(('_jes', '_jer')):#self.sysWeightList):
    
                    #if not sys.startswith('_nom'): #hardcoded, will remove; just in case someone tries to oaccidentally add a systematic in the nominal mode
                    selRecoJets[sys] = selRecoJets['_nom']
                    selRecoAK4bhadjets[sys] = selRecoAK4bhadjets['_nom']
                    selRecodR_AK8_hadAK4b[sys] = selRecodR_AK8_hadAK4b['_nom']
                    iRecoSel[sys] = iRecoSel['_nom']
                    passRecoSel[sys] = passRecoSel['_nom']

                    # PDF sets for most/all(to be checked if true) RunIISummer20UL samples seem to be NNPDF31_nnlo_as_0.118_mc_hessian_pdfas (pdfid=325300), definitely true for QCD_HT signal MC for dijets
                    # structure of the pdf set's array of 103 members a la: https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_as_0118_mc_hessian_pdfas/NNPDF31_nnlo_as_0118_mc_hessian_pdfas.info
                    # general structure of pdf weights = w_i,var/w_nom, where i runeth over 0-102, where [0]=>nominal from avg. over replicas, 
                    # [1-100]=> PDF eigenvectors of the covariance matrix in the parameter space, 
                    # [101,102]=> central value for (forced positive definite) a_s=[0.116,0.120]
                    if self.isSigMC:
                        
                        ############## b-tag weight variations for nominal+wtUnc runs
                        tempWtbtagup=1.
                        tempWtbtagdn=1.
                        bTagSFsUp = [x.btagSF_deepjet_M_up for x in selRecoAK4bjets]
                        bTagSFsDown = [x.btagSF_deepjet_M_down for x in selRecoAK4bjets]
                        for i,j in zip(bTagSFsUp,bTagSFsDown):
                            tempWtbtagup *= i
                            tempWtbtagdn *= j
                        self.btaggingWeightUp=tempWtbtagup
                        self.btaggingWeightDown=tempWtbtagdn 

                        ############## PDF and gen weights for nominal+wtUnc runs
                        pdfWeights =  getattr( event, 'LHEPdfWeight' ) #convert to a simple numpy array to extract replica weights and thereby the variations
                        self.pdfWeightAll = np.array([pdfWeights[i] for i in range(pdfWeights.GetSize())], dtype=np.float32)
                        self.pdfWeight = pdfWeights[0]
                        
                        # 'Naively' defining pdf uncertainty symmetrically by taking the envelope defined as per Eq. 15 of https://link.springer.com/article/10.1140/epjc/s10052-015-3318-8#Sec29 
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
                        self.l1PreFireWeightUp = event.L1PreFiringWeight_Up
                        self.l1PreFireWeightDown = event.L1PreFiringWeight_Dn


                        ############# Naive all sources up/down variations of lepton SF  for lepton weight calc.
                        self.leptonWeightAllUp = np.prod([ i[1] for i in leptonSFs]) if not self.leptonWeight==0 else 0
                        self.leptonWeightAllDown = np.prod([ i[2] for i in leptonSFs]) if not self.leptonWeight==0 else 0

                        ############# Separated, all sources up/down variations of lepton SF for lepton weight variation calc.
                        #lepton SFs=[[trig],[ID],[ISO],[RecoEff]], in each SF, order=nom,up,down

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
                    self.totalRecoWeight = self.evtGenWeight*self.puWeight*self.leptonWeight*self.btaggingWeight*self.l1PreFireWeight
            
            genJet = OrderedDict()
            recoJet = OrderedDict()
            #tmpRecoJet = OrderedDict()
            #tmpGenJet = OrderedDict()

            if self.isMC and (not passRecoSel[sys]): 
                #fill in 'missing' gen jets, ie those going into RM underflow since event didn't pass the nominal reco seleciton

                if passGenSel:
                    #### Gen Misses
                    self.miss=self.miss+1
                
                    # inp. to createNsubBasis: ( AK8 jet p4, evt. pointer, particle collection, isGen)
                    genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )

                    # inp. to fillBranches: ( evt. pointer, obj name/type, gen/reco jet, isDummy, sys name)
                    # if not a dummy fill set an isDummy flag (arg name is actually dummy)=False
                    if sys.startswith('_nom'): 
                
                        # In this case this is a selected gen jet, so fill dummy=False
                        self.fillAK8Branches( event, 'selGenJets'+sys, genJet, False, sys )
                        self.fillOtherBranches( event, 'selGenAK4bjetHadHem'+sys, selGenAK4bhadjets[0], 1)#, False, sys )
                        self.fillOtherBranches( event, 'selGenMu'+sys, selGenMuons[0], 1 )#, False, sys )
                        self.fillOtherBranches( event, 'selGenLeptW'+sys, selGenLeptW[0], 1)#, False, sys )
                        self.fillOtherBranches( event, 'selGenMET'+sys, selGenMET, 1)#, False, sys )
                
                    # filling for all jer/jes systematics separately from nominal case
                    # since in those cases the selected gen jet and smeared selected reco jet might not pass the deltaR match
                    # In this case this is not an accepted (for filling the centre of the RM) gen jet, so fill dummy=True to track misses
                    self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, True, sys ) #Will use dummies in accepgen arrays as masks to selGen to extract fakes' info
            
            if passRecoSel[sys]:
                self.recoLevel = self.recoLevel+1
                #tmpRecoJets[sys] = {}
                #tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                
                recoJet['Jet'] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )#tmpRecoJets[sys][0] 
                self.fillAK8Branches( event, 'selRecoJets'+sys, recoJet, False, sys ) #not a dummy fill so dummy=False
                self.fillOtherBranches( event, 'selRecoAK4bjetHadHem'+sys, selRecoAK4bhadjets[sys][0], 1)#, False, sys )
                if sys.startswith('_nom'):
                    self.fillOtherBranches( event, 'selRecoMu'+sys, selRecoMuons[0], 1)#, False, sys )
                    self.fillOtherBranches( event, 'selRecoLeptW'+sys, selRecoLeptW[0], 1)#, False, sys )
                    self.fillOtherBranches( event, 'selRecoMET'+sys, selRecoMET, 1)#, False, sys )
        
                if self.isMC: 
                    
                    deltaRmatch = False

                    if not passGenSel:  ##### fake reco
                        self.fakes = self.fakes+1
                        self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, True, sys ) #True=dummy fill, to track fake reco, similar strategy as with gen

                    else:
                        
                        genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                        
                        # consider event as missing gen if not deltaR matched 
                        if sys.startswith('_nom'): 
                            self.fillAK8Branches( event, 'selGenJets'+sys, genJet, False, sys )
                            self.fillOtherBranches( event, 'selGenAK4bjetHadHem'+sys, selGenAK4bhadjets[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenMu'+sys, selGenMuons[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenLeptW'+sys, selGenLeptW[0], 1)#, False, sys )
                            self.fillOtherBranches( event, 'selGenMET'+sys, selGenMET, 1)#, False, sys )
                    

                        #need to delta R match for true reco and accepted gen
                        if self.DrRapPhi( recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() ) < 0.4 :
                            deltaRmatch = True
                            self.response= self.response+1
                        
                        self.alsoDeltaRMatchedB=1 if self.DrRapPhi( selRecoAK4bhadjets[sys][0].p4(),selGenAK4bhadjets[0].p4())<0.2 else 0

                        
                        if deltaRmatch:
                            # fill only if deltaR matched, for eventual response matrix filling
                            self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, False, sys )    
                            self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, False, sys )
                        else:
                            self.fillAK8Branches( event, 'accepGenJets'+sys, genJet, True, sys )    #fill dummy=True to indicate fakes and misses 
                            self.fillAK8Branches( event, 'trueRecoJets'+sys, recoJet, True, sys )
                            self.fakes = self.fakes+1
                            self.miss = self.miss+1

            recoleptHemBjets=[]                    
            recoleptHemBjets = [x for x in selRecoAK4bjets if not x in selRecoAK4bhadjets[sys]]
            self.recoEventCategory=len(selRecoAK4bjets) if ((len(selRecoAK4bjets)>0 and len(selRecoAK4bhadjets[sys])>0)) else -1
            
            if sys.startswith('_nom') and self.isMC:
                genleptHemBjets=[]
                genleptHemBjets = [x for x in selGenAK4bjets if not x in selGenAK4bhadjets] 
                self.genEventCategory=len(selGenAK4bjets) if ((len(selGenAK4bjets)>0 and len(selGenAK4bhadjets)>0)) else -1
            
            if not self.isMC:
                self.out.fillBranch( 'recoEventCategory'+sys, self.recoEventCategory if passRecoSel[sys] else -1  )
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoBtags'+sys, len(selRecoAK4bjets) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoHadBtags'+sys, len(selRecoAK4bhadjets[sys]) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoLepBtags'+sys, len(recoleptHemBjets) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'selRecoHadHemDeltaR'+sys, selRecodR_AK8_hadAK4b[sys] if passRecoSel[sys] else 929.)
                self.out.fillBranch( 'FlagRecoLeptHemBjet'+sys, 1 if (passRecoSel[sys] and len(recoleptHemBjets)>0) else 0)
                self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                

            else:
                self.out.fillBranch( 'recoEventCategory'+sys, self.recoEventCategory if passRecoSel[sys] else -1)
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoBtags'+sys, len(selRecoAK4bjets) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoHadBtags'+sys, len(selRecoAK4bhadjets[sys]) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'nRecoLepBtags'+sys, len(recoleptHemBjets) if passRecoSel[sys] else 0)
                self.out.fillBranch( 'selRecoHadHemDeltaR'+sys, selRecodR_AK8_hadAK4b[sys] if passRecoSel[sys] else 929.)            
                self.out.fillBranch( 'FlagRecoLeptHemBjet'+sys, 1 if (passRecoSel[sys] and len(recoleptHemBjets)>0) else 0)
                self.out.fillBranch( 'FlagDeltaRMatchedBjets'+sys, self.alsoDeltaRMatchedB if (passRecoSel[sys] and passGenSel) else 0)
                self.out.fillBranch( 'totalRecoWeight'+sys, self.totalRecoWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'puWeightNom'+sys, self.puWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'l1prefiringWeightNom'+sys, self.l1PreFireWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'btagWeightNom'+sys, self.btaggingWeight if passRecoSel[sys] else 0.)
                self.out.fillBranch( 'leptonWeightNom'+sys, self.leptonWeight if passRecoSel[sys] else 0.)
                
                if sys.startswith('_nom'): 
                    self.out.fillBranch( 'genEventCategory'+sys, self.genEventCategory if passGenSel else -1)
                    self.out.fillBranch( 'nGenBtags'+sys, len(selGenAK4bjets) if passGenSel else 0)
                    self.out.fillBranch( 'nGenHadBtags'+sys, len(selRecoAK4bhadjets) if passGenSel else 0)
                    self.out.fillBranch( 'nGenLepBtags'+sys, len(genleptHemBjets) if passGenSel else 0)
                    self.out.fillBranch( 'selGenHadHemDeltaR'+sys, selGendR_AK8_hadAK4b if passGenSel else 929.)
                    self.out.fillBranch( 'FlagGenLeptHemBjet'+sys, 1 if (passGenSel and len(genleptHemBjets)>0) else 0)
                    self.out.fillBranch( 'evtGenWeight'+sys, self.evtGenWeight if passGenSel else 0.) 

                if self.isSigMC and not self.onlyUnc:
                    self.out.fillBranch( 'pdfWeightNom'+sys, self.pdfWeight if passGenSel else 0.)
                    self.out.fillBranch( 'pdfWeightAll'+sys, self.pdfWeightAll if passGenSel else np.zeros((103,),dtype=np.float32))
                    self.out.fillBranch( 'pdfWeightUp'+sys, self.pdfWeightUp if passGenSel else 0.)
                    self.out.fillBranch( 'pdfWeightDown'+sys, self.pdfWeightDown if passGenSel else 0.)
                    
                    self.out.fillBranch( 'isrWeightUp'+sys, self.isrWeightUp if passGenSel else 0.)
                    self.out.fillBranch( 'isrWeightDown'+sys, self.isrWeightDown if passGenSel else 0.)
                    
                    self.out.fillBranch( 'fsrWeightUp'+sys, self.fsrWeightUp if passGenSel else 0.)
                    self.out.fillBranch( 'fsrWeightDown'+sys, self.fsrWeightDown if passGenSel else 0.)
                    
                    self.out.fillBranch( 'puWeightUp'+sys, self.puWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'puWeightDown'+sys, self.puWeightDown if passRecoSel[sys] else 0.)

                    self.out.fillBranch( 'l1prefiringWeightUp'+sys, self.l1PreFireWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'l1prefiringWeightDown'+sys, self.l1PreFireWeightDown if passRecoSel[sys] else 0.)

                    self.out.fillBranch( 'btagWeightUp'+sys, self.puWeightUp if passRecoSel[sys] else 0.)
                    self.out.fillBranch( 'btagWeightDown'+sys, self.puWeightDown if passRecoSel[sys] else 0.)

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
        met = Object(event, 'MET')    
        
        ########### Lepton selection ###############
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minTightMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.tightId and abs(x.dxy)<0.2 and abs(x.dz)<0.5 and x.isGlobal and x.isTracker and x.highPtId and x.tkRelIso<0.3] #  applying tight selection on muons already here since we only veto for loose muons and loose electrons in event
        # descriptors for id/iso selectors used provided in: 
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#HighPt_Muon
        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons)

        ############################################

        ################### MET  #######################
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiM(met.pt, 0., met.phi, 0.)#met.sumEt)
        ################################################


        ############ Basic AK4 selection ###############
        recoAK4jets = [ x for x in jets if x.pt > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and (x.jetId >= 2)]# and x.btagDeepFlavB > self.minBDisc]
        recoAK4jets.sort(key=lambda x:x.pt,reverse=True)

        recoAK4bjets = [ x for x in recoAK4jets if x.btagDeepFlavB > self.minBDisc]
        recoAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ################################################
        
        if not len(AK8jets)==0:
            for ijets in AK8jets: 
                ijets.rapidity = ijets.p4().Rapidity()#self.etaToRapidity(ijets)
        recoAK8jets = {}
        passSel = {}
        iSel = {}
        recodR_AK8_hadAK4b = {}
        recoAK4bhadjets = {}
        for sys in self.sysSource:
             
            if sys.startswith(self.sysWeightList): sys = '_nom'

            ################### Basic AK8 jet selection ############################
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs(x.rapidity) < self.maxJetAK8Rap and (x.jetId >= 2) ]
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)
            AK8HT = sum( [ getattr( x, 'pt'+sys ) for x in recoAK8jets[sys] ] )
            ########################################################################
            
            ################################## Applying selection ##########################################
            
            # Select the event as generically passing W/top selection (ie, a boosted semileptonic ttbar selection):
            # >0 boosted heavy AK8, >=1 btag (w/ exactly one required btag on the
            # leading hadr. hem. AK4),==1 high pT mu, no further leptons passing criteria, 
            # large MET, large lept. W pT

            passSel[sys], iSel[sys], recodR_AK8_hadAK4b[sys], recoAK4bhadjets[sys] = self.WtopSelection( False, event, recoMuons, recoElectrons, recoAK4jets, recoAK4bjets, recoAK8jets[sys], MET, sys)
            
            #print (recodR_AK8_hadAK4b,recoAK4bhadjets)
            # fill in relevant arrays like had hem AK4 b's, objects like 
            # leptWpT, dR(AK8[sys], lead had. hem. AK4), etc., 
            # here, if general selection satisfied, else return 
            # do whole rigmarole of checking for nmuon, etc., as above commented lines do
            # in the WtopSelection method
            ################################################################################################################################
            
            


        ######################### Weights ############################

        # all these weights are to plot nominal, control histos below
        # weights for array masks are done in analyzer
        weight=1.
        
        if self.isMC:  
            #### b-tagging Weights #####
            bTagSFs=[]
            w=0.
            if len(recoAK4bjets)>0: 
                w=1.
                bTagSFs = [x.btagSF_deepjet_M for x in recoAK4bjets]
                for i in bTagSFs:
                    w *= i  
                self.btaggingWeight = w 
            else:
                self.btaggingWeight = 0.
            ####################################################

            #### Lepton Weights (also reconstruct leptW for control histo plots)####
            #if self.isMC:
            if len(recoMuons)>0:
                recoleptW = [recoMuons[0].p4()+MET ]
                leptonSFs = self.getleptonSF( "muon", recoMuons[0] )
                self.leptonWeight = np.prod([i[0] for i in leptonSFs])
            else: 
                recoleptW = []
                leptonSFs = [[0, 0, 0,],[0, 0, 0,],[0, 0, 0,],[0, 0, 0,]]
                self.leptonWeight=0.
            ####################################################


            self.puWeight = event.puWeight
            self.l1PreFireWeight = event.L1PreFiringWeight_Nom
            self.evtGenWeight = event.genWeight 

            #### Applying ALL remaining object-related weights ####
            weight =  self.evtGenWeight  * self.puWeight * self.leptonWeight * self.btaggingWeight *  self.l1PreFireWeight #btagweights['_nom'] # * self.topweight            
            #self.totalRecoWeight = weight

        else:
            leptonSFs = [1, 1, 1, 1.]
            recoleptW = [recoMuons[0].p4()+MET ] if len(recoMuons)>0 else []
            #self.leptonWeight = np.prod(leptonSFs)
            #self.btaggingWeight = 1.
            #self.puWeight = 1. 
            #self.totalRecoWeight = 1.#weight
            weight = 1.
            ##############################################################

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
                getattr( self, 'AK8jets_msoftdrop_noSelnoWeight' ).Fill( ijet.msoftdrop_nom)
            
            getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4jets) )
            getattr( self, 'nAK4bjets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
            getattr( self, 'nAK4bhadjets_noSelnoWeight' ).Fill( len(recoAK4bhadjets['_nom']) )
            getattr( self, 'dR_AK8_AK4bhad_noSelnoWeight' ).Fill( recodR_AK8_hadAK4b['_nom'] )

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

            getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )
            for lW in recoleptW:
                getattr( self, 'leptonicW_pt_noSelnoWeight' ).Fill( lW.Pt() )
                getattr( self, 'leptonicW_eta_noSelnoWeight' ).Fill( lW.Eta() )
                getattr( self, 'leptonicW_y_noSelnoWeight' ).Fill( lW.Rapidity() )
                getattr( self, 'leptonicW_phi_noSelnoWeight' ).Fill( lW.Phi() )
                getattr( self, 'leptonicW_mass_noSelnoWeight' ).Fill( lW.M() )
                

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
                getattr( self, 'AK8jets_msoftdrop_noSel' ).Fill( ijet.msoftdrop_nom, weight)

            getattr( self, 'nAK4jets_noSel' ).Fill( len(recoAK4jets), weight )
            getattr( self, 'nAK4bjets_noSel' ).Fill( len(recoAK4bjets) , weight)
            getattr( self, 'nAK4bhadjets_noSelnoWeight' ).Fill( len(recoAK4bhadjets['_nom']) , weight)
            getattr( self, 'dR_AK8_AK4bhad_noSelnoWeight' ).Fill( recodR_AK8_hadAK4b['_nom'] , weight)

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

            getattr( self, 'METPt_noSel' ).Fill( MET.Pt() , weight)
            
            for lW in recoleptW:
                getattr( self, 'leptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                getattr( self, 'leptonicW_eta_noSel' ).Fill( lW.Eta() , weight)
                getattr( self, 'leptonicW_y_noSel' ).Fill( lW.Rapidity() , weight)
                getattr( self, 'leptonicW_phi_noSel' ).Fill( lW.Phi() , weight)
                getattr( self, 'leptonicW_mass_noSel' ).Fill( lW.M() , weight)
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
                    getattr( self, 'AK8jets_msoftdrop'+iSel['_nom'] ).Fill( ijet.msoftdrop_nom, weight)
                
                getattr( self, 'nAK4jets'+iSel['_nom'] ).Fill( len(recoAK4jets) , weight )
                getattr( self, 'nAK4bjets'+iSel['_nom'] ).Fill( len(recoAK4bjets)  , weight)
                getattr( self, 'nAK4bhadjets'+iSel['_nom'] ).Fill( len(recoAK4bhadjets['_nom']) , weight)
                getattr( self, 'dR_AK8_AK4bhad'+iSel['_nom']  ).Fill( recodR_AK8_hadAK4b['_nom'] , weight)

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

                getattr( self, 'METPt'+iSel['_nom'] ).Fill( MET.Pt() , weight)
                
                for lW in recoleptW:
                    getattr( self, 'leptonicW_pt'+iSel['_nom'] ).Fill( lW.Pt(), weight )
                    getattr( self, 'leptonicW_eta'+iSel['_nom'] ).Fill( lW.Eta() , weight)
                    getattr( self, 'leptonicW_y'+iSel['_nom'] ).Fill( lW.Rapidity() , weight)
                    getattr( self, 'leptonicW_phi'+iSel['_nom'] ).Fill( lW.Phi() , weight)
                    getattr( self, 'leptonicW_mass'+iSel['_nom'] ).Fill( lW.M() , weight)
                    #getattr( self, 'leptonicWMT'+iSel['_nom']).Fill( lW.Mt() , weight)


                #### Checking nominal W-cand selection with weights
                if recodR_AK8_hadAK4b['_nom']>=0.8 and recoAK8jets['_nom'][0].pt>self.minLeadAK8JetPtW and self.minSDMassW<recoAK8jets['_nom'][0].msoftdrop_nom<=self.maxSDMassW and self.evtSelection.startswith('_Wtop'):
                    #nominally the above is the cut I will be placing on _nom/_jes*/_jer* versions of these quantities when separating W vs. top during offline histogramming
                    # basic reco histos
                    getattr( self, 'nPVs'+'_WSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
                    getattr( self, 'nleps'+'_WSel' ).Fill( nleptons , weight)
                    
                    for imuon in recoMuons:
                        getattr( self, 'muons_pt'+'_WSel' ).Fill( imuon.pt , weight)
                        getattr( self, 'muons_eta'+'_WSel' ).Fill( imuon.eta , weight)
                        getattr( self, 'muons_y'+'_WSel' ).Fill( imuon.p4().Rapidity() , weight)
                        getattr( self, 'muons_phi'+'_WSel' ).Fill( imuon.phi , weight)

                    for iele in recoElectrons:
                        getattr( self, 'eles_pt'+'_WSel' ).Fill( iele.pt , weight)
                        getattr( self, 'eles_eta'+'_WSel' ).Fill( iele.eta , weight)
                        getattr( self, 'eles_y'+'_WSel' ).Fill( iele.p4().Rapidity() , weight)
                        getattr( self, 'eles_phi'+'_WSel' ).Fill( iele.phi , weight)

                    getattr( self, 'nAK8jets'+'_WSel' ).Fill( len(recoAK8jets['_nom']) , weight)
                    getattr( self, 'HT'+'_WSel' ).Fill( AK8HT , weight)

                    for ijet in recoAK8jets['_nom']:
                        getattr( self, 'AK8jets_pt'+'_WSel' ).Fill( ijet.pt_nom , weight)
                        getattr( self, 'AK8jets_eta'+'_WSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK8jets_y'+'_WSel' ).Fill( ijet.rapidity , weight)
                        getattr( self, 'AK8jets_phi'+'_WSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK8jets_mass'+'_WSel' ).Fill( ijet.mass_nom, weight)
                        getattr( self, 'AK8jets_msoftdrop'+'_WSel' ).Fill( ijet.msoftdrop_nom, weight)
                    
                    getattr( self, 'nAK4jets'+'_WSel').Fill( len(recoAK4jets) , weight )
                    getattr( self, 'nAK4bjets'+'_WSel' ).Fill( len(recoAK4bjets)  , weight)
                    getattr( self, 'nAK4bhadjets'+'_WSel' ).Fill( len(recoAK4bhadjets['_nom']) , weight)
                    getattr( self, 'dR_AK8_AK4bhad'+'_WSel'  ).Fill( recodR_AK8_hadAK4b['_nom'] , weight)

                    for ijet in recoAK4jets:
                        getattr( self, 'AK4jets_pt'+'_WSel' ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4jets_eta'+'_WSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4jets_y'+'_WSel' ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4jets_phi'+'_WSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4jets_mass'+'_WSel' ).Fill( ijet.mass , weight)
                        
                    for ijet in recoAK4bjets:
                        getattr( self, 'AK4btaggedjets_pt'+'_WSel' ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4btaggedjets_eta'+'_WSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4btaggedjets_y'+'_WSel' ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4btaggedjets_phi'+'_WSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4btaggedjets_mass'+'_WSel' ).Fill( ijet.mass , weight)

                    getattr( self, 'METPt'+'_WSel' ).Fill( MET.Pt() , weight)
                    
                    for lW in recoleptW:
                        getattr( self, 'leptonicW_pt'+'_WSel' ).Fill( lW.Pt(), weight )
                        getattr( self, 'leptonicW_eta'+'_WSel' ).Fill( lW.Eta() , weight)
                        getattr( self, 'leptonicW_y'+'_WSel' ).Fill( lW.Rapidity() , weight)
                        getattr( self, 'leptonicW_phi'+'_WSel' ).Fill( lW.Phi() , weight)
                        getattr( self, 'leptonicW_mass'+'_WSel' ).Fill( lW.M() , weight)
                        #getattr( self, 'leptonicWMT'+'_WSel').Fill( lW.Mt() , weight)


                #### Checking nominal top-cand selection with weights                
                elif recodR_AK8_hadAK4b['_nom']<0.8 and recoAK8jets['_nom'][0].pt_nom>self.minLeadAK8JetPtTop and self.minSDMassTop<(recoAK8jets['_nom'][0].msoftdrop_nom/recoAK8jets['_nom'][0].msoftdrop_corr_PUPPI) and self.evtSelection.startswith('_Wtop'):#<=self.maxSDMassTop leaving this upper mass cut open-ended
                    # basic reco histos
                    getattr( self, 'nPVs'+'_topSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
                    getattr( self, 'nleps'+'_topSel' ).Fill( nleptons , weight)
                    
                    for imuon in recoMuons:
                        getattr( self, 'muons_pt'+'_topSel' ).Fill( imuon.pt , weight)
                        getattr( self, 'muons_eta'+'_topSel' ).Fill( imuon.eta , weight)
                        getattr( self, 'muons_y'+'_topSel' ).Fill( imuon.p4().Rapidity() , weight)
                        getattr( self, 'muons_phi'+'_topSel' ).Fill( imuon.phi , weight)

                    for iele in recoElectrons:
                        getattr( self, 'eles_pt'+'_topSel' ).Fill( iele.pt , weight)
                        getattr( self, 'eles_eta'+'_topSel' ).Fill( iele.eta , weight)
                        getattr( self, 'eles_y'+'_topSel' ).Fill( iele.p4().Rapidity() , weight)
                        getattr( self, 'eles_phi'+'_topSel' ).Fill( iele.phi , weight)

                    getattr( self, 'nAK8jets'+'_topSel' ).Fill( len(recoAK8jets['_nom']) , weight)
                    getattr( self, 'HT'+'_topSel' ).Fill( AK8HT , weight)

                    for ijet in recoAK8jets['_nom']:
                        getattr( self, 'AK8jets_pt'+'_topSel' ).Fill( ijet.pt_nom , weight)
                        getattr( self, 'AK8jets_eta'+'_topSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK8jets_y'+'_topSel' ).Fill( ijet.rapidity , weight)
                        getattr( self, 'AK8jets_phi'+'_topSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK8jets_mass'+'_topSel' ).Fill( ijet.mass_nom, weight)
                        getattr( self, 'AK8jets_msoftdrop'+'_topSel' ).Fill( ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI, weight) #since as per Alejandro's advice, the m_SD PUPPI corrections are only derived for the W
                    
                    getattr( self, 'nAK4jets'+'_topSel'  ).Fill( len(recoAK4jets) , weight )
                    getattr( self, 'nAK4bjets'+'_topSel'  ).Fill( len(recoAK4bjets)  , weight)
                    getattr( self, 'nAK4bhadjets'+'_topSel' ).Fill( len(recoAK4bhadjets['_nom']) , weight)
                    getattr( self, 'dR_AK8_AK4bhad'+'_topSel'  ).Fill( recodR_AK8_hadAK4b['_nom'] , weight)

                    for ijet in recoAK4jets:
                        getattr( self, 'AK4jets_pt'+'_topSel' ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4jets_eta'+'_topSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4jets_y'+'_topSel' ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4jets_phi'+'_topSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4jets_mass'+'_topSel' ).Fill( ijet.mass , weight)
                        
                    for ijet in recoAK4bjets:
                        getattr( self, 'AK4btaggedjets_pt'+'_topSel' ).Fill( ijet.pt , weight)
                        getattr( self, 'AK4btaggedjets_eta'+'_topSel' ).Fill( ijet.eta , weight)
                        getattr( self, 'AK4btaggedjets_y'+'_topSel' ).Fill( ijet.p4().Rapidity() , weight)
                        getattr( self, 'AK4btaggedjets_phi'+'_topSel' ).Fill( ijet.phi , weight)
                        getattr( self, 'AK4btaggedjets_mass'+'_topSel' ).Fill( ijet.mass , weight)

                    getattr( self, 'METPt'+'_topSel' ).Fill( MET.Pt() , weight)
                    
                    for lW in recoleptW:
                        getattr( self, 'leptonicW_pt'+'_topSel' ).Fill( lW.Pt(), weight )
                        getattr( self, 'leptonicW_eta'+'_topSel' ).Fill( lW.Eta() , weight)
                        getattr( self, 'leptonicW_y'+'_topSel' ).Fill( lW.Rapidity() , weight)
                        getattr( self, 'leptonicW_phi'+'_topSel' ).Fill( lW.Phi() , weight)
                        getattr( self, 'leptonicW_mass'+'_topSel' ).Fill( lW.M() , weight)
                        #getattr( self, 'leptonicWMT'+'_topSel').Fill( lW.Mt() , weight)


        return passSel, iSel, recodR_AK8_hadAK4b, recoMuons, recoAK4bjets, recoAK4bhadjets, recoleptW, recoAK8jets, MET

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
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        
        ngenLeptons = len(genMuons)+len(genElectrons)
        ##################################################

        #### MET (not sure if needed)
        genMET = ROOT.TLorentzVector()
        genMET.SetPtEtaPhiM(genmet.pt, 0., genmet.phi, 0)
        ##################################################

        #### Basic AK4 jet selection
        genAK4jets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta]# and abs(x.hadronFlavour)==5 ] 
        genAK4jets.sort(key=lambda x:x.pt,reverse=True)
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
        
        #### reconstruct leptonic W if possible
        if len(genMuons)>0:
            genleptW = [genMuons[0].p4()+genMET]
        else:
            genleptW = []#genMuons[0].p4()+genMET]
        
        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        passSel, iSel, gendR_AK8_hadAK4b, genAK4bhadjets = self.WtopSelection( True, event, genMuons, genElectrons, genAK4jets, genAK4bjets, genAK8jets, genMET, '' )

        #### Weight
        weight = event.genWeight
        #self.totalGenWeight = weight


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
            
            getattr( self, 'ngenAK4jets_noSel' ).Fill( len(genAK4jets), weight )
            getattr( self, 'ngenAK4bjets_noSel' ).Fill( len(genAK4bjets), weight )
            getattr( self, 'ngenAK4bhadjets_noSel' ).Fill( len(genAK4bhadjets), weight )
            getattr( self, 'dR_genAK8_genAK4bhad_noSel' ).Fill( gendR_AK8_hadAK4b, weight )
            
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
            
            getattr( self, 'genMETPt_noSel' ).Fill( genMET.Pt(), weight )

            for lW in genleptW:
                getattr( self, 'genleptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                getattr( self, 'genleptonicW_eta_noSel' ).Fill( lW.Eta(), weight )
                getattr( self, 'genleptonicW_y_noSel' ).Fill( lW.Rapidity(), weight )
                getattr( self, 'genleptonicW_phi_noSel' ).Fill( lW.Phi(), weight )
                getattr( self, 'genleptonicW_mass_noSel' ).Fill( lW.M(), weight )
                #getattr( self, 'genleptonicWMT_noSel' ).Fill( lW.Mt(), weight )
            
            
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
                getattr( self, 'ngenAK4bhadjets'+iSel ).Fill( len(genAK4bhadjets), weight )
                getattr( self, 'dR_genAK8_genAK4bhad'+iSel).Fill( gendR_AK8_hadAK4b, weight )

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

                getattr( self, 'genMETPt'+iSel ).Fill( genMET.Pt(), weight )
                
                for lW in genleptW:
                    getattr( self, 'genleptonicW_pt'+iSel  ).Fill( lW.Pt(), weight )
                    getattr( self, 'genleptonicW_eta'+iSel  ).Fill( lW.Eta(), weight )
                    getattr( self, 'genleptonicW_y'+iSel  ).Fill( lW.Rapidity(), weight )
                    getattr( self, 'genleptonicW_phi'+iSel  ).Fill( lW.Phi(), weight )
                    getattr( self, 'genleptonicW_mass'+iSel  ).Fill( lW.M(), weight )
                    #getattr( self, 'genleptonicWMT'+iSel ).Fill( lW.Mt(), weight )
                
                #### Checking nominal W-cand selection with weights
                if gendR_AK8_hadAK4b>=0.8 and genAK8jets[0].pt>self.minLeadAK8JetPtW and self.minSDMassW<genAK8jets[0].msoftdrop<=self.maxSDMassW and self.evtSelection.startswith('_Wtop'):

                    getattr( self, 'ngenleps'+'_WSel' ).Fill( len(genMuons)+len(genElectrons), weight )
                    for imuon in genMuons:
                        getattr( self, 'genmuons_pt'+'_WSel' ).Fill( imuon.pt, weight )
                        getattr( self, 'genmuons_eta'+'_WSel' ).Fill( imuon.eta, weight )
                        getattr( self, 'genmuons_y'+'_WSel' ).Fill( imuon.p4().Rapidity(), weight )
                        getattr( self, 'genmuons_phi'+'_WSel' ).Fill( imuon.phi, weight )

                    for iele in genElectrons:
                        getattr( self, 'geneles_pt'+'_WSel' ).Fill( iele.pt, weight )
                        getattr( self, 'geneles_eta'+'_WSel' ).Fill( iele.eta, weight )
                        getattr( self, 'geneles_y'+'_WSel' ).Fill( iele.p4().Rapidity(), weight )
                        getattr( self, 'geneles_phi'+'_WSel' ).Fill( iele.phi, weight )

                    getattr( self, 'ngenAK8jets'+'_WSel' ).Fill( len(genAK8jets), weight )
                    getattr( self, 'genHT'+'_WSel' ).Fill( genAK8HT, weight )

                    for ijet in genAK8jets:
                        getattr( self, 'AK8genjets_pt'+'_WSel' ).Fill( ijet.pt, weight )
                        getattr( self, 'AK8genjets_eta'+'_WSel' ).Fill( ijet.eta, weight )
                        getattr( self, 'AK8genjets_y'+'_WSel' ).Fill( ijet.rapidity, weight )
                        getattr( self, 'AK8genjets_phi'+'_WSel' ).Fill( ijet.phi, weight )
                        getattr( self, 'AK8genjets_mass'+'_WSel' ).Fill( ijet.mass, weight )
                        getattr( self, 'AK8genjets_msoftdrop'+'_WSel' ).Fill( ijet.msoftdrop, weight )

                    getattr( self, 'ngenAK4jets'+'_WSel'  ).Fill( len(genAK4jets), weight )
                    getattr( self, 'ngenAK4bjets'+'_WSel' ).Fill( len(genAK4bjets), weight )
                    getattr( self, 'ngenAK4bhadjets'+'_WSel' ).Fill( len(genAK4bhadjets), weight )
                    getattr( self, 'dR_genAK8_genAK4bhad'+'_WSel').Fill( gendR_AK8_hadAK4b, weight )

                    for ijet in genAK4jets:
                        getattr( self, 'AK4genjets_pt'+'_WSel'  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4genjets_eta'+'_WSel'  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4genjets_y'+'_WSel'  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4genjets_phi'+'_WSel'  ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4genjets_mass'+'_WSel'  ).Fill( ijet.mass, weight )

                    for ijet in genAK4bjets:
                        getattr( self, 'AK4bmatchedgenjets_pt'+'_WSel'  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4bmatchedgenjets_eta'+'_WSel'  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4bmatchedgenjets_y'+'_WSel'  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4bmatchedgenjets_phi'+'_WSel' ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4bmatchedgenjets_mass'+'_WSel' ).Fill( ijet.mass, weight )

                    getattr( self, 'genMETPt'+'_WSel' ).Fill( genMET.Pt(), weight )
                    
                    for lW in genleptW:
                        getattr( self, 'genleptonicW_pt'+'_WSel'  ).Fill( lW.Pt(), weight )
                        getattr( self, 'genleptonicW_eta'+'_WSel'  ).Fill( lW.Eta(), weight )
                        getattr( self, 'genleptonicW_y'+'_WSel'  ).Fill( lW.Rapidity(), weight )
                        getattr( self, 'genleptonicW_phi'+'_WSel'  ).Fill( lW.Phi(), weight )
                        getattr( self, 'genleptonicW_mass'+'_WSel'  ).Fill( lW.M(), weight )
                        #getattr( self, 'genleptonicWMT'+'_WSel' ).Fill( lW.Mt(), weight )
                    
                #### Checking nominal top-cand selection with weights
                elif gendR_AK8_hadAK4b<0.8 and genAK8jets[0].pt>self.minLeadAK8JetPtTop and self.minSDMassTop<genAK8jets[0].msoftdrop and self.evtSelection.startswith('_Wtop'):#<=self.maxSDMassTop:

                    getattr( self, 'ngenleps'+'_topSel' ).Fill( len(genMuons)+len(genElectrons), weight )
                    for imuon in genMuons:
                        getattr( self, 'genmuons_pt'+'_topSel' ).Fill( imuon.pt, weight )
                        getattr( self, 'genmuons_eta'+'_topSel' ).Fill( imuon.eta, weight )
                        getattr( self, 'genmuons_y'+'_topSel' ).Fill( imuon.p4().Rapidity(), weight )
                        getattr( self, 'genmuons_phi'+'_topSel' ).Fill( imuon.phi, weight )

                    for iele in genElectrons:
                        getattr( self, 'geneles_pt'+'_topSel' ).Fill( iele.pt, weight )
                        getattr( self, 'geneles_eta'+'_topSel' ).Fill( iele.eta, weight )
                        getattr( self, 'geneles_y'+'_topSel' ).Fill( iele.p4().Rapidity(), weight )
                        getattr( self, 'geneles_phi'+'_topSel' ).Fill( iele.phi, weight )

                    getattr( self, 'ngenAK8jets'+'_topSel' ).Fill( len(genAK8jets), weight )
                    getattr( self, 'genHT'+'_topSel' ).Fill( genAK8HT, weight )

                    for ijet in genAK8jets:
                        getattr( self, 'AK8genjets_pt'+'_topSel' ).Fill( ijet.pt, weight )
                        getattr( self, 'AK8genjets_eta'+'_topSel' ).Fill( ijet.eta, weight )
                        getattr( self, 'AK8genjets_y'+'_topSel' ).Fill( ijet.rapidity, weight )
                        getattr( self, 'AK8genjets_phi'+'_topSel' ).Fill( ijet.phi, weight )
                        getattr( self, 'AK8genjets_mass'+'_topSel' ).Fill( ijet.mass, weight )
                        getattr( self, 'AK8genjets_msoftdrop'+'_topSel' ).Fill( ijet.msoftdrop, weight )

                    getattr( self, 'ngenAK4jets'+'_topSel' ).Fill( len(genAK4jets), weight )
                    getattr( self, 'ngenAK4bjets'+'_topSel').Fill( len(genAK4bjets), weight )
                    getattr( self, 'ngenAK4bhadjets'+'_topSel' ).Fill( len(genAK4bhadjets), weight )
                    getattr( self, 'dR_genAK8_genAK4bhad'+'_topSel').Fill( gendR_AK8_hadAK4b, weight )

                    for ijet in genAK4jets:
                        getattr( self, 'AK4genjets_pt'+'_topSel'  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4genjets_eta'+'_topSel'  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4genjets_y'+'_topSel'  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4genjets_phi'+'_topSel'  ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4genjets_mass'+'_topSel'  ).Fill( ijet.mass, weight )

                    for ijet in genAK4bjets:
                        getattr( self, 'AK4bmatchedgenjets_pt'+'_topSel'  ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4bmatchedgenjets_eta'+'_topSel'  ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4bmatchedgenjets_y'+'_topSel'  ).Fill( ijet.p4().Rapidity(), weight )
                        getattr( self, 'AK4bmatchedgenjets_phi'+'_topSel' ).Fill( ijet.phi, weight )
                        getattr( self, 'AK4bmatchedgenjets_mass'+'_topSel' ).Fill( ijet.mass, weight )

                    getattr( self, 'genMETPt'+'_topSel' ).Fill( genMET.Pt(), weight )
                    
                    for lW in genleptW:
                        getattr( self, 'genleptonicW_pt'+'_topSel'  ).Fill( lW.Pt(), weight )
                        getattr( self, 'genleptonicW_eta'+'_topSel'  ).Fill( lW.Eta(), weight )
                        getattr( self, 'genleptonicW_y'+'_topSel'  ).Fill( lW.Rapidity(), weight )
                        getattr( self, 'genleptonicW_phi'+'_topSel'  ).Fill( lW.Phi(), weight )
                        getattr( self, 'genleptonicW_mass'+'_topSel'  ).Fill( lW.M(), weight )
                        #getattr( self, 'genleptonicWMT'+'_topSel' ).Fill( lW.Mt(), weight )
                    


        return passSel, iSel, gendR_AK8_hadAK4b, genMuons, genAK4bjets, genAK4bhadjets, genleptW, genAK8jets, genMET

    #############################################################################
    def WtopSelection( self, isGen, event, muons, electrons, AK4jets, AK4bjets, AK8jets, MET, ptLabel):
    
        if (len(muons)==1) and (len(electrons)==0) and (len(AK8jets)>=1) and (len(AK4bjets)>=1) and (len(AK4jets)>3) and MET.Pt()>self.METCutWtop:                
            
            # reconstruct leptonic hemisphere objects for selection here 
            # safe, since the contributing objects are checked to exist above

            leptW=muons[0].p4()+MET
            
            AK4lepjets = [x for x in AK4jets if (self.DrRapPhi(x.p4(), muons[0].p4() )>0.4 and self.DrRapPhi(x.p4(),muons[0].p4())<1.6)]# and x.btagDeepFlavB > self.minBDisc] 
            AK4lepjets.sort(key=lambda x:x.pt,reverse=True) 

            if not isGen: 
                AK4lepbjets = [x for x in AK4lepjets if x.btagDeepFlavB > self.minBDisc] 
                AK4lepbjets.sort(key=lambda x:x.pt,reverse=True) 
            else: 
                AK4lepbjets = [x for x in AK4lepjets if abs(x.hadronFlavour)==5] 
                AK4lepbjets.sort(key=lambda x:x.pt,reverse=True) 

            ######### check if possible, if AK8jets[0] is really the hadronic top/w or the leptonic top... 

            if len(AK4lepjets)>0 and leptW.Pt()>self.minLeptonicWPt and AK8jets[0].p4().DeltaPhi(muons[0].p4())>2. and self.DrRapPhi(AK8jets[0].p4(),AK4lepjets[0].p4())>1.6:
                # to rconstruct the hadronic hemisphere: ensure separation of lead AK8 jet in the event from the leptonic hemisphere mu and
                # the leading leptonic hemisphere AK4 (not requiring this to be a b-tagged AK4 currently),  
                # additionally, ensure hadronic hemisphere's topology is boosted, with the leptonic W pT cut 
                # (ie, the leptonic top's decay products )
                leadJetpT = getattr( AK8jets[0], 'pt'+ptLabel )
                leadJetmass = AK8jets[0].msoftdrop if isGen else AK8jets[0].msoftdrop_nom #should this be: getattr( AK8jets[0], 'mosftdrop'+ptLabel ) ??? 

                # find hadronic hemisphere AK4's within dR<1.6 of AK8 and btagged
                # if exactly one b-tagged AK4 close to AK8 keep event
                if not isGen: 
                    AK4hadbjets = [x for x in AK4jets if self.DrRapPhi(x.p4(),AK8jets[0].p4())<1.6 and x.btagDeepFlavB > self.minBDisc]
                    AK4hadbjets.sort(key=lambda x:x.pt,reverse=True) 
                else: 
                    AK4hadbjets = [x for x in AK4jets if self.DrRapPhi(x.p4(),AK8jets[0].p4())<1.6 and abs(x.hadronFlavour)==5]                    
                    AK4hadbjets.sort(key=lambda x:x.pt,reverse=True) 

                if len(AK4hadbjets)>0: 
                    dR_leadJet_leadhadAK4 = self.DrRapPhi(AK4hadbjets[0].p4(),AK8jets[0].p4()) 
                else: 
                    dR_leadJet_leadhadAK4 = -929
                    return False, None, dR_leadJet_leadhadAK4, []
                #print(dR_leadJet_leadhadAK4)
                if self.evtSelection.startswith('_Wtop'):
                    if leadJetmass>self.minleadJetMass and leadJetpT>self.minLeadAK8JetPtW and (len(AK4hadbjets)==1):
                        # apply a lowish inv. mass and pT cut to get rid of extra events, 
                        # require exactly 1-btagged had hem AK4 and keep info relevant to both W and top
                        # store deltaR of leading (b-tagged) hadronic hemisphere AK4, if any 
                        # keeping the chosen had hem b-jets to return to do plots with them
                        
                        return True, '_WtopSel', dR_leadJet_leadhadAK4, AK4hadbjets
                    else: 
                        return False, None, dR_leadJet_leadhadAK4, AK4hadbjets
                elif 'WSel' in self.evtSelection: 
                    if (((leadJetmass<self.maxSDMassW) and (leadJetmass>=self.minSDMassW)) and (leadJetpT>self.minLeadAK8JetPtW)) and (len(AK4hadbjets)==1) and (dR_leadJet_leadhadAK4>0.8 and dR_leadJet_leadhadAK4<1.6): 
                        #print (leptWpT.Mt(),isGen)                           
                        return True, '_WSel', dR_leadJet_leadhadAK4, AK4hadbjets 
                    else: return False, None, dR_leadJet_leadhadAK4, AK4hadbjets

                elif '_topSel' in self.evtSelection:
                    if ((leadJetmass/ (1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) >= self.minSDMassTop) and (leadJetmass/(1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop) and (leadJetpT> self.minLeadAK8JetPtTop)) and (len(AK4hadbjets)==1) and (dR_leadJet_leadhadAK4<=0.8): # 
                        return True, '_topSel', dR_leadJet_leadhadAK4, AK4hadbjets
                    else: return False, None, dR_leadJet_leadhadAK4, AK4hadbjets   
                
            #returning no had hem b jet info in these following two return statements since as per the requirements in these if blocks, 
            #the hadronic hemisphere is not entirely reconstructable as per the constraints above
            else: return False, None, -929, [] 
        else: return False, None, -929, []
        

       

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
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            tp = tp * p.puppiWeight if not isGen else tp
            #except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        #print ("pushing back candidates")
        for x in CandsPUPPIweightedVec:
            if self.DrRapPhi( AK8jet.p4(), x ) < 0.8: constituents.push_back(x)
        #print ("pushed back candidates")

        #### Computing n-subjetiness basis from PF PUPPI constituents
        nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
        nsub1 = self.nSub1.getTau( self.maxTau, constituents )
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
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]

        try: ak8jet['tau21_exkT'] = nsub1[1]/nsub1[0]
        except ZeroDivisionError: ak8jet['tau21_exkT'] = -1
        try: ak8jet['tau32_exkT'] = nsub1[2]/nsub1[1]
        except ZeroDivisionError: ak8jet['tau32_exkT'] = -1

        if isGen: #to add in softdrop mass as a variable to the selected (accep)gen jet branches
            sd_AK8jets = self.sd.result( constituents)
            if len(sd_AK8jets)>0: #stupidly, in some rare cases, this will not be true (with ptmin>=170) leading to errors, hence, switching to ptmin=0 in function calls (but, since not sure this is error, free I use this if-else block)
                ak8jet['msoftdrop'] = sd_AK8jets[0].m()
            else: ak8jet['msoftdrop'] = -1. 
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
    def fillAK8Branches( self, event, jetLabel, jetInfo, dummy=False, sys='_nom' ): 
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
            
            self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
            c=0
            for i,iJ in jetInfo.items():
                if c==0:
                    if 'reco' in jetLabel.lower():

                        self.out.fillBranch(jetLabel+"_msoftdrop",  getattr(iJ['jet'], 'msoftdrop'+sys) )
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI",  getattr(iJ['jet'], 'msoftdrop_corr_PUPPI') )
                        self.out.fillBranch(jetLabel+"_pt",  getattr(iJ['jet'], 'pt'+sys)  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass",  getattr(iJ['jet'], 'mass'+sys)  )
                        
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
                        self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN), iJ['0p5'+str(tauN)]  )
                        self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN), iJ['1'+str(tauN)]  )
                        self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN), iJ['2'+str(tauN)]  )
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
                        self.out.fillBranch(jetLabel+"_msoftdrop", dummyFill  )
                        self.out.fillBranch(jetLabel+"_msoftdrop_corr_PUPPI", dummyFill  )
                        self.out.fillBranch(jetLabel+"_pt", dummyFill  )#, 'pt'+jetLabel.split('Jets')[1] )  )
                        self.out.fillBranch(jetLabel+"_mass", dummyFill  )

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
                        self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN),  dummyFill  )
                        self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN),  dummyFill  )
                        self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN),  dummyFill  )
                c+=1


    #############################################################################
    def fillOtherBranches( self, event, objectLabel, objectInfo, lenVar): 
        #use dummy to fill dummy values into branches for pass/nonpass (gen/reco)sel to maintain correspondence between events in different branches
        # in principle since we don't do true/accep reco/gen discrimination between these reco/gen objects, no need to use dummy=True ever when calling this
        # functionality removed, hence
        #dummyFill=-929
        #### Filling branch with passAK8jet info after selection
        #print (objectLabel+"_pt", [ getattr( iJ['jet'], 'pt'+objectLabel.split('Jets')[1] ) )
        
            
        self.out.fillBranch( 'n'+objectLabel, lenVar )
    
        if 'ak4' in objectLabel.lower() or 'mu' in objectLabel.lower():
            self.out.fillBranch(objectLabel+"_pt",  objectInfo.pt  )
            self.out.fillBranch(objectLabel+"_eta",  objectInfo.eta  )
            self.out.fillBranch(objectLabel+"_y",  objectInfo.p4().Rapidity()  )
            self.out.fillBranch(objectLabel+"_phi",  objectInfo.phi  )
            self.out.fillBranch(objectLabel+"_mass",  objectInfo.mass  )
        else:
            #for MET and leptonic W
            self.out.fillBranch(objectLabel+"_pt",  objectInfo.Pt()  )
            self.out.fillBranch(objectLabel+"_phi",  objectInfo.Phi()  )
            if not 'met' in objectLabel.lower(): #for leptonic W
                self.out.fillBranch(objectLabel+"_eta",  objectInfo.Eta()  )
                self.out.fillBranch(objectLabel+"_y",  objectInfo.Rapidity()  )
                self.out.fillBranch(objectLabel+"_mass",  objectInfo.M()  )
            

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

    ##################### Helpful functions
    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        elif hasattr( c, "Perp"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        else:
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.perp(), c.eta(), c.phi(), c.m() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )

    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret
    
    #############################################################################
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)

    #############################################################################
    '''
    
    def mindRandIndex(self, recoAK4Coll, muon):
        mindR=999
        c=0
        for x in recoAK4Coll:
            #if x.p4().DeltaR(muon)>0.4:
            if x.p4().DeltaR(muon)<mindR:
                mindR=x.p4().DeltaR(muon)
            #else:
            #    if (x.p4()-muon).DeltaR(muon)<mindR:
            #        mindR=(x.p4()-muon).DeltaR(muon)
            c+=1
        return mindR, c-1

    def getpTRel(self, event, recoAK4, muon): 
        corrCollection=list(Collection(event, 'CorrT1METJet' ))

        
        if muon.DeltaR(recoAK4)<=0.4:
            recoAK4=recoAK4-muon

        p_mu = ROOT.TVector3(muon.Px(),muon.Py(),muon.Pz())#muon.Vect()#
        p_j = ROOT.TVector3(recoAK4.Px(),recoAK4.Py(),recoAK4.Pz())#recoAK4.Vect()#
        
        pTrel = ((p_mu.Cross(p_j)).Mag()/p_j.Mag()) if not p_j.Mag()<=0. else 0.#p_mu.Dot(p_mu) - (p_mu.Dot(p_j)/p_j.Mag())**2.
        
        if pTrel<0.: return 0.
        else: return pTrel#ROOT.TMath.Sqrt(pTrel2)

    def getBTagWeight(self, nBTagged=0, jet_SFs=[0]): # Now redundant to allow for more than just two btagged jets per event
        bTagWeight=0
        if len(jet_SFs)>2 or nBTagged>2:
            print "Error, only leading and subleading AK4 jets are considered: # of btagged jets cannot exceed 2"
        elif nBTagged>len(jet_SFs):
            print "Number of b-tagged jets cannot be greater than number of them for which SFs are provided!"
            return 0
        elif nBTagged==0 and len(jet_SFs)==0: return 1

        if len(jet_SFs)==1:
            SF = jet_SFs[0]

            for i in range(0,2):
                if i!=nBTagged: continue
                bTagWeight+=pow(SF,i)*pow(1-SF,1-i)

        elif len(jet_SFs)==2:
            SF1, SF2 = jet_SFs[0], jet_SFs[1]
            for i in range(0,2):
                for j in range(0,2):
                    if (i+j)!=nBTagged: continue
                    bTagWeight+=pow(SF1,i)*pow(1-SF1,1-i)*pow(SF2,j)*pow(1-SF2,1-j)

        return bTagWeight

    '''


    #############################################################################
    '''
    def getBTagWeight(self, nBTagged=0, jet_SFs=[0]): # Now deprecated, to allow for more than just two btagged jets per event
        bTagWeight=0
        if len(jet_SFs)>2 or nBTagged>2:
            print "Error, only leading and subleading AK4 jets are considered: # of btagged jets cannot exceed 2"
        elif nBTagged>len(jet_SFs):
            print "Number of b-tagged jets cannot be greater than number of them for which SFs are provided!"
            return 0
        elif nBTagged==0 and len(jet_SFs)==0: return 1

        if len(jet_SFs)==1:
            SF = jet_SFs[0]

            for i in range(0,2):
                if i!=nBTagged: continue
                bTagWeight+=pow(SF,i)*pow(1-SF,1-i)

        elif len(jet_SFs)==2:
            SF1, SF2 = jet_SFs[0], jet_SFs[1]
            for i in range(0,2):
                for j in range(0,2):
                    if (i+j)!=nBTagged: continue
                    bTagWeight+=pow(SF1,i)*pow(1-SF1,1-i)*pow(SF2,j)*pow(1-SF2,1-j)

        return bTagWeight

    '''
    ###############
    # useful for the top pT reweighting that's no longer being used
    '''    
    def hasFlags(self, topcand, *flags):
        """Check if one or more status flags are set
        Parameters
        ----------
            flags : str or list
                A list of flags that are required to be set true. If the first argument
                is a list, it is expanded and subsequent arguments ignored.
                Possible flags are enumerated in the `FLAGS` attribute
        Returns a boolean array
        """
        if not len(flags):
            raise ValueError("No flags specified")
        elif isinstance(flags[0], list):
            flags = flags[0]
        mask = 0
        for flag in flags:
            mask |= 1 << self.FLAGS.index(flag)
        return (topcand.statusFlags & mask) == mask
    '''

###################### Different gen functions
def getDaughters(GenParticle,gp):
    ret = []
    tmpListGenParticles = list(GenParticle)
    for part in GenParticle:
        if part != gp:
            if part.genPartIdxMother == tmpListGenParticles.index(gp):
                ret.append(part)
    return ret

def GenBquarkFromTop(GenParticle):
    ret = []
    for gP in GenParticle:
        if (abs(gP.pdgId) == 5) and gP.genPartIdxMother>0:
            mom = GenParticle[gP.genPartIdxMother]
            if abs(mom.pdgId)  == 6:
                dauids = [abs(dau.pdgId) for dau in getDaughters(GenParticle,mom)]
                if 6 not in dauids:
                    ret.append(gP)
    return ret

def GenQuarkFromW(GenParticle):
    ret = []
    for gP in GenParticle:
        if gP.genPartIdxMother>0:
            if (abs(gP.pdgId) <= 5) and (abs(GenParticle[gP.genPartIdxMother].pdgId) == 24):
                ret.append(gP)
        if (abs(gP.pdgId) <= 5) and gP.genPartIdxMother>0:
            mom = GenParticle[gP.genPartIdxMother]
            if abs(mom.pdgId)  == 24:
                dauids = [abs(dau.pdgId) for dau in getDaughters(GenParticle,mom)]
                if 24 not in dauids:
                    ret.append(gP)



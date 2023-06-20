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

    def __init__(self, sysSource=[], year='2017', isMC=True, onlyUnc='', onlyTrees=False, isSigMC=False):
        self.writeHistFile=True
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.isSigMC = isSigMC
        self.onlyTrees = onlyTrees
        self.runSDVariables = False

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

        ### Kinematics Cuts AK8Jets ###
        self.minLeadAK8JetPtDijet = 200.
        self.minAK8JetPt = 170.  ### this is the basic minimum, not the final
        #self.maxJetAK8Eta = 2.4
        self.maxLeadAK8JetRap = 1.7
        self.diffPt = [ 200, 10000 ] # 350, 500, 750, 1000,

        ### Kinenatic Cuts Muons ###
        self.minLooseMuonPt = 30.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 20.
        self.maxElectronEta = 2.5
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
        self.triggerTable[ 'AK8PFJet80' ] = {    #### from the list below, only the first two numbers (trigger turn on) are used. The others were a test - Alejandro
                    '2016_preVFP' : [ 200,   242 ],
                    '2016' : [ 200,   245 ],
                    '2017' : [ 200,   257, 50510.97, 6727.43, 6854.43, 23381.99, 38371.17 ],
                    '2018' : [ 200,   267, 17046.53, 34000.07, 33988.86, 33998.78 ],
                    }
        self.triggerTable[ 'AK8PFJet140' ] = {
                    '2016_preVFP' : [ 242,   308 ],
                    '2016' : [ 245,   310 ],
                    '2017' : [ 257,   323, 672.78, 1644.76, 1255.72, 2027.79, 2315.59 ],
                    '2018' : [ 267,   332, 1601.95, 1207.43, 1220.84, 1184.09 ],
                    }
        self.triggerTable[ 'AK8PFJet200' ] = {
                    '2016_preVFP' : [ 308,   373 ],
                    '2016' : [ 310,   375 ],
                    '2017' : [ 323,   389, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 332,   398, 368.33, 281.73, 284.86, 276.30 ],
                    }
        self.triggerTable[ 'AK8PFJet260' ] = {
                    '2016_preVFP' : [ 373,   439 ],
                    '2016' : [ 375,   440 ],
                    '2017' : [ 389,   455, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 398,   464, 132.00, 128.00, 122.91, 127.42 ],
                    }
        self.triggerTable[ 'AK8PFJet320' ] = {
                    '2016_preVFP' : [ 439,   526 ],
                    '2016' : [ 440,   525 ],
                    '2017' : [ 455,   543, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 464,   551, 49.29, 48.00, 47.28, 47.92 ],
                    }
        self.triggerTable[ 'AK8PFJet400' ] = {
                    '2016_preVFP' : [ 526,   580 ],
                    '2016' : [ 525,   580 ],
                    '2017' : [ 543,   598, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 551,   606, 16.39, 16.00, 15.92, 15.99 ],
                    }
        self.triggerTable[ 'AK8PFJet450' ] = {
                    '2016_preVFP' : [ 580,   635 ],
                    '2016' : [ 580,   635 ],
                    '2017' : [ 598,   653, 54.41, 382.72, 263.79, 380.97, 410.70  ],
                    '2018' : [ 606,   661, 8.43, 8.00, 7.98, 8.00 ],
                    }
        self.triggerTable[ 'AK8PFJet500' ] = {
                    '2016_preVFP' : [ 635,   10000000. ],
                    '2016' : [ 635,   1000000. ],
                    '2017' : [ 653,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                    '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
                    }
        #self.triggerTable[ 'AK8PFJet550' ] = {
        #            '2016_preVFP' : [ 635,   10000000. ],
        #            '2016' : [ 635,   1000000. ],
        #            '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
        #            '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
        #            }

        
        ### Defining nsubjetiness basis
        self.maxTau = 5
        #list of numbers in third column not used anymore, remove later
        self.nSub_labels = {
                        "_tau_0p5_1": [ 0., 1., 1100  ],
                        "_tau_0p5_2": [ 0., 0.9, 900  ],
                        "_tau_0p5_3": [ 0., 0.8, 800  ],
                        "_tau_0p5_4": [ 0., 0.7, 700  ],
                        "_tau_0p5_5": [ 0., 0.7, 700  ],
                        "_tau_1_1": [ 0., 1., 1100  ],
                        "_tau_1_2": [ 0., 0.7, 700  ],
                        "_tau_1_3": [ 0., 0.5, 500  ],
                        "_tau_1_4": [ 0., 0.5, 500  ],
                        "_tau_1_5": [ 0., 0.5, 500  ],
                        "_tau_2_1": [ 0., 0.5, 500  ],
                        "_tau_2_2": [ 0., 0.3, 600  ],
                        "_tau_2_3": [ 0., 0.2, 400  ],
                        "_tau_2_4": [ 0., 0.2, 400  ],
                        "_tau_2_5": [ 0., 0.2, 400  ]
                }
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
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

        ### Helpers
        #self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
        self.nJet = [ 'Jet'] #1', 'Jet2' ]
        if self.runSDVariables: self.nJet = [ 'sdJet']#, 'Jet2', 'sdJet1', 'sdJet2' ]

        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
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
        if not self.onlyUnc:
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   12, 0, 12) ) 
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            #self.addObject( ROOT.TH1F('nEvtBadPFCands',   ';Events with bad PFCands',   20, 0, 2) )
            if not self.isMC:
                for isel in [ '_only'+x+'_dijetSel' for x in self.triggerTable  ]:
                    self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                    self.addP4Hists( 'LeadPtJetAK8', isel )
                    self.addP4Hists( 'SubleadPtJetAK8', isel )
                    self.addP4Hists( 'CentralJetAK8', isel )
                    self.addP4Hists( 'ForwardJetAK8', isel )

            #### general selection
            for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
                self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
                self.addP4Hists( 'muons', isel )
                self.addP4Hists( 'eles', isel )
                self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
                self.addP4Hists( 'AK8jets', isel )
                self.addP4Hists( 'LeadingPtAK8jet', isel )
                self.addP4Hists( 'SubleadingPtAK8jet', isel )
                self.addP4Hists( 'CentralAK8jet', isel )
                self.addP4Hists( 'ForwardAK8jet', isel )
                self.addObject( ROOT.TH1F('recoPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('recoDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )
                self.addObject( ROOT.TH1F('HT'+isel,   ';HT (GeV)',   200, 0, 2000) )

                self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

            if self.isMC:
                for isel in [ '_noSel' ] + selList :
                    self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   20, 0, 20) )
                    self.addP4Hists( 'genmuons', isel )
                    self.addP4Hists( 'geneles', isel )
                    self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   20, 0, 20) )
                    self.addP4Hists( 'AK8genjets', isel )
                    self.addP4Hists( 'LeadingPtAK8genjet', isel )
                    self.addP4Hists( 'SubleadingPtAK8genjet', isel )
                    self.addP4Hists( 'CentralAK8genjet', isel )
                    self.addP4Hists( 'ForwardAK8genjet', isel )
                    self.addObject( ROOT.TH1F('genPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                    self.addObject( ROOT.TH1F('genDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )
                    self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )
                    self.addObject(ROOT.TH1F('genJetAK8_partonFlavour'+isel,   ';genpartonFlav (MCpns)',   50, -25., 25.) )

        #for isel in selList:
        #    for iJ in self.nJet:
        #        #self.addP4Hists( 'uforeco'+iJ, isel )
        #        #self.addP4Hists( 'ufogen'+iJ, isel )
        #        self.addP4Hists( 'reco'+iJ, isel ) #+'_sortedPt'
        #        if self.isMC: self.addP4Hists( 'gen'+iJ, isel )#+'_sortedPt'

    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  ';p_{T} (GeV)',   600, 0, 3600) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, ';#eta', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_y'+t, ';y', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, ';#phi', 100, -3.5, 3.5) )
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

        #tmp = inputFile.Get('Runs')
        #for i in tmp.LHEPdfSumw: print(i)

        #if self.onlyTrees:
        self.out = wrappedOutputTree
        tmplist=[]
        if not self.isMC:

            self.out.branch( 'eventCategory_nom',  "I")
            self.out.branch( 'totalRecoWeight_nom', "F" )
            self.out.branch( 'passRecoSel_nom', "I" )
            self.out.branch( 'recoSelectedEventNumber_nom', "L" )
            self.out.branch( 'good_nPVs_nom', "F" )
            #self.out.branch( 'eventWithBadPFCands', "L" )
            tmplist.append( 'selRecoJets_nom')
            tmplist.append( 'selRecoJetsF_nom')
            for x in self.triggerTable.keys():
                self.out.branch( 'passHLT_'+x, "I" )
            if not('2016' in self.year): self.out.branch( 'passHLT_AK8PFJet550',"I" )# 1 if getattr(event, 'HLT_AK8PFJet550')==1 and passRecoSel[sys] else 0)  

        if self.isMC:
            for sys in self.sysSource: 
                self.out.branch( 'eventCategory'+sys,  "I")

                #if sys.startswith('_nom'): self.out.branch( 'eventWithBadPFCands', "L" )

                self.out.branch( 'totalRecoWeight'+sys, "F" )
                self.out.branch( 'passRecoSel'+sys, "I" )
                self.out.branch( 'puWeightNom'+sys, "F" )
                self.out.branch( 'l1prefiringWeightNom'+sys, "F")

                self.out.branch( 'good_nPVs'+sys, "F" )

                if sys.startswith('_nom'): 
                    self.out.branch( 'evtGenWeight'+sys, "F" )
                    self.out.branch( 'passGenSel'+sys, "I" )

                self.out.branch( 'recoSelectedEventNumber'+sys, "L" )
                if sys.startswith('_nom'): self.out.branch( 'genSelectedEventNumber'+sys, "L" )
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
            tmplist_gen = ['selGenJets'+sys for sys in self.sysSource if sys.startswith('_nom')] #check if you really need this for the jes or otherwise just the nom branch=>checked, not needed
            tmplist_accepgen = ['accepGenJets'+sys for sys in self.sysSource]
            for x in tmplist_reco+tmplist_truereco+tmplist_gen+tmplist_accepgen: tmplist.append(x)

            tmplist_reco = [ 'selRecoJetsF'+sys for sys in self.sysSource ]
            tmplist_truereco = ['trueRecoJetsF'+sys for sys in self.sysSource]
            tmplist_gen = ['selGenJetsF'+sys for sys in self.sysSource if sys.startswith('_nom')] #check if you really need this for the jes or otherwise just the nom branch=>checked, not needed
            tmplist_accepgen = ['accepGenJetsF'+sys for sys in self.sysSource]
            for x in self.triggerTable.keys():
                self.out.branch( 'passHLT_'+x, "I" )
            if not('2016' in self.year): self.out.branch( 'passHLT_AK8PFJet550',"I" )

            for x in tmplist_reco+tmplist_truereco+tmplist_gen+tmplist_accepgen: tmplist.append(x)

        print ("Stored jet branches:", tmplist)
        for iJ in tmplist:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_y',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_msoftdrop',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_WTA',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21_exkT',  'F')#, lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32_exkT',  'F')#, lenVar='n'+iJ)
            
            for x in self.nSub_labels:
                self.out.branch(iJ+x, 'F')#, lenVar='n'+iJ )
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        if not self.onlyTrees:
            if self.isMC and not self.onlyUnc:
                self.genLevel = self.response+self.miss

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
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''
        self.isMC = event.run == 1
        #print(event.event)
        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenJets = self.genSelection( event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoJets = self.recoSelection( event )
        #print(event.event)


        if not( self.isMC) and not(passRecoSel['_nom']): return False



        if (self.isMC) and not(passGenSel) and not(passRecoSel['_nom']): 
            #print (passRecoSel[sys])
            return False
        
        #self.FlagBadPFCands=False

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
            
            genJet = OrderedDict()
            recoJet = OrderedDict()
            genJetF = OrderedDict()
            recoJetF = OrderedDict()
            tmpRecoJets = OrderedDict()
            tmpGenJets = OrderedDict()

            if self.isMC and (not passRecoSel[sys]):# and not self.onlyUnc.startswith(self.sysWeightList):

                if passGenSel:# and not self.onlyTrees:

                    #### Misses
                    self.miss = self.miss+1
                    self.eventCategory = 3
                    tmpGenJets={}
                    tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                    tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )

                    #if self.FlagBadPFCands==True:
                    #    self.FlagBadPFCandsCounter+=1
                    #    #print ("Found an event with bad PF Cands in gen will not retain this event",event.event,self.FlagBadPFCandsCounter,self.evtCounter,self.FlagBadPFCandsCounter/self.evtCounter*100.)
                    #    self.out.fillBranch('eventWithBadPFCands', event.event)
                    #    return False
                    #else: self.out.fillBranch('eventWithBadPFCands',-3)                         

                    if abs(tmpGenJets[0]['jet'].rapidity) > abs(tmpGenJets[1]['jet'].rapidity):
                        genJet['Jet'] = tmpGenJets[1]#self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                        genJetF['Jet'] = tmpGenJets[0]
                    else: 
                        genJet['Jet'] = tmpGenJets[0]#self.createNsubBasis( selGenJets[0], event, 'GenCands', True  )
                        genJetF['Jet'] = tmpGenJets[1]


                    if sys.startswith('_nom'): 
                        self.fillBranches( event, 'selGenJets'+sys, genJet, False, sys ) 
                        self.fillBranches( event, 'selGenJetsF'+sys, genJetF, False, sys ) 
                    self.fillBranches( event, 'accepGenJets'+sys, genJet, True, sys )
                    self.fillBranches( event, 'accepGenJetsF'+sys, genJetF, True, sys )
                else: self.evtGenWeight=0.
                self.totalRecoWeight=0.

            if passRecoSel[sys]:  #### Detector level dist.
                
                self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel
                self.eventCategory = 1

                tmpRecoJets[sys]={}
                tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                tmpRecoJets[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )
                #if self.FlagBadPFCands==True:
                    #self.FlagBadPFCandsCounter+=1
                    #print ("Found an event with bad PF Cands in reco, will not retain this event",event.event,self.FlagBadPFCandsCounter,self.evtCounter,self.FlagBadPFCandsCounter/self.evtCounter*100.)
                    #self.out.fillBranch('eventWithBadPFCands', event.event)
                    #return False

                


                #######################Testing block for pad pfnanos will remove#####################
                # prevent events with negative mass pfcands in the branches from 
                # being considered in the jets passing our selection to see if the differences are massive

                if passGenSel:
                    tmpGenJets={}
                    tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                    tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                    #if self.FlagBadPFCands==True:
                        #self.FlagBadPFCandsCounter+=1
                        #print ("Found an event with bad PF Cands in gen will not retain this event",event.event,self.FlagBadPFCandsCounter,self.evtCounter,self.FlagBadPFCandsCounter/self.evtCounter*100.)
                        #self.out.fillBranch('eventWithBadPFCands', event.event)                        
                        #return False
                    #else: self.out.fillBranch('eventWithBadPFCands',-3)                         
                #######################Testing block for pad pfnanos will remove#####################





                if abs(tmpRecoJets[sys][0]['jet'].rapidity) > abs(tmpRecoJets[sys][1]['jet'].rapidity):
                    recoJet['Jet'] = tmpRecoJets[sys][1]
                    recoJetF['Jet'] = tmpRecoJets[sys][0]#self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )
                else: 
                    recoJet['Jet'] = tmpRecoJets[sys][0]#self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )                    
                    recoJetF['Jet'] = tmpRecoJets[sys][1]

                self.fillBranches( event, 'selRecoJets'+sys, recoJet, False, sys ) #not a dummy fill so dummy=False
                self.fillBranches( event, 'selRecoJetsF'+sys, recoJetF, False, sys ) #not a dummy fill so dummy=False

            
                if self.isMC:
                    deltaRmatch=False
                    if not passGenSel:  ##### fakes
                        self.fakes = self.fakes+1
                        self.eventCategory = 2
                        self.fillBranches( event, 'trueRecoJets'+sys, recoJet, True, sys ) #True=dummyFill
                        self.fillBranches( event, 'trueRecoJetsF'+sys, recoJetF, True, sys ) #True=dummyFill
                        self.evtGenWeight=0.

                    else:        ##### go to matrix
                        
                        tmpGenJets={}
                        tmpGenJets[0] = self.createNsubBasis( selGenJets[0], event, 'GenCands', True )
                        tmpGenJets[1] = self.createNsubBasis( selGenJets[1], event, 'GenCands', True )

                        if abs(tmpGenJets[0]['jet'].rapidity) > abs(tmpGenJets[1]['jet'].rapidity):
                            genJet['Jet'] = tmpGenJets[1]#self.createNsubBasis( selGenJets[1], event, 'GenCands', True )
                            genJetF['Jet'] = tmpGenJets[0]
                        else: 
                            genJet['Jet'] = tmpGenJets[0]#self.createNsubBasis( selGenJets[0], event, 'GenCands', True  )
                            genJetF['Jet'] = tmpGenJets[1]

                        if sys.startswith('_nom'): 
                            self.fillBranches( event, 'selGenJets'+sys, genJet, False, sys ) 
                            self.fillBranches( event, 'selGenJetsF'+sys, genJetF, False, sys ) 

                        # small redundancy here since reco/genJet['Jet'] is already one of the two tmpReco/GenJets:
                        # basically, checking that the measurement jets as stored in the recoJet and genJet objects are in whack with one another, 
                        # while also checking that both of the pT ordered tmp jets in gen/reco are deltaR matched, since if not it's not safe to continue
                        # this block is to ensure that the measurement jet is pristinely selected in an event where the entire dijet system is consistent b/w gen- and reco-level 
                        # doing all this so that it doesn't need to be rechecked later via storing the subleading jet too (too much computational overhead for QCD samples)
                        if ( self.DrRapPhi( recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() ) < 0.4 ) and ( self.DrRapPhi( recoJetF['Jet']['jet'].p4(), genJetF['Jet']['jet'].p4()) < 0.4 ):#self.DrRapPhi( tmpRecoJets[sys][0]['jet'].p4(), tmpGenJets[0]['jet'].p4() ) < 0.4 and self.DrRapPhi( tmpRecoJets[sys][1]['jet'].p4(), tmpGenJets[1]['jet'].p4() ) < 0.4):
                            deltaRmatch = True
                            self.response= self.response+1
                            self.eventCategory = 4
                        
                        if deltaRmatch:
                            self.fillBranches( event, 'accepGenJets'+sys, genJet, False, sys )    
                            self.fillBranches( event, 'accepGenJetsF'+sys, genJetF, False, sys )    
                            self.fillBranches( event, 'trueRecoJets'+sys, recoJet, False, sys )
                            self.fillBranches( event, 'trueRecoJetsF'+sys, recoJetF, False, sys )
                        
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
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') if passRecoSel[sys] else 0.)
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
                self.out.fillBranch( 'good_nPVs'+sys, getattr( event, 'PV_npvsGood') if passRecoSel[sys] else 0.)
                if sys.startswith('_nom'): 
                    self.out.fillBranch( 'evtGenWeight'+sys, self.evtGenWeight if passGenSel or passRecoSel[sys] else 0.) 
                    self.out.fillBranch( 'passGenSel'+sys, 1 if passGenSel else 0) 
                    
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

        return True


    #############################################################################
    def recoSelection( self, event, sysUnc=[] ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))

        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and abs(x.p4().Eta()) < self.maxElectronEta and x.convVeto] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minLooseMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.mediumId and x.miniPFRelIso_all<0.2 ] #and x.pfIsoId>=2 and abs(x.dxy)<0.2 and abs(x.dz)<0.5

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons) 

        ##################################################
        if not len(AK8jets)==0:
            for ijets in AK8jets: 
                ijets.rapidity = ijets.p4().Rapidity()
        recoAK8jets = {}
        passSel = {}
        ptAsym = {}
        deltaPhi = {}
        iSel = {}
        for sys in self.sysSource:
            #if sys.startswith(self.sysWeightList): sys = '_nom'
            #### Basic AK8 jet selection
            #JetID (https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD): For UL16-18 samples jetId==2 means: pass tight ID, fail tightLepVeto, jetId==6 means: pass tight and tightLepVeto ID, so we use >2 (ie, 6) 
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minLeadAK8JetPtDijet and abs( x.rapidity ) < self.maxLeadAK8JetRap and (x.jetId > 2)] #self.maxJetAK8Eta
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ),reverse=True)
            ##################################################

            ##### Applying selection
            #print ("in reco sel", event.event, recoMuons, recoElectrons, recoAK8jets[sys], sys)
            #print (self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets[sys], sys ))
            passSel[sys], ptAsym[sys], deltaPhi[sys] = self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets[sys], sys )
            iSel[sys] = '_dijetSel' if passSel[sys] else None
        
        ##### Trigger weights
        if not self.isMC:
            if passSel['_nom']:
               
                for itrigger, itrValues in self.triggerTable.items():
                    if ( getattr(event, 'HLT_'+itrigger)==1 ):

                        getattr( self, 'nPVs_only'+itrigger+'_dijetSel' ).Fill( getattr( event, 'PV_npvsGood') )
                        if len(recoAK8jets['_nom'])>0:#and not self.onlyTrees:
                            getattr( self, 'LeadPtJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].pt_nom )   ### pt_nom here to ensure data process
                            getattr( self, 'LeadPtJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].eta )
                            getattr( self, 'LeadPtJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].rapidity )
                            getattr( self, 'LeadPtJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].phi )
                            getattr( self, 'LeadPtJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                            getattr( self, 'LeadPtJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )
                            
                            getattr( self, 'SubleadPtJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].pt_nom )   ### pt_nom here to ensure data process
                            getattr( self, 'SubleadPtJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].eta )
                            getattr( self, 'SubleadPtJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].rapidity )
                            getattr( self, 'SubleadPtJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].phi )
                            getattr( self, 'SubleadPtJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                            getattr( self, 'SubleadPtJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )
                            
                            
                            if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                                getattr( self, 'CentralJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].pt_nom )   ### pt_nom here to ensure data process
                                getattr( self, 'CentralJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].eta )
                                getattr( self, 'CentralJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].rapidity )
                                getattr( self, 'CentralJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].phi )
                                getattr( self, 'CentralJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                                getattr( self, 'CentralJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )

                                getattr( self, 'ForwardJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].pt_nom )   ### pt_nom here to ensure data process
                                getattr( self, 'ForwardJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].eta )
                                getattr( self, 'ForwardJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].rapidity )
                                getattr( self, 'ForwardJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].phi )
                                getattr( self, 'ForwardJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                                getattr( self, 'ForwardJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )
                            else:
                                getattr( self, 'CentralJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].pt_nom )   ### pt_nom here to ensure data process
                                getattr( self, 'CentralJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].eta )
                                getattr( self, 'CentralJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].rapidity )
                                getattr( self, 'CentralJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].phi )
                                getattr( self, 'CentralJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                                getattr( self, 'CentralJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )

                                getattr( self, 'ForwardJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].pt_nom )   ### pt_nom here to ensure data process
                                getattr( self, 'ForwardJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].eta )
                                getattr( self, 'ForwardJetAK8_y_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].rapidity )
                                getattr( self, 'ForwardJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].phi )
                                getattr( self, 'ForwardJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                                getattr( self, 'ForwardJetAK8_msoftdrop_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )

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
                getattr( self, 'AK8jets_msoftdrop_noSelnoWeight' ).Fill( ijet.msoftdrop_nom )
            if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                getattr( self, 'LeadingPtAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                getattr( self, 'LeadingPtAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                getattr( self, 'LeadingPtAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                getattr( self, 'LeadingPtAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                getattr( self, 'LeadingPtAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )
                if  len(recoAK8jets['_nom'])>1:
                    getattr( self, 'SubleadingPtAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                    getattr( self, 'SubleadingPtAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                    getattr( self, 'SubleadingPtAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                    getattr( self, 'SubleadingPtAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                    getattr( self, 'SubleadingPtAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                    getattr( self, 'SubleadingPtAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )

                    if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                        getattr( self, 'CentralAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                        getattr( self, 'CentralAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                        getattr( self, 'CentralAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                        getattr( self, 'CentralAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                        getattr( self, 'CentralAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )

                        getattr( self, 'ForwardAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                        getattr( self, 'ForwardAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                        getattr( self, 'ForwardAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                        getattr( self, 'ForwardAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                        getattr( self, 'ForwardAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )
                    else: 
                        getattr( self, 'CentralAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].pt_nom )
                        getattr( self, 'CentralAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].eta )
                        getattr( self, 'CentralAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].rapidity )
                        getattr( self, 'CentralAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].phi )
                        getattr( self, 'CentralAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].mass_nom )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom )

                        getattr( self, 'ForwardAK8jet_pt_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].pt_nom )
                        getattr( self, 'ForwardAK8jet_eta_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].eta )
                        getattr( self, 'ForwardAK8jet_y_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].rapidity )
                        getattr( self, 'ForwardAK8jet_phi_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].phi )
                        getattr( self, 'ForwardAK8jet_mass_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].mass_nom )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSelnoWeight' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )

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
                getattr( self, 'AK8jets_msoftdrop_noSel' ).Fill( ijet.msoftdrop_nom, weight )

            if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                getattr( self, 'LeadingPtAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                getattr( self, 'LeadingPtAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                getattr( self, 'LeadingPtAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                getattr( self, 'LeadingPtAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                getattr( self, 'LeadingPtAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )

                if  len(recoAK8jets['_nom'])>1:
                    getattr( self, 'SubleadingPtAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'SubleadingPtAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'SubleadingPtAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity , weight )
                    getattr( self, 'SubleadingPtAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi , weight )
                    getattr( self, 'SubleadingPtAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'SubleadingPtAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )

                    if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                        getattr( self, 'CentralAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                        getattr( self, 'CentralAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                        getattr( self, 'CentralAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                        getattr( self, 'CentralAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                        getattr( self, 'CentralAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )

                        getattr( self, 'ForwardAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                        getattr( self, 'ForwardAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                        getattr( self, 'ForwardAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                        getattr( self, 'ForwardAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi, weight )
                        getattr( self, 'ForwardAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )
                    else: 
                        getattr( self, 'CentralAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                        getattr( self, 'CentralAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][1].eta, weight )
                        getattr( self, 'CentralAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                        getattr( self, 'CentralAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][1].phi, weight )
                        getattr( self, 'CentralAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                        getattr( self, 'CentralAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )

                        getattr( self, 'ForwardAK8jet_pt_noSel' ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                        getattr( self, 'ForwardAK8jet_eta_noSel' ).Fill( recoAK8jets['_nom'][0].eta, weight )
                        getattr( self, 'ForwardAK8jet_y_noSel' ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                        getattr( self, 'ForwardAK8jet_phi_noSel' ).Fill( recoAK8jets['_nom'][0].phi, weight )
                        getattr( self, 'ForwardAK8jet_mass_noSel' ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                        getattr( self, 'ForwardAK8jet_msoftdrop_noSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )
                        
            getattr( self, 'recoPtAsym_noSel' ).Fill( ptAsym['_nom'], weight )
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
                    getattr( self, 'AK8jets_msoftdrop'+iSel['_nom'] ).Fill( ijet.msoftdrop_nom, weight )
                #if not len(recoAK8jets['_nom'])==0:
                getattr( self, 'LeadingPtAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                getattr( self, 'LeadingPtAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                getattr( self, 'LeadingPtAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                getattr( self, 'LeadingPtAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                getattr( self, 'LeadingPtAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                getattr( self, 'LeadingPtAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )

                getattr( self, 'SubleadingPtAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                getattr( self, 'SubleadingPtAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                getattr( self, 'SubleadingPtAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity , weight )
                getattr( self, 'SubleadingPtAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi , weight )
                getattr( self, 'SubleadingPtAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                getattr( self, 'SubleadingPtAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )

                if abs(recoAK8jets['_nom'][0].rapidity)<abs(recoAK8jets['_nom'][1].rapidity):
                    getattr( self, 'CentralAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                    getattr( self, 'CentralAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                    getattr( self, 'CentralAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                    getattr( self, 'CentralAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                    getattr( self, 'CentralAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                    getattr( self, 'CentralAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )

                    getattr( self, 'ForwardAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'ForwardAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'ForwardAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                    getattr( self, 'ForwardAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi, weight )
                    getattr( self, 'ForwardAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'ForwardAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )
                else: 
                    getattr( self, 'CentralAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].pt_nom, weight )
                    getattr( self, 'CentralAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].eta, weight )
                    getattr( self, 'CentralAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].rapidity, weight )
                    getattr( self, 'CentralAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].phi, weight )
                    getattr( self, 'CentralAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].mass_nom, weight )
                    getattr( self, 'CentralAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][1].msoftdrop_nom, weight )

                    getattr( self, 'ForwardAK8jet_pt'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].pt_nom, weight )
                    getattr( self, 'ForwardAK8jet_eta'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].eta, weight )
                    getattr( self, 'ForwardAK8jet_y'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].rapidity, weight )
                    getattr( self, 'ForwardAK8jet_phi'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].phi, weight )
                    getattr( self, 'ForwardAK8jet_mass'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].mass_nom, weight )
                    getattr( self, 'ForwardAK8jet_msoftdrop'+iSel['_nom'] ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom, weight )

                getattr( self, 'recoPtAsym'+iSel['_nom'] ).Fill( ptAsym['_nom'], weight )
                getattr( self, 'recoDeltaPhi'+iSel['_nom'] ).Fill( deltaPhi['_nom'], weight )
        return passSel, iSel, recoMuons, recoElectrons, recoAK8jets

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing reco information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genParticles = Collection(event, 'GenPart')

        ### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt > self.minLooseElectronPt and abs(x.eta)<self.maxElectronEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minLooseMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        #print ("in gensel", event.event, genMuons, genElectrons, genJetsAK8, '' )

        ##################################################

        #### Basic AK8 jet selection
        if not len(genJetsAK8)==0:
            for ijets in genJetsAK8: 
                ijets.rapidity = ijets.p4().Rapidity()
                ijets.msoftdrop = self.getGenJetAK8softdropmass( AK8jet=ijets, event=event, PFCollection='GenCands', isGen=True)
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minLeadAK8JetPtDijet and abs(x.rapidity) < self.maxLeadAK8JetRap]#maxJetAK8Eta 
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        #print ("in gensel", event.event, genMuons, genElectrons, genAK8jets, '' )
        passSel, ptAsym, deltaPhi = self.dijetSelection( event, genMuons, genElectrons, genAK8jets, '' )
        iSel = '_dijetSel' if passSel else None
        #### Weight
        weight = event.genWeight

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

            getattr( self, 'genPtAsym_noSel' ).Fill( ptAsym, weight )
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

                getattr( self, 'genPtAsym'+iSel ).Fill( ptAsym, weight )
                getattr( self, 'genDeltaPhi'+iSel ).Fill( deltaPhi, weight )
                getattr(self, 'genJetAK8_partonFlavour'+iSel).Fill(genAK8jets[0].partonFlavour, weight)

        return passSel, iSel, genMuons, genElectrons, genAK8jets

    #############################################################################
    def dijetSelection( self, event, muons, electrons, AK8jets, ptLabel ):

        #print ("in main dijetSel", event.event, muons, electrons, AK8jets, ptLabel)
        if (len(muons)+len(electrons))==0 and len(AK8jets)>1:# and abs(AK8jets[0].eta)<self.maxLeadAK8JetEta):#
            jet1Pt = getattr( AK8jets[0], 'pt'+ptLabel )
            jet2Pt = getattr( AK8jets[1], 'pt'+ptLabel )
            #print ("in main dijetSel", event.event, muons, electrons, AK8jets, ptLabel, jet1Pt, jet2Pt)

            tmpJet1 = ROOT.TLorentzVector()
            tmpJet1.SetPtEtaPhiM( getattr( AK8jets[0], 'pt'+ptLabel ), AK8jets[0].eta, AK8jets[0].phi, getattr( AK8jets[0], 'mass'+ptLabel )  )
            tmpJet2 = ROOT.TLorentzVector()
            tmpJet2.SetPtEtaPhiM( getattr( AK8jets[1], 'pt'+ptLabel ), AK8jets[1].eta, AK8jets[1].phi, getattr( AK8jets[1], 'mass'+ptLabel )  )
            
            deltaR_flag = True

            #print("DijetSel Loop flag",len(AK8jets))
            for i in range(1,len(AK8jets)):
                tmpJet = ROOT.TLorentzVector()
                tmpJet.SetPtEtaPhiM( getattr( AK8jets[i], 'pt'+ptLabel ), AK8jets[i].eta, AK8jets[i].phi, getattr( AK8jets[i], 'mass'+ptLabel )  )
                if self.DrRapPhi(tmpJet1,tmpJet)<1.6:# or abs(tmpJet1.DeltaPhi(tmpJet)<2.): #2*jet radius separation + deltaphi sep. between leading AK8 jet and any other fatjets in the event
                    deltaR_flag=False
                    break
            #print("Post-dijetSel Loop flag",len(AK8jets))

            if not deltaR_flag: return False, -999, -999

            else:
                #do calcs only when necessary
                ptAsym = ( jet1Pt - jet2Pt ) / (jet1Pt + jet2Pt)
                deltaPhi = tmpJet1.DeltaPhi( tmpJet2 )

                if (jet1Pt>self.minLeadAK8JetPtDijet) and (abs(deltaPhi)>2.):# and (ptAsym<0.3)  and
                    return True, ptAsym, deltaPhi
                else: return False, -999, -999
        
        else: return False, -999, -999

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
                ak8jet['msoftdrop'] = sd_AK8jets[0].m()
            else: ak8jet['msoftdrop'] = -1. 
        return ak8jet['msoftdrop']


    #############################################################################
    def fillBranches( self, event, jetLabel, jetInfo, dummy=False, sys='_nom' ): 
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
                    self.out.fillBranch(jetLabel+"_tau21_exkT", iJ['tau21_exkT']  )
                    self.out.fillBranch(jetLabel+"_tau32_exkT", iJ['tau32_exkT']  )
            
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
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)

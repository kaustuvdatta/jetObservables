import ROOT
import math, os, sys
import numpy as np
import pandas as pd
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nSubProd(Module):

    def __init__(self, sysSource=[], leptonSF={}, year='2017', triggerFile='../data/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017B.pkl', isMC=True):
        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.year = year
        self.isMC = isMC

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
        self.minAK8JetPt = 170  ### this is the basic minimum, not the final
        self.maxJetAK8Eta = 2.4

        ### Kinematics Cuts Jets ###
        self.minJetPt = 25
        self.maxJetEta = 2.4
        self.minBDisc = 0.3040  ### L: 0.0532, M: 0.3040, T: 0.7476, for DeepJet (ie, DeepFlavB)

        ### Kinenatic Cuts Muons ###
        self.minLooseMuonPt = 40.
        self.minTightMuonPt = 55.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 40.
        self.minTightElectronPt = 120.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.4]

        self.totalWeight = 1
        self.triggerWeight = 1

        #self.triggerInfo = pd.read_pickle( triggerFile )
        self.runTables = {
                '2017' : {
                    'low' : [ 297046, 299368, 302030, 303824, 305040  ],
                    'high' : [ 299329, 302029, 303434, 304797, 306462  ]
                    },
                '2018' : {
                    'low' : [ 315252, 317080, 319337, 320673  ],
                    'high' : [ 316995, 319310, 320065, 325175 ]
                    }
                }
        self.triggerTable = OrderedDict()
        self.triggerTable[ 'AK8PFJet80' ] = {
                    '2017' : [ 200,   238, 50510.97, 6727.43, 6854.43, 23381.99, 38371.17 ],
                    '2018' : [ 200,   246, 17046.53, 34000.07, 33988.86, 33998.78 ],
                    }
        self.triggerTable[ 'AK8PFJet140' ] = {
                    '2017' : [ 238,   304, 672.78, 1644.76, 1255.72, 2027.79, 2315.59 ],
                    '2018' : [ 246,   311, 1601.95, 1207.43, 1220.84, 1184.09 ],
                    }
        self.triggerTable[ 'AK8PFJet200' ] = {
                    '2017' : [ 304,   370, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 311,   376, 368.33, 281.73, 284.86, 276.30 ],
                    }
        self.triggerTable[ 'AK8PFJet260' ] = {
                    '2017' : [ 370,   436, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 376,   441, 132.00, 128.00, 122.91, 127.42 ],
                    }
        self.triggerTable[ 'AK8PFJet320' ] = {
                    '2017' : [ 436,   524, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 441,   528, 49.29, 48.00, 47.28, 47.92 ],
                    }
        self.triggerTable[ 'AK8PFJet400' ] = {
                    '2017' : [ 524,   579, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 528,   583, 16.39, 16.00, 15.92, 15.99 ],
                    }
        self.triggerTable[ 'AK8PFJet450' ] = {
                    '2017' : [ 579,   634, 54.41, 382.72, 263.79, 380.97, 410.70  ],
                    '2018' : [ 583,   637, 8.43, 8.00, 7.98, 8.00 ],
                    }
        self.triggerTable[ 'AK8PFJet500' ] = {
                    '2017' : [ 634.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                    '2018' : [ 637,   1000000., 1, 1, 1, 1 ],
                    }

        ### Defining nsubjetiness basis
        self.maxTau = 5
        self.nSub_labels = {
                        "_tau_0p5_1": [ 1, 200  ],
                        "_tau_0p5_2": [ 0.8, 800  ],
                        "_tau_0p5_3": [ 0.6, 600  ],
                        "_tau_0p5_4": [ 0.6, 600  ],
                        "_tau_0p5_5": [ 0.6, 600  ],
                        "_tau_1_1": [ 1, 200  ],
                        "_tau_1_2": [ 0.6, 600  ],
                        "_tau_1_3": [ 0.4, 800  ],
                        "_tau_1_4": [ 0.4, 800  ],
                        "_tau_1_5": [ 0.4, 800  ],
                        "_tau_2_1": [ 1, 200  ],
                        "_tau_2_2": [ 0.4, 800  ],
                        "_tau_2_3": [ 0.4, 800  ],
                        "_tau_2_4": [ 0.4, 800  ],
                        "_tau_2_5": [ 0.4, 800  ]
                }
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )
        self.nSub1_OP_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 ) ##### needed for genjet tau21 or tau32

        ### Softdrop quantities
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, self.minAK8JetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")

        ### Helpers
        self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
        self.nJet = [ 'Jet1', 'Jet2' ]

        ### Uncerstinties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        ## JES from https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
        self.JESLabels = [ "Total" ]
        #self.JESLabels = [ "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpMuZero", "PileUpEnvelope", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "Total", "TotalNoFlavor", "TotalNoTime", "TotalNoFlavorNoTime" ]

    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        tauBins = 100

        ### Booking histograms
        self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   12, 0, 12) )
        self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
        self.addObject( ROOT.TH1F('Btagweight',   ';BtagWeight',   25, 0, 2) )
        #### general selection
        selList = ([ '_'+x+'_dijetSel' for x in self.triggerTable  ] + [ '_'+x+'_weight_dijetSel' for x in self.triggerTable  ]) if not self.isMC else [ '_dijetSel' ]
        for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
            self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
            self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
            self.addP4Hists( 'muons', isel )
            self.addP4Hists( 'eles', isel )
            self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
            self.addP4Hists( 'AK8jets', isel )
            self.addObject( ROOT.TH1F('recoPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
            self.addObject( ROOT.TH1F('recoDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )

            self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

        if self.isMC:
            for isel in [ '_noSel' ] + selList :
                self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   20, 0, 20) )
                self.addP4Hists( 'genmuons', isel )
                self.addP4Hists( 'geneles', isel )
                self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   20, 0, 20) )
                self.addP4Hists( 'AK8genjets', isel )
                self.addObject( ROOT.TH1F('genPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                self.addObject( ROOT.TH1F('genDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )


        for isel in selList:
            for iJ in self.nJet:
                #self.addP4Hists( 'uforeco'+iJ, isel )
                #self.addP4Hists( 'ufogen'+iJ, isel )
                self.addP4Hists( 'reco'+iJ+'_sortedPt', isel )
                if self.isMC: self.addP4Hists( 'gen'+iJ+'_sortedPt', isel )
                for itype in ( [ 'gen', 'missgen', 'accepgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] ):

                    if itype.endswith('reco'):

                        for sysUnc in self.sysSource:
                            self.addP4Hists( itype+iJ, sysUnc+isel )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{21}', tauBins, 0, 1) )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{32}', tauBins, 0, 1) )

                            if itype.startswith('truereco') and self.isMC:
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau21'+sysUnc+isel, ';AK8 accepgen jet #tau_{21};AK8 truereco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau32'+sysUnc+isel, ';AK8 accepgen jet #tau_{32};AK8 truereco jet #tau_{32}', tauBins, 0, 1, tauBins, 0, 1) )
                                if sysUnc.endswith('nom'):
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_pt'+isel, ';AK8 reco/gen jet pt', 100, 0, 2) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_sdmass'+isel, ';AK8 reco/gen jet SD mass', 100, 0, 2) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau21'+isel, ';AK8 reco/gen jet #tau_{21}', 100, 0, 2) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau32'+isel, ';AK8 reco/gen jet #tau_{32}', 100, 0, 2) )

                            for x, y in self.nSub_labels.items():
                                self.addObject( ROOT.TH1F(itype+iJ+x+sysUnc+isel, ';AK8 '+itype+' jet #tau', y[1], 0, y[0] ) )
                                if itype.startswith('reco') and self.isMC:
                                    self.addObject( ROOT.TH2F('resp'+iJ+x+sysUnc+isel, ';AK8 gen jet '+x+';AK8 reco jet '+x, y[1], 0, y[0], y[1], 0, y[0] ) )
                                    if sysUnc.endswith('nom'): self.addObject( ROOT.TH1F('resol'+iJ+x+isel, ';AK8 reco/gen jet '+x, 200, 0, 2) )
                    else:

                        self.addP4Hists( itype+iJ, isel )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+isel, ';AK8 '+itype+' jet #tau_{21}', tauBins, 0, 1) )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+isel, ';AK8 '+itype+' jet #tau_{32}', tauBins, 0, 1) )

                        for x, y in self.nSub_labels.items():
                            self.addObject( ROOT.TH1F(itype+iJ+x+isel, ';AK8 '+itype+' jet '+x, y[1], 0, y[0] ) )


    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -4.0, 4.0 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14159) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 1000) )


    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree
        self.out.branch('triggerWeight',  "F")
        self.out.branch('eventCategory',  "I")

        for iJ in ([ 'Reco', 'Gen' ] if self.isMC else ['Reco'] ):
            self.out.branch('nsel'+iJ+'jets',  'I')  ### dummy for nanoAOD Tools
            self.out.branch('sel'+iJ+'jets_pt',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_eta',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_phi',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_mass',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_Tau21',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_Tau32',  'F', lenVar='nsel'+iJ+'jets')

            for x in self.nSub_labels:
                self.out.branch('sel'+iJ+'jets'+x, 'F', lenVar='nsel'+iJ+'jets' )
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        if self.isMC:
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

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenJets = self.genSelection(event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoJets = self.recoSelection( event )
        #print '3', passRecoSel, iRecoSel
        if self.isMC:
            if (not (iRecoSel or iGenSel)) and (not(passGenSel or passRecoSel)): return False
        if not self.isMC and not passRecoSel: return False

        #print iRecoSel, passRecoSel

        genJet = OrderedDict()
        recoJet = OrderedDict()
        if passRecoSel:  #### Detector level dist.

            tmpRecoJet1 = self.createNsubBasis( selRecoJets[0], event, 'PFCands' )
            tmpRecoJet2 = self.createNsubBasis( selRecoJets[1], event, 'PFCands' )
            getattr( self, 'recoJet1_sortedPt_pt'+iRecoSel ).Fill( getattr(tmpRecoJet1['jet'], 'pt_nom' ), self.totalWeight )
            getattr( self, 'recoJet1_sortedPt_eta'+iRecoSel ).Fill( getattr(tmpRecoJet1['jet'], 'eta' ), self.totalWeight )
            getattr( self, 'recoJet1_sortedPt_phi'+iRecoSel ).Fill( getattr(tmpRecoJet1['jet'], 'phi' ), self.totalWeight )
            getattr( self, 'recoJet1_sortedPt_mass'+iRecoSel ).Fill( getattr(tmpRecoJet1['jet'], 'msoftdrop_nom' ), self.totalWeight )
            getattr( self, 'recoJet2_sortedPt_pt'+iRecoSel ).Fill( getattr(tmpRecoJet2['jet'], 'pt_nom' ), self.totalWeight )
            getattr( self, 'recoJet2_sortedPt_eta'+iRecoSel ).Fill( getattr(tmpRecoJet2['jet'], 'eta' ), self.totalWeight )
            getattr( self, 'recoJet2_sortedPt_phi'+iRecoSel ).Fill( getattr(tmpRecoJet2['jet'], 'phi' ), self.totalWeight )
            getattr( self, 'recoJet2_sortedPt_mass'+iRecoSel ).Fill( getattr(tmpRecoJet2['jet'], 'msoftdrop_nom' ), self.totalWeight )
            if not self.isMC:
                tmpSel = iRecoSel.replace('_dijet', '_weight_dijet')
                getattr( self, 'recoJet1_sortedPt_pt'+tmpSel ).Fill( getattr(tmpRecoJet1['jet'], 'pt_nom' ), self.triggerWeight )
                getattr( self, 'recoJet1_sortedPt_eta'+tmpSel ).Fill( getattr(tmpRecoJet1['jet'], 'eta' ), self.triggerWeight )
                getattr( self, 'recoJet1_sortedPt_phi'+tmpSel ).Fill( getattr(tmpRecoJet1['jet'], 'phi' ), self.triggerWeight )
                getattr( self, 'recoJet1_sortedPt_mass'+tmpSel ).Fill( getattr(tmpRecoJet1['jet'], 'msoftdrop_nom' ), self.triggerWeight )
                getattr( self, 'recoJet2_sortedPt_pt'+tmpSel ).Fill( getattr(tmpRecoJet2['jet'], 'pt_nom' ), self.triggerWeight )
                getattr( self, 'recoJet2_sortedPt_eta'+tmpSel ).Fill( getattr(tmpRecoJet2['jet'], 'eta' ), self.triggerWeight )
                getattr( self, 'recoJet2_sortedPt_phi'+tmpSel ).Fill( getattr(tmpRecoJet2['jet'], 'phi' ), self.triggerWeight )
                getattr( self, 'recoJet2_sortedPt_mass'+tmpSel ).Fill( getattr(tmpRecoJet2['jet'], 'msoftdrop_nom' ), self.triggerWeight )


            if abs(tmpRecoJet1['jet'].eta) > abs(tmpRecoJet2['jet'].eta):
                recoJet['Jet1'] = tmpRecoJet1
                recoJet['Jet2'] = tmpRecoJet2
            else:
                recoJet['Jet1'] = tmpRecoJet2
                recoJet['Jet2'] = tmpRecoJet1

            self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel

            for iRJ,ireco in recoJet.items():

                WEIGHT =  self.totalWeight
                for sysUnc in self.sysSource:
                    if ( getattr(ireco['jet'], 'pt'+(sysUnc if not sysUnc.startswith('_pu') else '_nom')) > self.minLeadAK8JetPtDijet ):

                        if sysUnc.startswith('_pu'):
                            WEIGHT = event.genWeight * getattr( event, 'puWeight'+sysUnc.split('pu')[1] )
                            getattr( self, 'reco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                        else:
                            getattr( self, 'reco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt'+sysUnc ), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sysUnc ), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_eta'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_phi'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau21'+sysUnc+iRecoSel ).Fill( ireco['tau21'], WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau32'+sysUnc+iRecoSel ).Fill( ireco['tau32'], WEIGHT )
                        for tauN in range(1, self.maxTau+1):
                            getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['2'+str(tauN)], WEIGHT )

                        if not self.isMC:
                            tmpSel = iRecoSel.replace('_dijet', '_weight_dijet')
                            getattr( self, 'reco'+iRJ+'_pt'+sysUnc+tmpSel ).Fill( getattr(ireco['jet'], 'pt'+sysUnc ), self.triggerWeight )
                            getattr( self, 'reco'+iRJ+'_mass'+sysUnc+tmpSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sysUnc ), self.triggerWeight )
                            getattr( self, 'reco'+iRJ+'_eta'+sysUnc+tmpSel ).Fill( getattr(ireco['jet'], 'eta'), self.triggerWeight )
                            getattr( self, 'reco'+iRJ+'_phi'+sysUnc+tmpSel ).Fill( getattr(ireco['jet'], 'phi'), self.triggerWeight )
                            getattr( self, 'reco'+iRJ+'_tau21'+sysUnc+tmpSel ).Fill( ireco['tau21'], self.triggerWeight )
                            getattr( self, 'reco'+iRJ+'_tau32'+sysUnc+tmpSel ).Fill( ireco['tau32'], self.triggerWeight )
                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+tmpSel ).Fill( ireco['0p5'+str(tauN)], self.triggerWeight )
                                getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sysUnc+tmpSel ).Fill( ireco['1'+str(tauN)], self.triggerWeight )
                                getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sysUnc+tmpSel ).Fill( ireco['2'+str(tauN)], self.triggerWeight )

            if self.isMC:
                if passGenSel:  ##### go to matrix
                    self.response= self.response+1

                    tmpGenJet1 = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
                    tmpGenJet2 = self.createNsubBasis( selGenJets[1], event, 'GenCands' )

                    if abs(tmpGenJet1['jet'].eta) > abs(tmpGenJet2['jet'].eta):
                        genJet['Jet1'] = tmpGenJet1
                        genJet['Jet2'] = tmpGenJet2
                    else:
                        genJet['Jet1'] = tmpGenJet2
                        genJet['Jet2'] = tmpGenJet1

                    if ( iGenSel==iRecoSel ):

                        for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                            getattr( self, 'accepgen'+iGJ+'_pt'+iRecoSel ).Fill( igen['jet'].pt, event.genWeight )
                            getattr( self, 'accepgen'+iGJ+'_eta'+iRecoSel ).Fill( igen['jet'].eta, event.genWeight )
                            getattr( self, 'accepgen'+iGJ+'_mass'+iRecoSel ).Fill( igen['jet'].mass, event.genWeight )
                            getattr( self, 'accepgen'+iGJ+'_phi'+iRecoSel ).Fill( igen['jet'].phi, event.genWeight )
                            getattr( self, 'accepgen'+iGJ+'_tau21'+iRecoSel ).Fill( igen['tau21'], event.genWeight )
                            getattr( self, 'accepgen'+iGJ+'_tau32'+iRecoSel ).Fill( igen['tau32'], event.genWeight )

                            getattr( self, 'gen'+iGJ+'_pt'+iRecoSel ).Fill( igen['jet'].pt, event.genWeight )
                            getattr( self, 'gen'+iGJ+'_eta'+iRecoSel ).Fill( igen['jet'].eta, event.genWeight )
                            getattr( self, 'gen'+iGJ+'_mass'+iRecoSel ).Fill( igen['jet'].mass, event.genWeight )
                            getattr( self, 'gen'+iGJ+'_phi'+iRecoSel ).Fill( igen['jet'].phi, event.genWeight )
                            getattr( self, 'gen'+iGJ+'_tau21'+iRecoSel ).Fill( igen['tau21'], event.genWeight )
                            getattr( self, 'gen'+iGJ+'_tau32'+iRecoSel ).Fill( igen['tau32'], event.genWeight )

                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'accepgen'+iGJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( igen['2'+str(tauN)], event.genWeight )

                                getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                                getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                                getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( igen['2'+str(tauN)], event.genWeight )

                            WEIGHT =  self.totalWeight
                            for sysUnc in self.sysSource:
                                if ( getattr(ireco['jet'], 'pt'+(sysUnc if not sysUnc.startswith('_pu') else '_nom')) > self.minLeadAK8JetPtDijet ):

                                    if sysUnc.startswith('_pu'):
                                        WEIGHT = event.genWeight * getattr( event, 'puWeight'+sysUnc.split('pu')[1] )
                                        getattr( self, 'truereco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                    else:
                                        getattr( self, 'truereco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt'+sysUnc ), WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sysUnc ), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_eta'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_phi'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau21'+sysUnc+iRecoSel ).Fill( ireco['tau21'], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau32'+sysUnc+iRecoSel ).Fill( ireco['tau32'], WEIGHT )
                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'truereco'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['2'+str(tauN)], WEIGHT )

                                #### FIlling response matrices
                                if self.isMC:
                                    if sysUnc.endswith('nom'):
                                        getattr( self, 'resol'+iRJ+'_pt'+iRecoSel ).Fill( getattr(ireco['jet'], 'pt')/getattr(igen['jet'], 'pt'), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_sdmass'+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop')/getattr(igen['jet'], 'mass'), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau21'+iRecoSel ).Fill( ireco['tau21']/igen['tau21'], WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau32'+iRecoSel ).Fill( ireco['tau32']/igen['tau32'], WEIGHT )

                                    getattr( self, 'resp'+iRJ+'_tau21'+sysUnc+iRecoSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau32'+sysUnc+iRecoSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                                        getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )
                                        if sysUnc.endswith('nom'):
                                            getattr( self, 'resol'+iRJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( ( ireco['0p5'+str(tauN)]/igen['0p5'+str(tauN)] if igen['0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resol'+iRJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( ( ireco['1'+str(tauN)]/igen['1'+str(tauN)] if igen['1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resol'+iRJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( ( ireco['2'+str(tauN)]/igen['2'+str(tauN)] if igen['2'+str(tauN)]!=0 else -999 ), WEIGHT )


                    else: self.ufoResponse = self.ufoResponse+1

                else:           ##### fakes
                    self.fakes = self.fakes+1

                    for iRJ,ireco in recoJet.items():

                        WEIGHT =  self.totalWeight
                        for sysUnc in self.sysSource:
                            if ( getattr(ireco['jet'], 'pt'+(sysUnc if not sysUnc.startswith('_pu') else '_nom')) > self.minLeadAK8JetPtDijet ):

                                if sysUnc.startswith('_pu'):
                                    WEIGHT = event.genWeight * getattr( event, 'puWeight'+sysUnc.split('pu')[1] )
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                else:
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt'+sysUnc ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sysUnc ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_eta'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_phi'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau21'+sysUnc+iRecoSel ).Fill( ireco['tau21'], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau32'+sysUnc+iRecoSel ).Fill( ireco['tau32'], WEIGHT )
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['2'+str(tauN)], WEIGHT )

        else:  #### Misses
            self.miss = self.miss+1

            tmpGenJet1 = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
            tmpGenJet2 = self.createNsubBasis( selGenJets[1], event, 'GenCands' )

            if abs(tmpGenJet1['jet'].eta) > abs(tmpGenJet2['jet'].eta):
                genJet['Jet1'] = tmpGenJet1
                genJet['Jet2'] = tmpGenJet2
            else:
                genJet['Jet1'] = tmpGenJet2
                genJet['Jet2'] = tmpGenJet1

            if passGenSel:

                for iGJ,igen in genJet.items():

                    getattr( self, 'missgen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )

                    getattr( self, 'missgen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )

                    getattr( self, 'gen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )

                    getattr( self, 'gen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                    getattr( self, 'gen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )

                    for tauN in range(1, self.maxTau+1):
                        getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )

                        getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                        getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                        getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )

            else:
                self.ufo = self.ufo+1

        self.fillBranches( 'selRecojets', recoJet )
        if self.isMC: self.fillBranches( 'selGenjets', genJet )

        return True


    #############################################################################
    def recoSelection( self, event ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))

        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minTightMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.pfIsoId>=2 and x.tightId and abs(x.dxy)<0.2 and abs(x.dz)<0.5 and x.miniPFRelIso_all<0.10] # applying tight selection on muons already here since we only veto for loose muons and loose electrons in event

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons)

        ##################################################

        #### Basic AK8 jet selection
        recoAK8jets = [ x for x in AK8jets if x.pt_nom > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta and (x.jetId >= 2)]
        recoAK8jets.sort(key=lambda x:x.pt_nom,reverse=True)
        ##################################################

        ##### Applying selection
        passSel, ptAsym, deltaPhi = self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets )
        iSel = '_dijetSel' if passSel else None

        ##### Trigger weights
        if not self.isMC:
            if passSel:
                self.triggerWeight = 0
                triggerVersion = -9999
                for i in range( len(self.runTables[ self.year ]['low'] ) ):
                    if ( self.runTables[self.year]['low'][i] <= event.run ) and ( event.run <= self.runTables[self.year]['high'][i] ):
                        triggerVersion = 2+i

                for itrigger, itrValues in self.triggerTable.items():
                    if ( getattr(event, 'HLT_'+itrigger)==1 ):

                        if ( itrValues[self.year][0] < recoAK8jets[0].pt < itrValues[self.year][1] ) :
                            iSel = '_'+itrigger+'_dijetSel'
                            passSel = True
                            self.triggerWeight = itrValues[self.year][triggerVersion]
                            break
                        else:
                            passSel = False
                            iSel = None
                            self.triggerWeight = 0

#
#                        #print event.run, event.luminosityBlock, event.event
#                        tmpWeight = itrValues[self.year][triggerVersion]
#                        #triggerPrescale = self.triggerInfo.loc[ ( (self.triggerInfo.run==event.run) & (self.triggerInfo.lumi==event.luminosityBlock) & (abs(self.triggerInfo.event)==event.event ) ) ]
#                        #print self.triggerInfo.loc[ ( self.triggerInfo.run==event.run ) & (self.triggerInfo.lumi==event.luminosityBlock) ]
#                        #try: tmpWeight = 1 #triggerPrescale.iloc[0][ itrigger.split('_')[1] ]
#                        #except IndexError: break
#                        if (tmpWeight==1):
#                            triggerWeight = tmpWeight
#                            #print( 'passOne', event.event, itrigger, getattr(event, itrigger), recoAK8jets[0].pt, itrValues[self.year], triggerWeight  )
#                            break
#                        elif ( itrValues[self.year][0] < recoAK8jets[1].pt < itrValues[self.year][1] ) :
#                            triggerWeight = tmpWeight # triggerPrescale.iloc[0][ itrigger.split('_')[1] ]
#                            #print( 'pass', event.event, itrigger, getattr(event, itrigger), recoAK8jets[0].pt, itrValues[self.year], triggerWeight  )
#                            break
#                        else:
#                            #print( 'fail', event.event, itrigger, getattr(event, itrigger), itrValues[self.year][0], recoAK8jets[0].pt, itrValues[self.year][1], itrValues[self.year], triggerWeight  )
#                            pass



                    else: passSel = False
            else:
                self.triggerWeight = 0
                passSel = False
                iSel = None
        #print '2', iSel
        self.out.fillBranch( 'triggerWeight', self.triggerWeight )
        ##################################################


        #### Weight #########
        if self.isMC:
            weight = event.puWeight * event.genWeight
            getattr( self, 'PUweight' ).Fill( event.puWeight )
        else:
            weight = 1 #self.triggerWeight
        self.totalWeight = weight
        ##################################################

        #### Checking no selection without weights
        getattr( self, 'nPVs_noSelnoWeight' ).Fill( getattr( event, 'PV_npvsGood') )
        getattr( self, 'nleps_noSelnoWeight' ).Fill( nleptons )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSelnoWeight' ).Fill( imuon.pt )
            getattr( self, 'muons_eta_noSelnoWeight' ).Fill( imuon.eta )
            getattr( self, 'muons_phi_noSelnoWeight' ).Fill( imuon.phi )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSelnoWeight' ).Fill( iele.pt )
            getattr( self, 'eles_eta_noSelnoWeight' ).Fill( iele.eta )
            getattr( self, 'eles_phi_noSelnoWeight' ).Fill( iele.phi )
        getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recoAK8jets) )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
            getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop )

        #### Checking no selection with weights
        getattr( self, 'nPVs_noSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
        getattr( self, 'nleps_noSel' ).Fill( nleptons, weight )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt, weight )
            getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta, weight )
            getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi, weight )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSel' ).Fill( iele.pt, weight )
            getattr( self, 'eles_eta_noSel' ).Fill( iele.eta, weight )
            getattr( self, 'eles_phi_noSel' ).Fill( iele.phi, weight )
        getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets), weight )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop, weight )
        getattr( self, 'recoPtAsym_noSel' ).Fill( ptAsym, weight )
        getattr( self, 'recoDeltaPhi_noSel' ).Fill( deltaPhi, weight )


        reweight = self.totalWeight
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        if passSel and iSel:

            #### Basic reco histos
            getattr( self, 'nPVs'+iSel ).Fill( getattr( event, 'PV_npvsGood'), reweight )
            getattr( self, 'nleps'+iSel ).Fill( len(recoMuons)+len(recoElectrons), reweight )
            for imuon in recoMuons:
                getattr( self, 'muons_pt'+iSel ).Fill( imuon.pt, reweight )
                getattr( self, 'muons_eta'+iSel ).Fill( imuon.eta, reweight )
                getattr( self, 'muons_phi'+iSel ).Fill( imuon.phi, reweight )
            for iele in recoElectrons:
                getattr( self, 'eles_pt'+iSel ).Fill( iele.pt, reweight )
                getattr( self, 'eles_eta'+iSel ).Fill( iele.eta, reweight )
                getattr( self, 'eles_phi'+iSel ).Fill( iele.phi, reweight )
            getattr( self, 'nAK8jets'+iSel ).Fill( len(recoAK8jets), reweight )
            for ijet in recoAK8jets:
                getattr( self, 'AK8jets_pt'+iSel ).Fill( ijet.pt, reweight )
                getattr( self, 'AK8jets_eta'+iSel ).Fill( ijet.eta, reweight )
                getattr( self, 'AK8jets_phi'+iSel ).Fill( ijet.phi, reweight )
                getattr( self, 'AK8jets_mass'+iSel ).Fill( ijet.msoftdrop, reweight )
            getattr( self, 'recoPtAsym'+iSel ).Fill( ptAsym, reweight )
            getattr( self, 'recoDeltaPhi'+iSel ).Fill( deltaPhi, reweight )

        return passSel, iSel, recoMuons, recoElectrons, recoAK8jets

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing reco information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genParticles = Collection(event, 'GenPart')

        ### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        ##################################################

        #### Basic AK8 jet selection
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta ]
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        passSel, ptAsym, deltaPhi = self.dijetSelection( event, genMuons, genElectrons, genAK8jets )
        iSel = '_dijetSel' if passSel else None
        #### Weight
        weight = event.genWeight

        #### Checking no selection
        getattr( self, 'ngenleps_noSel' ).Fill( len(genMuons)+len(genElectrons), weight )
        for imuon in genMuons:
            getattr( self, 'genmuons_pt_noSel' ).Fill( imuon.pt, weight )
            getattr( self, 'genmuons_eta_noSel' ).Fill( imuon.eta, weight )
            getattr( self, 'genmuons_phi_noSel' ).Fill( imuon.phi, weight )
        for iele in genElectrons:
            getattr( self, 'geneles_pt_noSel' ).Fill( iele.pt, weight )
            getattr( self, 'geneles_eta_noSel' ).Fill( iele.eta, weight )
            getattr( self, 'geneles_phi_noSel' ).Fill( iele.phi, weight )
        getattr( self, 'ngenAK8jets_noSel' ).Fill( len(genAK8jets), weight )
        for ijet in genAK8jets:
            getattr( self, 'AK8genjets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK8genjets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK8genjets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'AK8genjets_mass_noSel' ).Fill( ijet.mass, weight )
        getattr( self, 'genPtAsym_noSel' ).Fill( ptAsym, weight )
        getattr( self, 'genDeltaPhi_noSel' ).Fill( deltaPhi, weight )


        ##### Filling histograms
        if passSel and iSel:

            #### Checking nominal selection with weights
            getattr( self, 'ngenleps'+iSel ).Fill( len(genMuons)+len(genElectrons), weight )
            for imuon in genMuons:
                getattr( self, 'genmuons_pt'+iSel ).Fill( imuon.pt, weight )
                getattr( self, 'genmuons_eta'+iSel ).Fill( imuon.eta, weight )
                getattr( self, 'genmuons_phi'+iSel ).Fill( imuon.phi, weight )
            for iele in genElectrons:
                getattr( self, 'geneles_pt'+iSel ).Fill( iele.pt, weight )
                getattr( self, 'geneles_eta'+iSel ).Fill( iele.eta, weight )
                getattr( self, 'geneles_phi'+iSel ).Fill( iele.phi, weight )
            getattr( self, 'ngenAK8jets'+iSel ).Fill( len(genAK8jets), weight )
            for ijet in genAK8jets:
                getattr( self, 'AK8genjets_pt'+iSel ).Fill( ijet.pt, weight )
                getattr( self, 'AK8genjets_eta'+iSel ).Fill( ijet.eta, weight )
                getattr( self, 'AK8genjets_phi'+iSel ).Fill( ijet.phi, weight )
                getattr( self, 'AK8genjets_mass'+iSel ).Fill( ijet.mass, weight )
            getattr( self, 'genPtAsym'+iSel ).Fill( ptAsym, weight )
            getattr( self, 'genDeltaPhi'+iSel ).Fill( deltaPhi, weight )

        return passSel, iSel, genMuons, genElectrons, genAK8jets

    #############################################################################
    def dijetSelection( self, event, muons, electrons, AK8jets ):

        if (len(muons)+len(electrons)==0) and (len(AK8jets)>1):
            ptAsym = ( AK8jets[0].pt - AK8jets[1].pt ) / (AK8jets[0].pt + AK8jets[1].pt)
            deltaPhi = AK8jets[0].p4().DeltaPhi( AK8jets[1].p4() )
            if (ptAsym<0.3) and (deltaPhi>2) and (AK8jets[1].pt>self.minLeadAK8JetPtDijet):         ### PT REQUIREMENT NEEDS TO BE REVISIT
                return True, ptAsym, deltaPhi
            else: return False, -999, -999
        else: return False, -999, -999

    #############################################################################
    def createNsubBasis(self, AK8jet, event, PFCollection ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        #### Run calculations of NSub bases and store for ungroomed AK8jets (default in CMS)

        #### Applying PUPPI weights to the PF candidates
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()
        for p in pfCands :
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            try: tp = tp * p.puppiWeight
            except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        for x in CandsPUPPIweightedVec:
            if abs(AK8jet.p4().DeltaR( x )) < 0.8: constituents.push_back(x)

        #### Computing n-subjetiness basis from PF PUPPI constituents
        nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
        nsub1 = self.nSub1.getTau( self.maxTau, constituents )
        nsub2 = self.nSub2.getTau( self.maxTau, constituents )
        nsub1_OP_kT = self.nSub1_OP_kT.getTau( 3, constituents )    ### needed for genjet tau21 tau32


        ### default in CMS OP_KT https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        try: ak8jet['tau21'] = nsub1_OP_kT[1]/nsub1_OP_kT[0]
        except ZeroDivisionError: ak8jet['tau21'] = -1
        try: ak8jet['tau32'] = nsub1_OP_kT[2]/nsub1_OP_kT[1]
        except ZeroDivisionError: ak8jet['tau32'] = -1


        #### filling histos and branches with nsub basis
        for tauN in range(self.maxTau):
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]

        return ak8jet

    #############################################################################
    def fillBranches( self, jetLabel, jetInfo ):

        #### Filling branch with passAK8jet info after selection
        self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
        self.out.fillBranch(jetLabel+"_mass", [ iJ['jet'].mass for i,iJ in jetInfo.items() ] )
        self.out.fillBranch(jetLabel+"_pt", [ iJ['jet'].pt for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_eta", [ iJ['jet'].eta for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_phi", [ iJ['jet'].phi for i,iJ in jetInfo.items()  ] )

        for tauN in range(1, self.maxTau+1):
            self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN),  [ iJ['0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN),  [ iJ['1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN),  [ iJ['2'+str(tauN)] for i,iJ in jetInfo.items() ] )



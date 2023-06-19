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

    def __init__(self, sysSource=[], year='2017', isMC=True, onlyUnc='', onlyTrees=False):
        self.writeHistFile=True
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.onlyTrees = onlyTrees
        self.runSDVariables = False

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

        ### Kinematics Cuts Jets ###
        self.minJetPt = 30.
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

        self.totalRecoWeight = 1.
        self.totalGenWeight = 1.
        #self.totalWeight = 1.
        self.triggerWeight = 1.
        self.puWeight = 1.

        self.pdfWeightUp = 1.
        self.puWeightUp = 1.
        self.isrWeightUp = 1.
        self.fsrWeightUp = 1.
        self.pdfWeightDown = 1.
        self.puWeightDown = 1.
        self.isrWeightDown = 1.
        self.fsrWeightDown = 1.

        self.eventCategory = -1
        self.dummy = 0

        #self.runTables = {
        #        '2017' : {
        #            'low' : [ 297046, 299368, 302030, 303824, 305040  ],
        #            'high' : [ 299329, 302029, 303434, 304797, 306462  ]
        #            },
        #        '2018' : {
        #            'low' : [ 315252, 317080, 319337, 320673  ],
        #            'high' : [ 316995, 319310, 320065, 325175 ]
        #            }
        #        }

        self.triggerTable = OrderedDict()
        self.triggerTable[ 'AK8PFJet80' ] = {    #### from the list below, only the first two numbers (trigger turn on) are used. The others were a test.
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
                    '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                    '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
                    }

        '''
                                self.triggerTable = OrderedDict()
                                self.triggerTable[ 'AK8PFJet80' ] = {
                                            '2017' : [ 200,   257, 50510.97, 6727.43, 6854.43, 23381.99, 38371.17 ],
                                            '2018' : [ 200,   267, 17046.53, 34000.07, 33988.86, 33998.78 ],
                                            }
                                self.triggerTable[ 'AK8PFJet140' ] = {
                                            '2017' : [ 257,   323, 672.78, 1644.76, 1255.72, 2027.79, 2315.59 ],
                                            '2018' : [ 267,   332, 1601.95, 1207.43, 1220.84, 1184.09 ],
                                            }
                                self.triggerTable[ 'AK8PFJet200' ] = {
                                            '2017' : [ 323,   389, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                                            '2018' : [ 332,   398, 368.33, 281.73, 284.86, 276.30 ],
                                            }
                                self.triggerTable[ 'AK8PFJet260' ] = {
                                            '2017' : [ 389,   455, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                                            '2018' : [ 398,   464, 132.00, 128.00, 122.91, 127.42 ],
                                            }
                                self.triggerTable[ 'AK8PFJet320' ] = {
                                            '2017' : [ 455,   543, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                                            '2018' : [ 464,   551, 49.29, 48.00, 47.28, 47.92 ],
                                            }
                                self.triggerTable[ 'AK8PFJet400' ] = {
                                            '2017' : [ 543,   598, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                                            '2018' : [ 551,   606, 16.39, 16.00, 15.92, 15.99 ],
                                            }
                                self.triggerTable[ 'AK8PFJet450' ] = {
                                            '2017' : [ 598,   653, 54.41, 382.72, 263.79, 380.97, 410.70  ],
                                            '2018' : [ 606,   661, 8.43, 8.00, 7.98, 8.00 ],
                                            }
                                self.triggerTable[ 'AK8PFJet500' ] = {
                                            '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                                            '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
                                            }
        '''
        ### Defining nsubjetiness basis
        self.maxTau = 4
        self.nSub_labels = {
                        "_tau_0p5_1": [ 0., 1.1, 1100  ],
                        "_tau_0p5_2": [ 0., 0.9, 900  ],
                        "_tau_0p5_3": [ 0., 0.8, 800  ],
                        "_tau_0p5_4": [ 0., 0.7, 700  ],
                        "_tau_0p5_5": [ 0., 0.7, 700  ],
                        "_tau_1_1": [ 0., 1.1, 1100  ],
                        "_tau_1_2": [ 0., 0.7, 700  ],
                        "_tau_1_3": [ 0., 0.5, 500  ],
                        "_tau_1_4": [ 0., 0.5, 500  ],
                        "_tau_1_5": [ 0., 0.4, 400  ],
                        "_tau_2_1": [ 0., 0.5, 500  ],
                        "_tau_2_2": [ 0., 0.3, 600  ],
                        "_tau_2_3": [ 0., 0.2, 400  ],
                        "_tau_2_4": [ 0., 0.2, 400  ],
                        "_tau_2_5": [ 0., 0.4, 800  ]
                }
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )
        self.nSub1_WTA_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 ) ##### needed for genjet tau21 or tau32

        ### Softdrop quantities
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, self.minAK8JetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")

        ### Helpers
        self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
        self.nJet = [ 'Jet'] #1', 'Jet2' ]
        if self.runSDVariables: self.nJet = [ 'Jet']#, 'Jet2', 'sdJet1', 'sdJet2' ]

        ### Uncerstinties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        if onlyUnc: self.sysSource = [ onlyUnc+i for i in [ 'Up', 'Down' ] ]
        self.sysWeightList = ( '_pu', '_pdf', '_ps', '_isr', '_fsr' )

    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        tauBins = 100

        #if not self.onlyTrees:
        ### Booking histograms
        selList = [ '_'+x+'_dijetSel' for x in self.triggerTable  ] if not self.isMC else [ '_dijetSel' ]
        #selList = ([ '_'+x+'_dijetSel' for x in self.triggerTable  ] + [ '_'+x+'_weight_dijetSel' for x in self.triggerTable  ]) if not self.isMC else [ '_dijetSel' ]
        if not self.onlyUnc:
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   12, 0, 12) )
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            self.addObject( ROOT.TH1F('Btagweight',   ';BtagWeight',   25, 0, 2) )
            if not self.isMC:
                for isel in [ '_only'+x+'_dijetSel' for x in self.triggerTable  ]:
                    self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                    self.addP4Hists( 'LeadJetAK8', isel )

            #### general selection
            for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
                self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
                self.addP4Hists( 'muons', isel )
                self.addP4Hists( 'eles', isel )
                self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
                self.addP4Hists( 'AK8jets', isel )
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
                    self.addObject( ROOT.TH1F('genPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                    self.addObject( ROOT.TH1F('genDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )
                    self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )
                    self.addObject(ROOT.TH1F('genJetAK8_partonFlavour'+isel,   ';genpartonFlav (MCpns)',   50, -25., 25.) )

        for isel in selList:
            for iJ in self.nJet:
                #self.addP4Hists( 'uforeco'+iJ, isel )
                #self.addP4Hists( 'ufogen'+iJ, isel )
                self.addP4Hists( 'reco'+iJ, isel ) #+'_sortedPt'
                #if self.isMC: self.addP4Hists( 'gen'+iJ+'_sortedPt', isel )
                for itype in ( [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] ):

                    if itype.endswith(('reco', 'gen')):

                        for sysUnc in self.sysSource:
                            self.addP4Hists( itype+iJ, sysUnc+isel )
                            binning = np.array([(i/2000) for i in np.arange(0*2000, 2.*2000)])
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{21}', len(binning)-1, binning) )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{32}', len(binning)-1, binning) )

                            if itype.startswith('truereco') and self.isMC:
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau21'+sysUnc+isel, ';AK8 accepgen jet #tau_{21};AK8 truereco jet #tau_{21}', len(binning)-1, binning, len(binning)-1, binning) )
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau32'+sysUnc+isel, ';AK8 accepgen jet #tau_{32};AK8 truereco jet #tau_{32}', len(binning)-1, binning, len(binning)-1, binning) )
                                if sysUnc.endswith('nom'):
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_pt'+isel, ';AK8 reco/gen jet pt', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_sdmass'+isel, ';AK8 reco m_{SD}/gen jet m', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_mass'+isel, ';AK8 reco inv. m/gen jet inv. m', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau21'+isel, ';AK8 reco/gen jet #tau_{21}', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau32'+isel, ';AK8 reco/gen jet #tau_{32}', 500, 0, 5) )

                            for x, y in self.nSub_labels.items():
                                binning = np.array([(i/1000) for i in np.arange(y[0]*1000, y[1]*1000)])
                                self.addObject( ROOT.TH1F(itype+iJ+x+sysUnc+isel, ';AK8 '+itype+' jet #tau', len(binning)-1, binning ) )
                                if itype.startswith('reco') and self.isMC:
                                    self.addObject( ROOT.TH2F('resp'+iJ+x+sysUnc+isel, ';AK8 gen jet '+x+';AK8 reco jet '+x, len(binning)-1, binning, len(binning)-1, binning ) )
                                    if sysUnc.endswith('nom'): self.addObject( ROOT.TH1F('resol'+iJ+x+isel, ';AK8 reco/gen jet '+x, 500, 0, 5) )
                    else:

                        self.addP4Hists( itype+iJ, isel )
                        binning = np.array([(i/2000) for i in np.arange(0*2000, 2.*2000)])
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+isel, ';AK8 '+itype+' jet #tau_{21}', len(binning)-1, binning) )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+isel, ';AK8 '+itype+' jet #tau_{32}', len(binning)-1, binning) )

                        for x, y in self.nSub_labels.items():
                            binning = np.array([(i/1000) for i in np.arange(y[0]*1000, y[1]*1000)])
                            self.addObject( ROOT.TH1F(itype+iJ+x+isel, ';AK8 '+itype+' jet '+x, len(binning)-1, binning ) )
        '''
        if self.isMC:
            extSelList = ['_'+str(x)+'_dijetSel' for x in self.diffPt if x<10000 ]
        else:
            extSelList = [ '_'+str(y)+'_'+x+'_dijetSel' for x in self.triggerTable  for y in self.diffPt if y<10000 ]  #+ [ '_'+str(y)+'_'+x+'_weight_dijetSel' for x in self.triggerTable for y in self.diffPt if y<10000 ]
        for isel in extSelList:
            for iJ in self.nJet:
                for itype in ( [ 'gen', 'reco' ] if self.isMC else [ 'reco' ] ):

                    if itype.endswith('reco'):

                        for sysUnc in self.sysSource:
                            self.addObject( ROOT.TH1F(itype+iJ+'_pt'+sysUnc+isel, ';AK8 '+itype+' jet pt', 200, 0, 2000) )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{21}', tauBins, 0, 1) )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{32}', tauBins, 0, 1) )

                            if self.isMC:
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau21'+sysUnc+isel, ';AK8 accepgen jet #tau_{21};AK8 truereco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau32'+sysUnc+isel, ';AK8 accepgen jet #tau_{32};AK8 truereco jet #tau_{32}', tauBins, 0, 1, tauBins, 0, 1) )

                            for x, y in self.nSub_labels.items():
                                self.addObject( ROOT.TH1F(itype+iJ+x+sysUnc+isel, ';AK8 '+itype+' jet #tau', y[1], 0, y[0] ) )
                                if itype.startswith('reco') and self.isMC:
                                    self.addObject( ROOT.TH2F('resp'+iJ+x+sysUnc+isel, ';AK8 gen jet '+x+';AK8 reco jet '+x, y[1], 0, y[0], y[1], 0, y[0] ) )
                    else:

                        self.addObject( ROOT.TH1F(itype+iJ+'_pt'+isel, ';AK8 '+itype+' jet pt', 200, 0, 2000) )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+isel, ';AK8 '+itype+' jet #tau_{21}', tauBins, 0, 1) )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+isel, ';AK8 '+itype+' jet #tau_{32}', tauBins, 0, 1) )

                        for x, y in self.nSub_labels.items():
                            self.addObject( ROOT.TH1F(itype+iJ+x+isel, ';AK8 '+itype+' jet '+x, y[1], 0, y[0] ) )
        '''
    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  ';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, ';#eta', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_y'+t, ';y', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, ';#phi', 100, -3.5, 3.5) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,';mass (GeV)', 250, 0, 1000) )
        if 'reco' in s.lower(): self.addObject( ROOT.TH1F(s+'_msoftdrop'+t,';Softdrop mass (GeV)', 100, 0, 1000) )


    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #tmp = inputFile.Get('Runs')
        #for i in tmp.LHEPdfSumw: print(i)

        #if self.onlyTrees:
        self.out = wrappedOutputTree
        self.out.branch('triggerWeight',  "F")
        self.out.branch('eventCategory',  "I")
        self.out.branch( 'totalRecoWeight', "F" )
        self.out.branch( 'totalGenWeight', "F" )
        self.out.branch('puWeight',  "F")

        if self.onlyUnc:
            if 'pdf' in self.onlyUnc:
                self.out.branch( 'pdfWeightUp', "F" )
                self.out.branch( 'pdfWeightDown', "F" )

            if 'isr' in self.onlyUnc:
                self.out.branch( 'isrWeightUp', "F" )
                self.out.branch( 'isrWeightDown', "F" )

            if 'fsr' in self.onlyUnc:
                self.out.branch( 'fsrWeightUp', "F" )
                self.out.branch( 'fsrWeightDown', "F" )

            if 'puWeight' in self.onlyUnc:
                self.out.branch( 'puWeightUp', "F" )
                self.out.branch( 'puWeightDown', "F" )

        tmplist = [ 'selRecoJets'+sys for sys in self.sysSource ]
        #print ('selrecojetlist:', tmplist)

        if self.isMC: 
            X = ['trueRecoJets'+sys for sys in self.sysSource]
            

            #print ("appended MC reco/genjets to nano output file")
        
            if not self.onlyUnc: 
                tmplist.append(X[0])
                tmplist.append( 'selGenJets_nom' )
                tmplist.append( 'accepGenJets_nom' )
            else:
                tmpGenList = [ 'selGenJets'+sys for sys in self.sysSource ]
                tmpAccepGenList = [ 'accepGenJets'+sys for sys in self.sysSource ] 
                for x in X:
                    tmplist.append(x)
                for x in tmpGenList:
                    tmplist.append(x)
                for x in tmpAccepGenList:
                    tmplist.append(x)
        print ("Stored tree branches:", tmplist)
        for iJ in tmplist:
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_y',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_msoftdrop',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_tau21',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_tau32',  'F', lenVar='n'+iJ)
            if self.runSDVariables:
                self.out.branch(iJ+'_SD_tau21',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_SD_tau32',  'F', lenVar='n'+iJ)

            for x in self.nSub_labels:
                self.out.branch(iJ+x, 'F', lenVar='n'+iJ )
                if self.runSDVariables:
                    self.out.branch(iJ+'SD'+x, 'F', lenVar='n'+iJ )
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

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenJets = self.genSelection(event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoJets = self.recoSelection( event )
        if not self.isMC and not passRecoSel['_nom']: return False

        for sys in self.sysSource:

            if sys.startswith(self.sysWeightList):
                selRecoJets[sys] = selRecoJets['_nom']
                iRecoSel[sys] = iRecoSel['_nom']
                passRecoSel[sys] = passRecoSel['_nom']

                if sys.startswith('_pdfWeight'):
                    arr =  getattr( event, 'LHEPdfWeight' )
                    size = arr.GetSize() # force loading of array contents for this entry
                    buf = arr.GetAddress() # still not sure this is always equivalent to `&arr[0]`
                    buf.SetSize(size)
                    npLHEPdfWeight = np.frombuffer(buf, count=size, dtype='float32')
                    pdfWeight = np.std(npLHEPdfWeight)

                    if sys.endswith('Up'):
                        self.pdfWeightUp = 1 + pdfWeight
                        self.totalGenWeight = event.genWeight * self.pdfWeightUp
                        WEIGHT = self.totalRecoWeight * self.pdfWeightUp 
                    else:
                        self.pdfWeightDown = 1 - pdfWeight
                        self.totalGenWeight = event.genWeight * self.pdfWeightDown
                        WEIGHT = self.totalRecoWeight * self.pdfWeightDown 

                elif sys.endswith('isrWeightDown'):
                    self.isrWeightDown = getattr( event, 'PSWeight' )[0]
                    self.totalGenWeight = event.genWeight * self.isrWeightDown
                    WEIGHT = self.totalRecoWeight * self.isrWeightDown 

                elif sys.endswith('fsrWeightDown'):
                    self.fsrWeightDown = getattr( event, 'PSWeight' )[1]
                    self.totalGenWeight = event.genWeight * self.fsrWeightDown
                    WEIGHT = self.totalRecoWeight * self.fsrWeightDown 

                elif sys.endswith('isrWeightUp'):
                    self.isrWeightUp = getattr( event, 'PSWeight' )[2]
                    self.totalGenWeight = event.genWeight * self.isrWeightUp
                    WEIGHT = self.totalRecoWeight * self.isrWeightUp 
                elif sys.endswith('fsrWeightUp'):
                    self.fsrWeightUp = getattr( event, 'PSWeight' )[3]
                    self.totalGenWeight = event.genWeight * self.fsrWeightUp
                    WEIGHT = self.totalRecoWeight * self.fsrWeightUp 

                elif sys.endswith('puWeightUp'):
                    self.puWeightUp = event.puWeightUp
                    self.totalGenWeight = event.genWeight 
                    WEIGHT = self.totalRecoWeight * self.puWeightUp / self.puWeight #divide out puweight nominal in this case, since that's multiplied into the weight in recoSelection() for the nominal case and other non-PU variations
                elif sys.endswith('puWeightDown'):
                    self.puWeightDown = event.puWeightDown
                    self.totalGenWeight = event.genWeight 
                    WEIGHT = self.totalRecoWeight * self.puWeightDown / self.puWeight
                #elif sys.startswith('_btag'):
                #    self.totalGenWeight = event.genWeight
                #    WEIGHT = self.totalRecoWeight # the up down modulation of btagweights is controlled in the WtopSel function
            else: 
                
                if self.isMC: 
                    self.totalGenWeight = event.genWeight#*self.triggerWeight
                    #self.totalRecoWeight = self.totalRecoWeight#*self.triggerWeight
                    WEIGHT =  self.totalRecoWeight
                else:
                    WEIGHT =  self.totalRecoWeight




            genJet = OrderedDict()
            recoJet = OrderedDict()
            tmpRecoJets = OrderedDict()

            if self.isMC and not passGenSel and not passRecoSel[sys]: return False

            if self.isMC: 
                deltaRmatch = False


            if passRecoSel[sys]:  #### Detector level dist.
                #print ( self.totalRecoWeight,  event.puWeight, event.genWeight, event.genWeight*event.puWeight, event.genWeight*self.puWeight, self.puWeight, self.triggerWeight)

                tmpRecoJets[sys] = {}
                if sys.startswith('_nom'):

                    tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                    #tmpRecoJets[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )

                    #if not self.onlyTrees:
                    getattr( self, 'recoJet_pt'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'pt_nom' ), self.totalRecoWeight )
                    getattr( self, 'recoJet_eta'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'eta' ), self.totalRecoWeight)
                    getattr( self, 'recoJet_y'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'y' ), self.totalRecoWeight)
                    getattr( self, 'recoJet_phi'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'phi' ), self.totalRecoWeight )
                    getattr( self, 'recoJet_mass'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'mass' ), self.totalRecoWeight )
                    getattr( self, 'recoJet_msoftdrop'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'msoftdrop_nom' ), self.totalRecoWeight )
                    #getattr( self, 'recoJet2_sortedPt_pt'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'pt_nom' ), self.totalWeight )
                    #getattr( self, 'recoJet2_sortedPt_eta'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'eta' ), self.totalWeight )
                    #getattr( self, 'recoJet2_sortedPt_phi'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'phi' ), self.totalWeight )
                    #getattr( self, 'recoJet2_sortedPt_mass'+iRecoSel['_nom'] ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'msoftdrop_nom' ), self.totalWeight )
#                        if not self.isMC:
#                            tmpSel = iRecoSel['_nom'].replace('_dijet', '_weight_dijet')
#                            getattr( self, 'recoJet1_sortedPt_pt'+tmpSel ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'pt_nom' ), self.triggerWeight )
#                            getattr( self, 'recoJet1_sortedPt_eta'+tmpSel ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'eta' ), self.triggerWeight )
#                            getattr( self, 'recoJet1_sortedPt_phi'+tmpSel ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'phi' ), self.triggerWeight )
#                            getattr( self, 'recoJet1_sortedPt_mass'+tmpSel ).Fill( getattr(tmpRecoJets[sys][0]['jet'], 'msoftdrop_nom' ), self.triggerWeight )
#                            getattr( self, 'recoJet2_sortedPt_pt'+tmpSel ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'pt_nom' ), self.triggerWeight )
#                            getattr( self, 'recoJet2_sortedPt_eta'+tmpSel ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'eta' ), self.triggerWeight )
#                            getattr( self, 'recoJet2_sortedPt_phi'+tmpSel ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'phi' ), self.triggerWeight )
#                            getattr( self, 'recoJet2_sortedPt_mass'+tmpSel ).Fill( getattr(tmpRecoJets[sys][1]['jet'], 'msoftdrop_nom' ), self.triggerWeight )
                else:

                    if '_nom' in tmpRecoJets.keys():
                        tmpRecoJets[sys][0] = tmpRecoJets['_nom'][0]
                        #tmpRecoJets[sys][1] = tmpRecoJets['_nom'][1]
                    else:
                        tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                        #tmpRecoJets[sys][1] = self.createNsubBasis( selRecoJets[sys][1], event, 'PFCands' )


                #if abs(tmpRecoJets[sys][0]['jet'].rapidity) > abs(tmpRecoJets[sys][1]['jet'].rapidity):
                #    recoJet['Jet1'] = tmpRecoJets[sys][0]
                #    recoJet['Jet2'] = tmpRecoJets[sys][1]
                #else:
                #    recoJet['Jet1'] = tmpRecoJets[sys][1]
                #    recoJet['Jet2'] = tmpRecoJets[sys][0]

                recoJet['Jet'] = tmpRecoJets[sys][0]
                self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel
                self.eventCategory = 1


                
                #if not self.onlyTrees:
                for iRJ,ireco in recoJet.items():

                    if ( getattr(ireco['jet'], 'pt'+( '_nom' if sys.startswith(self.sysWeightList) else sys ) ) < self.minLeadAK8JetPtDijet ): continue


                    tmpSel = iRecoSel[sys] #.replace('_dijet', '_weight_dijet') if not self.isMC else iRecoSel[sys]
                    #tmpWeight = WEIGHT if not self.isMC else self.triggerWeight
                    if sys.startswith(self.sysWeightList):
                        getattr( self, 'reco'+iRJ+'_pt'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_mass'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                    else:
                        getattr( self, 'reco'+iRJ+'_pt'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_mass'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_eta'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_phi'+sys+tmpSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_tau21'+sys+tmpSel ).Fill( ireco['tau21'], WEIGHT )
                    getattr( self, 'reco'+iRJ+'_tau32'+sys+tmpSel ).Fill( ireco['tau32'], WEIGHT )
                    if self.runSDVariables:
                        getattr( self, 'recosd'+iRJ+'_tau21'+sys+tmpSel ).Fill( ireco['sdtau21'], WEIGHT )
                        getattr( self, 'recosd'+iRJ+'_tau32'+sys+tmpSel ).Fill( ireco['sdtau32'], WEIGHT )
                    for tauN in range(1, self.maxTau+1):
                        getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sys+tmpSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sys+tmpSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sys+tmpSel ).Fill( ireco['2'+str(tauN)], WEIGHT )
                        if self.runSDVariables:
                            getattr( self, 'recosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+tmpSel ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                            getattr( self, 'recosd'+iRJ+'_tau_1_'+str(tauN)+sys+tmpSel ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                            getattr( self, 'recosd'+iRJ+'_tau_2_'+str(tauN)+sys+tmpSel ).Fill( ireco['sd2'+str(tauN)], WEIGHT )

                    ##### differential pt
                    #tmpRecoPt = getattr(ireco['jet'], 'pt'+( '_nom' if sys.startswith(self.sysWeightList) else sys ))
                    #for i in range( 0, len(self.diffPt)-1):
                    #    if ( tmpRecoPt > self.diffPt[i] ) & ( tmpRecoPt < self.diffPt[i+1] ):
                    #        getattr( self, 'reco'+iRJ+'_pt'+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( tmpRecoPt, tmpWeight )
                    #        getattr( self, 'reco'+iRJ+'_tau21'+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['tau21'], tmpWeight )
                    #        getattr( self, 'reco'+iRJ+'_tau32'+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['tau32'], tmpWeight )
                    #        if self.runSDVariables:
                    #            getattr( self, 'recosd'+iRJ+'_tau21'+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['sdtau21'], tmpWeight )
                    #            getattr( self, 'recosd'+iRJ+'_tau32'+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['sdtau32'], tmpWeight )
                    #        for tauN in range(1, self.maxTau+1):
                    #            getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['0p5'+str(tauN)], tmpWeight )
                    #            getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['1'+str(tauN)], tmpWeight )
                    #            getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sys+'_'+str(self.diffPt[i])+tmpSel ).Fill( ireco['2'+str(tauN)], tmpWeight )

                self.fillBranches( 'selRecoJets'+sys, recoJet )

                if self.isMC:
                    if not passGenSel:  ##### fakes
                        self.fakes = self.fakes+1
                        self.eventCategory = 2

                        #if not self.onlyTrees:
                        for iRJ,ireco in recoJet.items():

                            if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < self.minLeadAK8JetPtDijet ): continue

                            if sys.startswith(self.sysWeightList):
                                getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                            else:
                                getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                            if self.runSDVariables:
                                getattr( self, 'fakerecosd'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['sdtau21'], WEIGHT )
                                getattr( self, 'fakerecosd'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['sdtau32'], WEIGHT )
                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                if self.runSDVariables:
                                    getattr( self, 'fakerecosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakerecosd'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakerecosd'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd2'+str(tauN)], WEIGHT )


                    else:        ##### go to matrix
                        
                        #tmpGenJet1 = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
                        #tmpGenJet2 = self.createNsubBasis( selGenJets[1], event, 'GenCands' )
                        
                        genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )

                        #if abs(tmpGenJet1['jet'].rapidity) > abs(tmpGenJet2['jet'].rapidity):
                        #    genJet['Jet1'] = tmpGenJet1
                        #    genJet['Jet2'] = tmpGenJet2
                        #else:
                        #    genJet['Jet1'] = tmpGenJet2
                        #    genJet['Jet2'] = tmpGenJet1

                        for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                            #if sys.startswith('_nom'):
                            getattr( self, 'gen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                            if self.runSDVariables:
                                getattr( self, 'gensd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], self.totalGenWeight )
                                getattr( self, 'gensd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], self.totalGenWeight )

                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                if self.runSDVariables:
                                    getattr( self, 'gensd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['sd0p5'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'gensd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['sd1'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'gensd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['sd2'+str(tauN)], self.totalGenWeight )

                        #    tmpGenPt = getattr(igen['jet'], 'pt')
                        #    for i in range( 0, len(self.diffPt)-1):
                        #        if ( tmpGenPt > self.diffPt[i] ) & ( tmpGenPt < self.diffPt[i+1] ):
                        #            getattr( self, 'gen'+iGJ+'_pt'+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                        #            getattr( self, 'gen'+iGJ+'_tau21'+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                        #            getattr( self, 'gen'+iGJ+'_tau32'+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                        #            for tauN in range(1, self.maxTau+1):
                        #                getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                        #                getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                        #                getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                        #print (genJet)
                        self.fillBranches( 'selGenJets'+sys, genJet )


                        if ( genJet['Jet']['jet'].p4().DeltaR( recoJet['Jet']['jet'].p4() ) < 0.8 ):# and ( genJet['Jet2']['jet'].p4().DeltaR( recoJet['Jet2']['jet'].p4() ) > 0.8 ):
                            #iGenSel=''
                            #passGenSel=False
                            deltaRmatch = True
                            self.response= self.response+1
                            self.eventCategory = 4
                        
                        #if deltaR matched, fill hists for accepted gen jet and truereco jets, as well as the response matrices, of course!
                        if deltaRmatch:#( iGenSel==iRecoSel[sys] ) and not self.onlyTrees:

                            for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                                #if sys.startswith('_nom'):

                                if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < self.minLeadAK8JetPtDijet or getattr(igen['jet'],'pt')<self.minLeadAK8JetPtDijet): 
                                    continue

                                getattr( self, 'accepgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                                if self.runSDVariables:
                                    getattr( self, 'accepgensd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], self.totalGenWeight )
                                    getattr( self, 'accepgensd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], self.totalGenWeight )

                                
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'accepgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                    if self.runSDVariables:
                                        getattr( self, 'accepgensd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['sd0p5'+str(tauN)], self.totalGenWeight )
                                        getattr( self, 'accepgensd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['sd1'+str(tauN)], self.totalGenWeight )
                                        getattr( self, 'accepgensd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['sd2'+str(tauN)], self.totalGenWeight )

                                    
                                
                                if sys.startswith(self.sysWeightList):
                                    getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                else:
                                    getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                
                                getattr( self, 'truereco'+iRJ+'_eta'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_phi'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_tau21'+sys+iGenSel ).Fill( ireco['tau21'], WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_tau32'+sys+iGenSel ).Fill( ireco['tau32'], WEIGHT )
                                if self.runSDVariables:
                                    getattr( self, 'truerecosd'+iRJ+'_tau21'+sys+iGenSel ).Fill( ireco['sdtau21'], WEIGHT )
                                    getattr( self, 'truerecosd'+iRJ+'_tau32'+sys+iGenSel ).Fill( ireco['sdtau32'], WEIGHT )
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'truereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                    if self.runSDVariables:
                                        getattr( self, 'truerecosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'truerecosd'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                        getattr( self, 'truerecosd'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd2'+str(tauN)], WEIGHT )

                                #### Filling response matrices and resolution plots
                            
                                if sys.endswith('nom'):
                                    getattr( self, 'resol'+iRJ+'_pt'+iGenSel ).Fill( getattr(ireco['jet'], 'pt')/getattr(igen['jet'], 'pt'))#, WEIGHT )
                                    getattr( self, 'resol'+iRJ+'_sdmass'+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop')/getattr(igen['jet'], 'mass'))#, WEIGHT )
                                    try:
                                        getattr( self, 'resol'+iRJ+'_tau21'+iGenSel ).Fill( ireco['tau21']/igen['tau21'])#, WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau32'+iGenSel ).Fill( ireco['tau32']/igen['tau32'])#, WEIGHT )
                                    except ZeroDivisionError:
                                        getattr( self, 'resol'+iRJ+'_tau21'+iGenSel ).Fill( -1.)#, WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau32'+iGenSel ).Fill( -1.)#, WEIGHT )

                                    #getattr( self, 'resolsd'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['sdjet'], 'pt')/getattr(igen['sdjet'], 'pt'), WEIGHT )
                                    #getattr( self, 'resolsd'+iRJ+'_sdmass'+sys+iGenSel ).Fill( getattr(ireco['sdjet'], 'mass')/getattr(igen['sdjet'], 'mass'), WEIGHT )
                                    if self.runSDVariables:
                                        try:
                                            getattr( self, 'resolsd'+iRJ+'_tau21'+iGenSel ).Fill( ireco['sdtau21']/igen['sdtau21'])#, WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau32'+iGenSel ).Fill( ireco['sdtau32']/igen['sdtau32'])#, WEIGHT )
                                        except ZeroDivisionError:
                                            getattr( self, 'resolsd'+iRJ+'_tau21'+iGenSel ).Fill( -1)#, WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau32'+iGenSel ).Fill( -1)#, WEIGHT )

                                getattr( self, 'resp'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                                getattr( self, 'resp'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                                getattr( self, 'resp'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight - WEIGHT )
                                getattr( self, 'resp'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight - WEIGHT )

                                if self.runSDVariables:
                                    getattr( self, 'respsd'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], ireco['sdtau21'], WEIGHT )
                                    getattr( self, 'respsd'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], ireco['sdtau32'], WEIGHT )
                                    getattr( self, 'respsd'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], -1., self.totalGenWeight - WEIGHT )
                                    getattr( self, 'respsd'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], -1., self.totalGenWeight - WEIGHT )

                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight - WEIGHT )


                                    if self.runSDVariables:
                                        getattr( self, 'respsd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['sd0p5'+str(tauN)], ireco['sd0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'respsd'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['sd1'+str(tauN)], ireco['sd1'+str(tauN)], WEIGHT )
                                        getattr( self, 'respsd'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['sd2'+str(tauN)], ireco['sd2'+str(tauN)], WEIGHT )
                                        getattr( self, 'respsd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['sd0p5'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                        getattr( self, 'respsd'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['sd1'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                        getattr( self, 'respsd'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['sd2'+str(tauN)], -1., self.totalGenWeight - WEIGHT )

                                    if sys.endswith('nom'):
                                        getattr( self, 'resol'+iRJ+'_tau_0p5_'+str(tauN)+iGenSel).Fill( ( ireco['0p5'+str(tauN)]/igen['0p5'+str(tauN)] if igen['0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau_1_'+str(tauN)+iGenSel).Fill( ( ireco['1'+str(tauN)]/igen['1'+str(tauN)] if igen['1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau_2_'+str(tauN)+iGenSel).Fill( ( ireco['2'+str(tauN)]/igen['2'+str(tauN)] if igen['2'+str(tauN)]!=0 else -999 ), WEIGHT )
                                        if self.runSDVariables:
                                            getattr( self, 'resolsd'+iRJ+'_tau_0p5_'+str(tauN)+iGenSel).Fill( ( ireco['sd0p5'+str(tauN)]/igen['sd0p5'+str(tauN)] if igen['sd0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau_1_'+str(tauN)+iGenSel).Fill( ( ireco['sd1'+str(tauN)]/igen['sd1'+str(tauN)] if igen['sd1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau_2_'+str(tauN)+iGenSel).Fill( ( ireco['sd2'+str(tauN)]/igen['sd2'+str(tauN)] if igen['sd2'+str(tauN)]!=0 else -999 ), WEIGHT )

                                ##### differential pt
                                tmpRecoPt = getattr(ireco['jet'], 'pt'+( '_nom' if sys.startswith(self.sysWeightList) else sys ))
                                tmpGenPt = getattr(igen['jet'], 'pt')
                                #for i in range( 0, len(self.diffPt)-1):
                                #    tmpReco = ( tmpRecoPt > self.diffPt[i] ) & ( tmpRecoPt < self.diffPt[i+1] )
                                #    tmpGen = ( tmpGenPt > self.diffPt[i] ) & ( tmpGenPt < self.diffPt[i+1] )
                                #    if tmpReco & tmpGen:
                                #        getattr( self, 'resp'+iRJ+'_tau21'+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                                #        getattr( self, 'resp'+iRJ+'_tau32'+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                                #        getattr( self, 'resp'+iRJ+'_tau21'+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight - WEIGHT )
                                #        getattr( self, 'resp'+iRJ+'_tau32'+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight - WEIGHT )

                                #        for tauN in range(1, self.maxTau+1):
                                #        getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                                #            getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                                #            getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )

                                #            getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                #            getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                                #            getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sys+'_'+str(self.diffPt[i])+iGenSel ).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight - WEIGHT )
                            self.fillBranches( 'accepGenJets'+sys, genJet )    
                            self.fillBranches( 'trueRecoJets'+sys, recoJet )
                                                
                        #else: self.ufoResponse = self.ufoResponse+1   
                        
                        #if not deltaR matched, fill hists for misreconstructed gen jet and fakereco jets
                        else:

                            self.miss = self.miss+1
                            self.eventCategory = 3

                            for iGJ,igen in genJet.items():

                                if getattr(igen['jet'], 'pt') < self.minLeadAK8JetPtDijet: continue

                                getattr( self, 'missgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
            
                                getattr( self, 'missgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                                
                                getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                                getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )

                                if self.runSDVariables:
                                    getattr( self, 'missgensd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], self.totalGenWeight )
                                    getattr( self, 'missgensd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], self.totalGenWeight )
                                    getattr( self, 'respsd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                                    getattr( self, 'respsd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )

                                
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                    
                                    getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight )
                                    getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight )
                                    getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight )
                                
                                    if self.runSDVariables:
                                        getattr( self, 'missgensd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['sd0p5'+str(tauN)], self.totalGenWeight )
                                        getattr( self, 'missgensd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['sd1'+str(tauN)], self.totalGenWeight )
                                        getattr( self, 'missgensd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['sd2'+str(tauN)], self.totalGenWeight )

                                        getattr( self, 'respsd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight )
                                        getattr( self, 'respsd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight )
                                        getattr( self, 'respsd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight )

                            self.fakes = self.fakes+1
                            self.eventCategory = 2

                            #if not self.onlyTrees:
                            for iRJ,ireco in recoJet.items():

                                if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < self.minLeadAK8JetPtDijet ): continue

                                if sys.startswith(self.sysWeightList):
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                else:
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                                if self.runSDVariables:
                                    getattr( self, 'fakerecosd'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['sdtau21'], WEIGHT )
                                    getattr( self, 'fakerecosd'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['sdtau32'], WEIGHT )
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                    if self.runSDVariables:
                                        getattr( self, 'fakerecosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'fakerecosd'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                        getattr( self, 'fakerecosd'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd2'+str(tauN)], WEIGHT )
                                    

                #if self.onlyTrees and not sys.startswith(self.sysWeightList): self.fillBranches( 'selRecoJets'+sys, recoJet )

        #if passGenSel and ('_nom' in passRecoSel.keys()) and not passRecoSel['_nom'] :
        if self.isMC and not passRecoSel[sys]:# and not self.onlyUnc.startswith(self.sysWeightList):
            if passGenSel:# and not self.onlyTrees:

                #### Misses
                self.miss = self.miss+1
                self.eventCategory = 3

                #tmpGenJet1 = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
                #tmpGenJet2 = self.createNsubBasis( selGenJets[1], event, 'GenCands' )
                genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )

                #if abs(tmpGenJet1['jet'].eta) > abs(tmpGenJet2['jet'].eta):
                #    genJet['Jet1'] = tmpGenJet1
                #    genJet['Jet2'] = tmpGenJet2
                #else:
                #    genJet['Jet1'] = tmpGenJet2
                #    genJet['Jet2'] = tmpGenJet1

                for iGJ,igen in genJet.items():
                    if getattr(igen['jet'], 'pt') < self.minLeadAK8JetPtDijet: continue

                    getattr( self, 'gen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight)
                    getattr( self, 'gen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                    getattr( self, 'gen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                    getattr( self, 'gen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )

                    getattr( self, 'gen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                    getattr( self, 'gen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )

                    
                    if self.runSDVariables:
                        getattr( self, 'gensd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], event.genWeight )
                        getattr( self, 'gensd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], event.genWeight )

                    getattr( self, 'missgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                    getattr( self, 'missgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                    getattr( self, 'missgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                    getattr( self, 'missgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )

                    getattr( self, 'missgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                    getattr( self, 'missgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                    getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                    getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )

                    if self.runSDVariables:
                        getattr( self, 'missgensd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], self.totalGenWeight )
                        getattr( self, 'missgensd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], self.totalGenWeight )
                        getattr( self, 'respsd'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                        getattr( self, 'respsd'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )


                    for tauN in range(1, self.maxTau+1):
                        
                        getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], self.totalGenWeight)
                        getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                        if self.runSDVariables:
                            getattr( self, 'gensd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['sd0p5'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'gensd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['sd1'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'gensd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['sd2'+str(tauN)], self.totalGenWeight )

                        getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                        
                        getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight )
                        getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight )
                        getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight )
                    
                        if self.runSDVariables:
                            getattr( self, 'missgensd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['sd0p5'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'missgensd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['sd1'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'missgensd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['sd2'+str(tauN)], self.totalGenWeight )

                            getattr( self, 'respsd'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight )
                            getattr( self, 'respsd'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight )
                            getattr( self, 'respsd'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight )
                self.fillBranches( 'selGenJets'+sys, genJet )

            #else:
            #    self.ufo = self.ufo+1
            #    self.eventCategory = 10

            #if self.onlyTrees: self.fillBranches( 'selGenJets', genJet )

        #if self.onlyTrees:
        self.out.fillBranch( 'eventCategory', self.eventCategory )
        self.out.fillBranch( 'totalRecoWeight', self.totalRecoWeight if self.isMC else 1)
        self.out.fillBranch( 'totalGenWeight', self.totalGenWeight if self.isMC else -1)
        if not ('puWeight' in self.onlyUnc):
           self.out.fillBranch( 'puWeight', self.puWeight if self.isMC else 1)

        if self.onlyUnc:
            if 'pdf' in self.onlyUnc:
                self.out.fillBranch( 'pdfWeightUp', self.pdfWeightUp if self.isMC else 1)
                self.out.fillBranch( 'pdfWeightDown', self.pdfWeightDown if self.isMC else 1)
            if 'isr' in self.onlyUnc:
                self.out.fillBranch( 'isrWeightUp', self.isrWeightUp if self.isMC else 1)
                self.out.fillBranch( 'isrWeightDown', self.isrWeightDown if self.isMC else 1)
            if 'fsr' in self.onlyUnc:
                self.out.fillBranch( 'fsrWeightUp', self.fsrWeightUp if self.isMC else 1)
                self.out.fillBranch( 'fsrWeightDown', self.fsrWeightDown if self.isMC else 1)
            if 'puWeight' in self.onlyUnc:
                self.out.fillBranch( 'puWeightUp', event.puWeightUp if self.isMC else 1 )
                self.out.fillBranch( 'puWeightDown', event.puWeightDown if self.isMC else 1 )


        return True


    #############################################################################
    def recoSelection( self, event, sysUnc=[] ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))

        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minLooseMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.pfIsoId>=2 and x.mediumId and abs(x.dxy)<0.2 and abs(x.dz)<0.5 and x.miniPFRelIso_all<0.20] # applying tight selection on muons already here since we only veto for loose muons and loose electrons in event

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons) 
        ##################################################

        recoAK8jets = {}
        passSel = {}
        ptAsym = {}
        deltaPhi = {}
        iSel = {}
        for sys in self.sysSource:
            if sys.startswith(self.sysWeightList): sys = '_nom'
            #### Basic AK8 jet selection
            for ijets in AK8jets: ijets.rapidity = self.etaToRapidity(ijets)
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs( x.rapidity ) < self.maxLeadAK8JetRap and (x.jetId >= 2)] #self.maxJetAK8Eta
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ),reverse=True)
            ##################################################

            ##### Applying selection
            passSel[sys], ptAsym[sys], deltaPhi[sys] = self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets[sys], sys )
            iSel[sys] = '_dijetSel' if passSel[sys] else None
        '''
        ##### Trigger weights
        if not self.isMC:
            if passSel['_nom']:
                self.triggerWeight = 0
                triggerVersion = -9999
                for i in range( len(self.runTables[ self.year ]['low'] ) ):
                    if ( self.runTables[self.year]['low'][i] <= event.run ) and ( event.run <= self.runTables[self.year]['high'][i] ):
                        triggerVersion = 2+i

                            for itrigger, itrValues in self.triggerTable.items():
                                if ( getattr(event, 'HLT_'+itrigger)==1 ):

                                    getattr( self, 'nPVs_only'+itrigger+'_dijetSel' ).Fill( getattr( event, 'PV_npvsGood') )
                                    if len(recoAK8jets['_nom'])>0 and not self.onlyTrees:
                                        getattr( self, 'LeadJetAK8_pt_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].pt_nom )   ### pt_nom here to ensure data process
                                        getattr( self, 'LeadJetAK8_eta_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].eta )
                                        getattr( self, 'LeadJetAK8_phi_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].phi )
                                        getattr( self, 'LeadJetAK8_mass_only'+itrigger+'_dijetSel' ).Fill( recoAK8jets['_nom'][0].msoftdrop_nom )
                                    if ( itrValues[self.year][0] < recoAK8jets['_nom'][0].pt_nom < itrValues[self.year][1] ) :
                                        iSel['_nom'] = '_'+itrigger+'_dijetSel'
                                        passSel['_nom'] = True
                                        self.triggerWeight = 1 #itrValues[self.year][triggerVersion]
                                        break
                                    else:
                                        passSel['_nom'] = False
                                        iSel['_nom'] = None
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



            #                    else: passSel['_nom'] = False
            else:
                self.triggerWeight = 0
                passSel['_nom'] = False
                iSel['_nom'] = None
        #self.out.fillBranch( 'triggerWeight', self.triggerWeight )
        ##################################################
        '''

        #### Weight #########
        weight=1.

        if self.isMC:
            weight = event.puWeight * event.genWeight
            #if not self.onlyTrees and not self.onlyUnc: getattr( self, 'PUweight' ).Fill( event.puWeight )
            self.puWeight = event.puWeight

        else:
            weight = 1 #self.triggerWeight
        self.totalRecoWeight = weight

        ##################################################

        if not self.onlyUnc:
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
            getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recoAK8jets['_nom']) )
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt_nom )
                getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop_nom )

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
            getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets['_nom']), weight )
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt_nom, weight )
                getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop_nom, weight )
            getattr( self, 'recoPtAsym_noSel' ).Fill( ptAsym['_nom'], weight )
            getattr( self, 'recoDeltaPhi_noSel' ).Fill( deltaPhi['_nom'], weight )

            reweight = self.totalRecoWeight
            #### Creating Nsub basis, filling histos and creating branches IF passSel
            if ('_nom' in iSel.keys()) and passSel['_nom'] and iSel['_nom']:

                #### Basic reco histos
                getattr( self, 'nPVs'+iSel['_nom'] ).Fill( getattr( event, 'PV_npvsGood'), reweight )
                getattr( self, 'nleps'+iSel['_nom'] ).Fill( len(recoMuons)+len(recoElectrons), reweight )
                for imuon in recoMuons:
                    getattr( self, 'muons_pt'+iSel['_nom'] ).Fill( imuon.pt, reweight )
                    getattr( self, 'muons_eta'+iSel['_nom'] ).Fill( imuon.eta, reweight )
                    getattr( self, 'muons_phi'+iSel['_nom'] ).Fill( imuon.phi, reweight )
                for iele in recoElectrons:
                    getattr( self, 'eles_pt'+iSel['_nom'] ).Fill( iele.pt, reweight )
                    getattr( self, 'eles_eta'+iSel['_nom'] ).Fill( iele.eta, reweight )
                    getattr( self, 'eles_phi'+iSel['_nom'] ).Fill( iele.phi, reweight )
                getattr( self, 'nAK8jets'+iSel['_nom'] ).Fill( len(recoAK8jets['_nom']), reweight )
                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt'+iSel['_nom'] ).Fill( ijet.pt_nom, reweight )
                    getattr( self, 'AK8jets_eta'+iSel['_nom'] ).Fill( ijet.eta, reweight )
                    getattr( self, 'AK8jets_phi'+iSel['_nom'] ).Fill( ijet.phi, reweight )
                    getattr( self, 'AK8jets_mass'+iSel['_nom'] ).Fill( ijet.msoftdrop_nom, reweight )
                getattr( self, 'recoPtAsym'+iSel['_nom'] ).Fill( ptAsym['_nom'], reweight )
                getattr( self, 'recoDeltaPhi'+iSel['_nom'] ).Fill( deltaPhi['_nom'], reweight )
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

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minLooseMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        ##################################################

        #### Basic AK8 jet selection
        for ijets in genJetsAK8: ijets.rapidity = self.etaToRapidity(ijets)
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.rapidity) < self.maxLeadAK8JetRap]#maxJetAK8Eta 
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
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
            if len(genAK8jets)>0: getattr(self, 'genJetAK8_partonFlavour_noSel').Fill(genAK8jets[0].partonFlavour)


            ##### Filling histograms
            if passSel and iSel:
                #print (genAK8jets[0].partonFlavour)
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
                getattr(self, 'genJetAK8_partonFlavour'+iSel).Fill(genAK8jets[0].partonFlavour)

        return passSel, iSel, genMuons, genElectrons, genAK8jets

    #############################################################################
    def dijetSelection( self, event, muons, electrons, AK8jets, ptLabel ):

        if  (len(muons)+len(electrons)==0) and (len(AK8jets)>1):# and abs(AK8jets[0].eta)<self.maxLeadAK8JetEta):#

            jet1Pt = getattr( AK8jets[0], 'pt'+ptLabel )
            jet2Pt = getattr( AK8jets[1], 'pt'+ptLabel )
            ptAsym = ( jet1Pt - jet2Pt ) / (jet1Pt + jet2Pt)

            tmpJet1 = ROOT.TLorentzVector()
            tmpJet1.SetPtEtaPhiM( getattr( AK8jets[0], 'pt'+ptLabel ), AK8jets[0].eta, AK8jets[0].phi, getattr( AK8jets[0], 'mass'+ptLabel )  )
            tmpJet2 = ROOT.TLorentzVector()
            tmpJet2.SetPtEtaPhiM( getattr( AK8jets[1], 'pt'+ptLabel ), AK8jets[1].eta, AK8jets[1].phi, getattr( AK8jets[1], 'mass'+ptLabel )  )
            deltaPhi = tmpJet1.DeltaPhi( tmpJet2 )
            
            deltaR_flag = True
            for i in range(1,len(AK8jets)-1):
                tmpJet = ROOT.TLorentzVector()
                tmpJet.SetPtEtaPhiM( getattr( AK8jets[i], 'pt'+ptLabel ), AK8jets[i].eta, AK8jets[i].phi, getattr( AK8jets[i], 'mass'+ptLabel )  )
                if tmpJet1.DeltaR(tmpJet)<1.6:# or abs(tmpJet1.DeltaPhi(tmpJet)<2.): #2*jet radius separation + deltaphi sep. between leading AK8 jet and any other fatjets in the event
                    deltaR_flag=False
                    break

            
            

            if (jet1Pt>self.minLeadAK8JetPtDijet) and (abs(deltaPhi)>2.) and deltaR_flag:# and (ptAsym<0.3)  and
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
        nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, constituents )    ### needed for genjet tau21 tau32


        ### default in CMS WTA_kT https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        try: ak8jet['tau21'] = nsub1_WTA_kT[1]/nsub1_WTA_kT[0]
        except ZeroDivisionError: ak8jet['tau21'] = -1
        try: ak8jet['tau32'] = nsub1_WTA_kT[2]/nsub1_WTA_kT[1]
        except ZeroDivisionError: ak8jet['tau32'] = -1


        #### filling histos and branches with nsub basis
        for tauN in range(self.maxTau):
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]


        #### Computing Softdrop jets
        if self.runSDVariables:
            sdAK8jets = self.sd.result( constituents ) #CandsPUPPIweightedVec )

            ### Storing good jet as list for later use
            if len(sdAK8jets)>0:

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
                sd_nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, sd_constituents )

                try: ak8jet['sdtau21'] = sd_nsub1_WTA_kT[1]/sd_nsub1_WTA_kT[0]
                except ZeroDivisionError: ak8jet['sdtau21'] = -1
                try: ak8jet['sdtau32'] = sd_nsub1_WTA_kT[2]/sd_nsub1_WTA_kT[1]
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

        return ak8jet

    #############################################################################
    def fillBranches( self, jetLabel, jetInfo ):

        #### Filling branch with passAK8jet info after selection
        #print (jetLabel+"_pt", [ getattr( iJ['jet'], 'pt'+jetLabel.split('Jets')[1] ) for i,iJ in jetInfo.items()  ])
        self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
        self.out.fillBranch(jetLabel+"_mass", [ iJ['jet'].mass for i,iJ in jetInfo.items() ] )
        self.out.fillBranch(jetLabel+"_pt", [ iJ['jet'].pt for i,iJ in jetInfo.items()  ] )#, 'pt'+jetLabel.split('Jets')[1] ) for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_eta", [ iJ['jet'].eta for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_phi", [ iJ['jet'].phi for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau21", [ iJ['tau21'] for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau32", [ iJ['tau32'] for i,iJ in jetInfo.items()  ] )
        if self.runSDVariables:
            self.out.fillBranch(jetLabel+"_SD_Tau21", [ iJ['sdtau21'] for i,iJ in jetInfo.items()  ] )
            self.out.fillBranch(jetLabel+"_SD_Tau32", [ iJ['sdtau32'] for i,iJ in jetInfo.items()  ] )

        for tauN in range(1, self.maxTau+1):
            self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN),  [ iJ['0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN),  [ iJ['1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN),  [ iJ['2'+str(tauN)] for i,iJ in jetInfo.items() ] )
            if self.runSDVariables:
                self.out.fillBranch(jetLabel+"SD_tau_0p5_"+str(tauN),  [ iJ['sd0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
                self.out.fillBranch(jetLabel+"SD_tau_1_"+str(tauN),  [ iJ['sd1'+str(tauN)] for i,iJ in jetInfo.items() ] )
                self.out.fillBranch(jetLabel+"SD_tau_2_"+str(tauN),  [ iJ['sd2'+str(tauN)] for i,iJ in jetInfo.items() ] )

    #############################################################################
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)

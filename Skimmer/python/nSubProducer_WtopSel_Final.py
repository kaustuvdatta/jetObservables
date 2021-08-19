import ROOT
import math, os, sys
import numpy as np
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nSubProd(Module):

    def __init__(self, sysSource=[], leptonSF={}, year='2017', isMC=True, topreweight=False, onlyUnc='', onlyTrees=False):
        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.onlyTrees = onlyTrees
        
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
        self.maxJetAK8Eta = 2.4

        ### Cuts for selections
        self.minLeadAK8JetPtW = 200.
        self.minSDMassW = 65.#60.#   ### looser to pick low bins
        self.maxSDMassW = 115. #120.#  ### looser to pick higher bins
        self.minLeadAK8JetPtTop= 350.
        self.minSDMassTop = 140.
        self.maxSDMassTop = 220.
        self.METCutWtop = 50.
        self.minLeptonicWPt = 150.

        ### Kinematics Cuts Jets ###
        self.minJetPt = 30.
        self.maxJetEta = 2.4
        self.minBDisc = 0.3040  ### L: 0.0532, M: 0.3040, T: 0.7476, for DeepJet (ie, DeepFlavB)

        ### Kinenatic Cuts Muons ###
        self.minLooseMuonPt = 40.
        self.minTightMuonPt = 55.
        self.maxMuonEta = 2.4
        self.minMuonMETPt = 50.

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 40.
        self.minTightElectronPt = 120.
        self.minElectronMETPt = 80.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.4]

        self.totalWeight = 1.
        self.btaggingWeight = 1.
        self.topreweight = topreweight 
        self.pdfWeightUp = 0
        self.isrWeightUp = 0
        self.fsrWeightUp = 0
        self.pdfWeightDown = 0
        self.isrWeightDown = 0
        self.fsrWeightDown = 0
        self.eventCategory = -1
        self.leptonWeight=1.
        self.dummy = 0

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
        self.nJet = [ 'Jet']#, 'sdJet' ]

        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        if onlyUnc: self.sysSource = [ onlyUnc+i for i in [ 'Up', 'Down' ] ]
        print ("Sys sources:", self.sysSource)
        self.sysWeightList = ( '_pu',  '_ps', '_isr', '_fsr' )#'_pdf', 
        print ("Sys weight listL", self.sysWeightList)
    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        tauBins = 100
        if not self.onlyTrees:
            ### Booking histograms
            selList = [ '_WSel', '_topSel' ]
            if not self.onlyUnc:
                self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   25, 0, 25) )
                self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
                self.addObject( ROOT.TH1F('Lepweight',   ';LepWeight',   20, 0, 2) )
                self.addObject( ROOT.TH1F('Btagweight',   ';BtagWeight',   25, 0, 2) )
                self.addObject( ROOT.TH1F('Topweight',   ';BtagWeight',   25, 0, 2) )

                #### general selection
                for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
                    self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
                    self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
                    self.addP4Hists( 'muons', isel )
                    self.addP4Hists( 'eles', isel )
                    self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
                    self.addP4Hists( 'AK8jets', isel )
                    self.addObject( ROOT.TH1F('nAK4jets'+isel,   ';number of AK4 jets',   20, 0, 20) )
                    self.addP4Hists( 'AK4jets', isel )
                    self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )
                    self.addObject( ROOT.TH1F('HT'+isel,   ';HT (GeV)',   200, 0, 2000) )
                    #self.addObject( ROOT.TH1F('recoPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                    #self.addObject( ROOT.TH1F('recoDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )

                    self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

                if self.isMC:
                    for isel in [ '_noSel' ] + selList :
                        self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   20, 0, 20) )
                        self.addP4Hists( 'genmuons', isel )
                        self.addP4Hists( 'geneles', isel )
                        self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   20, 0, 20) )
                        self.addP4Hists( 'AK8genjets', isel )
                        self.addObject( ROOT.TH1F('ngenAK4jets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                        self.addP4Hists( 'AK4genjets', isel )
                        self.addObject( ROOT.TH1F('genMETPt'+isel,   ';gen MET (GeV)',   200, 0, 2000) )
                        self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )
                        #self.addObject( ROOT.TH1F('genPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
                        #self.addObject( ROOT.TH1F('genDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )


            for isel in selList:
                for iJ in self.nJet:
                    #self.addP4Hists( 'uforeco'+iJ, isel )
                    #self.addP4Hists( 'ufogen'+iJ, isel )
                    
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
                            #print ('Creating gen objects:', itype, iJ, isel)
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

        #tmp = inputFile.Get('Runs')
        #for i in tmp.LHEPdfSumw: print(i)
        if self.onlyTrees:
            self.out = wrappedOutputTree
            #self.out.branch('triggerWeight',  "F")
            self.out.branch('eventCategory',  "I")
            self.out.branch( 'totalWeight', "F" )
            self.out.branch('leptonWeight',  "F")
            self.out.branch('btagWeight',  "F")
            self.out.branch('top_pTreweight_SF',  "F")
        
            if not self.onlyUnc:
                #self.out.branch( 'pdfWeightUp', "F" )
                self.out.branch( 'isrWeightUp', "F" )
                self.out.branch( 'fsrWeightUp', "F" )
                self.out.branch( 'puWeightUp', "F" )
                #self.out.branch( 'pdfWeightDown', "F" )
                self.out.branch( 'isrWeightDown', "F" )
                self.out.branch( 'fsrWeightDown', "F" )
                self.out.branch( 'puWeightDown', "F" )
        
            tmplist = [ 'selRecoJets'+sys for sys in self.sysSource if not sys.startswith(self.sysWeightList) ]
            print ('selrecojetlist:', tmplist)
            if self.isMC: 
                tmplist.append( 'selGenJets' )
                #print ("appended genjets to nano output file")
            for iJ in tmplist:
                self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
                self.out.branch(iJ+'_pt',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_eta',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_phi',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_mass',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_Tau21',  'F', lenVar='n'+iJ)
                self.out.branch(iJ+'_Tau32',  'F', lenVar='n'+iJ)
                #self.out.branch(iJ+'_SD_Tau21',  'F', lenVar='n'+iJ)
                #self.out.branch(iJ+'_SD_Tau32',  'F', lenVar='n'+iJ)

                for x in self.nSub_labels:
                    self.out.branch(iJ+x, 'F', lenVar='n'+iJ )
                    #self.out.branch(iJ+'SD'+x, 'F', lenVar='n'+iJ )
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        if not self.onlyTrees:
            if self.isMC and not self.onlyUnc:
                self.genLevel = self.response+self.miss
                self.genLevelW = self.responseW+self.missW
                self.genLeveltop = self.responsetop+self.misstop

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

                getattr( self, 'cutflow_test' ).SetBinContent( 7, self.recoLevelW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 7, "WReco" )
                getattr( self, 'cutflow_test' ).SetBinContent( 8, self.responseW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 8, "WtrueRecoAccepGen" )
                getattr( self, 'cutflow_test' ).SetBinContent( 9, self.fakesW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 9, "WFakes" )
                getattr( self, 'cutflow_test' ).SetBinContent( 10, self.missW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 10, "WMiss" )
                getattr( self, 'cutflow_test' ).SetBinContent( 11, self.genLevelW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 11, "WGen" )

                getattr( self, 'cutflow_test' ).SetBinContent( 13, self.recoLeveltop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 13, "topReco" )
                getattr( self, 'cutflow_test' ).SetBinContent( 14, self.responsetop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 14, "toptrueRecoAccepGen" )
                getattr( self, 'cutflow_test' ).SetBinContent( 15, self.fakestop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 15, "topFakes" )
                getattr( self, 'cutflow_test' ).SetBinContent( 16, self.misstop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 16, "topMiss" )
                getattr( self, 'cutflow_test' ).SetBinContent( 17, self.genLeveltop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 17, "topGen" )

                getattr( self, 'cutflow_test' ).SetBinContent( 19, self.ufo )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 19, "allUfo" )
                getattr( self, 'cutflow_test' ).SetBinContent( 20, self.ufoResponse )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 20, "ufoResponse" )
                getattr( self, 'cutflow_test' ).SetBinContent( 21, self.ufoFake )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 21, "ufoFake" )
                getattr( self, 'cutflow_test' ).SetBinContent( 22, self.ufoMiss )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 22, "ufoMiss" )
                getattr( self, 'cutflow_test' ).SetBinContent( 23, self.ufoW )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 23, "ufoW" )
                getattr( self, 'cutflow_test' ).SetBinContent( 24, self.ufotop )
                getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 24, "ufotop" )

    #############################################################################
    def leptonSF(self, lepton, leptonP4 ):

        leptonP4eta = abs(leptonP4.eta)
        leptonP = ROOT.TMath.Sqrt(leptonP4.p4().Px()**2 + leptonP4.p4().Py()**2 + leptonP4.p4().Pz()**2)

        SFFileTrigger = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['Trigger'][0] )
        histoSFTrigger = SFFileTrigger.Get( self.leptonSFhelper[lepton]['Trigger'][1] )
        SFTrigger = histoSFTrigger.GetBinContent( histoSFTrigger.GetXaxis().FindBin( leptonP4eta ), histoSFTrigger.GetYaxis().FindBin( leptonP4.pt ) )

        SFFileID = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ID'][0] )
        histoSFID = SFFileID.Get( self.leptonSFhelper[lepton]['ID'][1] )
        histoSFID_X = histoSFID.GetXaxis().FindBin( leptonP4eta)
        histoSFID_Y = histoSFID.GetYaxis().FindBin( leptonP4.pt )
        SFID = histoSFID.GetBinContent( histoSFID_X, histoSFID_Y )
        SFID = SFID if SFID>0 else 1

        SFFileISO = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ISO'][0] )
        histoSFISO = SFFileISO.Get( self.leptonSFhelper[lepton]['ISO'][1] )
        histoSFISO_X = histoSFISO.GetXaxis().FindBin( leptonP4eta )
        histoSFISO_Y = histoSFISO.GetYaxis().FindBin( leptonP4.pt )
        SFISO = histoSFISO.GetBinContent( histoSFISO_X, histoSFISO_Y )
        SFISO = SFISO if SFISO>0 else 1
        
        SFFileRecoEff = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['RecoEff'][0] )
        histoSFRecoEff = SFFileRecoEff.Get( self.leptonSFhelper[lepton]['RecoEff'][1] )
        histoSFRecoEff_X = histoSFRecoEff.GetXaxis().FindBin( leptonP4eta )
        histoSFRecoEff_Y = histoSFRecoEff.GetYaxis().FindBin( leptonP )
        SFRecoEff = histoSFRecoEff.GetBinContent( histoSFRecoEff_X, histoSFRecoEff_Y )
        SFRecoEff = SFRecoEff if SFRecoEff>0 else 1

        #print (SFTrigger * SFID * SFISO), SFTrigger , SFID , SFISO, leptonP4.pt, leptonP4.eta
        return [SFTrigger , SFID , SFISO, SFRecoEff]


    #############################################################################
    def getBTagWeight(self, nBTagged=0, jet_SFs=[0]):
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


    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''
        self.isMC = event.run == 1

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenAK4bjets, selGenJets, selGenMET = self.genSelection(event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoAK4bjets, selRecoJets, selRecoMET = self.recoSelection( event )

        if not self.isMC and not passRecoSel['_nom']: return False
            
        for sys in self.sysSource:

            if sys.startswith(self.sysWeightList):
    
                selRecoJets[sys] = selRecoJets['_nom']
                iRecoSel[sys] = iRecoSel['_nom']
                passRecoSel[sys] = passRecoSel['_nom']

            genJet = OrderedDict()
            recoJet = OrderedDict()
            tmpRecoJets = OrderedDict()

            if passRecoSel[sys]:  #### Detector level dist.
            
                self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel
                if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                    self.recoLevelW = self.recoLevelW+1
                    #minLeadAK8JetPt = self.minLeadAK8JetPtW
                elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                    self.recoLeveltop = self.recoLeveltop+1
                    #minLeadAK8JetPt = self.minLeadAK8JetPtTop
                else: print('Weird reco', iRecoSel[sys])

                tmpRecoJets[sys] = {}

                if sys.startswith('_nom'):
                    tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                else:
                    if '_nom' in tmpRecoJets.keys():
                        tmpRecoJets[sys][0] = tmpRecoJets['_nom'][0]
                    else:
                        tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )

                recoJet['Jet'] = tmpRecoJets[sys][0] #self.createNsubBasis( tmpRecoJets[sys][0], event, 'PFCands' )#tmpRecoJet[sys][0]
                self.eventCategory=1

                if sys.startswith(self.sysWeightList):
                    ##### PDF does not work
                    #if sys.endswith('pdfWeightUp'):
                    #    self.pdfWeightUp = getattr( event, 'LHEPdfWeight' )[0]
                    #    WEIGHT = event.genWeight * self.pdfWeightUp
                    #elif sys.endswith('pdfWeightDown'):
                    #    self.pdfWeightDown = (getattr( event, 'LHEPdfWeight' )[0] )
                    #    WEIGHT = event.genWeight * self.pdfWeightDown
                    if sys.endswith('isrWeightDown'):
                        self.isrWeightDown = getattr( event, 'PSWeight' )[0]
                        WEIGHT = event.genWeight * self.isrWeightDown
                    elif sys.endswith('fsrWeightDown'):
                        self.fsrWeightDown = getattr( event, 'PSWeight' )[1]
                        WEIGHT = event.genWeight * self.isrWeightDown
                    elif sys.endswith('isrWeightUp'):
                        self.isrWeightUp = getattr( event, 'PSWeight' )[2]
                        WEIGHT = event.genWeight * self.isrWeightUp
                    elif sys.endswith('fsrWeightUp'):
                        self.fsrWeightUp = getattr( event, 'PSWeight' )[3]
                        WEIGHT = event.genWeight * self.fsrWeightUp
                    else:
                        WEIGHT = event.genWeight * getattr( event, sys.split('_')[1] )
                else: WEIGHT =  self.totalWeight
                #if passRecoSel[sys]:
                #    print ("#####################################################################################")
                #    print ("FINAL WEIGHT=",WEIGHT)
                #    print ("#####################################################################################")

                if not self.onlyTrees:
                    for iRJ,ireco in recoJet.items():

                        if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                            minLeadAK8JetPt = self.minLeadAK8JetPtW
                        elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                            minLeadAK8JetPt = self.minLeadAK8JetPtTop

                        if ( getattr(ireco['jet'], 'pt'+( '_nom' if sys.startswith(self.sysWeightList) else sys ) ) > minLeadAK8JetPt ):
                            #if passRecoSel[sys]: print ("filling reco histos, and sys is:",sys)
                            if sys.startswith(self.sysWeightList):
                                #print ("Storing recojet variables for sysweightlist quantities", sys)                            
                                
                                getattr( self, 'reco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                if iRecoSel[sys].startswith('_W'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                elif iRecoSel[sys].startswith('_top'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                                
                            else:
                                #print ("Storing recojet variables for nom, jer, jestot quantities",sys)
                                getattr( self, 'reco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                if iRecoSel[sys].startswith('_W'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                elif iRecoSel[sys].startswith('_top'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop'+sys )), WEIGHT) #??/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )

                            getattr( self, 'reco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                            #getattr( self, 'recosd'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['sdtau21'], WEIGHT )
                            #getattr( self, 'recosd'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['sdtau32'], WEIGHT )
                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                #getattr( self, 'recosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                #getattr( self, 'recosd'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                #getattr( self, 'recosd'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd2'+str(tauN)], WEIGHT )
                            
                            
                if self.isMC:
                    if passGenSel:  ##### go to matrix
                        #print ("yay passed gensel")
                        self.response = self.response+1

                        if iGenSel.startswith('_W'): self.responseW = self.responseW+1
                        elif iGenSel.startswith('_top'): self.responsetop = self.responsetop+1
                        else: print('Weird response', iRecoSel[sys])
                         
                        self.eventCategory = 4
                        
                        genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
                        #print ("yayx2", iGenSel, passGenSel)                       

                        if ( genJet['Jet']['jet'].p4().DeltaR( recoJet['Jet']['jet'].p4()) > 0.8 ): # we should tighten this?
                            print ("booo, DeltaR match failed", iGenSel, passGenSel, recoJet['Jet']['jet'].p4(), genJet['Jet']['jet'].p4() )
                            iGenSel=''
                            passGenSel=False


                        if passGenSel and iGenSel==iRecoSel[sys] and not self.onlyTrees:
                        #    print ("yayx3!", iGenSel, iRecoSel[sys] , sys )

                            for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                                if sys.startswith('_nom'):
                                #    print ("yayx4!", iGenSel, iRecoSel[sys] , sys, passGenSel, iGJ, igen['jet'], iRJ, ireco )

                                    getattr( self, 'accepgen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                                    getattr( self, 'accepgen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                                    getattr( self, 'accepgen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                                    getattr( self, 'accepgen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )
                                    #getattr( self, 'accepgensd'+iGJ+'_tau21'+iGenSel ).Fill( igen['sdtau21'], event.genWeight )
                                    #getattr( self, 'accepgensd'+iGJ+'_tau32'+iGenSel ).Fill( igen['sdtau32'], event.genWeight )

                                    getattr( self, 'gen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                                    getattr( self, 'gen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                                    getattr( self, 'gen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                                    getattr( self, 'gen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )
                                    getattr( self, 'gen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                                    getattr( self, 'gen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )
                                    #getattr( self, 'gensd'+iGJ+'_tau21'+iGenSel ).Fill( igen['sdtau21'], event.genWeight )
                                    #getattr( self, 'gensd'+iGJ+'_tau32'+iGenSel ).Fill( igen['sdtau32'], event.genWeight )

                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'accepgen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                                        getattr( self, 'accepgen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                                        getattr( self, 'accepgen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )
                                        #getattr( self, 'accepgensd'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['sd0p5'+str(tauN)], event.genWeight )
                                        #getattr( self, 'accepgensd'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['sd1'+str(tauN)], event.genWeight )
                                        #getattr( self, 'accepgensd'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['sd2'+str(tauN)], event.genWeight )

                                        getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                                        getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                                        getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )
                                        #getattr( self, 'gensd'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['sd0p5'+str(tauN)], event.genWeight )
                                        #getattr( self, 'gensd'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['sd1'+str(tauN)], event.genWeight )
                                        #getattr( self, 'gensd'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['sd2'+str(tauN)], event.genWeight )

                                if iGenSel.startswith('_W'):           #### counting recoLevelW
                                    minLeadAK8JetPt = self.minLeadAK8JetPtW
                                elif iGenSel.startswith('_top'):       #### counting recoLeveltop
                                    minLeadAK8JetPt = self.minLeadAK8JetPtTop

                                if  getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) > minLeadAK8JetPt: 

                                    
                                    if sys.startswith(self.sysWeightList):
                                        getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                        if iGenSel.startswith('_W'): getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                        elif iGenSel.startswith('_top'): getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                                        
                                    else:
                                        getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                    
                                    getattr( self, 'truereco'+iRJ+'_eta'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_phi'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau21'+sys+iGenSel ).Fill( ireco['tau21'], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau32'+sys+iGenSel ).Fill( ireco['tau32'], WEIGHT )
                                    #getattr( self, 'truerecosd'+iRJ+'_tau21'+sys+iGenSel ).Fill( ireco['sdtau21'], WEIGHT )
                                    #getattr( self, 'truerecosd'+iRJ+'_tau32'+sys+iGenSel ).Fill( ireco['sdtau32'], WEIGHT )
                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'truereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                        getattr( self, 'truereco'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                        #getattr( self, 'truerecosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                        #getattr( self, 'truerecosd'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                        #getattr( self, 'truerecosd'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( ireco['sd2'+str(tauN)], WEIGHT )

                                #### FIlling response matrices
                                if self.isMC:
                                    if sys.endswith('nom'):
                                        getattr( self, 'resol'+iRJ+'_pt'+iGenSel ).Fill( getattr(ireco['jet'], 'pt')/getattr(igen['jet'], 'pt'), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_sdmass'+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop')/getattr(igen['jet'], 'mass'), WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau21'+iGenSel ).Fill( ireco['tau21']/igen['tau21'], WEIGHT )
                                        getattr( self, 'resol'+iRJ+'_tau32'+iGenSel ).Fill( ireco['tau32']/igen['tau32'], WEIGHT )
                                        #getattr( self, 'resolsd'+iRJ+'_pt'+iGenSel ).Fill( getattr(ireco['sdjet'], 'pt')/getattr(igen['sdjet'], 'pt'), WEIGHT )
                                        #getattr( self, 'resolsd'+iRJ+'_sdmass'+iGenSel ).Fill( getattr(ireco['sdjet'], 'mass')/getattr(igen['sdjet'], 'mass'), WEIGHT )
                                        '''
                                        try:
                                            getattr( self, 'resolsd'+iRJ+'_tau21'+iGenSel ).Fill( ireco['sdtau21']/igen['sdtau21'], WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau32'+iGenSel ).Fill( ireco['sdtau32']/igen['sdtau32'], WEIGHT )
                                        except ZeroDivisionError:
                                            getattr( self, 'resolsd'+iRJ+'_tau21'+iGenSel ).Fill( -1, WEIGHT )
                                            getattr( self, 'resolsd'+iRJ+'_tau32'+iGenSel ).Fill( -1, WEIGHT )
                                        '''
                                    getattr( self, 'resp'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                                    getattr( self, 'resp'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                                    #getattr( self, 'respsd'+iRJ+'_tau21'+sys+iGenSel ).Fill( igen['sdtau21'], ireco['sdtau21'], WEIGHT )
                                    #getattr( self, 'respsd'+iRJ+'_tau32'+sys+iGenSel ).Fill( igen['sdtau32'], ireco['sdtau32'], WEIGHT )
                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                                        getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )
                                        #getattr( self, 'respsd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['sd0p5'+str(tauN)], ireco['sd0p5'+str(tauN)], WEIGHT )
                                        #getattr( self, 'respsd'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['sd1'+str(tauN)], ireco['sd1'+str(tauN)], WEIGHT )
                                        #getattr( self, 'respsd'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['sd2'+str(tauN)], ireco['sd2'+str(tauN)], WEIGHT )
                                        if sys.endswith('nom'):
                                            getattr( self, 'resol'+iRJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( ( ireco['0p5'+str(tauN)]/igen['0p5'+str(tauN)] if igen['0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resol'+iRJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( ( ireco['1'+str(tauN)]/igen['1'+str(tauN)] if igen['1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            getattr( self, 'resol'+iRJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( ( ireco['2'+str(tauN)]/igen['2'+str(tauN)] if igen['2'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            #getattr( self, 'resolsd'+iRJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( ( ireco['sd0p5'+str(tauN)]/igen['sd0p5'+str(tauN)] if igen['sd0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            #getattr( self, 'resolsd'+iRJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( ( ireco['sd1'+str(tauN)]/igen['sd1'+str(tauN)] if igen['sd1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                            #getattr( self, 'resolsd'+iRJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( ( ireco['sd2'+str(tauN)]/igen['sd2'+str(tauN)] if igen['sd2'+str(tauN)]!=0 else -999 ), WEIGHT )


                        else: self.ufoResponse = self.ufoResponse+1

                    else:           ##### fakes
                        self.fakes = self.fakes+1
                        if iRecoSel[sys].startswith('_W'): self.fakesW = self.fakesW+1
                        elif iRecoSel[sys].startswith('_top'): self.fakestop = self.fakestop+1
                        else: self.ufoFake = self.ufoFake+1
                        self.eventCategory = 2
                        if not self.onlyTrees:
                            for iRJ,ireco in recoJet.items():
                                if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                                    minLeadAK8JetPt = self.minLeadAK8JetPtW
                                elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                                    minLeadAK8JetPt = self.minLeadAK8JetPtTop

                                if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) > minLeadAK8JetPt ):

                                    if sys.startswith(self.sysWeightList):
                                        getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                        if iRecoSel[sys].startswith('_W'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                        elif iRecoSel[sys].startswith('_top'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                                    else:
                                        getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                        getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                                    #getattr( self, 'fakerecosd'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['sdtau21'], WEIGHT )
                                    #getattr( self, 'fakerecosd'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['sdtau32'], WEIGHT )
                                    for tauN in range(1, self.maxTau+1):
                                        getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                        getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                        getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                        #getattr( self, 'fakerecosd'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd0p5'+str(tauN)], WEIGHT )
                                        #getattr( self, 'fakerecosd'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd1'+str(tauN)], WEIGHT )
                                        #getattr( self, 'fakerecosd'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['sd2'+str(tauN)], WEIGHT )

                if self.onlyTrees and not sys.startswith(self.sysWeightList): self.fillBranches( 'selRecoJets'+sys, recoJet )

        if passGenSel and not self.onlyUnc.startswith(self.sysWeightList):
            #### Misses
            self.miss = self.miss+1
            if iGenSel.startswith('_W'): self.missW = self.missW+1
            elif iGenSel.startswith('_top'): self.misstop = self.misstop+1
            else: self.ufoMiss = self.ufoMiss+1
            self.eventCategory = 3

            genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )
            
            if passGenSel and not self.onlyTrees:
                
                for iGJ,igen in genJet.items():

                    getattr( self, 'missgen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )

                    getattr( self, 'missgen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )
                    #getattr( self, 'missgensd'+iGJ+'_tau21'+iGenSel ).Fill( igen['sdtau21'], event.genWeight )
                    #getattr( self, 'missgensd'+iGJ+'_tau32'+iGenSel ).Fill( igen['sdtau32'], event.genWeight )

                    getattr( self, 'gen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                    getattr( self, 'gen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )

                    getattr( self, 'gen'+iGJ+'_tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                    getattr( self, 'gen'+iGJ+'_tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )
                    #getattr( self, 'gensd'+iGJ+'_tau21'+iGenSel ).Fill( igen['sdtau21'], event.genWeight )
                    #getattr( self, 'gensd'+iGJ+'_tau32'+iGenSel ).Fill( igen['sdtau32'], event.genWeight )

                    for tauN in range(1, self.maxTau+1):
                        getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )
                        #getattr( self, 'missgensd'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['sd0p5'+str(tauN)], event.genWeight )
                        #getattr( self, 'missgensd'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['sd1'+str(tauN)], event.genWeight )
                        #getattr( self, 'missgensd'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['sd2'+str(tauN)], event.genWeight )

                        getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                        getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                        getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )
                        #getattr( self, 'gensd'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['sd0p5'+str(tauN)], event.genWeight )
                        #getattr( self, 'gensd'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['sd1'+str(tauN)], event.genWeight )
                        #getattr( self, 'gensd'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['sd2'+str(tauN)], event.genWeight )

            else:
                self.ufo = self.ufo+1
                if iGenSel.startswith('_W'): self.ufosW = self.ufosW+1
                elif iGenSel.startswith('_top'): self.ufostop = self.ufostop+1
                else: print('Weird ufo', iRecoSel[sys])
                self.eventCategory = 10

            if self.onlyTrees: self.fillBranches( 'selGenJets', genJet )
        
        if self.onlyTrees:
            self.out.fillBranch( 'totalWeight', self.totalWeight )
            self.out.fillBranch( 'eventCategory', self.eventCategory )
            self.out.fillBranch( 'leptonWeight', self.leptonWeight )
            self.out.fillBranch("btagWeight", self.btaggingWeight)
            self.out.fillBranch( 'top_pTreweight_SF', self.topreweight )

            if not self.onlyUnc:
                #self.out.fillBranch( 'pdfWeightUp', self.pdfWeightUp )
                self.out.fillBranch( 'isrWeightUp', self.isrWeightUp )
                self.out.fillBranch( 'fsrWeightUp', self.fsrWeightUp )
                #self.out.fillBranch( 'pdfWeightDown', self.pdfWeightDown )
                self.out.fillBranch( 'isrWeightDown', self.isrWeightDown )
                self.out.fillBranch( 'fsrWeightDown', self.fsrWeightDown )
                self.out.fillBranch( 'puWeightDown', event.puWeightDown if self.isMC else 1 )
                self.out.fillBranch( 'puWeightUp', event.puWeightUp if self.isMC else 1 )
        return True


    #############################################################################
    def recoSelection( self, event, sysUnc=[]  ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))
        met = Object(event, 'MET')
        if self.isMC: genParticles = Collection(event, 'GenPart')
        
        
        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minTightMuonPt and x.looseId and abs(x.p4().Eta()) < self.maxMuonEta and x.pfIsoId>=2 and x.tightId and abs(x.dxy)<0.2 and abs(x.dz)<0.5 and x.miniPFRelIso_all<0.10] # applying tight selection on muons already here since we only veto for loose muons and loose electrons in event

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons)

        ##################################################

        #### MET (not sure if needed)
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)
        ##################################################

        
        #### Basic AK4 b-jet cand. selection
        recoAK4bjets = [ x for x in jets if x.pt > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and x.btagDeepFlavB > self.minBDisc and (x.jetId >= 2)]
        recoAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        
                
        ##################################################

        recoAK8jets = {}
        passSel = {}
        iSel = {}
        #print ("btag wt, before btag:", self.btaggingWeight)
        for sys in self.sysSource:
            #if sys.startswith(self.sysWeightList): continue #previously was 'continue', so some histos of tf the uncertainty related quantities in the last skims didn't really include some weights?
            if sys.startswith(self.sysWeightList): sys = '_nom'
            #### Basic AK8 jet selection
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta and (x.jetId >= 2) ]
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)
            AK8HT = sum( [ getattr( x, 'pt'+sys ) for x in recoAK8jets[sys] ] )
            ##################################################

            ##### Applying selection
            
            passSel[sys], iSel[sys] = self.WtopSelection( False, event, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets[sys], MET, sys)
            if passSel[sys]: print (sys, passSel[sys], iSel[sys])
            #############################
        
        #### Weight #########
        weight=1.
        
        #### b-tagging Weights #####
        if self.isMC: weight = weight*self.btaggingWeight
        #### ################# #####
        

        #### Lepton Weights #####
        if self.isMC:
            if len(recoMuons)>0: leptonWeights = self.leptonSF( "muon", recoMuons[0] )
            else: leptonWeights = [0, 0, 0, 0]
        else: leptonWeights = [1, 1, 1, 1.]
        self.leptonWeight = np.prod(leptonWeights)
        #self.out.fillBranch("leptonWeight", np.prod(self.leptonWeight) )  ### dummy for nanoAOD Tools
        #### ############## #####


        if self.isMC:
            #### Top pT reweighting #####
            itempSel = {}
            tops=[]
            antiTops=[]
            if passSel[sys] and self.topreweight and iSel[sys].startswith('_top'):
                tops =  [x for x in genParticles if x.pdgId==6 and x.statusFlags==14]
                tops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)
                antiTops =  [x for x in genParticles if x.pdgId==-6 and x.statusFlags==14]
                antiTops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)
                
                if (len(tops)>0 and len(antiTops)>0):
                    #print ("top reweighting going on", iSel['_nom'])
                    topSF = math.exp(0.0615 - 0.0005 * tops[0].pt) 
                    antitopSF = math.exp(0.0615 - 0.0005 * antiTops[0].pt)
                    self.topreweight = math.sqrt(topSF*antitopSF)
                else: 
                    self.topreweight = 1.#????math.sqrt(topSF*antitopSF)
            else: 
                self.topreweight = 1.#????math.sqrt(topSF*antitopSF)

            
            #### Applying ALL remaining object-related weights ####
            weight = event.puWeight * event.genWeight * self.leptonWeight * self.topreweight * weight # last term includes btag weight from above
            #print ("All weights:", self.totalWeight, event.puWeight , event.genWeight , self.btaggingWeight, self.leptonWeight , self.topreweight)
            
            if not self.onlyTrees and not self.onlyUnc:
                getattr( self, 'PUweight' ).Fill( event.puWeight )
                getattr( self, 'Lepweight' ).Fill( self.leptonWeight )
                getattr( self, 'Topweight' ).Fill( self.topreweight )
            
        else:
            weight = 1
            #self.leptonWeight = 1.
        self.totalWeight = weight
        #print ("wt after btag* lepton*genwt*puwt at least for MC",self.totalWeight, self.isMC)
        #### ############################################# ####

        if not self.onlyTrees and not self.onlyUnc:

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
            getattr( self, 'HT_noSelnoWeight' ).Fill( AK8HT )
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop_nom )
            getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
            for ijet in recoAK4bjets:
                getattr( self, 'AK4jets_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'AK4jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK4jets_phi_noSelnoWeight' ).Fill( ijet.phi )
            getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )

            #### Checking no selection with weights
            #reweight = self.totalWeight 
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
            getattr( self, 'HT_noSel' ).Fill( AK8HT, weight )
            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop_nom, weight )
            getattr( self, 'nAK4jets_noSel' ).Fill( len(recoAK4bjets), weight )
            for ijet in recoAK4bjets:
                getattr( self, 'AK4jets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK4jets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK4jets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'METPt_noSel' ).Fill( MET.Pt(), weight )
            
            reweight = self.totalWeight 
            
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
                getattr( self, 'HT'+iSel['_nom'] ).Fill( AK8HT, reweight )
                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt'+iSel['_nom'] ).Fill( ijet.pt, reweight )
                    getattr( self, 'AK8jets_eta'+iSel['_nom'] ).Fill( ijet.eta, reweight )
                    getattr( self, 'AK8jets_phi'+iSel['_nom'] ).Fill( ijet.phi, reweight )
                    getattr( self, 'AK8jets_mass'+iSel['_nom'] ).Fill( ijet.msoftdrop_nom, reweight )
                getattr( self, 'nAK4jets'+iSel['_nom'] ).Fill( len(recoAK4bjets), reweight )
                for ijet in recoAK4bjets:
                    getattr( self, 'AK4jets_pt'+iSel['_nom'] ).Fill( ijet.pt, reweight )
                    getattr( self, 'AK4jets_eta'+iSel['_nom'] ).Fill( ijet.eta, reweight )
                    getattr( self, 'AK4jets_phi'+iSel['_nom'] ).Fill( ijet.phi, reweight )
                getattr( self, 'METPt'+iSel['_nom'] ).Fill( getattr(event,'MET_pt'), reweight )
            

        return passSel, iSel, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets, MET

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

        ### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        ##################################################

        #### MET (not sure if needed)
        genMET = ROOT.TLorentzVector()
        genMET.SetPtEtaPhiE(genmet.pt, 0., genmet.phi, 0)
        ##################################################

        #### Basic AK8 jet selection
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta]
        genAK8HT = sum( [ x.pt for x in genAK8jets ] )
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        #### Basic AK4 bjet selection
        genAK4bjets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta and abs(x.hadronFlavour)==5 ]
        genAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        
        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        passSel, iSel = self.WtopSelection( True, event, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET,'' )
        #### Weight
        weight = event.genWeight

        if not self.onlyTrees and not self.onlyUnc: 
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
            getattr( self, 'genHT_noSel' ).Fill( genAK8HT, weight )
            for ijet in genAK8jets:
                getattr( self, 'AK8genjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK8genjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8genjets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8genjets_mass_noSel' ).Fill( ijet.mass, weight )
            getattr( self, 'ngenAK4jets_noSel' ).Fill( len(genAK4bjets), weight )
            for ijet in genAK4bjets:
                getattr( self, 'AK4genjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK4genjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK4genjets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'genMETPt_noSel' ).Fill( genMET.Pt(), weight )
            #ptAsym = ( genAK8jets[0].pt - genAK8jets[1].pt ) / (genAK8jets[0].pt + genAK8jets[1].pt) if len(genAK8jets)>1 else 0
            #deltaPhi = genAK8jets[0].p4().DeltaPhi( genAK8jets[1].p4() ) if len(genAK8jets)>1 else 0
            #getattr( self, 'genPtAsym_noSel' ).Fill( ptAsym, weight )
            #getattr( self, 'genDeltaPhi_noSel' ).Fill( deltaPhi, weight )

            
            ##### Filling histograms
            if passSel and iSel:
                #print ("Filling gen histos")
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
                getattr( self, 'genHT'+iSel ).Fill( genAK8HT, weight )
                for ijet in genAK8jets:
                    getattr( self, 'AK8genjets_pt'+iSel ).Fill( ijet.pt, weight )
                    getattr( self, 'AK8genjets_eta'+iSel ).Fill( ijet.eta, weight )
                    getattr( self, 'AK8genjets_phi'+iSel ).Fill( ijet.phi, weight )
                    getattr( self, 'AK8genjets_mass'+iSel ).Fill( ijet.mass, weight )
                getattr( self, 'ngenAK4jets'+iSel ).Fill( len(genAK4bjets), weight )
                for ijet in genAK4bjets:
                    getattr( self, 'AK4genjets_pt'+iSel ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4genjets_eta'+iSel ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4genjets_phi'+iSel ).Fill( ijet.phi, weight )
                getattr( self, 'genMETPt'+iSel ).Fill( genMET.Pt(), weight )
                #getattr( self, 'genPtAsym'+iSel ).Fill( ptAsym, weight )
                #getattr( self, 'genDeltaPhi'+iSel ).Fill( deltaPhi, weight )

        return passSel, iSel, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET

    #############################################################################
    def WtopSelection( self, isGen, event, muons, electrons, AK4bjets, AK8jets, MET, ptLabel):
    
    
    
        if (len(muons)==1) and (len(electrons)== 0) and (len(AK8jets)>0) and (len(AK4bjets)>=1) and (MET.Pt()>self.METCutWtop):                
            leadJetpT = getattr( AK8jets[0], 'pt'+ptLabel )
            leptWpT=muons[0].p4()+MET
            if leptWpT.Pt()>self.minLeptonicWPt and leadJetpT >=self.minLeadAK8JetPtW:
            
                
                AK4bjets = [x for x in AK4bjets if abs(x.p4().DeltaPhi(muons[0].p4()))<2. and AK8jets[0].p4().DeltaR( x.p4() )>0.8 and x.p4().DeltaR( muons[0].p4() )>0.4 and (muons[0].p4().DeltaR( x.p4() )<np.pi/2.)] 
                
                if (len(AK8jets)>0) and (len(AK4bjets)>=1) and (len(AK4bjets)<3) and abs(AK8jets[0].p4().DeltaPhi(muons[0].p4()))>2.:
                    if not isGen:
                        if self.isMC:
                            bTagSFs =  [x.btagSF_deepjet_M for x in AK4bjets]
                            self.btaggingWeight =  self.getBTagWeight(nBTagged=len(AK4bjets), jet_SFs=bTagSFs)
                        else: self.btaggingWeight = 1
                        
                    jetMass = AK8jets[0].mass if isGen else AK8jets[0].msoftdrop_nom
                    
                    if ((jetMass<=self.maxSDMassW) and (jetMass>self.minSDMassW)) and (leadJetpT>=self.minLeadAK8JetPtW):
                        #print ("Kinematics", jetMass,leadJetpT,'WSel', isGen)
                        return True, '_WSel' 

                    elif (jetMass/ (1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) > self.minSDMassTop) and (jetMass/ (1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop) and (leadJetpT>self.minLeadAK8JetPtTop): 
                        #print ("Kinematics", jetMass,leadJetpT,'topSel', isGen)
                        return True, '_topSel'
                    else: return False, None
                else: return False, None 
            else: return False, None 
        else: return False, None
        
        
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
        '''
        #### Computing Softdrop jets
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
    '''
    def dijetSelection( self, event, muons, electrons, AK8jets, ptAsym, deltaPhi ):

        if (len(muons)+len(electrons)==0) and (len(AK8jets)>1) and (ptAsym<0.3) and (deltaPhi>2) and (AK8jets[1].pt>self.minLeadAK8JetPtDijet):         ### PT REQUIREMENT NEEDS TO BE REVISIT
            return True
        else: return False

    '''
    #############################################################################
    def fillBranches( self, jetLabel, jetInfo ):

        #### Filling branch with passAK8jet info after selection
        self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
        self.out.fillBranch(jetLabel+"_mass", [ iJ['jet'].mass for i,iJ in jetInfo.items() ] )
        self.out.fillBranch(jetLabel+"_pt", [ iJ['jet'].pt for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_eta", [ iJ['jet'].eta for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_phi", [ iJ['jet'].phi for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau21", [ iJ['tau21'] for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau32", [ iJ['tau32'] for i,iJ in jetInfo.items()  ] )
        #self.out.fillBranch(jetLabel+"_SD_Tau21", [ iJ['sdtau21'] for i,iJ in jetInfo.items()  ] )
        #self.out.fillBranch(jetLabel+"_SD_Tau32", [ iJ['sdtau32'] for i,iJ in jetInfo.items()  ] )
        
        for tauN in range(1, self.maxTau+1):
            self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN),  [ iJ['0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN),  [ iJ['1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN),  [ iJ['2'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_0p5_"+str(tauN),  [ iJ['sd0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_1_"+str(tauN),  [ iJ['sd1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_2_"+str(tauN),  [ iJ['sd2'+str(tauN)] for i,iJ in jetInfo.items() ] )



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

###############

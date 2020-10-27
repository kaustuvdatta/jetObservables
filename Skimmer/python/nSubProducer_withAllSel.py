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

    def __init__(self, selection='dijet', sysSource=[], leptonSF={}, year='2016'):
        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.selection = '_'+selection
        self.year = year

        ### Kinematics Cuts AK8Jets ###
        self.minAK8JetPt = 170  ### this is the basic minimum, not the final
        self.maxJetAK8Eta = 2.4

        ### Cuts for selections
        self.minLeadAK8JetPtDijet = 450.
        self.minLeadAK8JetPtW = 200.
        self.minSDMassW = 60.   ### looser to pick low bins
        self.maxSDMassW = 120.  ### looser to pick higher bins
        self.minLeadAK8JetPtTop= 350.
        self.minSDMassTop = 140.
        self.METCutWtop = 40.

        ### Kinematics Cuts Jets ###
        self.minJetPt = 20
        self.maxJetEta = 2.4
        self.minBDisc = 0.3093  ### L: 0.0614, M: 0.3093, T: 07221

        ### Kinenatic Cuts Muons ###
        self.minMuonPt = 20. if self.selection.startswith('_dijet') else 55.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minElectronPt = 35.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]

        self.totalWeight = 1

        ### Defining nsubjetiness basis
        self.maxTau = 5
	self.nSub_labels = ["_tau_0p5_", "_tau_1_", "_tau_2_"]
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
        self.nJet = [ 'Jet1', 'Jet2' ] if self.selection.startswith('_dijet') else [ 'Jet' ]

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
        self.addObject( ROOT.TH1F('cutflow_reco',   ';cutflow_reco',   10, 0, 10) )
        self.addObject( ROOT.TH1F('cutflow_gen',   ';cutflow_gen',   10, 0, 10) )
        self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
        self.addObject( ROOT.TH1F('Lepweight',   ';LepWeight',   20, 0, 2) )
        #### general selection
        selList = (['_dijetSel' ] if self.selection.startswith('_dijet') else [ '_WSel', '_topSel' ] )
        for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
            self.addObject( ROOT.TH1F('cutflowComparison'+isel,   ';cutflowComparison',   10, 0, 10) )
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
            self.addObject( ROOT.TH1F('recoPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
            self.addObject( ROOT.TH1F('recoDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )

            self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

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
            self.addObject( ROOT.TH1F('genPtAsym'+isel,   '; pt Asymmetry', 100, 0, 1.) )
            self.addObject( ROOT.TH1F('genDeltaPhi'+isel,   '; #Delta #Phi(j,j)', 100, 0, 5) )


        for isel in selList:
            for iJ in self.nJet:
                self.addP4Hists( 'gen'+iJ, isel )
                self.addP4Hists( 'missgen'+iJ, isel )
                self.addObject( ROOT.TH1F('gen'+iJ+'Tau21'+isel, ';AK8 gen jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH1F('gen'+iJ+'Tau32'+isel, ';AK8 gen jet #tau_{32}', tauBins, 0, 1) )
                self.addObject( ROOT.TH1F('missgen'+iJ+'Tau21'+isel, ';AK8 missgen jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH1F('missgen'+iJ+'Tau32'+isel, ';AK8 missgen jet #tau_{32}', tauBins, 0, 1) )

                self.addP4Hists( 'fakereco'+iJ, isel )
                self.addObject( ROOT.TH1F('fakereco'+iJ+'Tau21'+isel, ';AK8 fakereco jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH1F('fakereco'+iJ+'Tau32'+isel, ';AK8 fakereco jet #tau_{32}', tauBins, 0, 1) )

                self.addP4Hists( 'uforeco'+iJ, isel )
                self.addP4Hists( 'ufogen'+iJ, isel )

                for tauN in range(1, self.maxTau+1):
                    for x in self.nSub_labels:
                        if tauN>1:
                            maxValueTau=0.4
                            tauBins=400
                        else: maxValueTau=1
                        self.addObject( ROOT.TH1F('gen'+iJ+x+str(tauN)+isel, ';AK8 genjet #tau', tauBins, 0, maxValueTau ) )
                        self.addObject( ROOT.TH1F('missgen'+iJ+x+str(tauN)+isel, ';AK8 missgenjet #tau', tauBins, 0, maxValueTau ) )
                        self.addObject( ROOT.TH1F('fakereco'+iJ+x+str(tauN)+isel, ';AK8 fakejet #tau', tauBins, 0, maxValueTau ) )

                for sysUnc in self.sysSource:
                    self.addP4Hists( 'reco'+iJ, sysUnc+isel )
                    self.addObject( ROOT.TH1F('reco'+iJ+'Tau21'+sysUnc+isel, ';AK8 reco jet #tau_{21}', tauBins, 0, 1) )
                    self.addObject( ROOT.TH2F('resp'+iJ+'Tau21'+sysUnc+isel, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )
                    self.addObject( ROOT.TH1F('reco'+iJ+'Tau32'+sysUnc+isel, ';AK8 reco jet #tau_{32}', tauBins, 0, 1) )
                    self.addObject( ROOT.TH2F('resp'+iJ+'Tau32'+sysUnc+isel, ';AK8 gen jet #tau_{32};AK8 reco jet #tau_{32}', tauBins, 0, 1, tauBins, 0, 1) )
                    if sysUnc.endswith('nom'):
                        self.addObject( ROOT.TH1F('resol'+iJ+'Pt'+isel, ';AK8 reco/gen jet pt', 100, 0, 2) )
                        self.addObject( ROOT.TH1F('resol'+iJ+'SDmass'+isel, ';AK8 reco/gen jet SD mass', 100, 0, 2) )
                        self.addObject( ROOT.TH1F('resol'+iJ+'Tau21'+isel, ';AK8 reco/gen jet #tau_{21}', 100, 0, 2) )
                        self.addObject( ROOT.TH1F('resol'+iJ+'Tau32'+isel, ';AK8 reco/gen jet #tau_{32}', 100, 0, 2) )

                    for tauN in range(1, self.maxTau+1):
                        for x in self.nSub_labels:
                            if tauN>1:
                                maxValueTau=0.4
                                tauBins=400
                            else: maxValueTau=1
                            self.addObject( ROOT.TH1F('reco'+iJ+x+str(tauN)+sysUnc+isel, ';AK8 jet #tau', tauBins, 0, maxValueTau ) )
                            self.addObject( ROOT.TH2F('resp'+iJ+x+str(tauN)+sysUnc+isel, ';AK8 gen jet #tau_{'+str(tauN)+'};AK8 reco jet #tau_{'+str(tauN)+'}', tauBins, 0, maxValueTau, tauBins, 0, maxValueTau) )
                            if sysUnc.endswith('nom'): self.addObject( ROOT.TH1F('resol'+iJ+x+str(tauN)+isel, ';AK8 reco/gen jet #tau_{'+str(tauN)+'}', 100, 0, 2) )


    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -4.0, 4.0 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14159) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 1000) )


    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        self.out = wrappedOutputTree
        self.out.branch('leptonWeight',  "F")
        self.out.branch('eventCategory',  "I")

        for iJ in [ 'Reco', 'Gen' ]:
            self.out.branch('nsel'+iJ+'jets',  'I')  ### dummy for nanoAOD Tools
            self.out.branch('sel'+iJ+'jets_pt',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_eta',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_phi',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_mass',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_Tau21',  'F', lenVar='nsel'+iJ+'jets')
            self.out.branch('sel'+iJ+'jets_Tau32',  'F', lenVar='nsel'+iJ+'jets')

            for tauN in range(1, self.maxTau+1):
                for x in self.nSub_labels:
                    self.out.branch('sel'+iJ+'jets'+x+str(tauN), 'F', lenVar='nsel'+iJ+'jets' )


        pass

    #############################################################################    
    def getBTagWeight(self, nBTagged=0, jet_SFs=[0]): 
        bTagWeight=0
        if len(jet_SFs)>2 or nBTagged>2: 
            print "Error, only leading and subleading AK4 jets are considered: # of btagged jets cannot exceed 2"
        if nBTagged>len(jet_SFs): 
            print "Number of b-tagged jets cannot be greater than number of them for which SFs are provided!"
            return 0
        if nBTagged==0 and len(jet_SFs)==0: return 1
        
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
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    #############################################################################
    def leptonSF(self, lepton, leptonP4 ):

        if lepton.startswith("muon"): leptonP4eta = abs(leptonP4.eta)
        else: leptonP4eta = leptonP4.eta

        SFFileTrigger = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/"+self.leptonSFhelper[lepton]['Trigger'][0] )
        histoSFTrigger = SFFileTrigger.Get( self.leptonSFhelper[lepton]['Trigger'][1] )
        SFTrigger = histoSFTrigger.GetBinContent( histoSFTrigger.GetXaxis().FindBin( leptonP4.pt ), histoSFTrigger.GetYaxis().FindBin( leptonP4eta ) )

        SFFileID = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/"+self.leptonSFhelper[lepton]['ID'][0] )
        histoSFID = SFFileID.Get( self.leptonSFhelper[lepton]['ID'][1] )
        histoSFID_X = histoSFID.GetXaxis().FindBin( leptonP4.pt if self.leptonSFhelper[lepton]['ID'][2] else leptonP4eta )
        histoSFID_Y = histoSFID.GetYaxis().FindBin( leptonP4eta if self.leptonSFhelper[lepton]['ID'][2] else leptonP4.pt )
        SFID = histoSFID.GetBinContent( histoSFID_X, histoSFID_Y )
        SFID = SFID if SFID>0 else 1

        if self.year.startswith('2016') and lepton.startswith("muon"): leptonP4eta = leptonP4.eta    #### stupid fix for the stupid SF file
        SFFileISO = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/"+self.leptonSFhelper[lepton]['ISO'][0] )
        histoSFISO = SFFileISO.Get( self.leptonSFhelper[lepton]['ISO'][1] )
        histoSFISO_X = histoSFISO.GetXaxis().FindBin( leptonP4.pt if self.leptonSFhelper[lepton]['ISO'][2] else leptonP4eta )
        histoSFISO_Y = histoSFISO.GetYaxis().FindBin( leptonP4eta if self.leptonSFhelper[lepton]['ISO'][2] else leptonP4.pt )
        SFISO = histoSFISO.GetBinContent( histoSFISO_X, histoSFISO_Y )
        SFISO = SFISO if SFISO>0 else 1

        #print (SFTrigger * SFID * SFISO), SFTrigger , SFID , SFISO, leptonP4.pt, leptonP4.eta
        return [SFTrigger , SFID , SFISO]


    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''

        self.isMC = event.run == 1

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else: passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenAK4bjets, selGenJets, selGenMET = self.genSelection(event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoAK4bjets, selRecoJets, selRecoMET = self.recoSelection( event )

        genJet = OrderedDict()
        recoJet = OrderedDict()
        if (passGenSel and passRecoSel) and ( iGenSel and iRecoSel ) and ( iGenSel==iRecoSel ):

            getattr( self, 'cutflowComparison'+iRecoSel ).Fill( 1 )
            self.out.fillBranch( 'eventCategory', 1 )

            if self.selection.startswith('_dijet'):
                tmpGenJet1 = selGenJets[0] if abs(selGenJets[0].eta) > abs(selGenJets[1].eta) else selGenJets[1]
                tmpGenJet2 = selGenJets[1] if abs(selGenJets[0].eta) > abs(selGenJets[1].eta) else selGenJets[0]
                genJet['Jet1'] = self.createNsubBasis( tmpGenJet1, event, 'GenPartAK8' )
                genJet['Jet2'] = self.createNsubBasis( tmpGenJet2, event, 'GenPartAK8' )
                tmpRecoJet1 = selRecoJets[0] if abs(selRecoJets[0].eta) > abs(selRecoJets[1].eta) else selRecoJets[1]
                tmpRecoJet2 = selRecoJets[1] if abs(selRecoJets[0].eta) > abs(selRecoJets[1].eta) else selRecoJets[0]
                recoJet['Jet1'] = self.createNsubBasis( tmpRecoJet1, event, 'PFCandsAK8' )
                recoJet['Jet2'] = self.createNsubBasis( tmpRecoJet2, event, 'PFCandsAK8' )
            else:
                genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenPartAK8' )
                recoJet['Jet'] = self.createNsubBasis( selRecoJets[0], event, 'PFCandsAK8' )

            for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                getattr( self, 'gen'+iGJ+'_pt'+iRecoSel ).Fill( igen['jet'].pt, event.genWeight )
                getattr( self, 'gen'+iGJ+'_eta'+iRecoSel ).Fill( igen['jet'].eta, event.genWeight )
                getattr( self, 'gen'+iGJ+'_mass'+iRecoSel ).Fill( igen['jet'].mass, event.genWeight )
                getattr( self, 'gen'+iGJ+'_phi'+iRecoSel ).Fill( igen['jet'].phi, event.genWeight )

                getattr( self, 'gen'+iGJ+'Tau21'+iRecoSel ).Fill( igen['tau21'], event.genWeight )
                getattr( self, 'gen'+iGJ+'Tau32'+iRecoSel ).Fill( igen['tau32'], event.genWeight )
                for tauN in range(1, self.maxTau+1):
                    getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                    getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                    getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( igen['2'+str(tauN)], event.genWeight )


                if iRecoSel.startswith('_dijet'): minLeadAK8JetPt = self.minLeadAK8JetPtDijet
                elif iRecoSel.startswith('_W'): minLeadAK8JetPt = self.minLeadAK8JetPtW
                elif iRecoSel.startswith('_top'): minLeadAK8JetPt = self.minLeadAK8JetPtTop

                WEIGHT =  self.totalWeight
                for sysUnc in self.sysSource:
                    if ( getattr(ireco['jet'], 'pt'+(sysUnc if not sysUnc.startswith('_pu') else '_nom')) > minLeadAK8JetPt ):

                        if sysUnc.startswith('_pu'):
                            WEIGHT = event.genWeight * self.leptonWeight * getattr( event, 'puWeight'+sysUnc.split('pu')[1] )
                            getattr( self, 'reco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                        else:
                            getattr( self, 'reco'+iRJ+'_pt'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'pt'+sysUnc ), WEIGHT )
                            getattr( self, 'reco'+iRJ+'_mass'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sysUnc ), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_eta'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                        getattr( self, 'reco'+iRJ+'_phi'+sysUnc+iRecoSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                        getattr( self, 'reco'+iRJ+'Tau21'+sysUnc+iRecoSel ).Fill( ireco['tau21'], WEIGHT )
                        getattr( self, 'reco'+iRJ+'Tau32'+sysUnc+iRecoSel ).Fill( ireco['tau32'], WEIGHT )
                        for tauN in range(1, self.maxTau+1):
                            getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                            getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( ireco['2'+str(tauN)], WEIGHT )

                    #### FIlling response matrices
                    if self.isMC:
                        if sysUnc.endswith('nom'):
                            getattr( self, 'resol'+iRJ+'Pt'+iRecoSel ).Fill( getattr(ireco['jet'], 'pt')/getattr(igen['jet'], 'pt'), WEIGHT )
                            getattr( self, 'resol'+iRJ+'SDmass'+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop')/getattr(igen['jet'], 'mass'), WEIGHT )
                            getattr( self, 'resol'+iRJ+'Tau21'+iRecoSel ).Fill( ireco['tau21']/igen['tau21'], WEIGHT )
                            getattr( self, 'resol'+iRJ+'Tau32'+iRecoSel ).Fill( ireco['tau32']/igen['tau32'], WEIGHT )

                        getattr( self, 'resp'+iRJ+'Tau21'+sysUnc+iRecoSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                        getattr( self, 'resp'+iRJ+'Tau32'+sysUnc+iRecoSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                        for tauN in range(1, self.maxTau+1):
                            getattr( self, 'resp'+iRJ+'_tau_0p5_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                            getattr( self, 'resp'+iRJ+'_tau_1_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                            getattr( self, 'resp'+iRJ+'_tau_2_'+str(tauN)+sysUnc+iRecoSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )
                            if sysUnc.endswith('nom'):
                                getattr( self, 'resol'+iRJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( ( ireco['0p5'+str(tauN)]/igen['0p5'+str(tauN)] if igen['0p5'+str(tauN)]!=0 else -999 ), WEIGHT )
                                getattr( self, 'resol'+iRJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( ( ireco['1'+str(tauN)]/igen['1'+str(tauN)] if igen['1'+str(tauN)]!=0 else -999 ), WEIGHT )
                                getattr( self, 'resol'+iRJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( ( ireco['2'+str(tauN)]/igen['2'+str(tauN)] if igen['2'+str(tauN)]!=0 else -999 ), WEIGHT )

        #elif passGenSel and not passRecoSel: #( (iGenSel!=iRecoSel) if passRecoSel else False ):
        else:

            self.out.fillBranch( 'eventCategory', 3 )

            if passGenSel and iGenSel:
                self.out.fillBranch( 'eventCategory', 5 )
                getattr( self, 'cutflowComparison'+iGenSel ).Fill( 2 )

                if iGenSel.startswith('_dijet'):
                    tmpGenJet1 = selGenJets[0] if abs(selGenJets[0].eta) > abs(selGenJets[1].eta) else selGenJets[1]
                    tmpGenJet2 = selGenJets[1] if abs(selGenJets[0].eta) > abs(selGenJets[1].eta) else selGenJets[0]
                    genJet['Jet1'] = self.createNsubBasis( tmpGenJet1, event, 'GenPartAK8' )
                    genJet['Jet2'] = self.createNsubBasis( tmpGenJet2, event, 'GenPartAK8' )
                else:
                    genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenPartAK8' )

                for iGJ,igen in genJet.items():

                    getattr( self, 'missgen'+iGJ+'_pt'+iGenSel ).Fill( igen['jet'].pt, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_eta'+iGenSel ).Fill( igen['jet'].eta, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_mass'+iGenSel ).Fill( igen['jet'].mass, event.genWeight )
                    getattr( self, 'missgen'+iGJ+'_phi'+iGenSel ).Fill( igen['jet'].phi, event.genWeight )

                    getattr( self, 'missgen'+iGJ+'Tau21'+iGenSel ).Fill( igen['tau21'], event.genWeight )
                    getattr( self, 'missgen'+iGJ+'Tau32'+iGenSel ).Fill( igen['tau32'], event.genWeight )
                    for tauN in range(1, self.maxTau+1):
                        getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( igen['0p5'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( igen['1'+str(tauN)], event.genWeight )
                        getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( igen['2'+str(tauN)], event.genWeight )

            elif iGenSel:
                print event.event, passGenSel, iGenSel
                print 'pass', passGenSel, iGenSel
                print 'reco', passRecoSel, iRecoSel
                getattr( self, 'cutflowComparison'+iGenSel ).Fill( 4 )
                self.out.fillBranch( 'eventCategory', 6 )
                if iGenSel.startswith('_dijet'):
                    genJet['Jet1'] = selGenJets[0]
                    genJet['Jet2'] = selGenJets[1]
                else:
                    genJet['Jet'] = selGenJets[0]
                for iGJ,igen in genJet.items():
                    getattr( self, 'ufogen'+iGJ+'_pt'+iGenSel ).Fill( igen.pt, event.genWeight )
                    getattr( self, 'ufogen'+iGJ+'_eta'+iGenSel ).Fill( igen.eta, event.genWeight )
                    getattr( self, 'ufogen'+iGJ+'_mass'+iGenSel ).Fill( igen.mass, event.genWeight )
                    getattr( self, 'ufogen'+iGJ+'_phi'+iGenSel ).Fill( igen.phi, event.genWeight )

            if passRecoSel and iRecoSel: #( (iGenSel!=iRecoSel) if passGenSel else False ):

                getattr( self, 'cutflowComparison'+iRecoSel ).Fill( 3 )
                self.out.fillBranch( 'eventCategory', 7 )

                if iRecoSel.startswith('_dijet'):
                    tmpRecoJet1 = selRecoJets[0] if abs(selRecoJets[0].eta) > abs(selRecoJets[1].eta) else selRecoJets[1]
                    tmpRecoJet2 = selRecoJets[1] if abs(selRecoJets[0].eta) > abs(selRecoJets[1].eta) else selRecoJets[0]
                    recoJet['Jet1'] = self.createNsubBasis( tmpRecoJet1, event, 'PFCandsAK8' )
                    recoJet['Jet2'] = self.createNsubBasis( tmpRecoJet2, event, 'PFCandsAK8' )
                else:
                    recoJet['Jet'] = self.createNsubBasis( selRecoJets[0], event, 'PFCandsAK8' )

                for iRJ,ireco in recoJet.items():

                    if iRecoSel.startswith('_dijet'): minLeadAK8JetPt = self.minLeadAK8JetPtDijet
                    elif iRecoSel.startswith('_W'): minLeadAK8JetPt = self.minLeadAK8JetPtW
                    elif iRecoSel.startswith('_top'): minLeadAK8JetPt = self.minLeadAK8JetPtTop

                    if ( getattr(ireco['jet'], 'pt_nom') > minLeadAK8JetPt ):

                        getattr( self, 'fakereco'+iRJ+'_pt'+iRecoSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), self.totalWeight )
                        getattr( self, 'fakereco'+iRJ+'_mass'+iRecoSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), self.totalWeight )
                        getattr( self, 'fakereco'+iRJ+'_eta'+iRecoSel ).Fill( getattr(ireco['jet'], 'eta'), self.totalWeight )
                        getattr( self, 'fakereco'+iRJ+'_phi'+iRecoSel ).Fill( getattr(ireco['jet'], 'phi'), self.totalWeight )
                        getattr( self, 'fakereco'+iRJ+'Tau21'+iRecoSel ).Fill( ireco['tau21'], self.totalWeight )
                        getattr( self, 'fakereco'+iRJ+'Tau32'+iRecoSel ).Fill( ireco['tau32'], self.totalWeight )
                        for tauN in range(1, self.maxTau+1):
                            getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+iRecoSel ).Fill( ireco['0p5'+str(tauN)], self.totalWeight )
                            getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+iRecoSel ).Fill( ireco['1'+str(tauN)], self.totalWeight )
                            getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+iRecoSel ).Fill( ireco['2'+str(tauN)], self.totalWeight )

            elif iRecoSel:

                self.out.fillBranch( 'eventCategory', 8 )
                if iRecoSel.startswith('_dijet'):
                    recoJet['Jet1'] = selRecoJets[0]
                    recoJet['Jet2'] = selRecoJets[1]
                else:
                    recoJet['Jet'] = selRecoJets[0]
                getattr( self, 'cutflowComparison'+iRecoSel ).Fill( 4 )
                for iRJ,ireco in recoJet.items():
                    getattr( self, 'uforeco'+iRJ+'_pt'+iRecoSel ).Fill( ireco.pt, self.totalWeight )
                    getattr( self, 'uforeco'+iRJ+'_mass'+iRecoSel ).Fill( ireco.msoftdrop, self.totalWeight )
                    getattr( self, 'uforeco'+iRJ+'_eta'+iRecoSel ).Fill( ireco.eta, self.totalWeight )
                    getattr( self, 'uforeco'+iRJ+'_phi'+iRecoSel ).Fill( ireco.phi, self.totalWeight )

        self.fillBranches( 'selRecojets', recoJet )
        self.fillBranches( 'selGenjets', genJet )

        return True


    #############################################################################
    def recoSelection( self, event ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))
        met = Object(event, 'MET')

        #### Lepton selection
        recoElectrons  = [x for x in electrons if x.pt>self.minElectronPt and x.cutBased_HEEP and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))]
        recoElectrons.sort(key=lambda x:x.pt, reverse=True)

        recoMuons = [ x for x in muons if x.pt > self.minMuonPt and x.highPtId > 1 and abs(x.p4().Eta()) < self.maxMuonEta and x.pfRelIso03_all < 0.1]
        recoMuons.sort(key=lambda x:x.pt, reverse=True)
        nleptons = len(recoMuons)+len(recoElectrons)
        ##################################################

        #### MET (not sure if needed)
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)
        ##################################################

        #### Basic AK8 jet selection
        recoAK8jets = [ x for x in AK8jets if x.pt_nom > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta ] # and x.jetId>4]
        AK8HT = sum( [ x.pt for x in recoAK8jets ] )
        recoAK8jets.sort(key=lambda x:x.pt_nom,reverse=True)
        ##################################################

        #### Basic AK4 bjet selection
        recoAK4bjets = [ x for x in jets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta and x.btagDeepB > self.minBDisc ] # and x.jetId>4 ]
        recoAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        #### Weight
        if self.isMC and self.selection.startswith('_Wtop'):
            if len(recoMuons)>0: leptonWeights= self.leptonSF( "muon", recoMuons[0] )
            else: leptonWeights = [0, 0, 0]
        else: leptonWeights = [1, 1, 1]

        if self.isMC:
            weight = event.puWeight * event.genWeight * np.prod(leptonWeights)
            getattr( self, 'PUweight' ).Fill( event.puWeight )
            getattr( self, 'Lepweight' ).Fill( np.prod(leptonWeights) )
            self.leptonWeight = np.prod(leptonWeights)
        else:
            weight = 1
            self.leptonWeight = 1
        self.out.fillBranch("leptonWeight", np.prod(leptonWeights) )  ### dummy for nanoAOD Tools
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
        getattr( self, 'HT_noSelnoWeight' ).Fill( AK8HT )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
            getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop )
        getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
        for ijet in recoAK4bjets:
            getattr( self, 'AK4jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK4jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK4jets_phi_noSelnoWeight' ).Fill( ijet.phi )
        getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )

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
        getattr( self, 'HT_noSel' ).Fill( AK8HT, weight )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop, weight )
        getattr( self, 'nAK4jets_noSel' ).Fill( len(recoAK4bjets), weight )
        for ijet in recoAK4bjets:
            getattr( self, 'AK4jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK4jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK4jets_phi_noSel' ).Fill( ijet.phi, weight )
        getattr( self, 'METPt_noSel' ).Fill( MET.Pt(), weight )
        ptAsym = ( recoAK8jets[0].pt - recoAK8jets[1].pt ) / (recoAK8jets[0].pt + recoAK8jets[1].pt) if len(recoAK8jets)>1 else 0
        deltaPhi = recoAK8jets[0].p4().DeltaPhi( recoAK8jets[1].p4() ) if len(recoAK8jets)>1 else 0
        getattr( self, 'recoPtAsym_noSel' ).Fill( ptAsym, weight )
        getattr( self, 'recoDeltaPhi_noSel' ).Fill( deltaPhi, weight )

        ##### Applying selection
        passSel = False
        iSel = None

        if self.selection.startswith('_dijet'):

            passSel = self.dijetSelection( event, recoMuons, recoElectrons, recoAK8jets, ptAsym, deltaPhi )
            iSel = '_dijetSel' if passSel else None

        elif self.selection.startswith('_Wtop'):

            passSel, iSel = self.WtopSelection( False, event, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets, MET )

        #### Creating Nsub basis, filling histos and creating branches IF passSel
        if passSel and iSel:

            #### Basic reco histos
            getattr( self, 'nPVs'+iSel ).Fill( getattr( event, 'PV_npvsGood'), weight )
            getattr( self, 'nleps'+iSel ).Fill( len(recoMuons)+len(recoElectrons), weight )
            for imuon in recoMuons:
                getattr( self, 'muons_pt'+iSel ).Fill( imuon.pt, weight )
                getattr( self, 'muons_eta'+iSel ).Fill( imuon.eta, weight )
                getattr( self, 'muons_phi'+iSel ).Fill( imuon.phi, weight )
            for iele in recoElectrons:
                getattr( self, 'eles_pt'+iSel ).Fill( iele.pt, weight )
                getattr( self, 'eles_eta'+iSel ).Fill( iele.eta, weight )
                getattr( self, 'eles_phi'+iSel ).Fill( iele.phi, weight )
            getattr( self, 'nAK8jets'+iSel ).Fill( len(recoAK8jets), weight )
            getattr( self, 'HT'+iSel ).Fill( AK8HT, weight )
            for ijet in recoAK8jets:
                getattr( self, 'AK8jets_pt'+iSel ).Fill( ijet.pt, weight )
                getattr( self, 'AK8jets_eta'+iSel ).Fill( ijet.eta, weight )
                getattr( self, 'AK8jets_phi'+iSel ).Fill( ijet.phi, weight )
                getattr( self, 'AK8jets_mass'+iSel ).Fill( ijet.msoftdrop, weight )
            getattr( self, 'nAK4jets'+iSel ).Fill( len(recoAK4bjets), weight )
            for ijet in recoAK4bjets:
                getattr( self, 'AK4jets_pt'+iSel ).Fill( ijet.pt, weight )
                getattr( self, 'AK4jets_eta'+iSel ).Fill( ijet.eta, weight )
                getattr( self, 'AK4jets_phi'+iSel ).Fill( ijet.phi, weight )
            getattr( self, 'METPt'+iSel ).Fill( getattr(event,'MET_pt'), weight )
            getattr( self, 'recoPtAsym'+iSel ).Fill( ptAsym, weight )
            getattr( self, 'recoDeltaPhi'+iSel ).Fill( deltaPhi, weight )

        return passSel, iSel, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets, MET

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing reco information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genJets = list(Collection( event, 'GenJet'))
        genParticles = Collection(event, 'GenPart')
        genmet = Object( event, 'GenMET')

        ### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minElectronPt and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        ##################################################

        #### MET (not sure if needed)
        genMET = ROOT.TLorentzVector()
        genMET.SetPtEtaPhiE(genmet.pt, 0., genmet.phi, 0)
        ##################################################

        #### Basic AK8 jet selection
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta ]
        genAK8HT = sum( [ x.pt for x in genAK8jets ] )
        genAK8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        #### Basic AK4 bjet selection
        genAK4bjets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta and abs(x.hadronFlavour)==5 ]
        genAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

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
        ptAsym = ( genAK8jets[0].pt - genAK8jets[1].pt ) / (genAK8jets[0].pt + genAK8jets[1].pt) if len(genAK8jets)>1 else 0
        deltaPhi = genAK8jets[0].p4().DeltaPhi( genAK8jets[1].p4() ) if len(genAK8jets)>1 else 0
        getattr( self, 'genPtAsym_noSel' ).Fill( ptAsym, weight )
        getattr( self, 'genDeltaPhi_noSel' ).Fill( deltaPhi, weight )

        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        passSel = False
        iSel = None

        if self.selection.startswith('_dijet'):

            passSel = self.dijetSelection( event, genMuons, genElectrons, genAK8jets, ptAsym, deltaPhi )
            iSel = '_dijetSel' if passSel else None

        elif self.selection.startswith('_Wtop'):

            passSel, iSel = self.WtopSelection( True, event, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET )

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
            getattr( self, 'genPtAsym'+iSel ).Fill( ptAsym, weight )
            getattr( self, 'genDeltaPhi'+iSel ).Fill( deltaPhi, weight )

        return passSel, iSel, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET


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
    def dijetSelection( self, event, muons, electrons, AK8jets, ptAsym, deltaPhi ):

        if (len(muons)+len(electrons)==0) and (len(AK8jets)>1) and (ptAsym<0.3) and (deltaPhi>2) and (AK8jets[1].pt>self.minLeadAK8JetPtDijet):         ### PT REQUIREMENT NEEDS TO BE REVISIT
            return True
        else: return False

    #############################################################################
    def WtopSelection( self, isGen, event, muons, electrons, AK4bjets, AK8jets, MET ):

        if (len(muons)==1) and (len(AK8jets)>0) and (len(AK4bjets)>1) and (MET.Pt()>self.METCutWtop):

            ### removing ak4 jets inside leadAK8 jet
            for bjet in AK4bjets:
                if AK8jets[0].p4().DeltaR( bjet.p4() )<0.8 : AK4bjets.remove(bjet)

            ### defining muon isolation and leptonic top
            muonIso = []
            leptonicTop = []
            for bjet in AK4bjets:
                if bjet.p4().DeltaR( muons[0].p4() )<0.4: muonIso.append( False )
                if (muons[0].p4().DeltaR( bjet.p4() )>0.4) and (muons[0].p4().DeltaR( bjet.p4() )<1.5): leptonicTop.append( True )

            if all(muonIso) and ((MET+muons[0].p4()).Pt()>200) and any(leptonicTop) and (AK8jets[0].p4().DeltaR(muons[0].p4())>0.8) and (len(AK4bjets)>1):
                jetMass = AK8jets[0].mass if isGen else AK8jets[0].msoftdrop
                if (jetMass<self.maxSDMassW) and (jetMass>self.minSDMassW) and (AK8jets[0].pt>self.minLeadAK8JetPtW): return True, '_WSel'
                elif (jetMass>self.minSDMassTop) and (AK8jets[0].pt>self.minLeadAK8JetPtTop): return True, '_topSel'
                else: return False, None
            else: return False, None

        else: return False, None

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


import ROOT
import math, os, sys
import numpy as np
#import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nsubjettinessProducer(Module):
    def __init__(self):
        self.writeHistFile=True
        self.verbose = False
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.ptmin = 200
        self.maxTau = 7

        
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)
        self.beta4_W_par = [ 0.49631107, -0.91138405, 0.11273356, -0.26872, 0.04777313, 0.06521296, 0.51138633, -0.03617261] # 4-body N-subjettiness product observable params obtained via linear regression method of arXiv:1902:07180 
        ### Observable params correspond to optimization over boosted W samples from W'->WZ(Z->vv) with pT_min=200 GeV, loose mass cut of [70,115] 


        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )

        self.dummy = 0;

        self.jetBranchName = "FatJet"
        self.AK4JetBranchName = "Jet"
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        self.metBranchName = "MET"
        self.muonBranchName = "Muon"
        self.electronBranchName = "Electron"
        self.subJetBranchName = "SubJet"
        self.puppiMETBranchName = "PuppiMET"
        self.rawMETBranchName = "RawMET"
        self.eventBranchName = "event"
        self.pfCandsBranchName = "PFCandsAK8"
        self.bname = "sbd0"

        self.genCandsBranchName = "GenPartAK8"
        self.genJetBranchName = "GenJetAK8"
        self.AK4GenJetBranchName = "GenJet"
        self.genPartBranchName = "GenPart"
        self.genMETBranchName = "GenMET"
        self.subGenJetAK8BranchName = "SubGenJetAK8"

        ### Kinematics Cuts Jets ###
        self.minJetPt = 200.
        self.minJetSDMass = 50.
        self.maxJetSDMass = 150.
        self.maxJetEta = 2.4

        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 53.
        self.maxTightMuonEta = 2.1
        self.minLooseMuonPt = 20.
        self.maxLooseMuonEta = 2.4
        self.minMuonMET = 40.

        ### Kinenatic Cuts Electrons ###
        self.minTightElectronPt = 120.
        self.minLooseElectronPt = 35.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]  
        self.minElectronMET = 80.

        ### Kinematic Cuts b-jets ###
        self.minbJetPt = 30.
        self.maxbJetEta = 2.4
        
        ### b-tag
        self.minBDisc = 0.8484
        ### Medium https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
        ### from https://github.com/thaarres/WTopScalefactorProducer/blob/master/Skimmer/python/TTSkimmer.py


        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")




    def endJob(self):
        Module.endJob(self)
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for ijet in ["0", "1"] :
            self.out.branch("goodrecojet" + ijet + "_pt",  "F")
            self.out.branch("goodrecojet" + ijet + "_eta",  "F")
            self.out.branch("goodrecojet" + ijet + "_phi",  "F")
            self.out.branch("goodrecojet" + ijet + "_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_softdrop_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_tau21",  "F")
            self.out.branch("goodrecojet" + ijet + "_N21",  "F")
            #self.out.branch("goodrecojet" + ijet + "_Beta_4_W",  "F")
            #print 'beginFIle', ijet
            for tauN in range(self.maxTau):
                self.out.branch("goodrecojet" + ijet + "_tau_0p5_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_1_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_2_"+str(tauN),  "F")
	    self.out.branch("MET1",  "F")
            self.out.branch("leptonicW_pT1",  "F")            
	    self.out.branch("MET2",  "F")
            self.out.branch("leptonicW_pT2",  "F")            
	    self.out.branch("lepton_pT",  "F")

        pass


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        #isMC = event.run == 1

        self.dummy+=1
        if (self.dummy > 10000): return False
        if self.verbose: print ('Event : ', event.event)
        if self.dummy%1000==0: print ("Analyzing events...", self.dummy)
            

        ### Get W->jj candidate ###
        
        jets = list(Collection(event, self.jetBranchName ))
        pfCands = list(Collection(event, self.pfCandsBranchName ))
        electrons = list(Collection(event, self.electronBranchName))
        muons = list(Collection(event, self.muonBranchName))
        bjets = list(Collection(event, self.AK4JetBranchName))
        met = Object(event, self.metBranchName)
        muonTrigger = event.HLT_Mu50
        #electronTrigger = event.HLT_Ele115_CaloIdVT_GsfTrkIdT
        
        ### Applying selections as per recommendations of AN2016_215_v3 (Aarestad, T. et al) ###
        
        # applying basic selections to loose leptons, to decide on veto or not (only using muons for now)

        
        recoLooseElectrons  = [x for x in electrons if x.pt>self.minLooseElectronPt and x.cutBased_HEEP] #electron events not kept in the end

        recoLooseMuons = [ x for x in muons if x.pt > self.minLooseMuonPt and abs(x.p4().Eta()) < self.maxLooseMuonEta]# and x.pfRelIso03_all < 0.1]   
        recoLooseElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoLooseMuons.sort(key=lambda x:x.pt, reverse=True)
        
        recoLepton=ROOT.TLorentzVector()
        lepflag = -1
        if len(recoLooseMuons)+len(recoLooseElectrons)==1:
            if len(recoLooseMuons)==1:                
                if muonTrigger==0: return False
                #if recoLooseMuons[0].pt<self.minTightMuonPt: return False
                lepflag=1
                recoLepton = recoLooseMuons[0].p4()
            if len(recoLooseElectrons)==1: #accepting only events with muons right now
                return False
        else:
            return False

        #if lepflag==1 and met.pt<self.minMuonMET: return False
        MET1 = ROOT.TLorentzVector()
        MET1.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)

	MET2 = ROOT.TLorentzVector()
        MET2.SetPtEtaPhiE(met.pt, 0., met.phi, 0)

        leptonicWCand1 = MET1+recoLepton
        leptonicWCand2 = MET2+recoLepton
	
        ### applying basic selection to jets if any, remaining AK8 jets are thus hadronic  W Candidates

        recojets = [ x for x in jets if x.p4().Perp() > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and x.msoftdrop>self.minJetSDMass]
        recojets.sort(key=lambda x:x.p4().Perp(),reverse=True)
	

        if len(recojets)<1: #exit if no AK8jets in event
            return False
        
        #recoAK8 = ROOT.TLorentzVector()
        #recoAK8.SetPtEtaPhiM(recojets[0].pt,recojets[0].eta,recojets[0].phi,recojets[0].mass)
        #dR_AK8Lepton = recoAK8.DeltaR(recoLepton)
	#if abs(dR_AK8Lepton)<1. : return False #take care of lepton overlap

        ### applying basic selection to ak4-jets if any, if passed these are b-jet candidates
        recoAK4jets = [ x for x in bjets if x.p4().Perp() > self.minbJetPt and abs(x.p4().Eta()) < self.maxbJetEta and x.btagCSVV2 > self.minBDisc]
	#for x in recoAK4jets:
	#    print x.btagCSVV2

	if len(recoAK4jets)<1: return False #exit if no AK4-jets in event
        

	#recoAK4 = ROOT.TLorentzVector()
	#recoAK4.SetPtEtaPhiM(recoAK4jets[0].pt, recoAK4jets[0].eta,recoAK4jets[0].phi,recoAK4jets[0].mass)
	#if abs(recoAK4.DeltaR(recoLepton))<0.3:
	#    return False
	#if abs(recoAK4.DeltaR(recoAK8))<0.8:
        #    return False


        pfCandsVec = ROOT.vector("TLorentzVector")()

        for p in pfCands :
            pfCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        

        for irecojet,recojet in enumerate(recojets):

            # Cluster only the particles near the appropriate jet to save time
            constituents = ROOT.vector("TLorentzVector")()

            for x in pfCandsVec:

                if recojet.p4().DeltaR( x ) < 0.8:
                    constituents.push_back(x)
                    nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
                    nsub1 = self.nSub1.getTau( self.maxTau, constituents )
                    nsub2 = self.nSub2.getTau( self.maxTau, constituents )
                    
            if (irecojet < 2 ): #to extract only the leading and sub-leading jet

                self.out.fillBranch("goodrecojet" + str(irecojet) + "_pt",  recojet.p4().Pt() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_eta",  recojet.p4().Eta() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_phi",  recojet.p4().Phi() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_mass",  recojet.p4().M() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_softdrop_mass", recojet.msoftdrop)
                self.out.fillBranch("leptonicW_pT1",leptonicWCand1.Perp())                	
                self.out.fillBranch("leptonicW_pT2",leptonicWCand2.Perp())                	
                self.out.fillBranch("lepton_pT",recoLepton.Perp())
                self.out.fillBranch("MET1", MET1.Perp())
                self.out.fillBranch("MET2", MET2.Perp())
                
                if recojet.tau1 > 0.: 
                    tau21 = recojet.tau2/recojet.tau1
                else:
                    tau21 = -1. 
                    
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau21", tau21)
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_N21",  recojet.n2b1)

                for tauN in range(self.maxTau):
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_1_"+str(tauN),  nsub1[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_2_"+str(tauN),  nsub2[tauN]  )
                
        return True

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

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#nSub0 = lambda : nsubjettinessProducer(beta=0.0, zcut=0.1, bname="sdb0")
#nSub1 = lambda : nsubjettinessProducer(beta=1.0, zcut=0.1, bname="sdb1")

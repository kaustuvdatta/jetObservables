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
        #self.beta4_W_par = [ 0.49631107, -0.91138405, 0.11273356, -0.26872, 0.04777313, 0.06521296, 0.51138633, -0.03617261] # 4-body N-subjettiness product observable params obtained via linear regression method of arXiv:1902:07180 
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
        #self.genPartBranchName = "GenPart"
        self.genMETBranchName = "GenMET"
        self.subGenJetAK8BranchName = "SubGenJetAK8"

        ### Kinematics Cuts Jets ###
        self.minJetPt = 200.
        self.minJetSDMass = 0.
        self.maxJetSDMass = 250.
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

        for ijet in ["0"] :
            self.out.branch("goodgenjet" + ijet + "_pt",  "F")
            self.out.branch("goodgenjet" + ijet + "_eta",  "F")
            self.out.branch("goodgenjet" + ijet + "_phi",  "F")
            self.out.branch("goodgenjet" + ijet + "_mass",  "F")
            self.out.branch("goodgenjet" + ijet + "_tau21", "F")
            self.out.branch("dR_gen_reco_AK8", "F")
            self.out.branch("goodrecojet" + ijet + "_pt",  "F")
            self.out.branch("goodrecojet" + ijet + "_eta",  "F")
            self.out.branch("goodrecojet" + ijet + "_phi",  "F")
            self.out.branch("goodrecojet" + ijet + "_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_softdrop_mass",  "F")
            self.out.branch("genEventNo_taus_are_0",  "I")
            self.out.branch("recoEventNo_taus_are_0",  "I")
            self.out.branch("goodrecojet" + ijet + "_tau21",  "F")
            self.out.branch("goodrecojet" + ijet + "_N21",  "F")
            #self.out.branch("goodrecojet" + ijet + "_Beta_4_W",  "F")
            #print 'beginFIle', ijet
            for tauN in range(self.maxTau):
                self.out.branch("goodrecojet" + ijet + "_tau_0p5_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_1_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_2_"+str(tauN),  "F")
		self.out.branch("goodgenjet" + ijet + "_tau_0p5_"+str(tauN),  "F")
                self.out.branch("goodgenjet" + ijet + "_tau_1_"+str(tauN),  "F")
                self.out.branch("goodgenjet" + ijet + "_tau_2_"+str(tauN),  "F")

        pass


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        isMC = event.run == 1

        self.dummy+=1
        #if (self.dummy > 50000): return False
        #if self.verbose: print 'Event : ', event.event
        if self.dummy%1000==0: print ("Analyzing events...", self.dummy)
            

        ### Get W->jj candidate ###
        
        jets = list(Collection(event, self.jetBranchName ))
        pfCands = list(Collection(event, self.pfCandsBranchName ))
        gen_ak8 = list(Collection(event, self.genJetBranchName ))
        #genJets = Collection(event, self.genJetBranchName )
        #genCands = Collection(event, self.genCandsBranchName )

        if isMC==True:
            genParticles = Collection(event, self.genCandsBranchName)


        recojets = [ x for x in jets if x.p4().Perp() > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and self.minJetSDMass<x.msoftdrop<self.maxJetSDMass]# and (x.p4().DeltaR(recoLepton)>=1.)]
        recojets.sort(key=lambda x:x.p4().Perp(),reverse=True)
        
        if len(recojets)==0: #exit if no AK8jets in event
            return False
        
        recoAK8 = ROOT.TLorentzVector()
        recoAK8.SetPtEtaPhiM(recojets[0].pt,recojets[0].eta,recojets[0].phi,recojets[0].mass)

        genjets = [x for x in gen_ak8]
        genjets.sort(key=lambda x:x.pt,reverse=True)

        dRmin=[0.2,0]

        for i in xrange(0,len(genjets)):
            genjet_4v = ROOT.TLorentzVector()
            genjet_4v.SetPtEtaPhiM(genjets[i].pt, genjets[i].eta, genjets[i].phi, genjets[i].mass)
            if abs(genjet_4v.DeltaR(recoAK8))<dRmin[0]:
                dRmin[0] = abs(genjet_4v.DeltaR(recoAK8))
                dRmin[1] = i

	
	goodgenjet = genjets[dRmin[1]]	
	genlevelAK8 = ROOT.TLorentzVector()
        genlevelAK8.SetPtEtaPhiM(goodgenjet.pt, goodgenjet.eta, goodgenjet.phi, goodgenjet.mass)	

        ### applying basic selection to ak4-jets if any, if pass these are b-jet candidates
        

        pfCandsVec = ROOT.vector("TLorentzVector")()

        for p in pfCands :
            pfCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )

        for irecojet,recojet in enumerate(recojets):
            ##recojetsGroomed[irecojet] = {}
            # Cluster only the particles near the appropriate jet to save time

            constituents = ROOT.vector("TLorentzVector")()

            for x in pfCandsVec:

                if recojet.p4().DeltaR( x ) < 0.8:
                    constituents.push_back(x)
                nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
                nsub1 = self.nSub1.getTau( self.maxTau, constituents )
                nsub2 = self.nSub2.getTau( self.maxTau, constituents )
                    
            if (irecojet < 1 ): #to extract only the leading and sub-leading jets 
                #print "irecojet", irecojet

                self.out.fillBranch("goodrecojet" + str(irecojet) + "_pt",  recojet.p4().Pt() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_eta",  recojet.p4().Eta() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_phi",  recojet.p4().Phi() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_mass",  recojet.p4().M() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_softdrop_mass", recojet.msoftdrop)
               
                
                if recojet.tau1 > 0.0: 
                    tau21 = recojet.tau2/recojet.tau1
                else:
                    tau21 = -1. 
                    
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau21", tau21)
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_N21",  recojet.n2b1)
                #self.out.fillbranch("goodrecojet" + str(irecojet) + "_Beta_4_W",  "F")
           

                for tauN in range(self.maxTau):
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_1_"+str(tauN),  nsub1[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_2_"+str(tauN),  nsub2[tauN]  )
                #continue
            

	genpart_CandsVec = ROOT.vector("TLorentzVector")()

        for p in genParticles :
            genpart_CandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        

        # Cluster only the particles near the appropriate jet to save time
        constituents = ROOT.vector("TLorentzVector")()

        for x in genpart_CandsVec:

            if abs(goodgenjet.p4().DeltaR( x )) < 0.8:
                constituents.push_back(x)
            nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
            nsub1 = self.nSub1.getTau( self.maxTau, constituents )
            nsub2 = self.nSub2.getTau( self.maxTau, constituents )
                    
	self.out.fillBranch("goodgenjet0_pt",  goodgenjet.pt)
        self.out.fillBranch("goodgenjet0_eta",  goodgenjet.p4().Eta() )
        self.out.fillBranch("goodgenjet0_phi",  goodgenjet.p4().Phi() )
        self.out.fillBranch("goodgenjet0_mass",  goodgenjet.p4().M() )
        self.out.fillBranch("dR_gen_reco_AK8",  dRmin[0])   
        tau11 = 0
        tau21 = 0 
        for tauN in range(self.maxTau):
            if tauN==0:
                tau11 = nsub1[tauN]
		if nsub0p5[tauN]==0. or nsub1[tauN]==0. or nsub2[tauN]==0.:
		    self.out.fillBranch("genEventNo_taus_are_0", event.event) 
                    #print nsub0p5, nsub1, nsub2
            if tauN==1:
                tau21 = nsub1[tauN]
		if nsub0p5[tauN]==0. or nsub1[tauN]==0. or nsub2[tauN]==0.:
		    self.out.fillBranch("genEventNo_taus_are_0", event.event) 
                    #print nsub0p5, nsub1, nsub2
	    
                
            self.out.fillBranch("goodgenjet0_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
            self.out.fillBranch("goodgenjet0_tau_1_"+str(tauN),  nsub1[tauN]  )
            self.out.fillBranch("goodgenjet0_tau_2_"+str(tauN),  nsub2[tauN]  )
        
        if tau11!=0.:
            self.out.fillBranch("goodgenjet0_tau21", tau21/tau11 )
        else:
            self.out.fillBranch("goodgenjet0_tau21", -1. )
        
		

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

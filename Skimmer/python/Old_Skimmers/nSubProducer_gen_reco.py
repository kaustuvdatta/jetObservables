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
        self.maxTau = 10

       
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)
  

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

        self.genCandsAK8BranchName = "GenPartAK8"
        self.genCandsBranchName = "GenPart"
        self.genJetBranchName = "GenJetAK8"
        self.AK4GenJetBranchName = "GenJet"
        #self.genPartBranchName = "GenPart"
        self.genMETBranchName = "GenMET"
        #self.subGenJetAK8BranchName = "SubGenJetAK8"

        
        ### Kinematics Cuts Jets ###
        self.minJetPt = 200.
        self.minJetSDMass = 40.
        self.maxJetSDMass = 250.
        self.maxJetEta = 2.4

        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 53.
        self.maxTightMuonEta = 2.1
        self.minMuonPt = 20.
        self.maxMuonEta = 2.4
        self.minMuonMET = 40.

        ### Kinenatic Cuts Electrons ###
        self.minTightElectronPt = 120.
        self.minElectronPt = 35.
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
            self.out.branch("goodrecojet" + ijet + "_pt",  "F")
            self.out.branch("goodrecojet" + ijet + "_eta",  "F")
            self.out.branch("goodrecojet" + ijet + "_phi",  "F")
            self.out.branch("goodrecojet" + ijet + "_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_softdrop_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_tau21",  "F")
            self.out.branch("goodrecojet" + ijet + "_N21",  "F")

            self.out.branch("goodgenjet" + ijet + "_pt",  "F")
            self.out.branch("goodgenjet" + ijet + "_eta",  "F")
            self.out.branch("goodgenjet" + ijet + "_phi",  "F")
            self.out.branch("goodgenjet" + ijet + "_mass",  "F")
            self.out.branch("goodgenjet" + ijet + "_tau21", "F")

            self.out.branch("dr_LepJet",  "F")
            self.out.branch("dphi_LepJet",  "F")
            self.out.branch("dphi_LepMet", "F") 
            self.out.branch("dphi_MetJet",  "F")
            self.out.branch("dphi_WJet"  ,  "F")
        
            #self.out.branch("goodrecojet" + ijet + "_Beta_4_W",  "F")
            #print 'beginFIle', ijet
            for tauN in range(self.maxTau):
                self.out.branch("goodrecojet" + ijet + "_tau_0p5_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_1_"+str(tauN),  "F")
                self.out.branch("goodrecojet" + ijet + "_tau_2_"+str(tauN),  "F")
                self.out.branch("goodgenjet" + ijet + "_tau_0p5_"+str(tauN),  "F")
                self.out.branch("goodgenjet" + ijet + "_tau_1_"+str(tauN),  "F")
                self.out.branch("goodgenjet" + ijet + "_tau_2_"+str(tauN),  "F")

            self.out.branch("dR_gen_reco_AK8", "F")
            self.out.branch("dR_genW_genAK8", "F")
            self.out.branch("MET",  "F")
            self.out.branch("leptonicW_pT",  "F")            
            self.out.branch("lepton_pT",  "F")
            self.out.branch("passedMETfilters",  "I")
            self.out.branch("genmatchedrecoAK8",  "I")
            self.out.branch("genmatchedgenAK8",  "I")
            self.out.branch("genEventNo_taus_are_0",  "I")
            self.out.branch("recoEventNo_taus_are_0",  "I")
            self.out.branch("PID_mW",  "I")
            self.out.branch("PID_umW",  "I")
            self.out.branch("PID_W",  "I")
            self.out.branch("GPIdx_W",  "I")
            self.out.branch("lepFlag", "I")
            self.out.branch("dR_AK8qi_mW",  "F")
            self.out.branch("dR_AK8motherW_mW",  "F")
            self.out.branch("dR_AK8qi_umW",  "F")
            self.out.branch("dR_AK8motherW_W",  "F")

        pass


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        isMC = event.run == 1

        self.dummy+=1
        #if (self.dummy > 10000): return False
        #if self.verbose: print ('Event : ', event.event)
        #if self.dummy%1000==0: print ("Analyzing events...", self.dummy)
            
    
        ### Get W->jj candidate ###
        gen_ak8 = list(Collection(event, self.genJetBranchName ))
        jets = list(Collection(event, self.jetBranchName ))
        pfCands = list(Collection(event, self.pfCandsBranchName ))
        electrons = list(Collection(event, self.electronBranchName))
        muons = list(Collection(event, self.muonBranchName))
        bjets = list(Collection(event, self.AK4JetBranchName))
        met = Object(event, self.metBranchName)

        #electronTrigger = event.HLT_Ele115_CaloIdVT_GsfTrkIdT
        
        if isMC==True:
            genParticles = Collection(event, self.genCandsAK8BranchName)
            genPartAll = Collection(event, self.genCandsBranchName)

        #list real W's, match gen-level hadronically daughting W's to selected W candidate jets at reco-level
        genWmoms = []
        genWdaughts = []
        
        
        Tmom = [x for x in genPartAll if x.pt>10 and abs(x.pdgId)==6]
        Tdaughts = [x for x in genPartAll if x.pt>1 and (abs(x.pdgId)==5 or abs(x.pdgId)==24)]

        Wdaughts = [x for x in genPartAll if x.pt>1. and 0<abs(x.pdgId)<9]
        
        Wmom = [x for x in genPartAll if x.pt>10. and abs(x.pdgId)==24]
   

        if len(Wdaughts)>0 and len(Wmom)>0 and len(Tdaughts)>0 and len(Tmom)>0 :
            for x in Wdaughts:
                for y in Wmom:
                    try:
                        if y==Wmom[x.genPartIdxMother] and (y in Tdaughts):
                            self.out.fillBranch("PID_W",y.pdgId)
                            #print y, "W from top decaying to qq"
                            #self.out.fillBranch("GPIdx_W",x.genPartIdxMother)
                            genWmoms.append(y)
                            genWdaughts.append(x)
                                        
                    except:
                        continue


        ### Applying selections
        
        # applying basic selections to loose leptons, to decide on veto or not (only using muons for now)

        if not event.HLT_Mu50: return False
        if not event.nMuon > 0: return False
        if not event.nFatJet > 0: return False  
        muonTrigger = event.HLT_Mu50

        recoElectrons  = [x for x in electrons if x.pt>self.minElectronPt and x.cutBased_HEEP and ((self.range1ElectronEta[0]<abs(x.p4().Eta())<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.p4().Eta())<self.range2ElectronEta[1]))]

        recoMuons = [ x for x in muons if x.pt > self.minMuonPt and x.highPtId > 1 and abs(x.p4().Eta()) < self.maxMuonEta and x.pfRelIso03_all < 0.1]  

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt, reverse=True)
        
        if not len(recoMuons) > 0 or not recoMuons[0].pt > self.minTightMuonPt or not abs(recoMuons[0].eta) < self.maxTightMuonEta: return False
        #if not recoMuons[0].highPtId >= 2 or not recoMuons[0].isPFcand or not recoMuons[0].pfIsoId >= 4: return False

        recoLepton=ROOT.TLorentzVector()
        lepflag = -1
        if len(recoMuons)+len(recoElectrons)>=1:
            if len(recoMuons)>=1:                
                if muonTrigger==0: return False
                lepflag=1
                recoLepton = recoMuons[0].p4()

            if len(recoElectrons)>=1:
                return False
        else:
            return False

        if lepflag==1 and met.pt<self.minMuonMET:
            return False
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)
    

        passedMETFilters = False
        try:
          if event.Flag_BadChargedCandidateFilter and event.Flag_BadPFMuonFilter and event.Flag_EcalDeadCellTriggerPrimitiveFilter and event.Flag_HBHENoiseFilter and event.Flag_HBHENoiseIsoFilter and event.Flag_METFilters and event.Flag_ecalBadCalibFilter and event.Flag_globalTightHalo2016Filter and event.Flag_goodVertices:
            passedMETFilters = True
        except:
           passedMETFilters = False


        leptonicWCand = MET+recoLepton
        #if leptonicWCand.Pt()<200: return False
        ### applying selection to fatjets, ak4jets, angular selections between lepton and ak4jet, lepton and fatjets, and ak4 and ak8 jets

        recojets = [ x for x in jets if x.pt_nom > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and x.msoftdrop_nom>self.minJetSDMass]
        recojets.sort(key=lambda x:x.pt_nom,reverse=True)


        if len(recojets)<1: #exit if no AK8jets in event
            return False

        recoAK8 = ROOT.TLorentzVector()
        recoAK8.SetPtEtaPhiM(recojets[0].pt_nom,recojets[0].eta,recojets[0].phi,recojets[0].mass_nom)
    

        dR_AK8Lepton = recoAK8.DeltaR(recoLepton)
        if abs(dR_AK8Lepton)<1. : return False #take care of lepton overlap

        ### applying basic selection to ak4-jets if any, if passed these are b-jet candidates
        recoAK4jets = [ x for x in bjets if x.pt > self.minbJetPt and abs(x.p4().Eta()) < self.maxbJetEta and x.btagCSVV2 > self.minBDisc and recoAK8.DeltaR(x.p4())>1.]
    
        if len(recoAK4jets)<1: return False #exit if no AK4-jets in event
            

        recoAK4 = ROOT.TLorentzVector()
        recoAK4.SetPtEtaPhiM(recoAK4jets[0].pt, recoAK4jets[0].eta,recoAK4jets[0].phi,recoAK4jets[0].mass)
        if abs(recoAK4.DeltaR(recoLepton))<0.3:
            return False
        if abs(recoAK4.DeltaR(recoAK8))<0.8:
                return False
        

        genjets = [x for x in gen_ak8 if x.pt > 0.8*self.minJetPt]
        genjets.sort(key=lambda x:x.pt,reverse=True)
        

        if len(genjets)<1: 
            print "#exit if no gen-jets in event", event.event
            return False #exit if no AK4-jets in event


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

        self.realW = 0
        ###Matching
        if isMC:
            for W in genWmoms:
                genW_4v = ROOT.TLorentzVector()
                genW_4v.SetPtEtaPhiM(W.pt,W.eta,W.phi,W.mass)
                
                if abs(recoAK8.DeltaR(genW_4v))<0.8:

                    ndec = 0
                    for Wdec in genWdaughts:

                        gendec_4v = ROOT.TLorentzVector()
                        gendec_4v.SetPtEtaPhiM(Wdec.pt,Wdec.eta,Wdec.phi,Wdec.mass)
                    
                        if abs(recoAK8.DeltaR(gendec_4v))<0.6:
                            ndec +=1

                    if ndec>1: 
                        self.out.fillBranch("dR_AK8motherW_mW", abs(recoAK8.DeltaR(genW_4v)))
                        self.out.fillBranch("dR_AK8qi_mW", abs(recoAK8.DeltaR(gendec_4v)))
                        self.out.fillBranch("PID_mW", W.pdgId)
                        self.realW = 1
                    else: 
                        self.out.fillBranch("PID_umW", W.pdgId)
                        self.realW = 0 

        self.isgenW=0      
        if isMC:
            for W in genWmoms:
                genW_4v = ROOT.TLorentzVector()
                genW_4v.SetPtEtaPhiM(W.pt,W.eta,W.phi,W.mass)
                
                if abs(genlevelAK8.DeltaR(genW_4v))<0.8:

                    ndec = 0
                    for Wdec in genWdaughts:
                        gendec_4v = ROOT.TLorentzVector()
                        gendec_4v.SetPtEtaPhiM(Wdec.pt,Wdec.eta,Wdec.phi,Wdec.mass)
                    
                        if abs(genlevelAK8.DeltaR(gendec_4v))<0.6:
                            ndec +=1

                    if ndec>1: 
                        #self.out.fillBranch("genmatchedgenAK8", 1)
                        self.isgenW=1
                        self.out.fillBranch("dR_genW_genAK8", abs(genlevelAK8.DeltaR(genW_4v)))
                    else:
                        #self.out.fillBranch("genmatchedgenAK8", 0)
                        self.isgenW=0
                        
                        

        pfCandsVec = ROOT.vector("TLorentzVector")()

        for p in pfCands :
            pfCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        

        for irecojet,recojet in enumerate(recojets):

            # Cluster only the particles near the appropriate jet to save time
            constituents = ROOT.vector("TLorentzVector")()

            for x in pfCandsVec:

                if abs(recojet.p4().DeltaR( x )) < 0.8:
                    constituents.push_back(x)
            nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
            nsub1 = self.nSub1.getTau( self.maxTau, constituents )
            nsub2 = self.nSub2.getTau( self.maxTau, constituents )
                    
            if (irecojet < 1 ): #to extract only the leading jet

                self.out.fillBranch("goodrecojet" + str(irecojet) + "_pt",  recojet.pt_nom )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_eta",  recojet.p4().Eta() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_phi",  recojet.p4().Phi() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_mass",  recojet.p4().M() )
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_softdrop_mass", recojet.msoftdrop_nom)
                self.out.fillBranch("leptonicW_pT",leptonicWCand.Pt())                  
                self.out.fillBranch("lepton_pT",recoLepton.Pt())
                self.out.fillBranch("MET", met.pt)
                self.out.fillBranch("dr_LepJet"  , abs(dR_AK8Lepton))
                self.out.fillBranch("dphi_LepJet", abs(recoAK8.DeltaPhi(recoLepton)))
                self.out.fillBranch("dphi_LepMet", abs(recoLepton.DeltaPhi(MET)))
                self.out.fillBranch("dphi_MetJet", abs(recoAK8.DeltaPhi(MET)))
                self.out.fillBranch("dphi_WJet" , abs(recoAK8.DeltaPhi(leptonicWCand)))
                self.out.fillBranch("genmatchedrecoAK8",  self.realW)   

                if recojet.tau1 > 0.: 
                    tau21 = recojet.tau2/recojet.tau1
                else:
                    tau21 = -1. 
                    
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau21", tau21)
                self.out.fillBranch("goodrecojet" + str(irecojet) + "_N21",  recojet.n2b1)
                self.out.fillBranch("passedMETfilters",passedMETFilters)

                for tauN in range(self.maxTau):
		    if nsub0p5[tauN]==0. or nsub1[tauN]==0. or nsub2[tauN]==0.:
  			self.out.fillBranch("recoEventNo_taus_are_0", event.event) 
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_1_"+str(tauN),  nsub1[tauN]  )
                    self.out.fillBranch("goodrecojet" + str(irecojet) + "_tau_2_"+str(tauN),  nsub2[tauN]  )
        
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
        self.out.fillBranch("genmatchedgenAK8", self.isgenW)
        #self.out.fillBranch("matched_genAK8_recoAK8",  self.match_gen_reco)   

        tau11 = 0
        tau21 = 0 
        for tauN in range(self.maxTau):
            if tauN==0:
                tau11 = nsub1[tauN]
            if tauN==1:
                tau21 = nsub1[tauN]
	    if nsub0p5[tauN]==0. or nsub1[tauN]==0. or nsub2[tauN]==0.:
		self.out.fillBranch("genEventNo_taus_are_0", event.event) 
                print nsub0p5
                print nsub1
                print nsub2

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


import os, sys, glob, gc
#import BTagWeightsFromCSV 
#from BTagWeightsFromCSV import *
#from BranchTools import *
import ROOT
from ROOT import *
from collections import OrderedDict, namedtuple
import operator
import functools
from array import array
import numba
from numba import *

class btaggingUtils():#AddCollectionsModule):
    # a la https://github.com/rappoccio/usercode/blob/bd74c352c7b5ec1e4bf0fe3470d9d05b436c7916/EDSHyFT/test/bTaggingEfficiency/makeBTaggingEfficiencyMapAK5PF.py#L976C5-L977C56
    
    '''
    Class to help in producing efficiency maps for b-tagging AK4 jets in boosted W/top-enriched selections in ttbar; produces the efficiency maps for b,c,udsg AK4s that are b-tagged (in the leptonic hemisphere) in events passing the selections. 
    For events passing the W/top selection , it is assumed that the leading AK4 in the leptonic hemisphere is the b-candidate of interest if the AK4 passes the medium WP of DeepJet. The events that pass the overall selection without an explicit b-tag on the leading leptonic hemisphere AK4 are used to define the denominators for the eff. maps.
    '''
    
    def __init__(self, effMap_btagged, effMap_denom, 
                 discCut=0.2783, events=None,
                 minAK4pt=30., AK4jetIDCut=2, AK4jetMaxEta=2.4, AK4hadronFlavour=5,
                 selAK8pt=200., selAK8mass_min=65., selAK8mass_max=300.,
                 selMETpt=30., selLeptWpt=100., selMupt=55., selLeptWmass=125.,
                 outBinsX=array('d',[]), outBinsY=array('d',[]), 
                 withDR=False, selDRMax=1.6, selDRMin=0.4,
                 verbose=False, verbosityN=1000, nAK8s_min=1, 
                 massType="msoftdrop_nom",withMETXYcorr = True
                ): 
        self.withMETXYcorr=withMETXYcorr
        self.massType = massType
        self.events=events 
        self.effMap_btagged=effMap_btagged
        self.effMap_denom=effMap_denom 


        self.discriminatorCut=discCut
        self.selAK8pt=selAK8pt 
        self.selAK8mass_min=selAK8mass_min 
        self.selAK8mass_max=selAK8mass_max
        self.nAK8s_min = nAK8s_min
        self.selMETpt=selMETpt
        self.selLeptWpt=selLeptWpt
        self.selLeptWmass=selLeptWmass
        self.selMupt=selMupt
        self.outBinsX=outBinsX
        self.outBinsY=outBinsY
        
        if len(self.outBinsX)==0 or len(self.outBinsY)==0:
            print("WARNING, output bin arrays not set, using input eff. TH2's binnings as default")
            self.outBinsX=[]
            self.outBinsY=[]
            for i in range(1,effMap_btagged.GetNbinsX()+2):
                self.outBinsX.append(effMap_btagged.GetXaxis().GetBinLowEdge(i))
            self.outBinsX=array('d',self.outBinsX)
            
            for i in range(1,effMap_btagged.GetNbinsY()+2):
                self.outBinsY.append(effMap_btagged.GetYaxis().GetBinLowEdge(i))
            self.outBinsY=array('d',self.outBinsY)
        self.verbose=verbose
        self.minAK4pt=minAK4pt
        self.AK4jetIDCut=AK4jetIDCut
        self.AK4jetMaxEta=AK4jetMaxEta 
        self.AK4hadronFlavour=AK4hadronFlavour
        self.verbosityN=verbosityN
        
        ROOT.TH1.SetDefaultSumw2()
        ROOT.TH2.SetDefaultSumw2()
        ROOT.TH1.StatOverflows()
        ROOT.TH2.StatOverflows()
        self.denominatorOut = TH2D('denominator_b' , '; AK4 p_T; AK4 #eta', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)
        self.numeratorOut = TH2D('numerator_b' , '; AK4 p_T; AK4 #eta', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)
        self.efficiencyOut = TH2D('efficiency_b' , '; AK4 p_T; AK4 #eta', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)
        
        self.pTmin=0.
        self.pTmax=3000.
        self.pTnbins=300
        self.etamin=-2.4
        self.etamax=2.4
        self.etanbins=48
        
        self.AK8JetMassHisto = TH1D('mass_measurementAK8', ' ; AK8 jet m_{inv.} [GeV/c^{2}]', 60, 20., 320.)#
        self.AK8JetPtHisto = TH1D('pt_measurementAK8', ' ; AK8 jet p_{T} [GeV/c^{2}]', 300, 170., 3170.)#
        self.DeltaRHisto = TH1D('DeltaR_AK4b_mu', ' ; #Delta R(#mu,AK4 lept. hem b)', 20, 0., 2.0)#
        self.AK4LeptHemBMassHisto = TH1D('mass_AK4b', ' ; AK4 lept. hem. b-cand. mass', 60, 20., 320.)#
        self.AK4LeptHemBPtHisto = TH1D('pt_AK4b', ' ; AK4 lept. hem. b-cand. p_{T}', 125, 20., 2520.)#
        
        self.denominatorOut.Sumw2()
        self.numeratorOut.Sumw2()
        self.efficiencyOut.Sumw2()
        self.withDR=withDR
        self.selDRMax=selDRMax
        self.selDRMin=selDRMin
            
    def getEffNumAndDenom(self):
        '''
        Produces TH2s binned in pT vs. eta of numerator and denominator for final b-tagging efficiency map.
        Respectively, these correspond to the number of b-tagged AK4 jets in the selection that are hadron flavour matched, 
        and all AK4's in the selection that are hadron flavour (ghost) matched as originating from B's. Maps for efficiency 
        of other flavours are also derived thusly.

        '''

        c=0
        nEntries = self.events.GetEntries()
        n=0

        for event in self.events:

            n+=1
            
            if not(event.passRecoSel_nom==1):
                continue

            bjets_truth=[]
            bjets_truthAndTagged=[]
            nBtagged = 0
            nAK4WithHadronFlavour = 0 
            #print(event.event)
            if event.selRecoJets_nom_pt>=self.selAK8pt and event.selRecoJets_nom_mass>self.selAK8mass_min and event.selRecoJets_nom_mass<=self.selAK8mass_max and event.selRecoMET_nom_pt>=self.selMETpt and event.selRecoLeptW_nom_pt>=self.selLeptWpt and event.selRecoMu_nom_pt>=self.selMupt and event.selRecoAK4bjetLeptHem_nom_pt>=self.minAK4pt and event.selRecoMu_nom_tkRelIso<=0.05: #and event.nselRecoJets_nom>=self.nAK8s_min event.selRecoLeptW_nom_mass<=self.selLeptWmass#/getattr(event,f"selRecoJets_nom_{massType}")_corr_PUPPI
                
                AK8vect = makeVect(event.selRecoJets_nom_pt,event.selRecoJets_nom_eta,event.selRecoJets_nom_phi,event.selRecoJets_nom_mass)
                if AK8vect==makeVect(0.,0.,0.,0):
                    print('No selected AK8')
                    continue 
                
                Muvect = makeVect(event.selRecoMu_nom_pt,event.selRecoMu_nom_eta,event.selRecoMu_nom_phi,event.selRecoMu_nom_mass)
                if Muvect==makeVect(0.,0.,0.,0):
                    print('No selected AK8')
                    continue 
          
                #store necessary variables/vectors to avoid second loop over jets
                leptHemBCandVect = None
                leptHemBCandFlav = -1
                leptHemBCand_btagScore = 0
                DeltaR = -1
                DeltaR_AK8 = -1

                #loop not necessary anymore since we require exactly 1 btag, on the closest lept hem AK4 to the muon; keeping in for future generalisation if necessary
                #for i in range(len(event.Jet_pt_nom)):


                #if not(event.selRecoAK4bjetLeptHem_nom_pt==event.Jet_pt_nom[i]): 
                #    #work only with the jet that matches the chosen lept. hem. b-cand AK4
                #    continue

                #leptHemBCandVect = makeVect(event.Jet_pt_nom[i],event.Jet_eta[i],event.Jet_phi[i],event.Jet_mass_nom[i])
                leptHemBCandVect = makeVect(event.selRecoAK4bjetLeptHem_nom_pt,event.selRecoAK4bjetLeptHem_nom_eta,event.selRecoAK4bjetLeptHem_nom_phi,event.selRecoAK4bjetLeptHem_nom_mass)
                
                DeltaR=DrRapPhi(Muvect,leptHemBCandVect)
                DeltaR_AK8=DrRapPhi(AK8vect,leptHemBCandVect)

                leptHemBCandFlav = event.selRecoAK4bjetLeptHem_nom_hadronFlavour
                leptHemBCand_btagScore = event.selRecoAK4bjetLeptHem_nom_btagDeepFlavB

                #if event.Jet_pt_nom[i]>=self.minAK4pt and abs(event.Jet_eta[i])<=self.AK4jetMaxEta and event.Jet_jetId[i]>=self.AK4jetIDCut and (DeltaR<=self.selDRMax and DeltaR>=self.selDRMin) and DeltaR_AK8>=1.2: 
                if event.selRecoAK4bjetLeptHem_nom_pt>=self.minAK4pt and abs(event.selRecoAK4bjetLeptHem_nom_eta)<=self.AK4jetMaxEta and event.selRecoAK4bjetLeptHem_nom_jetId>=self.AK4jetIDCut and (DeltaR<=self.selDRMax and DeltaR>=self.selDRMin) and DeltaR_AK8>=1.2: 
                    
                    #if event.selRecoAK4bjetLeptHem_nom_hadronFlavour==self.AK4hadronFlavour:
                    if leptHemBCandFlav==self.AK4hadronFlavour:
                        nAK4WithHadronFlavour+=1 #count number of truth matches
                        #if event.selRecoAK4bjetLeptHem_nom_btagDeepFlavB>=self.discriminatorCut:
                        if leptHemBCand_btagScore>=self.discriminatorCut:
                            nBtagged +=1 #count number of truth matches that are b-tagged
                
                if nAK4WithHadronFlavour==0: continue
                    
                #for i in range(len(event.Jet_pt_nom)):
                     
                #if not(event.selRecoAK4bjetLeptHem_nom_pt==event.Jet_pt_nom[i]): 
                #    continue

                #leptHemBCandVect = makeVect(event.Jet_pt_nom[i],event.Jet_eta[i],event.Jet_phi[i],event.Jet_mass_nom[i])
                #leptHemBCandVect = makeVect(event.selRecoAK4bjetLeptHem_nom_pt,event.selRecoAK4bjetLeptHem_nom_eta,event.selRecoAK4bjetLeptHem_nom_phi,event.selRecoAK4bjetLeptHem_nom_mass_nom)
                #DeltaR=DrRapPhi(Muvect,leptHemBCandVect)
                #DeltaR_AK8=DrRapPhi(AK8vect,leptHemBCandVect)

                #leptHemBCandFlav = event.Jet_hadronFlavour[i]
                #leptHemBCand_btagScore = event.Jet_btagDeepFlavB[i]


                self.AK8JetMassHisto.Fill(event.selRecoJets_nom_mass)
                self.AK8JetPtHisto.Fill(event.selRecoJets_nom_pt)


                #if event.Jet_pt_nom[i]>=self.minAK4pt and abs(event.Jet_eta[i])<=self.AK4jetMaxEta and event.Jet_jetId[i]>=self.AK4jetIDCut and (DeltaR<=self.selDRMax and DeltaR>self.selDRMin) and DeltaR_AK8>=1.2:
                    #print((event.Jet_pt_nom[i],event.Jet_eta[i],event.Jet_phi[i],event.Jet_mass_nom[i]), DrRapPhi(AK8vect,leptHemBCandVect), event.Jet_btagDeepFlavB[i], event.Jet_hadronFlavour[i], )
                if event.selRecoAK4bjetLeptHem_nom_pt>=self.minAK4pt and abs(event.selRecoAK4bjetLeptHem_nom_eta)<=self.AK4jetMaxEta and event.selRecoAK4bjetLeptHem_nom_jetId>=self.AK4jetIDCut and (DeltaR<=self.selDRMax and DeltaR>=self.selDRMin) and DeltaR_AK8>=1.2: 
                   
                    if leptHemBCandFlav==self.AK4hadronFlavour:# and nBtagged<=1:

                        self.DeltaRHisto.Fill(DeltaR)
                        #print(event.event,i)

                        self.effMap_denom.Fill(leptHemBCandVect.Pt(), leptHemBCandVect.Eta())
                        #self.effMap3DDen.Fill(event.Jet_pt_nom[i], event.Jet_eta[i], DeltaR)
                        bjets_truth.append(leptHemBCandVect)

                        #nAK4WithHadronFlavour+=1

                        if leptHemBCand_btagScore>=self.discriminatorCut:

                            self.AK4LeptHemBMassHisto.Fill(leptHemBCandVect.M())
                            self.AK4LeptHemBPtHisto.Fill(leptHemBCandVect.Pt())
                            #print(event.event,i)

                            self.effMap_btagged.Fill(leptHemBCandVect.Pt(), leptHemBCandVect.Eta())
                            #self.effMap3DNum.Fill(event.Jet_pt_nom[i], event.Jet_eta[i], DeltaR)
                            bjets_truthAndTagged.append(leptHemBCandVect)


            else:
                continue


            c+=1
            if self.verbose and n%self.verbosityN==0: 
                print(f"Processed {n}/{nEntries} events, i.e., {n/nEntries *100.}% of entries in tree ")
                print(f"Number of events that pass selection so far = {c}", '\n')
                    
        return 1

    def getBtaggingEfficiencyMap(self):
        '''
        Produces final 2D efficiency map for b-tagging in the provided event sub-selection (boosted W/top enriched as per settings for AK8 jet selections)
        '''
        
        #updateNumDen = self.getEffNumAndDenom() #if not self.withDR else self.getEffNumAndDenomWithDR()

        xShift = self.effMap_denom.GetXaxis().GetBinWidth(1)/2.
        yShift = self.effMap_denom.GetYaxis().GetBinWidth(1)/2.

        #self.denominatorOut = self.denOutTH2D('denominator_b' , '; AK4 #eta; AK4 p_T', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)
        #self.numeratorOut   = self.numOutTH2D('numerator_b' , '; AK4 #eta; AK4 p_T', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)
        #self.efficiencyOut  = self.effOutTH2D('efficiency_b' , '; AK4 #eta; AK4 p_T', (len(self.outBinsX)-1), self.outBinsX, (len(self.outBinsY)-1), self.outBinsY)

        for i in range(1, self.denominatorOut.GetXaxis().GetNbins()+1):
            for j in range(1, self.denominatorOut.GetYaxis().GetNbins()+1):

                binXMin = self.effMap_denom.GetXaxis().FindBin(self.denominatorOut.GetXaxis().GetBinLowEdge(i)+xShift)
                binXMax = self.effMap_denom.GetXaxis().FindBin(self.denominatorOut.GetXaxis().GetBinUpEdge(i)-xShift)
                binYMinPos = self.effMap_denom.GetYaxis().FindBin(self.denominatorOut.GetYaxis().GetBinLowEdge(j)+yShift)
                binYMaxPos = self.effMap_denom.GetYaxis().FindBin(self.denominatorOut.GetYaxis().GetBinUpEdge(j)-yShift)
                binYMinNeg = self.effMap_denom.GetYaxis().FindBin(-self.denominatorOut.GetYaxis().GetBinUpEdge(j)+yShift)
                binYMaxNeg = self.effMap_denom.GetYaxis().FindBin(-self.denominatorOut.GetYaxis().GetBinLowEdge(j)-yShift)

                denominator = self.effMap_denom.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
                denominator = denominator + self.effMap_denom.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)
                numerator = self.effMap_btagged.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
                numerator = numerator + self.effMap_btagged.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)

                if(i==self.denominatorOut.GetXaxis().GetNbins()): # also add overflow to the last bin in jet pT
                    denominator = denominator + self.effMap_denom.Integral(binXMax+1,self.effMap_denom.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
                    denominator = denominator + self.effMap_denom.Integral(binXMax+1,self.effMap_denom.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)
                    numerator = numerator + self.effMap_btagged.Integral(binXMax+1,self.effMap_btagged.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
                    numerator = numerator + self.effMap_btagged.Integral(binXMax+1,self.effMap_btagged.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)

                self.denominatorOut.SetBinContent(i,j,denominator)
                self.denominatorOut.SetBinError(i,j,ROOT.TMath.Sqrt(denominator))
                self.numeratorOut.SetBinContent(i,j,numerator)
                self.numeratorOut.SetBinError(i,j,ROOT.TMath.Sqrt(numerator))
                #if(denominator>0.): 
                eff = numerator/(1.*denominator) if denominator>0. else 0.
                err = ROOT.TMath.Sqrt(eff * (1. - eff) / denominator) if denominator > 0 else 0.
                self.efficiencyOut.SetBinContent(i,j,eff)
                self.efficiencyOut.SetBinError(i,j,err)

            ## check if there are any bins with 0 or 100% efficiency
            for i in range(1,self.denominatorOut.GetXaxis().GetNbins()+1):
                for j in range(1,self.denominatorOut.GetYaxis().GetNbins()+1):
            
                    efficiency = self.efficiencyOut.GetBinContent(i,j)
                    if(efficiency==0. or efficiency==1.):
                        print ('Warning! Bin(%i,%i) for %s jets has a b-tagging efficiency of %.3f'%(i,j,'b',efficiency))

            # set efficiencies in overflow bins
            for i in range(1,self.denominatorOut.GetXaxis().GetNbins()+1):
                self.efficiencyOut.SetBinContent(i, self.denominatorOut.GetYaxis().GetNbins()+1, self.efficiencyOut.GetBinContent(i, self.denominatorOut.GetYaxis().GetNbins()))

            for j in range(1,self.denominatorOut.GetYaxis().GetNbins()+2):
                self.efficiencyOut.SetBinContent(self.denominatorOut.GetXaxis().GetNbins()+1, j, self.efficiencyOut.GetBinContent(self.denominatorOut.GetXaxis().GetNbins(), j))
                
        return 1.#self.efficiencyOut, self.numeratorOut, self.denominatorOut

    def getBtaggingEfficiencyMap3D(self):
        '''
        Produces a 3D efficiency map for b-tagging as a function of AK4 jet pT, eta, and deltaR with AK8.
        '''
        
        xShift = self.effMap_denom.GetXaxis().GetBinWidth(1)/2.
        yShift = self.effMap_denom.GetYaxis().GetBinWidth(1)/2.
        zShift = self.effMap_denom.GetZaxis().GetBinWidth(1)/2.

        for i in range(1, self.efficiencyMap3DOut.GetNbinsX()+1):
            for j in range(1, self.efficiencyMap3DOut.GetNbinsY()+1):
                for k in range(1, self.efficiencyMap3DOut.GetNbinsZ()+1):
                    
                    # Find corresponding bin ranges in the input histograms
                    binXMin = self.effMap3DDen.GetXaxis().FindBin(self.efficiencyMap3DOut.GetXaxis().GetBinLowEdge(i)+xShift)
                    binXMax = self.effMap3DDen.GetXaxis().FindBin(self.efficiencyMap3DOut.GetXaxis().GetBinUpEdge(i)-xShift)
                    
                    binYMinPos = self.effMap3DDen.GetYaxis().FindBin(self.efficiencyMap3DOut.GetYaxis().GetBinLowEdge(j)+yShift)
                    binYMaxPos = self.effMap3DDen.GetYaxis().FindBin(self.efficiencyMap3DOut.GetYaxis().GetBinUpEdge(j)-yShift)
                    binYMinNeg = self.effMap3DDen.GetYaxis().FindBin(-self.efficiencyMap3DOut.GetYaxis().GetBinUpEdge(j)+yShift)
                    binYMaxNeg = self.effMap3DDen.GetYaxis().FindBin(-self.efficiencyMap3DOut.GetYaxis().GetBinLowEdge(j)-yShift)

                    binZMin = self.effMap3DDen.GetZaxis().FindBin(self.efficiencyMap3DOut.GetZaxis().GetBinLowEdge(k)+yShift)
                    binZMax = self.effMap3DDen.GetZaxis().FindBin(self.efficiencyMap3DOut.GetZaxis().GetBinUpEdge(k)-yShift)

                    n_btagged_pos = self.effMap3DNum.Integral(binXMin, binXMax, binYMinPos, binYMaxPos, binZMin, binZMax)
                    n_total_pos = self.effMap3DDen.Integral(binXMin, binXMax, binYMinPos, binYMaxPos, binZMin, binZMax)

                    n_btagged_neg = self.effMap3DNum.Integral(binXMin, binXMax, binYMinNeg, binYMaxNeg, binZMin, binZMax)
                    n_total_neg = self.effMap3DDen.Integral(binXMin, binXMax, binYMinNeg, binYMaxNeg, binZMin, binZMax)

                    n_btagged = n_btagged_pos + n_btagged_neg
                    n_total = n_total_pos + n_total_neg
                    
                    if(i==self.effMap3DDen.GetXaxis().GetNbins()): # also add overflow to the last bin in jet pT
                        
                        
                        n_total = n_total + self.effMap3DDen.Integral(binXMax+1, self.effMap3DDen.GetXaxis().GetNbins()+1, 
                                                                      binYMinPos, binYMaxPos, 
                                                                      binZMin, binZMax) 
                        n_total = n_total + self.effMap3DDen.Integral(binXMax+1,self.effMap3DDen.GetXaxis().GetNbins()+1,
                                                                       binYMinNeg,binYMaxNeg, 
                                                                       binZMin, binZMax) 
                        
                        n_btagged = n_btagged + self.effMap3DNum.Integral(binXMax+1, self.effMap3DNum.GetXaxis().GetNbins()+1, 
                                                                          binYMinPos, binYMaxPos, 
                                                                          binZMin, binZMax) 
                        n_btagged = n_btagged + self.effMap3DNum.Integral(binXMax+1,self.effMap3DNum.GetXaxis().GetNbins()+1,
                                                                          binYMinNeg,binYMaxNeg, 
                                                                          binZMin, binZMax)
                        
                    # Calculate efficiency and error
                    efficiency = n_btagged /(1.*n_total) if n_total > 0 else 0.
                    error = ROOT.TMath.Sqrt(efficiency * (1. - efficiency) / n_total) if n_total > 0 else 0.

                    # Set the efficiency and error for the output bin
                    self.efficiencyMap3DOut.SetBinContent(i, j, k, efficiency)
                    self.efficiencyMap3DOut.SetBinError(i, j, k, error)

        # Handle overflows
        self.set3DOverflow(self.efficiencyMap3DOut)

        return self.efficiencyMap3DOut

    def set3DOverflow(self, hist3D):
        '''
        Set overflow bins for a 3D histogram.
        '''
        #for i in range(1, hist3D.GetNbinsX()+1):
        for j in range(1, hist3D.GetNbinsY()+2):
            for k in range(1, hist3D.GetNbinsZ()+2):
                #if i == hist3D.GetNbinsX():  # Overflow for X
                hist3D.SetBinContent(hist3D.GetNbinsX()+1, j, k, hist3D.GetBinContent(hist3D.GetNbinsX(), j, k))
                hist3D.SetBinError(hist3D.GetNbinsX()+1, j, k, hist3D.GetBinError(hist3D.GetNbinsX(), j, k))
                
                
def DrRapPhi(va,vb ):
    dy = va.Rapidity()-vb.Rapidity()
    dphi = va.DeltaPhi(vb)
    return ROOT.TMath.Sqrt(dy*dy+dphi*dphi)


def makeVect(pt,eta,phi,m):
    vect = ROOT.TLorentzVector()
    vect.SetPtEtaPhiM(pt,eta,phi,m)
    return vect

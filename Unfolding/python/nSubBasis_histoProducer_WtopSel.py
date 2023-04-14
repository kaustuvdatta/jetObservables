'''
import uproot
from uproot import *
import coffea
'''


from coffea import processor
import coffea.processor
from coffea.processor import defaultdict_accumulator,dict_accumulator#,IterativeExecutor,FuturesExecutor
import time
import re
from datasets_WtopSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
from collections import OrderedDict
import numpy as np
import awkward as ak
import uproot
import hist
from hist import Hist

from collections import defaultdict
import numba


class nSubBasis_unfoldingHistoProd_WtopSel(processor.ProcessorABC):
    
    def __init__(self, sampleName, selection='_topSel', withLeptHemBtag=False, sysSource=[],year='2017', era='', 
                 isMC=True, isSigMC=True, onlyUnc='', wtUnc=False, verbose=True,
                 sampleDict=dictSamples, test=False, sysUnc=False):

        self.test=test
        self.year = year
        self.isMC = isMC
        self.isSigMC = isSigMC
        self.era = era
        self.withLeptHemBtag = withLeptHemBtag
        
        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc       
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        self.verbose=verbose

        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with while handling data' )
        
        
        ### Helpers
        self.listOfUnfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 

        self.listOfControlHistTypes =  [ 'gen',  'reco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        
        self.nJet = [ 'Jet']
        self.selObjects = ['Mu','LeptW','MET','AK4bjetHadHem']

        self.kinematic_labels=['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels=['good_nPVs','nRecoBtags','nRecoHadBtags','nRecoLepBtags', 'selRecoHadHemDeltaR']
        self.gen_only_labels=['nGenBtags','nGenHadBtags','nGenLepBtags', 'selGenHadHemDeltaR']

        self.dict_variables_kinematics_AK8 = {

                                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),
                                            "_msoftdrop": np.array([i for i in np.arange(0., 300, 5.)]),

                                        } 

        self.dict_variables_kinematics_Muon = {

                                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),

                                        } 

        self.dict_variables_kinematics_LeptW = {

                                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),

                                        } 

        self.dict_variables_kinematics_AK4Had = {

                                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),

                                        }

        self.dict_variables_kinematics_MET = {

                                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),

                                        } 

        self.dict_variables_reco = { 
                                        "good_nPVs": np.array([i for i in np.arange(0., 100, 1.)]),
                                        "nRecoBtags": np.array([i for i in np.arange(0., 4, 1.)]),
                                        "nRecoHadBtags": np.array([i for i in np.arange(0., 3, 1.)]),
                                        "nRecoLepBtags": np.array([i for i in np.arange(0., 3, 1.)]),
                                        "selRecoHadHemDeltaR": np.array([i for i in np.arange(0., 1., 0.05)]),

                                   }   
        
        self.dict_variables_gen = { 
                                        "nGenBtags": np.array([i for i in np.arange(0., 4., 1.)]),
                                        "nGenHadBtags": np.array([i for i in np.arange(0., 3., 1.)]),
                                        "nGenLepBtags": np.array([i for i in np.arange(0., 3., 1.)]),
                                        "selGenHadHemDeltaR": np.array([i for i in np.arange(0., 1., 0.05)]),

                                   }   
        self.dict_variables_toUnfold = {
                    
                    "_tau_0p5_1": np.array([(i/1000) for i in np.arange(0.*1000, 1.*1001)]),
                    "_tau_0p5_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.9*1001)]),
                    "_tau_0p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.8*1001)]),
                    "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.7*1001)]),
                    "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.7*1001)]),
                    "_tau_1_1": np.array([(i/1000) for i in np.arange(0.*1000, 0.9*1001)]),
                    "_tau_1_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.6*1001)]),
                    "_tau_1_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.4*1001)]),
                    "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1001)]),
                    "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1001)]),
                    "_tau_2_1": np.array([(i/1000) for i in np.arange(0.*1000, 0.5*1001)]),
                    "_tau_2_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1001)]),
                    "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.2*1001)]),
                    "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.2*1001)]),
                    "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.2*1001)]),
                    "_tau21": np.array([(i/1000) for i in np.arange(0.*1000, 1.6*1001)]),#for one-pass kT minimization as per CMS
                    "_tau32": np.array([(i/1000) for i in np.arange(0.*1000, 1.6*1001)]),#for one-pass kT minimization as per CMS
                    "_tau21_WTA": np.array([(i/1000) for i in np.arange(0.*1000, 1.4*1001)]),#for WTA-kT for comparison
                    "_tau32_WTA": np.array([(i/1000) for i in np.arange(0.*1000, 1.4*1001)]),#for WTA-kT for comparison
                    "_tau21_exkT": np.array([(i/1000) for i in np.arange(0.*1000, 2.*1001)]),#for excl.-kT and E-scheme as per basis
                    "_tau32_exkT": np.array([(i/1000) for i in np.arange(0.*1000, 2.*1001)]),#for excl.-kT and E-scheme as per basis
                    "_mass": np.array([(i/2.0) for i in np.arange(0.*2, 300*2.1)]),
                    "_msoftdrop": np.array([(i/2.0) for i in np.arange(0.*2, 200*2.1)]),
                    "_pt": np.array([(i/2.0) for i in np.arange(170.*2, 2500.*2.1)]),
                        }
        
        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not( isys.endswith('nom') or self.sysUnc) ]
        self.sysWeightList = ( '_pu', '_isr', '_pdf', '_fsr', '_l1prefiring', '_lepton', '_btag' ) #'_ps',
        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_btagWeight', '_l1prefiringWeight', '_leptonWeightAll']#, '_leptonWeightISO', '_leptonWeightID', '_leptonWeightTrig', '_leptonWeightRecoEff'] if self.wtUnc else [] 
        if self.onlyUnc: self.sysSource = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] #################### Not using for right now
        
        #puWeights used only to change reco event weight (up/down) without application to gen weight, others applied to modify the genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)]
        
        self.selList = [ selection ] 
        
                
        
        self._branchesToRead = self.getBranchesToRead(
                                                     dirname=self.inputDir, year=self.year,
                                                     kinematic_labels=self.kinematic_labels,
                                                     reco_only_labels=self.reco_only_labels,
                                                     nSub_labels=self.dict_variables_toUnfold,
                                                    )# isMC=True, isSigMC=False, wtUnc=False, sysUnc=False,
        #self._accumulator = self.buildDictOfHistograms()#processor.dict_accumulator()
        

        self.top_ptmin = 400.
        self.top_mSDmin = 140.
        self.top_mSDmax = 250.
        self.W_ptmin = 200.
        self.W_mSDmin = 60.
        self.W_mSDmax = 120.

    #@property
    #def accumulator(self):
    #    return self._accumulator    
    #def fill_KinematicHists(self, events, mask, recORgen='reco', obj_name='Mu', dict_kin_labels=self.dict_variables_kinematics_Muon):
    #    for i in dict_kin_labels:
    #        (events[f'selRecoJets{s}{varToFill}'][selRecoMask],weight=totalRecoWeight[selRecoMask],
    #                                                     threads=8)



    def process(self, events):
        '''fill Hist histograms to accumulate event stats, convert to root or whatever else after returned by processor'''

        output = self.buildDictOfHistograms()
        branches = self._branchesToRead
        events = events[branches]
        
        if self.verbose: print (f'systematic sources being considered:{self.sysSource} for nevents={len(events)}')
        for sys in self.sysSource:
            
            s='_nom' if (sys.startswith(self.sysWeightList)) else sys
            
            if self.verbose: 
                print('making individual W and topSel mask',sys,s)
                
            #if self.withLeptHemBtag==False:

            topRecoMask = (events[f'selRecoJets{s}_pt']>self.top_ptmin) & (events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop_corr_PUPPI']>self.top_mSDmin) & (events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop_corr_PUPPI']<self.top_mSDmax) & (events[f'selRecoHadHemDeltaR{s}']<0.8) & (events[f'totalRecoWeight{s}']!=0.)

            WRecoMask = (events[f'selRecoJets{s}_pt']>self.W_ptmin) & (events[f'selRecoJets{s}_mass']>self.W_mSDmin) & (events[f'selRecoJets{s}_mass']<self.W_mSDmax) & (events[f'selRecoHadHemDeltaR{s}']>0.8) & (events[f'selRecoHadHemDeltaR{s}']<1.6) & (events[f'totalRecoWeight{s}']!=0.)

            if self.isMC or self.isSigMC:                       
                topGenMask = (events[f'selGenJets{s}_pt']>self.top_ptmin) & (events[f'selGenJets{s}_msoftdrop']>self.top_mSDmin) & (events[f'selGenJets{s}_msoftdrop']<self.top_mSDmax) & (events[f'selGenHadHemDeltaR{s}']<0.8) & (events[f'evtGenWeight_nom'] !=0.)

            WGenMask = (events[f'selGenJets{s}_pt']>self.W_ptmin) & (events[f'selGenJets{s}_mass']>self.W_mSDmin) & (events[f'selGenJets{s}_mass']<self.W_mSDmax) & (events[f'selGenHadHemDeltaR{s}']>0.8) & (events[f'selGenHadHemDeltaR{s}']<1.6) & (events[f'evtGenWeight_nom'] !=0.)

            #else:

            #    topRecoMask = (events[f'selRecoJets{s}_pt']>self.top_ptmin) & (events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop_corr_PUPPI']>self.top_mSDmin) & (events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop']/events[f'selRecoJets{s}_msoftdrop_corr_PUPPI']<self.top_mSDmax) & (events[f'selRecoHadHemDeltaR{s}']<0.8) & (events[f'FlagRecoLeptHemBjet{s}']==1) & (events[f'totalRecoWeight{s}']!=0.)

            #    WRecoMask = (events[f'selRecoJets{s}_pt']>self.W_ptmin) & (events[f'selRecoJets{s}_mass']>self.W_mSDmin) & (events[f'selRecoJets{s}_mass']<self.W_mSDmax) & (events[f'selRecoHadHemDeltaR{s}']>0.8) & (events[f'selRecoHadHemDeltaR{s}']<1.6) & (events[f'FlagRecoLeptHemBjet{s}']==1) & (events[f'totalRecoWeight{s}']!=0.)

            #    if self.isMC or self.isSigMC: 

            #        topGenMask = (events[f'selGenJets{s}_pt']>self.top_ptmin) & (events[f'selGenJets{s}_msoftdrop']>self.top_mSDmin) & (events[f'selGenJets{s}_msoftdrop']<self.top_mSDmax) & (events[f'selGenHadHemDeltaR{s}']<0.8) & (events[f'FlagGenLeptHemBjet{s}']==1) & (events[f'evtGenWeight_nom'] !=0.)

            #        WGenMask = (events[f'selGenJets{s}_pt']>self.W_ptmin) & (events[f'selGenJets{s}_mass']>self.W_mSDmin) & (events[f'selGenJets{s}_mass']<self.W_mSDmax) & (events[f'selGenHadHemDeltaR{s}']>0.8) & (events[f'selGenHadHemDeltaR{s}']<1.6) & (events[f'FlagGenLeptHemBjet{s}']==1) & (events[f'evtGenWeight_nom'] !=0.)

            if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) & (~(topRecoMask|topGenMask))
                    if self.isMC or self.isSigMC:                        
                        WtopGenMask = (WGenMask) & (~(topRecoMask|topGenMask))

            elif '_topSel' in self.selList: 
                WtopRecoMask = (topRecoMask) &(~(WRecoMask|WGenMask))
                if self.isMC or self.isSigMC:                        
                    WtopGenMask = (topGenMask) &(~(WGenMask|WRecoMask))    

            
            
            if self.isMC and self.isSigMC:# or self.isaltSigMC:
                #event_mask = ((WtopRecoMask) | (WtopGenMask))
                #events = events[event_mask]
                
                #s='_nom' if (sys.startswith(self.sysWeightList)) else sys
                        
                selRecoMask = (events[f'totalRecoWeight{s}']!=0.) & (WtopRecoMask) #& (event_mask)
                selGenMask = (events[f'evtGenWeight_nom'] !=0.) & (WtopGenMask) #& (event_mask)

                trueRecoMask = (events[f'trueRecoJets{s}_pt']>0.) & (selRecoMask) & (selGenMask) 
                fakeRecoMask = (events[f'trueRecoJets{s}_pt']<0.) & (selRecoMask) & (~selGenMask) #((selRecoMask) ^ (trueRecoMask))#
                accepGenMask = (events[f'accepGenJets{s}_pt']>0.) & (selGenMask) & (selRecoMask)
                missGenMask =  (events[f'accepGenJets{s}_pt']<0.) & (selGenMask) & (~selRecoMask)#((selGenMask) ^ (accepGenMask))  #
                
                selRecoMask = (trueRecoMask) | (fakeRecoMask) #& (event_mask)
                selGenMask = (accepGenMask) | (missGenMask) #& (event_mask)

            elif self.isMC and (not self.isSigMC): #not so relevant for dijets but for background MC's in W/top
                
                #event_mask = (WtopRecoMask) | (WtopGenMask)
                #events = events[event_mask]
                
                if sys.endswith('nom'):
                    selRecoMask = (events[f'totalRecoWeight{sys}']!=0.) & (WtopRecoMask) #& (event_mask)
                    selGenMask = (events[f'evtGenWeight_nom'] !=0.) & (WtopGenMask) #& (event_mask)

                    trueRecoMask = (events[f'trueRecoJets{sys}_pt']>0.) & (selRecoMask)# & (selGenMask) 
                    fakeRecoMask = (events[f'trueRecoJets{sys}_pt']<0.) & (selRecoMask) #((selRecoMask) ^ (trueRecoMask))#
                    accepGenMask = (events[f'accepGenJets{sys}_pt']>0.) & (selGenMask)# & (selRecoMask)
                    missGenMask =  (events[f'accepGenJets{sys}_pt']<0.) & (selGenMask) #((selGenMask) ^ (accepGenMask))  #

                    selRecoMask = (trueRecoMask) | (fakeRecoMask) #& (event_mask)
                    selGenMask = (accepGenMask) | (missGenMask) #& (event_mask)
                    #selRecoMask = (events[f'totalRecoWeight{sys}']!=0.) & (WtopRecoMask)
                    #selGenMask = (events[f'evtGenWeight{sys}']!=0.) & (WtopGenMask)
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                #event_mask = (WtopRecoMask) 
                #events = events[event_mask]
                selRecoMask = (events[f'totalRecoWeight{s}']!=0.) & (WtopRecoMask)
                
            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            if self.verbose: 
                #print (selRecoMask[0:20],selRecoMask[-20:-1])
                #print (trueRecoMask[0:20],trueRecoMask[-20:-1])
                #print (selGenMask[0:20],selGenMask[-20:-1])
                #print (accepGenMask[0:20],accepGenMask[-20:-1])
                print("Going to fill histos")
                
            for isel in self.selList:
                
                
                listOfJetOutputHistosNames = [k for k,h in output.items() if ((('Jet' in k) and (sys in k)) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                if self.verbose: print (listOfJetOutputHistosNames)
                        
                #listOfMiscOutputHistosNames = [k for k,h in output.items() if ((not ('Jet' in k) and (sys in k)) or ('resol' in k and not(sys in k)))]
                #if self.verbose: print (listOfMiscOutputHistosNames)
                
                for k in listOfJetOutputHistosNames:
                    key=k
                    
                    #if self.verbose: print(key,sys)
                    ############### Safety checks ##################
                    if not(s in key) and not('resol' in key): 
                        continue

                    #if not( self.isMC) and not('AK8PF' in key):
                    #    continue

                    #if not( self.isMC) and 'AK8PF' in key:
                    #    if not( isel.split('_')[1] in key): 
                    #        continue #to ensure that for data, we only fill the histos relevant to a given trigger (stored in isel as '_AKPFJetXYZ_dijetSel') correctly

                    ############### ############ ##################
                    if self.verbose: print(key,sys,s)
                    # Decide on trigger variable labels
    
                    if not('_tau' in key):
                        whichKinVar = [var for var in self.dict_variables_kinematics_AK8.keys() if var[1:] in key]#[0]
                        varToFill = whichKinVar[0]

                    elif ('_tau' in key):# and 'nom' in key):
                        temp='WTA' if 'WTA' in key else 'tau' #hack, :(, to fix tau21_nonOPkT being histogrammed incorrectly
                        temp='exkT' if 'exkT' in key else 'tau'
                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                        varToFill = whichnSub[0]
                        

                    ################## Filling histos from accumulated event arrays ##################

                    # Decide what weights to use
                    if not (sys.startswith(self.sysWeightList)):
                        s=sys
                        totalRecoWeight = events[f'evtGenWeight_nom']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'btagWeightNom{s}']*events[f'leptonWeightNom{s}'] if self.isMC else events[f'totalRecoWeight_nom']
                        totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None
                    else:
                        s='_nom'
                        if 'pu' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'l1prefiringWeightNom{s}']*events[f'btagWeightNom{s}']*events[f'leptonWeightNom{s}']
                            
                        elif 'l1' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'btagWeightNom{s}']*events[f'leptonWeightNom{s}']

                        elif 'btag' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'leptonWeightNom{s}']

                        elif 'lepton' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'btagWeightNom{s}']
                            
                        elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                            totalGenWeight = events[f'evtGenWeight{s}']*events[f'{sys.split("_")[1]}{s}']
                            totalRecoWeight = totalGenWeight*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'btagWeightNom{s}']*events[f'leptonWeightNom{s}']

                    if (key.lower().startswith(('reco','good'))):
                        if 'nPV' in key:
                            #output[key].fill(events[f'good_nPVs{s}'][selRecoMask],weight=totalRecoWeight[selRecoMask],
                            #                 threads=8)
                            dummy=1
                        else: 
                            if self.verbose: print('all reco', key,s,sys, len(events[f'selRecoJets{s}{varToFill}'][selRecoMask]))
                            output[key].fill(events[f'selRecoJets{s}{varToFill}'][selRecoMask],weight=totalRecoWeight[selRecoMask],
                                             threads=8)

                    # Stuff below should have diff wts already available to them so no need to touch anything down here
                    if self.isMC and (not varToFill.startswith(tuple(self.reco_only_labels))):

                        if (key.lower().startswith('true')) and self.isSigMC :
                            if self.verbose: print('true', key,s,sys, len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask]))
                            output[key].fill(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                        elif (key.lower().startswith('fake')) and self.isSigMC:
                            if self.verbose: print('fake', key,s,sys,len(events[f'selRecoJets{s}{varToFill}'][fakeRecoMask]))
                            output[key].fill(events[f'selRecoJets{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                             threads=8)

                        if (key.lower().startswith('genjet')):
                            #events[f'trueRecoJets{s}{varToFill}'][trueRecoMask]
                            if sys.endswith('nom'):
                                output[key].fill(events[f'selGenJets_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)#=self._listofHistograms(histoName) #hist.Hist("Events", hist.Cat("branch", branch), hist.Bin("value", branch, 100, 0, 1000))    

                        elif (key.lower().startswith('accepgenjet')) and self.isSigMC:# and not key.startswith(('accep','miss')):
                            output[key].fill(events[f'accepGenJets{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                             threads=8)

                        elif (key.lower().startswith('missgenjet')) and self.isSigMC:
                            output[key].fill(events[f'selGenJets_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                             threads=8)
                        
                        elif ( 'resp' in key.lower() and not('miss' in key.lower())) and self.isSigMC:
                            
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)     
                        
                        elif ('respwithmiss' in key.lower()) and self.isSigMC:
                            
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{s}{varToFill}'][trueRecoMask], weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)
                            
                            #fill missgen weight
                            output[key].fill(gen=events[f'selGenJets_nom{varToFill}'][missGenMask], reco=-1.*np.ones(len(events[f'selGenJets_nom{varToFill}'][missGenMask])),
                                             weight=totalGenWeight[missGenMask],
                                             threads=8)
                            

                        elif ('resol' in key.lower()) and sys.endswith('nom') and self.isSigMC:
                            zeroMask=(events[f'accepGenJets{s}{varToFill}']!=0.)&(accepGenMask)

                            response = events[f'trueRecoJets{s}{varToFill}'][zeroMask]/events[f'accepGenJets{s}{varToFill}'][zeroMask]
                            response = np.nan_to_num(response,nan=-999.)
                            
                            output[key].fill(response, weight=totalRecoWeight[zeroMask])
                            if 'noWt_' in key: output[key].fill(response)
        if self.verbose: print("Histos filled!")
        l=[]
        for x,y in output.items(): #y.SetDirectory(0)
        
            if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss','true','fake' )) and '_nom' in x:
                l.append(x)
            #elif self.isMC and not(self.isSigMC) and x.startswith(('accepgen','miss','true','fake','resol' )) and '_nom' in x:
            #    l.append(x)
        for k in l:
            del(output[k])
        
        return output

    
    def postprocess(self, accumulator):
        pass
        
        
   
    def buildDictOfHistograms(self):
        '''build dictionary of Hist histograms, convert to root or whatever else after filled and returned by processor'''
        dictOfHists = OrderedDict()

       #self.selList = [ '_WSel' if self.selList.startswith(('_W')) else '_topSel' ]  
        if self.verbose: print(self.selList)

        #build from list of histo names for unfolding histograms

        for isel in self.selList:
            #if self.verbose:    print (self.listOfUnfHistTypes)

            for itype in self.listOfUnfHistTypes:
                iJ=self.nJet[0] 

                for sysUnc in self.sysSource:
                    
                    for x, y in self.dict_variables_kinematics_AK8.items():
                        binning = y
                        if sysUnc.endswith('nom') or self.onlyUnc:
                            if not x in tuple(self.reco_only_labels): 
                                dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())        
                                
                    if itype.startswith('truereco') and self.isMC and sysUnc.endswith('nom') and self.isSigMC:
                        dictOfHists['resol'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_pt'+isel, label='AK8 reco/gen jet pt').Weight())
                        dictOfHists['resol'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD}/gen jet m_{SD}').Weight())
                        dictOfHists['resol'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_mass'+isel, label='AK8 reco inv. m/gen jet inv. m').Weight())
                    
                    #print('building unfolding histos')     
                    for x, y in self.dict_variables_toUnfold.items():
                        binning = y
                        
                        #binning_coarse=np.array([binning[i] for i in range(len(binning)) if i%5==0])
                        
                        dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label='AK8 '+itype+' jet #tau', underflow=True,overflow=True).Weight())
                        
                        if itype.startswith('truereco') and self.isMC and self.isSigMC:

                            dictOfHists['resp'+iJ+x+sysUnc+isel] = (
                                                                    hist.Hist.new
                                                                    .Variable(binning,name='gen', label='AK8 gen jet'+x, underflow=True,overflow=True)
                                                                    .Variable(binning,name='reco',label='AK8 reco jet'+x, underflow=True,overflow=True)
                                                                    .Weight()
                                                                   )
                            
                            dictOfHists['respWithMiss'+iJ+x+sysUnc+isel] = (
                                                                            hist.Hist.new
                                                                            .Variable(binning,name='gen', label='AK8 gen jet'+x, underflow=True,overflow=True)
                                                                            .Variable(binning,name='reco',label='AK8 reco jet'+x, underflow=True,overflow=True)
                                                                            .Weight()
                                                                           )
                            if sysUnc.endswith('_nom') and self.isSigMC:
                                dictOfHists['resol'+iJ+x+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
                                dictOfHists['noWt_resol'+iJ+x+isel] = (hist.Hist.new.Regular(500, 0, 5, name='noWt_resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
        

        #build from list of histos for control plots
        if not(self.onlyUnc or self.sysUnc):
            for isel in self.selList:
                sysUnc= '_nom'

                #if self.verbose:    print (self.listOfControlHistTypes)

                for itype in self.listOfControlHistTypes:
                    for iO in self.selObjects:

                        #making control histos only for nominal case for now
                            
                        if 'muon' in iO.lower():
                            for x, y in self.dict_variables_kinematics_Muon.items():
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} muon {x}').Weight())        
                                
                        elif 'met' in iO.lower():
                            for x, y in self.dict_variables_kinematics_MET.items():
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} MET {x}').Weight())        
                                
                        elif 'leptw' in iO.lower():
                            for x, y in self.dict_variables_kinematics_LeptW.items():
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} leptonic W {x}').Weight())        
                                
                        elif 'ak4' in iO.lower():
                            for x, y in self.dict_variables_kinematics_AK4Had.items():
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} AK4 {x}').Weight())        

                if self.verbose: print (self.dict_variables_reco.keys())

                for itype in self.dict_variables_reco:
                    bins=self.dict_variables_reco[itype]
                    
                    if 'good_nPVs' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'Good nPVs').Weight())  
                    elif 'nRecoBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco').Weight())  
                    elif 'nRecoHadBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco (hadr. hem.)').Weight())  
                    elif 'nRecoLepBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco (lept. hem.)').Weight())  
                    elif 'selRecoHadHemDeltaR' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, AK4)').Weight())  


                if self.isMC or self.isSigMC:
                    if self.verbose: print (self.dict_variables_gen.keys())

                    for itype in self.dict_variables_gen:
                        bins=self.dict_variables_gen[itype]
                        
                        if 'nGenBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen').Weight())  
                        elif 'nGenHadBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen (hadr. hem.)').Weight())  
                        elif 'nGenLepBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen (lept. hem.)').Weight())  
                        elif 'selGenHadHemDeltaR' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, AK4)').Weight())  
        
        if self.verbose: print (dictOfHists)
        
        return dictOfHists
    
    
    
    def getBranchesToRead(self, dirname, year, kinematic_labels, 
                          reco_only_labels, nSub_labels, 
                          jesSources=[], jerSources=[]):#,'pdfWeightAll'
        reco_list = []
        gen_list = [] if self.isMC else None
        gen_reweights_list = [] if self.isSigMC else None
        reco_reweights_list = [] if self.isSigMC else None

        if not self.isSigMC: #assurances against running into issues accidentally when processing non-signal MC
            self.wtSources=[]
            wtUnc=False

        
        if self.wtUnc: 
            self.sysSource = ['_nom'] + [ iwt+i for i in [ 'Up', 'Down' ] for iwt in self.wtSources if not iwt.endswith(('nom','pdfWeightAll')) ]
            
            if 'pdfWeightAll' in self.wtSources: self.sysSources = self.sysSources+['pdfWeightAll'] 
        #else: self.sysSource = ['_nom']
            
        #print (sysSources)
        for sys in self.sysSource:
            #if not wtUnc and not sysUnc:
            for i in kinematic_labels+list(nSub_labels.keys()): 
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                    reco_list.append('selRecoJets'+sys+i)
                
                if not ('softdrop' in i) and not('tau' in i):
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                        if i in self.dict_variables_kinematics_Muon:
                            reco_list.append('selRecoMu'+sys+i)
                        if i in self.dict_variables_kinematics_MET:
                            reco_list.append('selRecoMET'+sys+i)
                        if i in self.dict_variables_kinematics_LeptW:
                            reco_list.append('selRecoLeptW'+sys+i)
                        if i in self.dict_variables_kinematics_AK4Had:
                            reco_list.append('selRecoAK4bjetHadHem'+sys+i)

                if (self.isMC and (sys.endswith('nom') or self.sysUnc)): 
                    reco_list.append('trueRecoJets'+sys+i)
                    if not i in reco_only_labels: gen_list.append('accepGenJets'+sys+i)
                        
                if self.isMC and sys.endswith('nom'): 
                    if not i in reco_only_labels: 
                        gen_list.append('selGenJets'+sys+i)
                        if not 'softdrop' in i and not('tau' in i):
                            if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                                if i in self.dict_variables_kinematics_Muon:
                                    gen_list.append('selGenMu'+sys+i)
                                if i in self.dict_variables_kinematics_MET:
                                    gen_list.append('selGenMET'+sys+i)
                                if i in self.dict_variables_kinematics_LeptW:
                                    gen_list.append('selGenLeptW'+sys+i)
                                if i in self.dict_variables_kinematics_AK4Had:
                                    gen_list.append('selGenAK4bjetHadHem'+sys+i)

            if self.isMC and (sys.endswith('nom') or (self.isSigMC and self.sysUnc)):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("btagWeightNom"+sys)
                #reco_list.append("leptonWeightISO"+sys)
                #reco_list.append("leptonWeightID"+sys)
                #reco_list.append("leptonWeightTrig"+sys)
                #reco_list.append("leptonWeightRecoEff"+sys)
                reco_list.append("leptonWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("pdfWeightNom"+sys)

                reco_list.append("totalRecoWeight"+sys)
                reco_list.append(f"selRecoJets{sys}_msoftdrop_corr_PUPPI")#+sys)
                    
                reco_list.append("FlagRecoLeptHemBjet"+sys)
                reco_list.append("selRecoHadHemDeltaR"+sys)

                gen_list.append("evtGenWeight"+sys)
                gen_list.append("FlagGenLeptHemBjet"+sys)
                gen_list.append("selGenHadHemDeltaR"+sys)
                
                
                for i in self.dict_variables_reco:
                    reco_list.append(i+sys)
                for i in self.dict_variables_gen:
                    gen_list.append(i+sys)
                    
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)               
                reco_list.append("FlagRecoLeptHemBjet"+sys)
                reco_list.append("selRecoHadHemDeltaR"+sys)    
                for i in self.dict_variables_reco:
                    reco_list.append(i+sys)
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or  self.sysUnc):
                if 'pu' in sys or 'l1' in sys or 'lepton' in sys or 'btag' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list
        elif (not self.isMC): branchesToRead=reco_list
        
        #if self.verbose: print(branchesToRead)

        return branchesToRead
    
    
'''
import uproot
from uproot import *
import coffea
'''
from os.path import exists

from coffea import processor
import coffea.processor
from coffea.processor import defaultdict_accumulator,dict_accumulator#,IterativeExecutor,FuturesExecutor
import time
import re
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
from collections import OrderedDict
import numpy as np
import awkward as ak
import uproot
import hist
from hist import Hist
import gc 
from collections import defaultdict
import numba



class nSubBasis_unfoldingHistoProd_Dijets(processor.ProcessorABC):
    
    def __init__(self, sampleName, sysSource=[],year='2017', era='', isMC=True, isSigMC=True, 
                 onlyUnc='', wtUnc=False, sampleDict=dictSamples,test=False, sysUnc=False, 
                 saveParquet=False, verbose=True, splitCount='0',jetType='Central'):#isaltSigMC=False,
        self.test=test
        self.jetType=jetType
        self.year = year
        self.isMC = isMC
        self.isSigMC = isSigMC
        self.era = era
        #self.isaltSigMC = isaltSigMC
        self.verbose=verbose
        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc
        self.splitCount=splitCount
        self.saveParquet=saveParquet
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        
        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with data' )
        
        
        ### Helpers
        self.listOfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        #self.outputDir = 'UL17and18_nano_tests/'
        
        self.nJet = [ 'Jet'] if 'Central' in self.jetType else ['JetF']
        self.jetFlag = self.nJet[0].split('Jet')[1]
        if self.verbose: print (f"Initialising for {self.nJet}/{jetType}/{self.jetFlag}")


        self.triggerTable = OrderedDict()
        self.triggerTable[ 'AK8PFJet80' ] = {# from the list below, only the first two numbers (trigger turn on) are used. The others were a test - Alejandro
                                                '2016_preVFP' : [ 200,   242 ],
                                                '2016' : [ 200,   245 ],
                                                '2017' : [ 200,   257 ],
                                                '2018' : [ 200,   267 ],
                                                }
        self.triggerTable[ 'AK8PFJet140' ] = {
                                                '2016_preVFP' : [ 242,   308 ],
                                                '2016' : [ 245,   310 ],
                                                '2017' : [ 257,   323 ],
                                                '2018' : [ 267,   332 ],
                                                }
        self.triggerTable[ 'AK8PFJet200' ] = {
                                                '2016_preVFP' : [ 308,   373 ],
                                                '2016' : [ 310,   375 ],
                                                '2017' : [ 323,   389 ],
                                                '2018' : [ 332,   398 ],
                                                }
        self.triggerTable[ 'AK8PFJet260' ] = {
                                                '2016_preVFP' : [ 373,   439 ],
                                                '2016' : [ 375,   440 ],
                                                '2017' : [ 389,   455 ],
                                                '2018' : [ 398,   464 ],
                                                }
        self.triggerTable[ 'AK8PFJet320' ] = {
                                                '2016_preVFP' : [ 439,   526 ],
                                                '2016' : [ 440,   525 ],
                                                '2017' : [ 455,   543 ],
                                                '2018' : [ 464,   551 ],
                                                }
        self.triggerTable[ 'AK8PFJet400' ] = {
                                                '2016_preVFP' : [ 526,   580 ],
                                                '2016' : [ 525,   580 ],
                                                '2017' : [ 543,   598 ],
                                                '2018' : [ 551,   606 ],
                                                }
        self.triggerTable[ 'AK8PFJet450' ] = {
                                                '2016_preVFP' : [ 580,   635 ],
                                                '2016' : [ 580,   635 ],
                                                '2017' : [ 598,   653 ],
                                                '2018' : [ 606,   661 ],
                                                }
        self.triggerTable[ 'AK8PFJet500' ] = {
                                                '2016_preVFP' : [ 635,   10000000. ],
                                                '2016' : [ 635,   1000000. ],
                                                '2017' : [ 653.,  1000000. ],
                                                '2018' : [ 661,   1000000. ],
                                                }
        
        
        self.dict_variables_toUnfold = {
                    
                                "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.*200, 1.*201)]),
                                "_tau_0p5_2": np.array([(i/200) for i in np.arange(0.*200, 0.9*201)]),
                                "_tau_0p5_3": np.array([(i/200) for i in np.arange(0.*200, 0.8*201)]),
                                "_tau_0p5_4": np.array([(i/200) for i in np.arange(0.*200, 0.7*201)]),
                                "_tau_0p5_5": np.array([(i/200) for i in np.arange(0.*200, 0.7*201)]),
                                "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.9*201)]),
                                "_tau_1_2": np.array([(i/200) for i in np.arange(0.*200, 0.6*201)]),
                                "_tau_1_3": np.array([(i/200) for i in np.arange(0.*200, 0.4*201)]),
                                "_tau_1_4": np.array([(i/200) for i in np.arange(0.*200, 0.3*201)]),
                                "_tau_1_5": np.array([(i/200) for i in np.arange(0.*200, 0.3*201)]),
                                "_tau_2_1": np.array([(i/200) for i in np.arange(0.*200, 0.5*201)]),
                                "_tau_2_2": np.array([(i/200) for i in np.arange(0.*200, 0.3*201)]),
                                "_tau_2_3": np.array([(i/200) for i in np.arange(0.*200, 0.2*201)]),
                                "_tau_2_4": np.array([(i/200) for i in np.arange(0.*200, 0.2*201)]),
                                "_tau_2_5": np.array([(i/200) for i in np.arange(0.*200, 0.2*201)]),
                                "_tau21": np.array([(i/200) for i in np.arange(0.*200, 1.6*201)]),#for one-pass kT minimization as per CMS
                                "_tau32": np.array([(i/200) for i in np.arange(0.*200, 1.6*201)]),#for one-pass kT minimization as per CMS
                                "_tau21_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.4*201)]),#for WTA-kT for comparison
                                "_tau32_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.4*201)]),#for WTA-kT for comparison
                                "_tau21_exkT": np.array([(i/200) for i in np.arange(0.*200, 2.*201)]),#for excl.-kT and E-scheme as per basis
                                "_tau32_exkT": np.array([(i/200) for i in np.arange(0.*200, 2.*201)]),#for excl.-kT and E-scheme as per basis

                               }
        self.kinematic_labels = ['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels = ['_good_nPVs']

        self.dict_variables_kinematics = {
                                            "_pt": np.array([i for i in np.arange(70., 2570., 50.)]),
                                            "_eta": np.array([i for i in np.arange(-2.2, 2.4, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.2, 2.4, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.2, 3.4, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 360., 10.)]),
                                            "_msoftdrop": np.array([i for i in np.arange(0., 360., 10.)]),
                                            "_good_nPVs": np.array([i for i in np.arange(0., 101., 1.)]),
                                         }
        
        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not( isys.endswith('nom'))]# or self.sysUnc) ]
        
        self.sysWeightList = ( '_pu', '_pdf', '_isr', '_fsr', '_l1prefiring' ) #'_ps',

        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight'] if self.wtUnc else [] 
        self.recoWtSources=['_pu', '_l1'] if self.isMC else [] 
        
        if self.onlyUnc: self.sysSource = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] 

        if self.verbose: print('In __init__',self.sysSource)

        self.selList = '_dijetSel' 
        
        self._branchesToRead = self.getBranchesToRead(
                                                         dirname=self.inputDir, year=self.year,
                                                         kinematic_labels=self.kinematic_labels,
                                                         reco_only_labels=self.reco_only_labels,
                                                         nSub_labels=self.dict_variables_toUnfold,
                                                     )
            
        
    def process(self, events):
        

        output = self.buildDictOfHistograms()
        branches = self._branchesToRead
        events = events[branches]
        
        if self.verbose: 
            print (f'Now, processing histograms; systematic sources being considered:{self.sysSource}') 
            
        for sys in self.sysSource:

            ##############################
            #### Building event masks ####
            ############################## 
            if self.isMC and self.isSigMC:# or self.isaltSigMC:
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys
                

                selRecoMask = events[f'totalRecoWeight{s}']!=0
                selGenMask = events[f'evtGenWeight_nom']!=0

                trueRecoMask = (selRecoMask) & (selGenMask) & ((events[f'trueRecoJets{s}_pt']>0.) & (events[f'trueRecoJetsF{s}_pt']>0.))
                
                fakeRecoMask = (selRecoMask) & ((events[f'trueRecoJets{s}_pt']<0.) | (events[f'trueRecoJetsF{s}_pt']<0.) ) #& (~selGenMask)
                
                accepGenMask = (selGenMask) & (selRecoMask) & ((events[f'accepGenJets{s}_pt']>0.) & (events[f'accepGenJetsF{s}_pt']>0.))
                
                missGenMask =  (selGenMask)  & ((events[f'accepGenJetsF{s}_pt']<0.) | (events[f'accepGenJets{s}_pt']<0.) ) #& (~selRecoMask)
                if (self.verbose and (sys.endswith('nom') or self.sysUnc)): 
                    print('#### Building event masks ####')
                    print(sys,s, 'masked array lengths for reco,true,fake,gen,accep,miss', len(events),
                          len(events[selRecoMask]),
                          len(events[trueRecoMask]),
                          len(events[fakeRecoMask]),
                          len(events[selGenMask]),
                          len(events[accepGenMask]),
                          len(events[missGenMask]),
                          self.splitCount
                         )

                if self.saveParquet: 
                    events['trueRecoMask']=trueRecoMask
                    events['selRecoMask']=selRecoMask
                    events['fakeRecoMask']=fakeRecoMask
                    events['selGenMask']=selGenMask
                    events['accepGenMask']=accepGenMask
                    events['missGenMask']=missGenMask
                    

                    if (sys.endswith('nom')): 
                        print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_V2_{self.splitCount}.parquet'}")
                        
                        ak.to_parquet(events,self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_V2_{self.splitCount}.parquet')
                        #return 1 # dummy 
                    elif self.sysUnc and self.onlyUnc: 
                        print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_{s}Unc_{self.jetType}Jet_V2_{self.splitCount}.parquet'}")
                        
                        ak.to_parquet(events, self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_{s}Unc_{self.jetType}Jet_V2_{self.splitCount}.parquet')
                        #return 1 # dummy 
                    
                    
                
            elif self.isMC and (not self.isSigMC): #not so relevant for dijets but for background MC's in W/top
                
                if sys.endswith('nom'):
                    selRecoMask = events[f'totalRecoWeight{sys}']!=0.
                    selGenMask = events[f'evtGenWeight{sys}']!=0.
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                selRecoMasks = OrderedDict()
                #masking now trigger dependent to create different histos for data events passing a given prescale trigger exclusively or an unprescaled trigger
                print("Building trigger masks")
                for itrigger, itrValues in self.triggerTable.items():    
                    triggerList=list(self.triggerTable.keys())
                    thistrigInd=triggerList.index(itrigger)
                    othertrigsInds=[i for i in range(0,len(triggerList)) if i!=thistrigInd]
                    
                    selRecoMasks[itrigger]=(events[f'totalRecoWeight{sys}']!=0.) & (events[f'passHLT_{itrigger}']==1 ) & ((events[f'selRecoJets{sys}_pt']>itrValues[self.year][0]) | (events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0])) & ((events[f'selRecoJets{sys}_pt']  < itrValues[self.year][1]) & (events[f'selRecoJetsF{sys}_pt']  < itrValues[self.year][1]))
                    

                    if self.saveParquet: 
                        
                        
                        for itrigger, itrValues in self.triggerTable.items():
                            events['selRecoMask_'+itrigger]=selRecoMasks[itrigger]

                        if sys.endswith('nom'): 
                            print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_nomWts_{self.jetType}Jet_V2_{self.splitCount}.parquet'}")
                            ak.to_parquet(events, (self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_nomWts_{self.jetType}Jet_V2_{self.splitCount}.parquet'))                                    
                    
                   

            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            
            for isel in self.selList:
                
                
                listOfOutputHistosNames = [k for k,h in output.items() if ((sys in k) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                #if self.verbose: 
                #    print(listOfOutputHistosNames)

                for k in listOfOutputHistosNames:
                    key=k
                    
                    ############### Safety checks ##################
                    if not(sys in key) and not('resol' in key): 
                        #if self.verbose: print(sys, key)
                        continue

                    if not( self.isMC) and not('AK8PF' in key):
                        continue

                    if not( self.isMC) and 'AK8PF' in key:
                        if not( isel.split('_')[1] in key): 
                            continue #to ensure that for data, we only fill the histos relevant to a given trigger (stored in isel as '_AKPFJetXYZ_dijetSel') correctly

                    ############### ############ ##################

                    # Decide on trigger (for data) and variable labels
                    if (not self.isMC) and 'AK8PF' in key:                                
                        whichTrigList = [trig for trig in self.triggerTable.keys() if trig in key]
                        whichTrig = whichTrigList[0]#f'AK8PF{key.split("AK8PF")[1]}'

                    if not('_tau' in key):
                        whichKinVar = [var for var in self.dict_variables_kinematics.keys() if var[1:] in key]#[0]
                        varToFill = whichKinVar[0]

                    elif ('_tau' in key or '__' in key):# and 'nom' in key):
                        temp='WTA' if 'WTA' in key else 'tau' #hack, :(, to fix tau21_nonOPkT being histogrammed incorrectly
                        temp='exkT' if 'exkT' in key else 'tau'
                        if not 'tau' in key: temp="__"
                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                        varToFill = whichnSub[0] if 'tau' in whichnSub[0] else whichnSub[0].replace('__','_')
                    
                    #if self.verbose and sys.endswith('nom'): print(varToFill)

                    ################## Filling histos from accumulated event arrays ##################

                    # Decide what weights to use
                    if not (sys.startswith(self.sysWeightList)):
                        s=sys
                        
                        totalRecoWeight = events[f'evtGenWeight_nom']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']  if self.isMC else events[f'totalRecoWeight_nom']
                        totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None

                    else:
                        s='_nom'
                        if 'pu' in sys:
                            if self.verbose: print(s,sys,key,f'evtGenWeight{s}',f'{sys.split("_")[1]}{s}')                            
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'l1prefiringWeightNom{s}']
                            
                        elif 'l1' in sys:
                            if self.verbose: print(s,sys,key,f'evtGenWeight{s}',f'{sys.split("_")[1]}{s}')                           
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']
                            
                        elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                            if self.verbose: print(s,sys,key,f'evtGenWeight{s}',f'{sys.split("_")[1]}{s}')
                            totalGenWeight = events[f'evtGenWeight{s}']*events[f'{sys.split("_")[1]}{s}']
                            totalRecoWeight = totalGenWeight*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']

                    if (key.lower().startswith(('reco','good'))):
                        if 'nPV' in key:
                            output[key].fill(events[f'good_nPVs{s}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)
                        else: 
                            if self.verbose: print(s,sys, key,f"filling from recoJet{self.nJet[0]}")#,events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask])
                            
                            output[key].fill(events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)

                    # Stuff below should have diff wts already available to them so no need to touch anything down here
                    if self.isMC and (not varToFill.startswith(tuple(self.reco_only_labels))):
                        if self.verbose: 
                            print(s, sys, key)
                            
                        if (key.lower().startswith('true')):
                            output[key].fill(events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                        elif (key.lower().startswith('fake')):
                            output[key].fill(events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                             threads=8)

                        if (key.lower().startswith('genjet')):
                            output[key].fill(events[f'selGenJets{self.jetFlag}_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)#=self._listofHistograms(histoName) #hist.Hist("Events", hist.Cat("branch", branch), hist.Bin("value", branch, 100, 0, 1000))    

                        elif (key.lower().startswith('accepgenjet')):# and not key.startswith(('accep','miss')):
                            output[key].fill(events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                             threads=8)

                        elif (key.lower().startswith('missgenjet')):
                            output[key].fill(events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                             threads=8)
                        
                        elif ( 'resp' in key.lower() and not('miss' in key.lower())):
                            #if self.verbose: print("filling resp")
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)     
                        
                        elif ('respwithmiss' in key.lower()):
                            #if self.verbose: print("filling resp with miss")
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask], weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)
                            
                            #fill missgen weight
                            output[key].fill(gen=events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask], reco=-1.*np.ones(len(events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask])),
                                             weight=totalGenWeight[missGenMask],
                                             threads=8)
                            

                        elif ('resol' in key.lower()) and sys.endswith('nom'):
                            zeroMask=(events[f'accepGenJets{self.jetFlag}{s}{varToFill}']!=0.)&(accepGenMask)

                            response = events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]/events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][zeroMask]
                            response = np.nan_to_num(response,nan=-999.)
                            
                            output[key].fill(response, weight=totalRecoWeight[zeroMask])
                            if 'noWt_' in key: output[key].fill(response)
        l=[]
        
        for x,y in output.items(): #y.SetDirectory(0)
        
            if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss','true','fake' )) and '_nom' in x:
                l.append(x)
        for k in l:
            del(output[k])
        gc.collect()
        return output

    def postprocess(self, accumulator):
        pass
        '''
        with uproot.recreate(f'UL17and18_nano_tests/{self.dictSamples[self.sampleName][self.year]["skimmerHisto"]}') as fout:#outputTest_{sample.split("_Tune")[0]}_{year}.root'
            #if i% 10==0:
            print (f'Starting to save output file:{fout}')

            for key,value in accumulator.items():
                if '_nom' in key and self.onlyUnc: continue
                fout[key]=processed_events[key]
            #if i% 10==0:print (f'Done with creating output file:{fout.file_path}')
            fout.close()
        
        return accumulator#rite(accumulator, "output.root", "myhist")
        '''
        
   
    def buildDictOfHistograms(self):
        '''build dictionary of Hist histograms, convert to root or whatever else after filled and returned by processor'''
        dictOfHists = OrderedDict()

        self.selList = [ '_'+x+'_dijetSel' for x in self.triggerTable  ] if not self.isMC else [ '_dijetSel' ]
        if not self.isMC: print(self.selList)
        for isel in self.selList:

            for itype in self.listOfHistTypes:
                iJ=self.nJet[0] #stupid hard-coding, to remove eventually

                for sysUnc in self.sysSource:
                    
                    for x, y in self.dict_variables_kinematics.items():
                        binning = y
                        if sysUnc.endswith('nom') or self.onlyUnc:
                            
                            if not x in tuple(self.reco_only_labels): 
                                dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())        
                            
                            else: 
                                #for nPVs
                                if sysUnc.endswith('nom') or sysUnc.startswith(tuple(self.recoWtSources)): 
                                    dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())  
                                
                    if itype.startswith('truereco') and self.isMC and sysUnc.endswith('nom'):
                            dictOfHists['resol'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_pt'+isel, label='AK8 reco/gen jet pt', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD}/gen jet m_{SD}', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_mass'+isel, label='AK8 reco inv. m/gen jet inv. m', underflow=True,overflow=True).Weight())
                    #if self.verbose and sysUnc.endswith('nom'): print('building unfolding histos from',self.dict_variables_toUnfold.keys())     
                    for x, y in self.dict_variables_toUnfold.items():
                        binning = y
                        
                        #binning_coarse=np.array([binning[i] for i in range(len(binning)) if i%5==0])
                        
                        dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label='AK8 '+itype+' jet #tau', underflow=True,overflow=True).Weight())
                        
                        if itype.startswith('truereco') and self.isMC:

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
                            if sysUnc.endswith('_nom'):
                                dictOfHists['resol'+iJ+x+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
                                dictOfHists['noWt_resol'+iJ+x+isel] = (hist.Hist.new.Regular(500, 0, 5, name='noWt_resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
        
        if self.verbose: print("Final dict of histograms keys:",dictOfHists.keys())      
        return dictOfHists
    
    
    
    def getBranchesToRead(self, dirname, year, kinematic_labels, 
                          reco_only_labels, nSub_labels): #jesSources=[], jerSources=[]):#,'pdfWeightAll'
                          
        reco_list = []
        gen_list = [] if self.isMC else None
        gen_reweights_list = [] if self.isSigMC else None
        reco_reweights_list = [] if self.isSigMC else None
        triggerBit_list = [] if not(self.isMC) else None

        if not self.isSigMC: #assurances against running into issues accidentally when processing non-signal MC
            self.wtSources=[]
            self.wtUnc=False

        
        if self.wtUnc: 
            self.sysSource = ['_nom'] + [ iwt+i for i in [ 'Up', 'Down' ] for iwt in self.wtSources if not iwt.endswith(('nom','pdfWeightAll')) ]
            
            #if 'pdfWeightAll' in self.wtSources: self.sysSources = self.sysSources+['pdfWeightAll'] 
            
        if self.verbose: print ("Preparing branches to read",self.sysSource)
            
        for sys in self.sysSource:
            #if not wtUnc and not sysUnc:
            
            if 'Central' in self.jetType:
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                    reco_list.append(f'selRecoJetsF{sys}_pt')

                if (self.isMC and not(self.isSigMC) and sys.endswith('nom')) or (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    if not ('selGenJetsF_nom_pt' in gen_list): gen_list.append('selGenJetsF_nom_pt')
                    gen_list.append(f'accepGenJetsF{sys}_pt')
                    reco_list.append(f'trueRecoJetsF{sys}_pt')
                
                for i in kinematic_labels+list(nSub_labels.keys()): 
                    if '__' in i: i=i.replace('__','_')
                    
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                        reco_list.append('selRecoJets'+sys+i)
                        
                    if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                        reco_list.append('trueRecoJets'+sys+i)
                        if not i in reco_only_labels: gen_list.append('accepGenJets'+sys+i)
                            
                    if self.isMC and sys.endswith('nom'): 
                        if not i in reco_only_labels: gen_list.append('selGenJets'+sys+i)

            elif 'Forward' in self.jetType: 
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                    reco_list.append(f'selRecoJets{sys}_pt')

                if (self.isMC and not(self.isSigMC) and sys.endswith('nom')) or (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    if not ('selGenJets_nom_pt' in gen_list): gen_list.append('selGenJets_nom_pt')
                    gen_list.append(f'accepGenJets{sys}_pt')
                    reco_list.append(f'trueRecoJets{sys}_pt')
                
                for i in kinematic_labels+list(nSub_labels.keys()): 
                    if '__' in i: i=i.replace('__','_')
                    
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                        reco_list.append('selRecoJetsF'+sys+i)
                        
                    if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                        reco_list.append('trueRecoJetsF'+sys+i)
                        if not i in reco_only_labels: gen_list.append('accepGenJetsF'+sys+i)
                            
                    if self.isMC and sys.endswith('nom'): 
                        if not i in reco_only_labels: gen_list.append('selGenJetsF'+sys+i)

            if (self.isMC and sys.endswith('nom')) or (self.isSigMC and self.sysUnc):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("totalRecoWeight"+sys)
                if not ("evtGenWeight_nom" in gen_list): gen_list.append("evtGenWeight_nom")#+sys)
                reco_list.append("good_nPVs"+sys)
                
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("good_nPVs"+sys)
                for itrigger in self.triggerTable.keys():    
                    if not(self.year.endswith('VFP') and self.era=='B'): triggerBit_list.append(f'HLT_{itrigger}')
                    triggerBit_list.append(f'passHLT_{itrigger}')
                
                if '2018' in self.year or '2017' in self.year: 
                    triggerBit_list.append(f'passHLT_AK8PFJet550')
                    triggerBit_list.append(f'HLT_AK8PFJet550')
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or self.sysUnc):
                if 'pu' in sys or 'l1' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list
        elif (not self.isMC): branchesToRead=reco_list+triggerBit_list
        if self.verbose: 
            print("Prepared branches to read",branchesToRead)        
        return branchesToRead
    


#################################### Helper for histo accumulation ##############################################

from itertools import islice
#import concurrent
#from concurrent import futures
#from coffea.processor import FuturesExecutor,Runner
from coffea.processor import accumulate
from uproot import ThreadPoolExecutor

#from coffea.processor import defaultdict_accumulator,dict_accumulator#,IterativeExecutor,FuturesExecutor
import time
import re
from rich.progress import track

########modify for data, potentially just put this in as the postprocessor in the histoproducer module
#{writeAccumulatedORhadd: 0=write accumulated hists in a separate file after saving or not saving constituent file chunks,
#1 = hadd a posteriori the chunks}
def histoMaker_BasicButSplit(myProcessor, sampleIdentifier='qcd_ht', y='2017', sampleDict_PFNano=OrderedDict(),
                             sampleDict_local=OrderedDict(), writeChunks=True, writeAccumulatedORhadd=1, verbose=False, saveParquet=False,
                             isSigMC=True, isMC=True, wtUnc=True, sysUnc=False, onlyUnc='', outputdir='processorTests/', sysSource=[],
                             era='', ext='_nomWts', splitchunks=5, nWorkers=10, stepSize="2048 MB",jetType='Central'):#for jetType: options='Central','Forward' (1 at a time)
    nchunk=copy.deepcopy(splitchunks)
    for sample in sampleDict_local.keys():
        yl=['2016_preVFP', '2016', '2017', '2018'] if y=='all' else [y]
        #if verbose: print (sample,sampleIdentifier,splitchunks)
        
        for year in yl:#['2017','2018','2016','2016_preVFP']:
            gc.collect()
            
            if not( sample.lower().startswith(sampleIdentifier.lower())):# or ('100to200' in sample.lower()) or ('200to300' in sample.lower()) or ('300to500' in sample.lower()) or ('500to700' in sample.lower()): 
                continue
            tstart0=time.time()

            print(f'Histogramming for {sample} in year: {year}')
            if verbose: print(f'using nanoskims from {sampleDict_local[sample][year]["t3_dirs"]}{chr(10)}','\n')
            
            
            if not os.path.exists(f'{outputdir}'): os.makedirs(f'{outputdir}')

            inpdir=sampleDict_local[sample][year]['t3_dirs']
            dirfnames=inpdir[0].split('/0000/')[0]+'/000*/jetObservables_nanoskim_*.root'
            #if verbose: print (dirfnames)
            SC='0'
            if verbose: print (sysSource)
            my_processor=myProcessor(sampleName=sample, sampleDict=sampleDict_PFNano,
                                     isMC=isMC, isSigMC=isSigMC,year=year,saveParquet=saveParquet,sysSource=sysSource,
                                     wtUnc=wtUnc,sysUnc=sysUnc,onlyUnc='' if 'jes' in onlyUnc else onlyUnc,era=era,
                                     verbose=False,splitCount=SC )
            if splitchunks>0:
                fl=[]
                n_tosplit = []
                c=0
                for x in range(len(inpdir)):
                    flist=os.listdir(inpdir[x])
                    flist = sorted([i for i in flist if 'nanoskim' in i], key=lambda s: int(re.search(r'\d+', s).group()))
                    
                    for i in flist:
                        if 'nanoskim' in i:

                            fl.append(inpdir[x]+i)
                            if c%splitchunks==0 and c!=0:
                                n_tosplit.append(splitchunks)
                            c=c+1
                
                if sum(n_tosplit)!=len(fl):
                    if len(fl)%splitchunks!=0: n_tosplit.append(len(fl)%splitchunks)
                    else: n_tosplit.append(splitchunks)

                if verbose: print(f"Splitting input filelist into the following sublists of file chunks:{n_tosplit}")
                ifl=iter(fl)
                ifls=[list(islice(ifl,x)) for x in n_tosplit]
                #if verbose: print(ifls)
                c=0
                hists=None
            else:
                ifls=[dirfnames]
                c=0
            if verbose: print(f"Loading events into arrays")# from the following file chunks:{n_tosplit}")    
            
            
            for i in track(ifls):#
                #print (c,i)
                if splitchunks>0:
                    ns = [int(x.split('jetObservables_nanoskim_')[1].split('.root')[0]) for x in i]
                    s = str(ns).split('[')[1].split(']')[0]
                    stringfnames=i[0].split('jetObservables_nanoskim_')[0]+"jetObservables_nanoskim_{"+f'{s}'+"}.root"
                else: stringfnames=i
                    
                stringfnames=stringfnames.replace(' ','')
                
                if verbose:
                    tstart1 = time.time()
                
                events = uproot.concatenate(stringfnames+':Events', my_processor._branchesToRead, 
                                            step_size=stepSize, library='ak', num_workers=nWorkers,
                                            #interpretation_executor=uproot.ThreadPoolExecutor(num_workers=nWorkers), 
                                            #decompression_executor=uproot.ThreadPoolExecutor(num_workers=nWorkers),
                                           )
                
                if verbose:
                    elapsed1 = time.time()-tstart1
                    print(f"Time taken to load {len(events)} events from {n_tosplit[c]} files in {sample} = {elapsed1}")
                
                if verbose: print(f'nEvents in this file chunk: {len(events)}')
                    
                if verbose: tstart2=time.time()

                if verbose: print("#####Processing events#####")
                    
                processed_events=my_processor.process(events)
                 
                if verbose:
                    elapsed2 = time.time()-tstart2
                    print(f"Time taken to build histos from {len(events)} events = {elapsed2}")
                
                my_processor.splitCount=chr(int(SC)+1)
                
                if c==0: 
                    hists=processed_events
                    if verbose: print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                elif c>0 and writeAccumulatedORhadd==0: 
                    hists=accumulate([hists,processed_events])
                    if verbose: print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                
                if not sysUnc: 
                    if 'pt' in sample.lower(): 
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]+sample.split("_Tune")[0].split("_")[2]}_UL{year}{ext}'
                    else:
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]}_UL{year}{ext}'

                    fnstem = f'{outputdir}/jetObservables_histograms_{string}'    
                    fn = f'{fnstem}_ForwardJet_{c}.root' if 'Forward' in jetType else f'{fnstem}_CentralJet_{c}.root'
                else: 
                    
                    fnstem = f'{outputdir}/{sampleDict_local[sample][year]["skimmerHisto"].split(".root")[0]}{ext}'
                    if 'jes' in onlyUnc: 
                        fnstem=f'{outputdir}/combinedJES/{fnstem.split(outputdir+"/")[1]}' 
                        if not os.path.exists(f'{outputdir}/combinedJES/'): os.makedirs(f'{outputdir}/combinedJES/')

                    fn = f'{fnstem}_ForwardJet_{c}.root' if 'Forward' in jetType else f'{fnstem}_CentralJet_{c}.root'
                    
                if splitchunks==0 or len(n_tosplit)==1: 
                    fn=f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root'
                if verbose: print (fn,fnstem)

                if writeChunks or writeAccumulatedORhadd==1:# or len(n_tosplit)==1:
                    if verbose: print(f"Current integral of pT nominal hists being saved: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")            
                    

                    #print (processed_events)
                    print ("Events proccesed, now writing output file(s)",fn)
                    
                    with uproot.recreate(fn) as fout:#outputTest_{sample.split("_Tune")[0]}_{year}.root'
                        for key in processed_events.keys():
                            fout[key]=processed_events[key]
                        print (f'Done with creating output file:{fn}')
                        fout.close()
                    
                c=c+1
                del(processed_events)
                del(events)
                gc.collect()
                
            if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):

                print(f"Integral of (reco) nominal hists being saved: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")            
                
                if not os.path.exists(f'{outputdir}/combinedJES/'): os.makedirs(f'{outputdir}/combinedJES/')
                
                print (f'Output file stem for accumulated histos: {fnstem}')
                with uproot.recreate(f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root') as fout_hadd:
                    for key in hists.keys():
                        fout_hadd[key]=hists[key]
                    print (f'Done with creating accumulated output file:{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root')
                    fout_hadd.close()
            #elif writeAccumulatedORhadd==1 and c>0 and splitchunks>0:
            #    if 'jes' in onlyUnc: 
            #        fnstem=f'{outputdir}/combinedJES/{fnstem.split(outputdir+"/")[1]}'
            #    input_haddfilelist = []
            #    for file in os.listdir(fnstem.split('jetObservables_')[0]):
            #        full_filepath=os.path.abspath(fnstem.split('jetObservables_')[0])+file
            #        if fnstem in full_filepath: 
            #            input_haddfilelist.append(file)
            #    print (f'Input files for hadded histos: {input_haddfilelist}')
            #    events = uproot.concatenate(input_haddfilelist)
            #    with uproot.recreate(f'{fnstem}.root') as fout_hadd:
            #        for key in events.fields():#keys():
            #            fout_hadd[key]=events[key]
            #        print (f'Done with creating accumulated output file:{fnstem}.root')
            #        fout_hadd.close()
            #else:
            #    print("separated output files written, the rest is in your hands! :)")
                
            elapsed0 = time.time() - tstart0

            print (f'Time for {sample} {year}:{elapsed0}')        

            del(my_processor)
            #del(events)
            #del(processed_events)
            del(hists)
            gc.collect()
            
            clear_output(wait=True)


    return 1


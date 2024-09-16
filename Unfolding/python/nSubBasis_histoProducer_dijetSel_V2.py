'''
import uproot
from uproot import *
import coffea
'''
import os,time,re,gc,copy
from os.path import exists

from coffea import processor
import coffea.processor
from coffea.processor import defaultdict_accumulator,dict_accumulator#,IterativeExecutor,FuturesExecutor
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
from collections import OrderedDict
import numpy as np
import awkward as ak
import uproot
import hist
from hist import Hist
from collections import defaultdict
from itertools import islice
from coffea.processor import accumulate
from uproot import ThreadPoolExecutor
from rich.progress import track



class nSubBasis_unfoldingHistoProd_Dijets():#processor.ProcessorABC
    
    def __init__(self, sampleName, sysSource=[],year='2017', era='', 
                 isMC=True, isSigMC=True, onlyUnc='', wtUnc=False, verbose=False, saveParquet=False, onlyParquet=False,
                 sampleDict=dictSamples,test=False, sysUnc=False, splitCount='0', jetType='Central', trigTest=False, 
                 trigUpDownVal=10, withLepVeto=False, parquetDir='/scratch/kadatta/dijetChecks/parquets/'):
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
        self.onlyParquet=onlyParquet   
        self.parquetDir=parquetDir
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        self.trigTest = trigTest
        self.trigUpDownVal = trigUpDownVal if self.trigTest else 0.
        self.withLepVeto = withLepVeto
            
        if (not self.isMC) and self.era=='': 
            
            print (f'Data-loading error: You need to specify what era if you want to work with data' )
        
        
        ### Helpers
        self.listOfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        #self.outputDir = 'UL17and18_nano_tests/'
        
        self.nJet = [ 'Jet'] if 'Central' in self.jetType else ['JetF']
        self.jetFlag = self.nJet[0].split('Jet')[1]
        if self.verbose: print (f"Initialising for {self.nJet}/{jetType}/{self.jetFlag}")


        self.triggerTable = OrderedDict()
        '''
        #new turn on (V2) values from better modelling of the plateau
        self.triggerTable[ 'AK8PFJet80' ] = {    
                    '2016_preVFP': [130.8, 203.52],
                    '2016': [128.81, 201.90],
                    '2017': [158.3, 229.03],
                    '2018': [165.83, 237.12],
                    }
        self.triggerTable[ 'AK8PFJet140' ] = {
                    '2016_preVFP': [203.52, 276.25],
                    '2016': [201.90, 274.99],
                    '2017': [229.03, 299.75],
                    '2018': [237.12, 308.41],
                    }
        self.triggerTable[ 'AK8PFJet200' ] = {
                    '2016_preVFP': [276.25, 348.98],
                    '2016': [274.99, 348.07],
                    '2017': [299.75, 370.48],
                    '2018': [308.41, 379.7],
                    }
        self.triggerTable[ 'AK8PFJet260' ] = {
                    '2016_preVFP': [348.98, 421.7],
                    '2016': [348.07, 421.16],
                    '2017': [370.48, 441.21],
                    '2018': [379.7, 450.99],
                    }
        self.triggerTable[ 'AK8PFJet320' ] = {
                    '2016_preVFP': [421.7, 518.67],
                    '2016': [421.16, 518.61],
                    '2017': [441.21, 535.51],
                    '2018': [450.99, 546.05],
                    }
        self.triggerTable[ 'AK8PFJet400' ] = {
                    '2016_preVFP': [518.67, 579.28],
                    '2016': [518.61, 579.52],
                    '2017': [535.51, 594.45],
                    '2018': [546.05, 605.46],  
                    }
        self.triggerTable[ 'AK8PFJet450' ] = {
                    '2016_preVFP': [579.28, 6500.],#639.88],
                    '2016': [579.52, 6500.],#640.42],
                    '2017': [594.45, 653.39],
                    '2018': [605.46, 664.86],
                    }
        self.triggerTable[ 'AK8PFJet500' ] = {
                    #'2016_preVFP': [639.88, 6500.],
                    #'2016': [640.42, 6500.],
                    '2017': [653.39, 6500.],#712.33],
                    '2018': [664.86, 6500.],#724.27],
                    }
        #if '2017' in self.year or '2018' in self.year:
        #    self.triggerTable[ 'AK8PFJet550' ] = {
        #                '2017' : [ 712.33, 6500.],
        #                '2018' : [ 724.27, 6500.],
        #                }
        '''   
        
        #new turn on (V3) values from better modelling of the plateau
        self.triggerTable['AK8PFJet80'] = {
                                            '2016': [117.649,191.368],
                                            '2016_preVFP': [121.938,192.909],
                                            '2017': [146.47,220.122],
                                            '2018': [141.586,224.871],
                                           }
        self.triggerTable['AK8PFJet140'] = {
                                            '2016': [191.368,263.426],
                                            '2016_preVFP': [192.909,266.569],
                                            '2017': [220.122,294.317],
                                            '2018': [224.871,299.834],
                                           }
        self.triggerTable['AK8PFJet200'] = {
                                            '2016': [263.426,335.765],
                                            '2016_preVFP': [266.569,338.288],
                                            '2017': [294.317,363.808],
                                            '2018': [299.834,372.496],
                                           }
        self.triggerTable['AK8PFJet260'] = {
                                            '2016': [335.765,409.866],
                                            '2016_preVFP': [338.288,411.402],
                                            '2017': [363.808,435.629],
                                            '2018': [372.496,443.459],
                                           }
        self.triggerTable['AK8PFJet320'] = {
                                            '2016': [409.866,509.961],
                                            '2016_preVFP': [411.402,510.262],
                                            '2017': [435.629,530.812],
                                            '2018': [443.459,540.081],
                                           }
        self.triggerTable['AK8PFJet400'] = {
                                            '2016': [509.961,560.765],
                                            '2016_preVFP': [510.262,568.293],
                                            '2017': [530.812,580.872],
                                            '2018': [540.081,585.585],
                                           }
        self.triggerTable['AK8PFJet450'] = {
                                            '2016': [560.765, 6500.],#634.024],
                                            '2016_preVFP': [568.293, 6500.],#631.693],
                                            '2017': [580.872,640.092],
                                            '2018': [585.585,645.01],
                                           }
        if '2017' in self.year or '2018' in self.year:
        
            self.triggerTable['AK8PFJet500'] = {
                                                #'2016': [634.024,6500.0],
                                                #'2016_preVFP': [631.693,6500.0],
                                                '2017': [640.092,6500.0],#,700.13],
                                                '2018': [645.01,6500.0],#,701.393],
                                               }
        #if '2017' in self.year or '2018' in self.year:
        #    self.triggerTable['AK8PFJet550'] = {
        #                                        '2017': [700.13,6500.0],
        #                                        '2018': [701.393,6500.0],
        #                                       }
        if self.trigTest:
            for it in list(self.triggerTable.keys()):
                if self.year in self.triggerTable[it].keys():
                    self.triggerTable[it][self.year][0] += trigUpDownVal
                    if not(self.triggerTable[it][self.year][1]==6500.):
                        self.triggerTable[it][self.year][1] += trigUpDownVal
            
            
        
        
        self.dict_variables_toUnfold = {
                    
                                "_tau_0p25_1": np.array([(i/200) for i in np.arange(0.*200, 1.005*200)]),
                                "_tau_0p25_2": np.array([(i/200) for i in np.arange(0.*200, 0.905*200)]),
                                "_tau_0p25_3": np.array([(i/500) for i in np.arange(0.*500, 0.852*500)]),
                                "_tau_0p25_4": np.array([(i/500) for i in np.arange(0.*500, 0.852*500)]),
                                "_tau_0p25_5": np.array([(i/500) for i in np.arange(0.*500, 0.802*500)]),
            
                                "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.*200, 1.005*200)]),
                                "_tau_0p5_2": np.array([(i/200) for i in np.arange(0.*200, 0.855*200)]),
                                "_tau_0p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.752*500)]),
                                "_tau_0p5_4": np.array([(i/500) for i in np.arange(0.*500, 0.702*500)]),
                                "_tau_0p5_5": np.array([(i/500) for i in np.arange(0.*500, 0.652*500)]),
            
                                "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.855*200)]),
                                "_tau_1_2": np.array([(i/200) for i in np.arange(0.*200, 0.505*200)]),
                                "_tau_1_3": np.array([(i/500) for i in np.arange(0.*500, 0.352*500)]),
                                "_tau_1_4": np.array([(i/500) for i in np.arange(0.*500, 0.302*500)]),
                                "_tau_1_5": np.array([(i/500) for i in np.arange(0.*500, 0.252*500)]),
                                
                                "_tau_1p5_1": np.array([(i/200) for i in np.arange(0.*200, 0.605*200)]),
                                "_tau_1p5_2": np.array([(i/200) for i in np.arange(0.*200, 0.325*200)]),
                                "_tau_1p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.252*500)]),
                                "_tau_1p5_4": np.array([(i/500) for i in np.arange(0.*500, 0.202*500)]),
                                "_tau_1p5_5": np.array([(i/500) for i in np.arange(0.*500, 0.182*500)]),
                                
                                
                                "_tau_2_1": np.array([(i/200) for i in np.arange(0.*200, 0.505*200)]),
                                "_tau_2_2": np.array([(i/200) for i in np.arange(0.*200, 0.225*200)]),
                                "_tau_2_3": np.array([(i/500) for i in np.arange(0.*500, 0.152*500)]),
                                "_tau_2_4": np.array([(i/500) for i in np.arange(0.*500, 0.102*500)]),
                                "_tau_2_5": np.array([(i/500) for i in np.arange(0.*500, 0.072*500)]),
                                
                                "_tau21": np.array([(i/500) for i in np.arange(0.*500, 1.202*500)]),#for one-pass kT minimization as per CMS
                                "_tau32": np.array([(i/500) for i in np.arange(0.*500, 1.202*500)]),#for one-pass kT minimization as per CMS
                                
                                "_tau21_WTA": np.array([(i/500) for i in np.arange(0.*500, 1.102*501)]),#for WTA-kT for comparison
                                "_tau32_WTA": np.array([(i/500) for i in np.arange(0.*500, 1.102*501)]),#for WTA-kT for comparison
                                
                                "_tau21_exkT": np.array([(i/500) for i in np.arange(0.*500, 1.602*500)]),#for excl.-kT and E-scheme as per basis
                                "_tau32_exkT": np.array([(i/500) for i in np.arange(0.*500, 1.602*500)]),#for excl.-kT and E-scheme as per basis

                               }
        self.kinematic_labels = ['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels = ['_good_nPVs']
        #self.extra_labels = ['deltaPhi', ]

        self.dict_variables_kinematics = {
                                            "_pt": np.array([i for i in np.arange(70., 3570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 455., 5.)]),
                                            "_msoftdrop": np.array([i for i in np.arange(0., 455., 5.)]),
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
        

        if not(self.onlyParquet): output = self.buildDictOfHistograms()
        branches = self._branchesToRead
        events = events[branches]
        
        if self.verbose and not(self.onlyParquet):
            print (f'Now, processing histograms; systematic sources being considered:{self.sysSource} for nevents={len(events)} (before masking)')
        elif self.verbose and (self.onlyParquet):
            print (f'Now, producing .parquet files; for nevents={len(events)} (before masking)')
            
        for sys in self.sysSource:

            ##############################
            #### Building event masks ####
            ############################## 
            if self.isMC and self.isSigMC:# or self.isaltSigMC:
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys
                #print(s)
                

                #if not self.isMC: selRecoMask = (events[f'totalRecoWeight{s}']!=0) & (events[f'passRecoSel{s}']!=0) & (events[f'recoSelectedEventNumber{s}']>=1) #all reco events that are flagged as non-negative should be accepted, >=1 to deal with the potential edge case that first event considered by skimmer isn't accepted in a given (set of) file(s) and stored as 0.
                #else:
                selRecoMask = (events[f'totalRecoWeight{s}']!=0) & (events[f'passRecoSel{s}']==1) 
                selGenMask = (events[f'evtGenWeight_nom']!=0)  & (events[f'passGenSel_nom']==1)
                                                                    
                if self.withLepVeto:
                    selRecoMask = (selRecoMask) & (events[f'nRecoLeptons_nom']==0)
                    selGenMask = (selGenMask) & (events[f'nGenLeptons_nom']==0)
                    
                trueRecoMask = (selRecoMask) & (selGenMask) & ((events[f'trueRecoJets{s}_pt']>0.) & (events[f'trueRecoJetsF{s}_pt']>0.))
                
                fakeRecoMask = (selRecoMask) & ((events[f'trueRecoJets{s}_pt']<0.) | (events[f'trueRecoJetsF{s}_pt']<0.) ) #& (~selGenMask)
                
                accepGenMask = (selGenMask) & (selRecoMask) & ((events[f'accepGenJets{s}_pt']>0.) & (events[f'accepGenJetsF{s}_pt']>0.))
                
                missGenMask =  (selGenMask)  & ((events[f'accepGenJetsF{s}_pt']<0.) | (events[f'accepGenJets{s}_pt']<0.) ) #& (~selRecoMask)
                if (self.verbose and (sys.endswith('nom') or self.sysUnc)) or self.onlyParquet: 
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
                        print(f"Saving .parquet files with file-stem: {self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco/gen/truereco/..._{self.splitCount}.parquet'}")#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'
                        
                        ak.to_parquet(events[selRecoMask],self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco_{self.splitCount}.parquet')#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/recoMasked/'+self.year+'/'
                        ak.to_parquet(events[selGenMask],self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_gen_{self.splitCount}.parquet')#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/genMasked/'+self.year+'/'
                        #ak.to_parquet(events[trueRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_truereco_{self.splitCount}.parquet')
                        #ak.to_parquet(events[accepGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_accepgen_{self.splitCount}.parquet')
                        #ak.to_parquet(events[fakeRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_fakereco_{self.splitCount}.parquet')
                        #ak.to_parquet(events[missGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_{self.jetType}Jet_missgen_{self.splitCount}.parquet')
                        
                        return 1 # dummy 
                    
                
            elif self.isMC and (not self.isSigMC) and sys.endswith('nom'): #not so relevant for dijet histogramming but for background MC's in W/top or simple parquet production in dijets
                
                
                selRecoMask = (events[f'totalRecoWeight{sys}']!=0.) & (events[f'passRecoSel{sys}']==1) 
                selGenMask = (events[f'evtGenWeight{sys}']!=0.) & (events[f'passGenSel{sys}']==1)
                
                if self.withLepVeto:
                    selRecoMask = (selRecoMask) & (events[f'nRecoLeptons_nom']==0)
                    selGenMask = (selGenMask) & (events[f'nGenLeptons_nom']==0)
                '''
                if self.saveParquet: 
                    #events['selRecoMask']=selRecoMask
                    #events['selGenMask']=selGenMask

                    if (sys.endswith('nom')): 
                        print(f"Saving .parquet files with file-stem: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco/gen/truereco/..._{self.splitCount}'}")
                        
                        ak.to_parquet(events[selRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco_{self.splitCount}.parquet')
                        ak.to_parquet(events[selGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_gen_{self.splitCount}.parquet')
                    
                    return 1 
                '''
                if self.saveParquet: 
                    
                    s= '_nom'
                    
                    trueRecoMask = (selRecoMask) & (selGenMask) & ((events[f'trueRecoJets{s}_pt']>0.) & (events[f'trueRecoJetsF{s}_pt']>0.))
                
                    fakeRecoMask = (selRecoMask) & ((events[f'trueRecoJets{s}_pt']<0.) | (events[f'trueRecoJetsF{s}_pt']<0.) ) #& (~selGenMask)

                    accepGenMask = (selGenMask) & (selRecoMask) & ((events[f'accepGenJets{s}_pt']>0.) & (events[f'accepGenJetsF{s}_pt']>0.))

                    missGenMask =  (selGenMask)  & ((events[f'accepGenJetsF{s}_pt']<0.) | (events[f'accepGenJets{s}_pt']<0.) ) #& (~selRecoMask)
                    if (self.verbose and (sys.endswith('nom') or self.sysUnc)) or self.onlyParquet: 
                        print('#### Building event masks ####')
                        print(sys,s, 'masked array lengths for all, reco,true,fake,gen,accep,miss', len(events),
                              len(events[selRecoMask]),
                              len(events[trueRecoMask]),
                              len(events[fakeRecoMask]),
                              len(events[selGenMask]),
                              len(events[accepGenMask]),
                              len(events[missGenMask]),
                              self.splitCount
                             )

                    events['trueRecoMask']=trueRecoMask
                    events['selRecoMask']=selRecoMask
                    events['fakeRecoMask']=fakeRecoMask
                    events['selGenMask']=selGenMask
                    events['accepGenMask']=accepGenMask
                    events['missGenMask']=missGenMask
                    

                    if (sys.endswith('nom')): 
                        print(f"Saving .parquet files with file-stem: {self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco/gen/truereco/..._{self.splitCount}.parquet'}")#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'
                        
                        ak.to_parquet(events[selRecoMask],self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco_{self.splitCount}.parquet')#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/recoMasked/'+self.year+'/'
                        ak.to_parquet(events[selGenMask],self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_gen_{self.splitCount}.parquet')#self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/genMasked/'+self.year+'/'
                        
                        return 1 # dummy 
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                selRecoMasks = OrderedDict()
                
                
                print("Building trigger masks using the following triggers:")
                
                print(self.triggerTable.keys())
                
                for itrigger, itrValues in self.triggerTable.items():    
                
                    triggerList=list(self.triggerTable.keys())
                    thistrigInd=triggerList.index(itrigger)
                    othertrigsInds=[i for i in range(0,len(triggerList)) if i!=thistrigInd]
                    
                    selRecoMasks[itrigger]=(events[f'totalRecoWeight{sys}']!=0.) & (events[f'passRecoSel{sys}']!=0) & (events[f'recoSelectedEventNumber{sys}']>=1) & (events[f'passHLT_{itrigger}']==1 ) & ((events[f'selRecoJets{sys}_pt']>itrValues[self.year][0]) | (events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0])) & ((events[f'selRecoJets{sys}_pt']  < itrValues[self.year][1]) & (events[f'selRecoJetsF{sys}_pt']  < itrValues[self.year][1]))
                    
                    if self.withLepVeto:
                        selRecoMasks[itrigger] = (selRecoMasks[itrigger]) & (events[f'nRecoLeptons{sys}']==0)
                        
                    

                    if self.saveParquet: 
                        
                        
                        for itrigger, itrValues in self.triggerTable.items():
                            events['selRecoMask_'+itrigger]=selRecoMasks[itrigger]


                        print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_nomWts_{self.jetType}Jet_{self.splitCount}.parquet'}")
                        ak.to_parquet(events, (self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Dijets_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_nomWts_{self.jetType}Jet_{self.splitCount}.parquet'))    
                        return 1
                    

            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            
            for isel in self.selList:
                
                
                listOfOutputHistosNames = [k for k,h in output.items() if ((sys in k) or ('residual' in k and not(sys in k)) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                
                for k in listOfOutputHistosNames:
                    key=k
                    
                    ############### Safety checks ##################
                    if not(sys in key) and not('residual' in key) and not('resol' in key): 
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
                            

                        elif ('residual' in key.lower() or 'resol' in key.lower()) and sys.endswith('nom') and self.isSigMC:
                            zeroMask=(events[f'accepGenJets{self.jetFlag}{s}{varToFill}']!=0.)&(accepGenMask)

                            response = events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]/events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][zeroMask]
                            response = np.nan_to_num(response,nan=-999.)
                            residual = events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]-events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][zeroMask]
                            
                            if 'noWt_' in key: output[key].fill(response)#, weight=totalRecoWeight[zeroMask])
                            elif 'residual' in key: output[key].fill(residual)
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
                iJ=self.nJet[0] 
                if self.verbose: print(iJ,self.jetFlag,itype)
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
                                
                    if itype.startswith('truereco') and self.isSigMC and sysUnc.endswith('nom'):
                        dictOfHists['residual'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_pt'+isel, label='AK8 reco pt - gen jet pt', underflow=True,overflow=True).Weight())
                        dictOfHists['residual'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD} - gen jet m_{SD}', underflow=True,overflow=True).Weight())
                        dictOfHists['residual'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_mass'+isel, label='AK8 reco inv. m - gen jet inv. m', underflow=True,overflow=True).Weight())
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
                            if sysUnc.endswith('_nom') and self.isSigMC:
                                dictOfHists['residual'+iJ+x+isel] = (hist.Hist.new.Regular(6000, -1.5, 1.5, name='residual'+iJ+x+isel, label='AK8 reco #tau - gen jet #tau', underflow=True,overflow=True).Weight())
                                dictOfHists['noWt_resol'+iJ+x+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
        
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
                    if not('passRecoSel'+sys in reco_list): reco_list.append('passRecoSel'+sys)
                
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
                    if not('passRecoSel'+sys in reco_list): reco_list.append('passRecoSel'+sys)
                
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
                if not('passGenSel'+sys in gen_list) and 'nom' in sys: gen_list.append('passGenSel'+sys)              
                if not('passRecoSel'+sys in reco_list): reco_list.append('passRecoSel'+sys)
                if sys.endswith('nom'): reco_list.append("good_nPVs"+sys)
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("good_nPVs"+sys)
                if not('passRecoSel'+sys in reco_list): reco_list.append('passRecoSel'+sys)
                for itrigger in self.triggerTable.keys():    
                    if not(self.year.endswith('VFP') and self.era=='B'): triggerBit_list.append(f'HLT_{itrigger}') #to remove
                    triggerBit_list.append(f'passHLT_{itrigger}')
                
                if '2018' in self.year or '2017' in self.year: 
                    triggerBit_list.append(f'passHLT_AK8PFJet550')
                    triggerBit_list.append(f'HLT_AK8PFJet550')
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or self.sysUnc):
                if 'pu' in sys or 'l1' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list+['nRecoLeptons_nom','nGenLeptons_nom']#+['recoSelectedEventNumber_nom']
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list+['nRecoLeptons_nom','nGenLeptons_nom']#+['recoSelectedEventNumber_nom']
        elif (not self.isMC): branchesToRead=reco_list+triggerBit_list+['recoSelectedEventNumber_nom','nRecoLeptons_nom']
        if self.verbose: 
            print("Prepared branches to read",branchesToRead)        
        return branchesToRead
    


def histoMaker(myProcessor, sampleIdentifier='qcd_ht', y='2017', sampleDict_PFNano=OrderedDict(),
               sampleDict_local=OrderedDict(), writeChunks=True, writeAccumulatedORhadd=1, verbose=False, 
               saveParquet=False,
               isSigMC=True, isMC=True, wtUnc=True, sysUnc=False, onlyUnc='', outputdir='processorTests/', 
               sysSource=[],
               era='', ext='_nomWts', splitchunks=10, nWorkers=40, stepSize="2048 MB",
               jetType='Central',forceProduction=False,onlyParquet=False):
    
    nchunk=copy.deepcopy(splitchunks)
    cz=0

    if not os.path.exists(f'{outputdir}'): os.makedirs(f'{outputdir}')
    
    for sample in sampleDict_local.keys():
        yl=['2016_preVFP', '2016', '2017', '2018'] if y=='all' else [y]
        
        
        for year in yl:
            gc.collect()
            pathExistsFlag=False
            if not( sample.lower().startswith(sampleIdentifier.lower())):# or ('170to300' in sample.lower() or jetType=='Central'):
                #print(sample,year)
                continue
                
            tstart0=time.time()

            print(f'Histogramming for {sample} in year: {year}')
            if verbose: print(f'using nanoskims from {sampleDict_local[sample][year]["t3_dirs"]}{chr(10)}','\n')
            
            

            inpdir=sampleDict_local[sample][year]['t3_dirs']
            dirfnames=inpdir[0].split('/0000/')[0]+'/000*/jetObservables_nanoskim_*.root'

            SC='0'
            if verbose: print (sysSource)
            my_processor=myProcessor(sampleName=sample, sampleDict=sampleDict_PFNano,
                                     isMC=isMC, isSigMC=isSigMC,year=year,
                                     saveParquet=saveParquet,onlyParquet=onlyParquet,
                                     sysSource=sysSource,
                                     wtUnc=wtUnc,sysUnc=sysUnc,
                                     onlyUnc='' if 'jes' in onlyUnc else onlyUnc,era=era,
                                     verbose=False,splitCount=SC,jetType=jetType )
            if cz==0:
                if verbose: print("Branches being read \n", my_processor._branchesToRead)
                cz=cz+1
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
            
            
            for i in ifls:#()track(ifls):#
                if c%100==0: print (c, sample, year, jetType)
                if 'jes' in onlyUnc and c%5==0: print (c, sample, year)
                if splitchunks>0:
                    ns = [int(x.split('jetObservables_nanoskim_')[1].split('.root')[0]) for x in i]
                    s = str(ns).split('[')[1].split(']')[0]
                    stringfnames=i[0].split('jetObservables_nanoskim_')[0]+"jetObservables_nanoskim_{"+f'{s}'+"}.root"
                else: stringfnames=i
                    
                stringfnames=stringfnames.replace(' ','')
                
                if verbose:
                    tstart1 = time.time()
                
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
                    
                if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):
                    fn=f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root'
                    
                if os.path.exists(fn) and not(forceProduction):
                    print("WARNING: file already  exists; recheck this sample")
                    #if not(forceProduction):
                    print("WARNING: skipping this sample")
                    pathExistsFlag=True
                    break
                    
                else:
                
                    events = uproot.concatenate(stringfnames+':Events', my_processor._branchesToRead, 
                                                step_size=stepSize, library='ak', num_workers=nWorkers,

                                               )
                    if not(isMC):
                        events_df = ak.to_pandas(copy.deepcopy(events))
                        print(len(events_df))
                        del(events)
                        events_df_nodup = events_df.drop_duplicates()#subset='recoSelectedEventNumber_nom',keep='first'----> not using this since it seems in this iteration the reco event number branch wasn't properly updated in skimmer logic
                        print(f"Dropped {len(events_df)-len(events_df_nodup)} duplicated events; left with {len(events_df_nodup)} events in this sample")

                        f = events_df_nodup.to_parquet(f'tempQCD_{year}_{era}.parquet')
                        events = ak.from_parquet(f'tempQCD_{year}_{era}.parquet')

                    if verbose:
                        elapsed1 = time.time()-tstart1
                        print(f"Time taken to load {len(events)} events from {n_tosplit[c]} files in {sample} = {elapsed1}")

                    if verbose: 
                        print(f'nEvents in this file chunk: {len(events)}')

                    if verbose: 
                        tstart2=time.time()

                    if verbose: print("#####Processing events#####")

                    processed_events=my_processor.process(events)

                    if verbose:
                        elapsed2 = time.time()-tstart2
                        print(f"Time taken to build histos from {len(events)} events, for {fnstem,jetType} = {elapsed2}")

                    my_processor.splitCount=chr(int(SC)+1)

                    if c==0: 
                        hists=processed_events
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                            elif jetType=="Forward":
                                print(f"Current integral of pT nominal hist: {hists['recoJetF_pt_nom_dijetSel'].sum(flow=True)}")
                    elif c>0 and writeAccumulatedORhadd==0: 
                        hists=accumulate([hists,processed_events])
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                            elif jetType=="Forward":
                                print(f"Current integral of pT nominal hist: {hists['recoJetF_pt_nom_dijetSel'].sum(flow=True)}")
                    
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
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                            elif jetType=="Forward":
                                print(f"Current integral of pT nominal hist: {hists['recoJetF_pt_nom_dijetSel'].sum(flow=True)}")

                        #print (processed_events)
                        print ("Events proccesed, now writing output file(s)",fn,c,len(i),len(ifls))

                        with uproot.recreate(fn) as fout:#outputTest_{sample.split("_Tune")[0]}_{year}.root'
                            for key in processed_events.keys():
                                fout[key]=processed_events[key]
                            print (f'Done with creating output file:{fn}')
                            fout.close()
                    elif writeAccumulatedORhadd==0:# or len(n_tosplit)==1:
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of gen pT nominal hist: {hists['genJet_pt_nom_dijetSel'].sum(flow=True)}")
                            elif jetType=="Forward":
                                print(f"Current integral of pT nominal hist: {hists['genJetF_pt_nom_dijetSel'].sum(flow=True)}")

                        #print (processed_events)
                        #print ("Events proccesed, now writing output file(s)",fn)   
                    c=c+1
                    del(processed_events)
                    del(events)
                    gc.collect()
            if not(pathExistsFlag):
                if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):

                    if verbose: 
                        if jetType=='Central':
                            print(f"Final integral of reco pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                            print(f"Final integral of gen pT nominal hist: {hists['genJet_pt_nom_dijetSel'].sum(flow=True)}")

                        elif jetType=="Forward":
                            print(f"Final integral of reco pT nominal hist: {hists['recoJetF_pt_nom_dijetSel'].sum(flow=True)}")
                            print(f"Final integral of gen pT nominal hist: {hists['genJetF_pt_nom_dijetSel'].sum(flow=True)}")

                    if not os.path.exists(f'{outputdir}/combinedJES/'): os.makedirs(f'{outputdir}/combinedJES/')

                    print (f'Output file stem for accumulated histos: {fnstem}')
                    with uproot.recreate(f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root') as fout_hadd:
                        for key in hists.keys():
                            fout_hadd[key]=hists[key]
                        print (f'Done with creating accumulated output file:{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root')
                        fout_hadd.close()


                elapsed0 = time.time() - tstart0

                print (f'Time for {sample} {year}:{elapsed0}')        

                del(my_processor)
                del(hists)
            gc.collect()
           
    return 1




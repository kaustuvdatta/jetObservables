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
                 sampleDict=dictSamples,test=False, sysUnc=False, splitCount='0', jetType='Central', 
                 vetoMap=None,
                 applyVetoMap=False,
                 trigTest=False, 
                 onlyRedonePDFandAlphaSWts=True,#False
                 withAllPdfVariations=False,#False
                 withAlphaSVariations=True,#False
                 trigUpDownVal=10, withLepVeto=False, parquetDir='/scratch/kadatta/dijetChecks/parquets/',
                 minLeadPt=200., minSubLeadPt=200., parquetExt='',onlyControlHistos=False
                ):
        self.onlyControlHistos = onlyControlHistos
        self.parquetExt = parquetExt
        self.test=test
        self.jetType=jetType
        self.year = year
        self.isMC = isMC
        self.isSigMC = isSigMC
        self.minLeadingJetPt=minLeadPt
        self.minSubLeadingJetPt=minSubLeadPt
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
        self.onlyRedonePDFandAlphaSWts = onlyRedonePDFandAlphaSWts
        self.withAlphaSVariations = withAlphaSVariations
        self.withAllPdfVariations = withAllPdfVariations
        self.events = None            
        if (not self.isMC) and self.era=='': 
            
            print (f'Data-loading error: You need to specify what era if you want to work with data' )
        
        
        ### Helpers
        if not(self.sysUnc):
            self.listOfHistTypes =  [ 'gen', 'accepgen', 'missgen', 'reco', 'fakereco', 'truereco',  ] if self.isMC else [ 'reco', 'recoLeading', 'recoSubleading' ] #'genLeading', 'genSubleading', 'genSub2leading', #'recoLeading', 'recoSubleading', 'recoSub2leading',
        elif self.sysUnc and self.isSigMC:
            self.listOfHistTypes =  [ 'gen', 'accepgen', 'missgen', 'reco', 'fakereco', 'truereco',  ] 
        
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        #self.outputDir = 'UL17and18_nano_tests/'
        
        self.nJet = [ 'Jet'] if 'Central' in self.jetType else ['JetF']
        self.jetFlag = self.nJet[0].split('Jet')[1]
        if self.verbose: print (f"Initialising for {self.nJet}/{jetType}/{self.jetFlag}")


        self.triggerTable = OrderedDict()
        
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
            
            
        self.applyVetoMap = applyVetoMap
        if self.applyVetoMap:
            
            if not(vetoMap==None):
                vetoMap_cset = correctionlib.CorrectionSet.from_file(vetoMap)
                if self.year.startswith('2016'):
                    mapName = 'Summer19UL16_V1'
                elif self.year=='2017':
                    mapName = 'Summer19UL17_V1'
                elif self.year=='2018':
                    mapName = 'Summer19UL18_V1'
                self.vetoMap = vetoMap_cset[mapName] 
            else:
                raise ValueError("jet veto map JSON is not provided!")

        self.recoMask = None
        self.recoWeights = None
        if self.isMC:# and self.onlyParquet:
            self.genMask = None
            self.accepGenMask = None
            self.trueRecoMask = None
            self.genWeights = None
            
        if not self.onlyControlHistos:
            self.dict_variables_toUnfold = {

                                    #"_pt": np.array([i for i in np.arange(70., 3570., 10.)]),
                                    #"_mass": np.array([i for i in np.arange(0., 455., 5.)]),

                                    "_tau_0p25_1": np.array([(i/200) for i in np.arange(0.*200, 1.005*200)]),
                                    "_tau_0p25_2": np.array([(i/500) for i in np.arange(0.*500, 0.952*500)]),
                                    "_tau_0p25_3": np.array([(i/500) for i in np.arange(0.*500, 0.902*500)]),
                                    "_tau_0p25_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.901*1000)]),
                                    "_tau_0p25_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.851*1000)]),

                                    "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.*200, 1.005*200)]),
                                    "_tau_0p5_2": np.array([(i/500) for i in np.arange(0.*500, 0.922*500)]),
                                    "_tau_0p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.852*500)]),
                                    "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.801*1000)]),
                                    "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.751*1000)]),

                                    "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.925*200)]),
                                    "_tau_1_2": np.array([(i/500) for i in np.arange(0.*500, 0.602*500)]),
                                    "_tau_1_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.451*1000)]),
                                    "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.401*1000)]),
                                    "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.351*1000)]),

                                    "_tau_1p5_1": np.array([(i/500) for i in np.arange(0.*500, 0.702*500)]),
                                    "_tau_1p5_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.421*1000)]),
                                    "_tau_1p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.351*1000)]),
                                    "_tau_1p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.301*1000)]),
                                    "_tau_1p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.281*1000)]),


                                    "_tau_2_1": np.array([(i/500) for i in np.arange(0.*500, 0.702*500)]),
                                    "_tau_2_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.421*1000)]),
                                    "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.251*1000)]),
                                    "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.201*1000)]),
                                    "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.141*1000)]),

                                    "_tau21": np.array([(i/500) for i in np.arange(0.*500, 1.202*500)]),#for one-pass kT minimization as per CMS
                                    "_tau32": np.array([(i/500) for i in np.arange(0.*500, 1.202*500)]),#for one-pass kT minimization as per CMS

                                    "_tau21_WTA": np.array([(i/500) for i in np.arange(0.*500, 1.102*501)]),#for WTA-kT for comparison
                                    "_tau32_WTA": np.array([(i/500) for i in np.arange(0.*500, 1.102*501)]),#for WTA-kT for comparison

                                    "_tau21_exkT": np.array([(i/500) for i in np.arange(0.*500, 1.602*500)]),#for excl.-kT and E-scheme as per basis
                                    "_tau32_exkT": np.array([(i/500) for i in np.arange(0.*500, 1.602*500)]),#for excl.-kT and E-scheme as per basis

                                   }
        else:
            self.dict_variables_toUnfold = {

                                   
                                    "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.*200, 1.005*200)]),
        
                                    }

        self.kinematic_labels = ['_pt','_eta', '_y', '_phi', '_mass']#, '_msoftdrop_new']
        self.reco_only_labels = ['_good_nPVs']#, '_JERfactor', '_JECfactor', '_pt_raw', ]#'_HT',
        #self.extra_labels = ['deltaPhi', ]

        self.dict_variables_kinematics = {
                                            "_pt": np.array([i for i in np.arange(70., 3570., 10.)]),
                                            #"_pt_raw": np.array([i for i in np.arange(70., 3570., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.42, 0.02)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.42, 0.02)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.4, 0.02)]),
                                            
                                            #"_JERfactor": np.array([i for i in np.arange(-2., 8.1, 0.01)]),
                                            #"_JECfactor": np.array([i for i in np.arange(-2., 8.1, 0.01)]),
            
                                            "_mass": np.array([i for i in np.arange(0., 555., 5.)]),
                                            #"_msoftdrop_new": np.array([i for i in np.arange(0., 555., 5.)]),
                                            "_good_nPVs": np.array([i for i in np.arange(0., 101., 1.)]),
                                            #"_HT" : np.array([i for i in np.arange(0., 3570., 10.)]),
                                            
                                         }
        
        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not( isys.endswith('nom'))]# or self.sysUnc) ]
        
        self.sysWeightList = ( '_pu', '_pdf', '_isr', '_fsr', '_l1prefiring' ) #'_ps',
        self.constJESList = ( '_constituentJES_neutral', '_constituentJES_charged', '_constituentJES_photon' )
        if not self.onlyRedonePDFandAlphaSWts:
            self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight'] if self.wtUnc else [] 
        else:
            self.wtSources=['_pdfWeight'] if self.wtUnc else [] 
            

        self.recoWtSources=['_pu', '_l1'] if self.isMC else [] 
        
        if self.onlyUnc!='' : self.sysSource = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] 

        if self.verbose: print('In __init__',self.sysSource)

        self.selList = '_dijetSel' 
        
        self._branchesToRead = self.getBranchesToRead(
                                                         dirname=self.inputDir, year=self.year,
                                                         kinematic_labels=self.kinematic_labels,
                                                         reco_only_labels=self.reco_only_labels,
                                                         nSub_labels=self.dict_variables_toUnfold,
                                                     )
    #def add_vetoMapOutput_to_recoWeights(self):
    #    #if self.apply
    #    tstart=time.time()
    #    events_jetVeto_fields = self.append_recoWeights_forJetVeto()
    #    if self.verbose: 
    #        print(f"Changing 'totalRecoWeight_nom/sysUpDown' reco weight in branches to easily add in info from jet veto maps to data/MC response matrices for unfolding in sig MC without changing genWeights in MC")
    #    
    #    for key, values in events_jetVeto_fields.items():
    #        events[key] = values
    #    elapsed = time.time() - tstart
    #    if self.verbose: print (f'Finished adding in jet veto info. Time taken:{elapsed}') 
    
    
    def add_new_PDF_and_AlphaS_wts(self):
        #using implementation by S. Rothman https://github.com/ssrothman/EECpostprocessing/blob/566c077de831613f1cbfeec9179acf9fcef00da5/selections/theorySF.py#L90C5-L90C18
        if self.onlyRedonePDFandAlphaSWts:
            pdf_weights = self.events['pdfWeightAll_nom']
            nevt = len(self.events['pdfWeightAll_nom'])
            pdf_nom = np.ones(len(self.events['pdfWeightUp_nom']))
            #self.wtSources[0] = 'pdfWeight2'
            # Hessian weights a la Eq. 21 in https://arxiv.org/pdf/1510.03865v1.pdf
            to_sum = pdf_weights[:,1:-2]-np.ones((nevt,100))
            summed_up = ak.sum(np.square(to_sum),axis=1)
            pdfWt_unc = np.sqrt( (1./99.) * summed_up )
            pdfWt_up = pdfWt_unc + pdf_nom
            pdfWt_dn = pdf_nom - pdfWt_unc
            #weights.add('wt_PDF', nom, pdf_up, pdf_dn)
            self.events['pdfWeightUp_nom'] = pdfWt_up
            self.events['pdfWeightDown_nom'] = pdfWt_dn
            
            
            if self.withAllPdfVariations:
                self.sysWeightList=( '_pdf', )#, '_aS', '_alphaSandPDF')#, '_isr', '_fsr', '_l1prefiring' )

                #don't do scale variations in same run
                for i in range(1,101):
                    self.events[f'{i}pdfWeight_nom'] =  pdf_weights[:,i]
                    #self.wtSources.append(f'_pdf{i}Weight_nom')
                    self.sysWeightList=tuple([wt for wt in list(self.sysWeightList)]+[f'_{i}pdfWeight'])
                    self.sysSource+=[f'_{i}pdfWeight']
                    
            if self.withAlphaSVariations:
                
                # alpha_S weights; Eq. 27 of https://arxiv.org/pdf/1510.03865v1.pdf
                alphaS_unc = 0.5*(pdf_weights[:,102] - pdf_weights[:,101])
                self.events['alphaSWeightUp_nom'] = pdf_nom+alphaS_unc
                self.events['alphaSWeightDown_nom'] = pdf_nom-alphaS_unc

                # PDF + alpha_S weights; Eq. 28
                pdf_and_alphaS_unc = np.sqrt( np.square(pdfWt_unc) + np.square(alphaS_unc) )
                self.events['aSandPDFWeightUp_nom'] = pdf_nom+pdf_and_alphaS_unc
                self.events['aSandPDFWeightDown_nom'] = pdf_nom-pdf_and_alphaS_unc

                self.wtSources+=['_alphaSWeight','_aSandPDFWeight']
                self.sysWeightList=(  '_pdf', '_alphaS', '_aSandPDF')# '_pu', '_isr', '_fsr', '_l1prefiring' )
                self.sysSource += [ iwt+i for i in [ 'Up', 'Down' ] for iwt in self.wtSources if not (iwt+i in self.sysSource) and not( iwt.endswith(('nom','pdfWeightAll'))) ]
            
        
    def process(self, events):
        '''fill Hist histograms to accumulate over chunks of processed files, convert to root or whatever else after returned by processor'''
        
        
        # Get branches to be read 
        branches = self._branchesToRead
        
        self.events=events               

        #hard-coding since only parquet prod or histo prod modes in use, so onlyParquet flag is redundant 
        self.onlyParquet = self.saveParquet

        ################## jet vetoes ##################
                
        #if self.applyVetoMap:
        #
        #    self.add_vetoMapOutput_to_recoWeights()#recoMask=selRecoMask)
        if self.isSigMC and self.wtUnc:
            if self.onlyRedonePDFandAlphaSWts or self.withAllPdfVariations:
                self.add_new_PDF_and_AlphaS_wts()
           
        
        if not(self.onlyParquet): 
            output = self.buildDictOfHistograms()
            if self.verbose:
                print (f'Now, processing histograms; systematic sources being considered:{self.sysSource} for nevents={len(events)} (before masking)')
            if self.verbose: 
                print (output.keys())

        elif self.verbose and (self.onlyParquet):
            print (f'Now, producing .parquet files; for nevents={len(events)} (before masking)')

            
        for sys in self.sysSource:
            
            if not(self.onlyRedonePDFandAlphaSWts):

                s='_nom' if (sys.startswith(self.sysWeightList) or 'const' in sys) or not(self.isSigMC) else sys
            
            else:
                
                s='_nom' if (sys.startswith(('_isr','_fsr','_pdf','_a')) or 'pdf' in sys or 'const' in sys) or not(self.isSigMC) else sys
                #print(sys,s)
                
            
           
            #########################################################################
            ################## Build/modify event weights ###########################
            #########################################################################        
            
            if self.isMC: 
                #print(sys,s,sum(self.events[f'evtGenWeight{s}']),len(self.events[f'evtGenWeight{s}'][self.events[f'evtGenWeight{s}']!=0]),sum(self.events[f'totalRecoWeight{s}']),len(self.events[f'totalRecoWeight{s}'][self.events[f'totalRecoWeight{s}']!=0]))
                #handles default nom weight and jes/jer reweightings also for sigMC
                self.genWeights = self.events[f'evtGenWeight{s}']
                self.recoWeights = self.genWeights * self.events[f'puWeightNom{s}'] * self.events[f'l1prefiringWeightNom{s}']
                
                #specifically handle weight variations on sigMC and modify above clas-level variables accordingly
                if self.isSigMC:
                    if sys.startswith(('_isr','_fsr','_pdf', '_a')) or 'pdf' in sys.lower():
                        extraWeights = self.events[f'{sys.split("_")[1]}{s}']
                        #print(f'{sys.split("_")[1]}{s}',sum(extraWeights))
                        self.genWeights = self.events[f'evtGenWeight{s}']*extraWeights          
                        self.recoWeights = self.genWeights * self.events[f'puWeightNom{s}'] * self.events[f'l1prefiringWeightNom{s}']
                    #print(f'{sys.split("_")[1]}{s}')
                    elif 'pu' in sys:
                        self.recoWeights = self.genWeights*self.events[f'{sys.split("_")[1]}{s}']*self.events[f'l1prefiringWeightNom{s}']
                    elif 'l1' in sys:                      
                        self.recoWeights = self.genWeights*self.events[f'{sys.split("_")[1]}{s}']*self.events[f'puWeightNom{s}']
                    elif 'const' in sys:
                        self.recoWeights = self.genWeights*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']     
                
                #modify total recoWeights branch to include effects of reweighitng on event yield for MC, ie, particularly for veto maps, by cutting on totalRecoweight{s}!=0; this in sigMC also propagates to accepgen counts naturally without further hard-coding below.
                self.events[f'totalRecoWeight{s}'] = self.recoWeights
                
            else:
                #data
                self.recoWeights = self.events[f'totalRecoWeight_nom']#*(self.events[f'new{self.sel}_eventJetVeto_nom'] if self.applyVetoMap else np.ones(len(self.events[f'passRecoSel_nom']))) 
              
            
            #if self.isMC: 
            #    #print(sys,s,sum(self.genWeights), len(self.genWeights[self.genWeights!=0]), sum(self.recoWeights), len(self.recoWeights[self.recoWeights!=0]))
            
            totalRecoWeight = self.recoWeights
            totalGenWeight = self.genWeights if self.isMC else None


            ##############################
            #### Building event masks ####
            ############################## 

            selRecoMask =  ( (self.events[f'passRecoSel{s}']==1) & 
                                 ( ((self.events[f'selRecoJets{s}_pt']>self.minLeadingJetPt) & (self.events[f'selRecoJetsF{s}_pt']>self.minSubLeadingJetPt)) |
                                   ((self.events[f'selRecoJetsF{s}_pt']>self.minLeadingJetPt) & (self.events[f'selRecoJets{s}_pt']>self.minSubLeadingJetPt)) 
                                 )
                                )  #(self.events[f'totalRecoWeight{s}']!=0.) &#include if jet vetoes are used

            if self.isMC:

                selGenMask = ( (self.events[f'passGenSel{s}']==1) & 
                               ( ((self.events[f'selGenJets_nom_pt']>self.minLeadingJetPt) & (self.events[f'selGenJetsF_nom_pt']>self.minSubLeadingJetPt)) |
                                 ((self.events[f'selGenJetsF_nom_pt']>self.minLeadingJetPt) & (self.events[f'selGenJets_nom_pt']>self.minSubLeadingJetPt)) 
                               )
                             )  #(self.events[f'evtGenWeight{s}']!=0.) &

            if self.withLepVeto:
                selRecoMask = (selRecoMask) & (self.events[f'nRecoLeptons_nom']==0)
                if self.isMC: selGenMask = (selGenMask) & (self.events[f'nGenLeptons_nom']==0)

            if self.isSigMC:
                
                #s='_nom' if (sys.startswith(self.sysWeightList)) or not(self.isSigMC) else sys
                if not(self.onlyRedonePDFandAlphaSWts):

                    s='_nom' if (sys.startswith(self.sysWeightList) or 'const' in sys) or not(self.isSigMC) else sys

                else:

                    s='_nom' if (sys.startswith(('_isr','_fsr','_pdf','_a')) or 'pdf' in sys or 'const' in sys) or not(self.isSigMC) else sys
                    #print(sys,s)
                    
                self.recoMask = selRecoMask
                self.genMask = selGenMask
                
                trueRecoMask = (selRecoMask) & (selGenMask) & (self.events[f'trueRecoJets{self.jetFlag}{s}_pt']>0.) #& (self.events[f'accepGenJets{self.jetFlag}{s}_pt']>0.)
                accepGenMask = (selGenMask) & (selRecoMask) & (self.events[f'accepGenJets{self.jetFlag}{s}_pt']>0.) #& (self.events[f'trueRecoJets{self.jetFlag}{s}_pt']>0.)
                
                self.accepGenMask = accepGenMask
                self.trueRecoMask = trueRecoMask

                fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) 
                
                missGenMask =  ((selGenMask) & (~accepGenMask)) 
                """
                if ('nom' in sys) or ('jer' in sys):
                    #print(sum(self.events['genWeight']))
                    print(sum(self.events['genWeight'][selRecoMask | selGenMask]), sum(self.events['genWeight'][selRecoMask]), sum(self.events['genWeight'][selGenMask]), )
                    print(sum(self.events['genWeight'][selRecoMask | selGenMask]), sum(self.events['genWeight'][trueRecoMask]), sum(self.events['genWeight'][accepGenMask]), )
                    print(sum(self.events['genWeight'][selRecoMask | selGenMask]), sum(self.events['genWeight'][fakeRecoMask]), sum(self.events['genWeight'][missGenMask]), )
                    
                    
                    print(sys,len(self.events[selRecoMask]),len(self.events[selGenMask]),sum(self.events[f'evtGenWeight{s}'][selRecoMask]),sum(self.events[f'evtGenWeight{s}'][selGenMask]),sum(self.events[f'totalRecoWeight{s}'][selRecoMask]),sum(self.events[f'totalRecoWeight{s}'][selGenMask]))
                    print(sys,len(self.events[trueRecoMask]),len(self.events[accepGenMask]),sum(self.events[f'evtGenWeight{s}'][trueRecoMask]),sum(self.events[f'evtGenWeight{s}'][accepGenMask]),sum(self.events[f'totalRecoWeight{s}'][trueRecoMask]),sum(self.events[f'totalRecoWeight{s}'][accepGenMask]))
                    print(sys,len(self.events[fakeRecoMask]),len(self.events[missGenMask]),sum(self.events[f'evtGenWeight{s}'][fakeRecoMask]),sum(self.events[f'evtGenWeight{s}'][missGenMask]),sum(self.events[f'totalRecoWeight{s}'][fakeRecoMask]),sum(self.events[f'totalRecoWeight{s}'][missGenMask]))
                """
                if (self.verbose and (sys.endswith('nom') or self.sysUnc)) or self.saveParquet: 
                    print('#### Building event masks ####')
                    print(sys,s, 'masked array lengths for reco,true,fake,gen,accep,miss', len(events),
                          len(self.events[selRecoMask]),
                          len(self.events[trueRecoMask]),
                          len(self.events[fakeRecoMask]),
                          len(self.events[selGenMask]),
                          len(self.events[accepGenMask]),
                          len(self.events[missGenMask]),
                          self.splitCount
                         )

                    
                if self.saveParquet and (sys.endswith('nom')): 
                        
                    self.events['trueRecoMask'] = self.trueRecoMask
                    self.events['selRecoMask'] = self.recoMask
                    self.events['fakeRecoMask'] = fakeRecoMask
                    self.events['selGenMask'] = self.genMask
                    self.events['accepGenMask'] = self.accepGenMask
                    self.events['missGenMask'] = missGenMask    

                    print(f"Saving .parquet files with file-stem: {self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_..._{self.splitCount}.parquet'}")
                    #ak.to_parquet(events,self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco_{self.splitCount}{self.parquetExt}.parquet')#[selRecoMask]

                    ak.to_parquet(events,self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_{self.splitCount}{self.parquetExt}.parquet')
                    


                    return 1 
                    
                
            elif self.isMC and not(self.isSigMC) and sys.endswith('nom'): #not so relevant for dijet histogramming but for background MC's in W/top or simple parquet production in dijets
                
                s= '_nom'
                
                self.recoMask = selRecoMask
                self.genMask = selGenMask
                
                #trueRecoMask = (selRecoMask) & (selGenMask) & (self.events[f'trueRecoJets{self.jetFlag}{s}_pt']>0.) #& (self.events[f'accepGenJets{self.jetFlag}{s}_pt']>0.)
                #accepGenMask = (selGenMask) & (selRecoMask) & (self.events[f'accepGenJets{self.jetFlag}{s}_pt']>0.) #& (self.events[f'trueRecoJets{self.jetFlag}{s}_pt']>0.)
                
                #self.accepGenMask = accepGenMask
                #self.trueRecoMask = trueRecoMask

                #fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) 
                
                #missGenMask =  ((selGenMask) & (~accepGenMask)) 
                #print(sys,len(self.events[trueRecoMask]),len(self.events[accepGenMask]))
                
                if self.verbose and self.onlyParquet: 
                    print('#### Building event masks ####')
                    print(sys,s, 'masked array lengths for all, reco,true,fake,gen,accep,miss', len(events),
                          len(self.events[selRecoMask]),
                          len(self.events[trueRecoMask]),
                          len(self.events[fakeRecoMask]),
                          len(self.events[selGenMask]),
                          len(self.events[accepGenMask]),
                          len(self.events[missGenMask]),
                          self.splitCount
                         )

                if self.saveParquet:
                        
                    self.events['trueRecoMask'] = self.trueRecoMask
                    self.events['selRecoMask'] = self.recoMask
                    self.events['fakeRecoMask'] = fakeRecoMask
                    self.events['selGenMask'] = self.genMask
                    self.events['accepGenMask'] = self.accepGenMask
                    self.events['missGenMask'] = missGenMask    
                    
                    print(f"Saving .parquet files with file-stem: {self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_..._{self.splitCount}.parquet'}")
                    #ak.to_parquet(self.events[selRecoMask],self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_reco_{self.splitCount}.parquet')

                    ak.to_parquet(self.events,self.parquetDir +self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts_dijetSel_{self.jetType}Jet_OC_{self.splitCount}{self.parquetExt}.parquet')

                    
                    return 1 
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                selRecoMasks = OrderedDict()
                self.recoMask = OrderedDict()
                
                print("Building trigger masks using the following triggers:")
                
                print(self.triggerTable.keys())
                
                for itrigger, itrValues in self.triggerTable.items():    
                
                    triggerList=list(self.triggerTable.keys())
                    thistrigInd=triggerList.index(itrigger)
                    othertrigsInds=[i for i in range(0,len(triggerList)) if i!=thistrigInd]
                    
                    selRecoMasks[itrigger]= ( (self.events[f'passRecoSel{sys}']!=0) & 
                                              (self.events[f'passHLT_{itrigger}']==1 ) & 
                                              ((self.events[f'selRecoJets{sys}_pt']>itrValues[self.year][0]) | 
                                               (self.events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0])
                                              ) & 
                                              ((self.events[f'selRecoJets{sys}_pt']  < itrValues[self.year][1]) & 
                                               (self.events[f'selRecoJetsF{sys}_pt']  < itrValues[self.year][1])
                                              ) & 
                                              ( ((self.events[f'selRecoJets{sys}_pt'] > self.minLeadingJetPt) & 
                                                 (self.events[f'selRecoJetsF{sys}_pt']  > self.minSubLeadingJetPt)
                                                ) | 
                                                ((self.events[f'selRecoJetsF{sys}_pt'] > self.minLeadingJetPt) & 
                                                 (self.events[f'selRecoJets{sys}_pt']  > self.minSubLeadingJetPt)
                                                ) 
                                              )
                                            )
                    
                    self.recoMask[itrigger]=selRecoMasks[itrigger]
                    
                    if self.withLepVeto:
                        selRecoMasks[itrigger] = (selRecoMasks[itrigger]) & (self.events[f'nRecoLeptons{sys}']==0)
                    

                if self.saveParquet: 


                    for itrigger, itrValues in self.triggerTable.items():
                        self.events['selRecoMask_'+itrigger]=self.recoMask[itrigger]#selRecoMasks[itrigger]


                    print(f"Saving the following .parquet file: {self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_{self.jetType}Jet_{self.splitCount}.parquet'}")
                    ak.to_parquet(events, (self.parquetDir + self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'_{self.jetType}Jet_{self.splitCount}{self.parquetExt}.parquet'))    
                    return 1
                    

            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            
            for isel in self.selList:
                
                
                listOfOutputHistosNames = [k for k,h in output.items() if ((sys in k) or ('residual' in k and not(sys in k))  or ('smear' in k.lower() and not(sys in k)) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                
                for k in listOfOutputHistosNames:
                    key=k
                    
                    ############### Safety checks ##################
                    if not(sys in key) and not('residual' in key) and not('resol' in key) and not('smear' in key.lower()): 
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
                        #whichKinVar = [var for var in self.dict_variables_kinematics.keys() if var[1:] in key]#[0]
                        #varToFill = whichKinVar[0]
                        if ('gen' in key.lower()) and ('softdrop' in key.lower()):
                            whichKinVar = ['_mSD'] 
                            #print(key)
                            #print(whichKinVar)
                        else:
                            whichKinVar = [var for var in self.dict_variables_kinematics.keys() if var[1:] in key]#[0]
                            
                        varToFill = whichKinVar[0]
                        

                    elif (('_tau' in key and ('21' in key or '32' in key)) or '__' in key):# and 'nom' in key):
                        if 'WTA' in key: temp='WTA' 
                        elif 'exkT' in key: temp='exkT'
                        else: temp='tau' 

                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]

                        if temp!='tau':
                            keep = whichnSub[0]
                        else:
                            for w in whichnSub:
                                if not('exkT' in w) and not('WTA'in w):
                                    keep = w
                        
                        varToFill = keep

                    elif (('_tau' in key and not('21' in key or '32' in key)) or '__' in key):# and 'nom' in key):
                        temp='tau'
                        #temp='exkT' if 'exkT' in key else 'tau'
                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                        varToFill = whichnSub[0]
                    

                    ################## Filling histos from accumulated event arrays ##################
                    
                    if (key.lower().startswith(('reco','good','ht','lhe'))):
                        if 'nPV' in key and 'nom' in sys:# .startswith(tuple(self.recoWtSources)):
                            output[key].fill(self.events[f'good_nPVs{s}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)
                        #elif 'HT' in key and self.isMC:
                        #    output[key].fill(self.events[f'LHE_genHT{s}'][selGenMask],
                        #                     weight=totalGenWeight[selGenMask],
                        #                     threads=8)
                        else: 
                            
                            if self.isSigMC and 'const' in sys:
                                s=sys

                            if self.verbose: print(s,sys, key,f"filling from recoJet{self.nJet[0]}")#,self.events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask])
                            
                            if not ('leading' in key.lower()):
                                output[key].fill(self.events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                 threads=8)
                            else:
                                if self.jetType=='Central' and not('sub2' in key.lower()):
                                    if 'sub' in key.lower():
                                        output[key].fill(self.events[f'selRecoSubleadingJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         threads=8)
                                    else:
                                        output[key].fill(self.events[f'selRecoLeadingJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         threads=8)
                                elif self.jetType=='Central' and ('sub2' in key.lower()):
                                    output[key].fill(self.events[f'selRecoSub2leadingJets{self.jetFlag}{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                                         threads=8)
                                    
                                    
                                

                    # Stuff below should have diff wts already available to them so no need to touch anything down here
                    if self.isMC and (not varToFill.startswith(tuple(self.reco_only_labels))):
                        if self.isSigMC and 'const' in sys:
                            s=sys

                        if self.verbose: 
                            print(s, sys, key)
                            #############################################################################
                            # self.listOfHistTypes =  [ 'gen', 'genLeading', 'genSubleading', 'accepgen', 'missgen', 'reco', 'recoLeading', 'recoSubleading', 'fakereco', 'truereco',  ] if self.isMC else [ 'reco', 'recoLeading', 'recoSubleading' ] 
                            #############################################################################
                        if self.isSigMC:
                            if (key.lower().startswith('true')):
                                output[key].fill(self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                            elif (key.lower().startswith('fake')):
                                output[key].fill(self.events[f'selRecoJets{self.jetFlag}{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                                 threads=8)

                            elif (key.lower().startswith('accepgenjet')):# and not key.startswith(('accep','miss')):
                                #if self.verbose: 
                                #    print('accepgen', key,s,sys, len(self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask]))
                                output[key].fill(self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                                 threads=8)

                            elif (key.lower().startswith('missgenjet')):
                                #if self.verbose: 
                                #    print('missgen', key,s,sys, len(self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask]))
                                output[key].fill(self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                                 threads=8)
                                
                            elif ( 'resp' in key.lower() and not('miss' in key.lower())):
                                #if self.verbose: print("filling resp")
                                #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                                output[key].fill(gen=self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                                output[key].fill(gen=self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask])),
                                                 weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                                 threads=8)     

                            elif ('respwithmiss' in key.lower()):
                                #if self.verbose: print("filling resp with miss")
                                #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                                output[key].fill(gen=self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask], weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                                output[key].fill(gen=self.events[f'accepGenJets{self.jetFlag}{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][trueRecoMask])),
                                                 weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill missgen weight
                                output[key].fill(gen=self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask], reco=-1.*np.ones(len(self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][missGenMask])),
                                                 weight=totalGenWeight[missGenMask],
                                                 threads=8)


                            elif ('residual' in key.lower() or 'resol' in key.lower() or 'smear' in key.lower()) and sys.endswith('nom') and self.isSigMC:
                                genVarToFill = '_mSD' if 'msoftdrop_new' in varToFill else varToFill

                                zeroMask=(self.events[f'accepGenJets{self.jetFlag}{s}{genVarToFill}']!=0.)&(accepGenMask)

                                response = self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]/self.events[f'accepGenJets{self.jetFlag}{s}{genVarToFill}'][zeroMask]
                                response = np.nan_to_num(response,nan=-999.)
                                residual = self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]-self.events[f'accepGenJets{self.jetFlag}{s}{genVarToFill}'][zeroMask]
                                relativeRes = residual/self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask]
                                #if 'pt' in key: print(key, residual[0:10],relativeRes[0:10],self.events[f'trueRecoJets{self.jetFlag}{s}{varToFill}'][zeroMask][0:10])
                                if 'noWt_' in key: output[key].fill(response)#, weight=totalRecoWeight[zeroMask])
                                elif 'mSmear' in key: output[key].fill(relativeRes)#, weight=totalRecoWeight[zeroMask])
                                elif 'residual' in key: output[key].fill(residual)
                            

                        if (key.lower().startswith('genjet')):
                            #if self.verbose or ('softdrop' in key): 
                            #    print('gen', key,s,sys, len(self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][selGenMask]))
                            output[key].fill(self.events[f'selGenJets{self.jetFlag}_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)#=self._listofHistograms(histoName) #hist.Hist("Events", hist.Cat("branch", branch), hist.Bin("value", branch, 100, 0, 1000))    
                        elif (key.lower().startswith('genleading')) and self.jetType=='Central':# and not key.startswith(('accep','miss')):
                            #if self.verbose: 
                            #    print('leading gen', key,s,sys, len(self.events[f'selGenLeadingJets{varToFill}'][selGenMask]))
                            output[key].fill(self.events[f'selGenLeadingJets{self.jetFlag}_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)
                        elif (key.lower().startswith('gensubleading')) and self.jetType=='Central':# and not key.startswith(('accep','miss')):
                            #if self.verbose: 
                            #    print('subleading gen', key,s,sys, len(self.events[f'selGenSubleadingJets{varToFill}'][selGenMask]))
                            output[key].fill(self.events[f'selGenSubleadingJets{self.jetFlag}_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)
                            
                        elif (key.lower().startswith('gensub2leading')) and self.jetType=='Central':# and not key.startswith(('accep','miss')):
                            #if self.verbose: 
                            #    print('subleading gen', key,s,sys, len(self.events[f'selGenSubleadingJets{varToFill}'][selGenMask]))
                            output[key].fill(self.events[f'selGenSub2leadingJets{self.jetFlag}_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)
                        
                        
        l=[]
        
        for x,y in output.items(): #y.SetDirectory(0)
        
            if self.sysUnc and self.onlyUnc!=''  and x.startswith(('accepgen','miss' )) and '_nom' in x:#,'true','fake'
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
                        if sysUnc.endswith('nom') or self.onlyUnc!='' :
                            
                            if not (x in tuple(self.reco_only_labels)): 
                                dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())        
                            
                            else: 
                                #for nPVs
                                if (sysUnc.endswith('nom') or sysUnc.startswith(tuple(self.recoWtSources))) and not('gen' in itype) :
                                    
                                    if 'npv' in x.lower():
                                    
                                        dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight()) 
                                        
                                    #elif 'ht' in x.lower():# in x.lower():
                                    #    #print("BRANCH:", x[1:]+sysUnc+isel)
                                    #    dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label='LHE H_{T}', underflow=True,overflow=True).Weight())
                                    else:
                                        dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())     
                                
                    if itype.startswith('truereco') and self.isSigMC and sysUnc.endswith('nom'):
                        dictOfHists['residual'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(1000, -50, 50, name='residual'+iJ+'_pt'+isel, label='AK8 reco pt - gen jet pt', underflow=False,overflow=False).Weight())
                        
                        dictOfHists['noWt_resol'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+'_pt'+isel, label='AK8 reco/gen jet ', underflow=False,overflow=False).Weight())
                        
                        dictOfHists['mSmearCheck'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(1000, -5, 5, name='mSmearCheck'+iJ+'_pt'+isel, label=f'(reco pt - gen pt) / reco pt ', underflow=False,overflow=False).Weight())
                        
                        if not(self.sysUnc) and ('_mSD' in self.dict_variables_kinematics.keys() or '_msoftdrop' in self.dict_variables_kinematics.keys() or '_msoftdrop_new' in self.dict_variables_kinematics.keys()):
                            
                            dictOfHists['residual'+iJ+'_msoftdrop_new'+isel] = (hist.Hist.new.Regular(500, -20, 20, name='residual'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD} - gen jet m_{SD}', underflow=False,overflow=False).Weight())

                            dictOfHists['residual'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(500, -20, 20, name='residual'+iJ+'_mass'+isel, label='AK8 reco inv. m - gen jet inv. m', underflow=False,overflow=False).Weight())



                            dictOfHists['noWt_resol'+iJ+'_msoftdrop_new'+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+'_msoftdrop'+isel, label='AK8 reco/gen jet ', underflow=False,overflow=False).Weight())

                            dictOfHists['noWt_resol'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+'_mass'+isel, label='AK8 reco/gen jet ', underflow=False,overflow=False).Weight())

                            dictOfHists['mSmearCheck'+iJ+'_msoftdrop_new'+isel] = (hist.Hist.new.Regular(1000, -5, 5, name='mSmearCheck'+iJ+'_msoftdrop'+isel, label=f'(reco mSD - gen mSD / reco mSD ', underflow=False,overflow=False).Weight())

                            dictOfHists['mSmearCheck'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(1000, -5, 5, name='mSmearCheck'+iJ+'_mass'+isel, label=f'(reco mass - gen mass) / reco mass ', underflow=False,overflow=False).Weight())
                        
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
                                bound = 20 if not('tau' in x) else 1.5
                                dictOfHists['residual'+iJ+x+isel] = (hist.Hist.new.Regular(1000, -1.*bound, bound, name='residual'+iJ+x+isel, label=f'AK8 reco {"#tau" if "tau" in x else x} - gen jet {"#tau" if "tau" in x else x}', underflow=True,overflow=True).Weight())
                                dictOfHists['noWt_resol'+iJ+x+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
                                dictOfHists['mSmearCheck'+iJ+x+isel] = (hist.Hist.new.Regular(1000, -5, 5, name='mSmearCheck'+iJ+x+isel, label=f'(reco {x} - gen {x}) / reco {x} ', underflow=True,overflow=True).Weight())
        
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
            
        if self.verbose: print ("Preparing branches to read",self.sysSource,self.wtSources)
            
        for sys in self.sysSource:
            #if not wtUnc and not sysUnc:
            
            if 'Central' in self.jetType:
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                    reco_list.append(f'selRecoJetsF{sys}_pt')
                    
                if (self.isMC and not(self.isSigMC) and sys.endswith('nom')) or (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    if not ('selGenJetsF_nom_pt' in gen_list): 
                        gen_list.append('selGenJetsF_nom_pt')
                        
                    gen_list.append(f'accepGenJetsF{sys}_pt')
                    reco_list.append(f'trueRecoJetsF{sys}_pt')
                
                for i in kinematic_labels+list(nSub_labels.keys()): 
                    if '__' in i: i=i.replace('__','_')
                    
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                        reco_list.append('selRecoJets'+sys+i)
                        if not(self.sysUnc):
                            reco_list.append('selRecoLeadingJets'+sys+i)
                            reco_list.append('selRecoSubleadingJets'+sys+i)
                            #reco_list.append('selRecoSub2leadingJets'+sys+i)
                        
                    if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                        reco_list.append('trueRecoJets'+sys+i)
                        if not i in reco_only_labels: 
                            gen_list.append('accepGenJets'+sys+i)
                            
                    if self.isMC and sys.endswith('nom'): 
                        if not i in reco_only_labels: 
                            gen_list.append('selGenJets'+sys+i)
                            if not(self.sysUnc):
                                gen_list.append('selGenLeadingJets'+sys+i)
                                gen_list.append('selGenSubleadingJets'+sys+i)
                                #gen_list.append('selGenSub2leadingJets'+sys+i)

            elif 'Forward' in self.jetType: 
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                    reco_list.append(f'selRecoJets{sys}_pt')
                                    
                if (self.isMC and not(self.isSigMC) and sys.endswith('nom')) or (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    if not ('selGenJets_nom_pt' in gen_list): 
                        gen_list.append('selGenJets_nom_pt')
                        
                    gen_list.append(f'accepGenJets{sys}_pt')
                    reco_list.append(f'trueRecoJets{sys}_pt')
                
                for i in kinematic_labels+list(nSub_labels.keys()): 
                    if '__' in i: i=i.replace('__','_')
                    
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                        reco_list.append('selRecoJetsF'+sys+i)
                        #reco_list.append('selRecoLeadingJets'+sys+i)
                        #reco_list.append('selRecoSubleadingJets'+sys+i)
                        
                    if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                        reco_list.append('trueRecoJetsF'+sys+i)
                        if not i in reco_only_labels: 
                            gen_list.append('accepGenJetsF'+sys+i)
                            
                    if self.isMC and sys.endswith('nom'): 
                        if not i in reco_only_labels: 
                            gen_list.append('selGenJetsF'+sys+i)
                            #gen_list.append('selGenLeadingJets'+sys+i)
                            #gen_list.append('selGenSubleadingJets'+sys+i)


            if (self.isMC and (sys.endswith('nom')) or (self.isSigMC and self.sysUnc and not('const' in sys))):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("totalRecoWeight"+sys)
                #if not ("evtGenWeight_nom" in gen_list): 
                gen_list.append("evtGenWeight"+sys)
                if not('passGenSel'+sys in gen_list):# and 'nom' in sys: 
                    gen_list.append('passGenSel'+sys)              
                    
                if not('passRecoSel'+sys in reco_list): 
                    reco_list.append('passRecoSel'+sys)
                if sys.endswith('nom'): 
                    reco_list.append("good_nPVs"+sys)
                                            
                    
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("good_nPVs"+sys)
                
                        
                if not('passRecoSel'+sys in reco_list): 
                    reco_list.append('passRecoSel'+sys)
                    
                    
                for itrigger in self.triggerTable.keys():    
                    if not(self.year.endswith('VFP') and self.era=='B'): triggerBit_list.append(f'HLT_{itrigger}') #to remove
                    triggerBit_list.append(f'passHLT_{itrigger}')
                
                if '2018' in self.year or '2017' in self.year: 
                    triggerBit_list.append(f'passHLT_AK8PFJet550')
                    triggerBit_list.append(f'HLT_AK8PFJet550')
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or self.sysUnc):
                if 'pu' in sys or 'l1' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: 
                    
                    if self.onlyRedonePDFandAlphaSWts:
                        if 'pdf' in sys.lower():
                            if not( 'pdfWeightAll_nom'  in gen_reweights_list): 
                                gen_reweights_list.append('pdfWeightAll_nom')
                            
                        gen_reweights_list.append(sys.split('_')[1]+'_nom')

                        #    gen_reweights_list.append((sys.split('_')[1]+'_nom').replace('pdfWeight', 'pdfWeight2'))
                        #    gen_reweights_list.append((sys.split('_')[1]+'_nom').replace('pdf', 'AlphaS'))
                    elif not(self.onlyRedonePDFandAlphaSWts):
                        gen_reweights_list.append(sys.split('_')[1]+'_nom')


        if self.isSigMC and not(self.sysUnc): 
            branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list+['nRecoLeptons_nom','nGenLeptons_nom', 'pt_asymm_nom','delta_phi_nom','delta_R_nom', 'gen_pt_asymm_nom','gen_delta_phi_nom','gen_delta_R_nom', 'genWeight']#+['recoSelectedEventNumber_nom']
            
        elif self.isSigMC and (self.sysUnc): 
            branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list+['genWeight']
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list+['nRecoLeptons_nom','nGenLeptons_nom', 'pt_asymm_nom','delta_phi_nom','delta_R_nom', 'gen_pt_asymm_nom','gen_delta_phi_nom','gen_delta_R_nom']#+['recoSelectedEventNumber_nom']
        elif (not self.isMC): branchesToRead=reco_list+triggerBit_list+['recoSelectedEventNumber_nom','nRecoLeptons_nom', 'pt_asymm_nom','delta_phi_nom','delta_R_nom']
            
        if self.isMC or self.isSigMC:
            for x in range(len(branchesToRead)): 
                if ('gen' in branchesToRead[x].lower()) and ('msoftdrop_new' in branchesToRead[x].lower()):
                #    print(branchesToRead[x])
                    branchesToRead[x] = branchesToRead[x].replace('msoftdrop_new','mSD')
                #    print(branchesToRead[x])
                    
        ##if self.verbose: 
        #print("Prepared branches to read",branchesToRead)        
        return branchesToRead
    
    def append_recoWeights_forJetVeto(self):
    
        #for non-signal MC, signal MC sys variations on jes/jer, or signal MC when run without exp/theory weight unc. variations
        #if (self.isMC and not(self.isSigMC)) or (self.isSigMC and not(self.wtUnc)): 

        #systematic_suffix_map = {
        #                         'central': 'Nom',
        #                        }

        if not(self.sysUnc) or ('const' in self.onlyUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        new_fields = OrderedDict()


        if not (self.sysUnc) or 'const' in self.onlyUnc: 
            new_fields[f'new{self.sel}_eventJetVeto_nom'] = np.ones(len(self.events['selRecoJets_nom_eta']))
            new_fields[f'new{self.sel}_CentralJetVeto_nom'] = np.ones(len(self.events['selRecoJets_nom_eta']))
            new_fields[f'new{self.sel}_ForwardJetVeto_nom'] = np.ones(len(self.events['selRecoJetsF_nom_eta']))
        else:
            for sys in self.sysSources:
                new_fields[f'new{self.sel}_eventJetVeto{sys}'] = np.ones(len(self.events[f'selRecoJets{sys}_eta']))
                new_fields[f'new{self.sel}_CentralJetVeto{sys}'] = np.ones(len(self.events[f'selRecoJets{sys}_eta']))
                new_fields[f'new{self.sel}_ForwardJetVeto{sys}'] = np.ones(len(self.events[f'selRecoJetsF{sys}_eta']))

        for sys in sysList:     
            print(sys)
            #CentralJet_pt = self.events[f'selRecoAK4bjetLeptHem{sys}_pt']
            CentralJet_eta = self.events[f'selRecoJets{sys}_eta'][self.recoMask]
            CentralJet_phi = self.events[f'selRecoJets{sys}_phi'][self.recoMask]
            #CentralJet_mass = self.events[f'selRecoAK4bjetLeptHem{sys}_mass']

            #ForwardJet_pt = self.events[f'selRecoJets{sys}_pt']
            ForwardJet_eta = self.events[f'selRecoJetsF{sys}_eta'][self.recoMask]
            ForwardJet_phi = self.events[f'selRecoJetsF{sys}_phi'][self.recoMask]
            #ForwardJet_mass = events[f'selRecoJets{sys}_mass']

            vetoForward = self.vetoMap.evaluate('jetvetomap',
                                            ForwardJet_eta,
                                            ForwardJet_phi                                      
                                           )
            vetoCentral = self.vetoMap.evaluate('jetvetomap',
                                            CentralJet_eta,
                                            CentralJet_phi                                      
                                           )

            #print(vetoForward,(vetoForward==0), len(vetoForward), len(vetoForward==0))
            #print(vetoCentral,(vetoCentral==0), len(vetoCentral), len(vetoCentral==0))

            new_fields[f'new{self.sel}_ForwardJetVeto{sys}'][self.recoMask] = (vetoForward==0)
            new_fields[f'new{self.sel}_CentralJetVeto{sys}'][self.recoMask] = (vetoCentral==0)


            #new_fields[f'new{self.sel}_eventJetVeto{sys}'][~self.recoMask] = ((vetoForward!=0) | (vetoCentral!=0))[~self.recoMask]
            new_fields[f'new{self.sel}_eventJetVeto{sys}'][self.recoMask] = ((vetoForward==0) & (vetoCentral==0))#[self.recoMask]

        return new_fields
    


def histoMaker(myProcessor, sampleIdentifier='qcd_ht', y='2017', sampleDict_PFNano=OrderedDict(),
               sampleDict_local=OrderedDict(), writeChunks=True, writeAccumulatedORhadd=1, verbose=False, 
               saveParquet=False,
               isSigMC=True, isMC=True, wtUnc=True, sysUnc=False, onlyUnc='', outputdir='processorTests/', 
               sysSource=[],
               era='', ext='_nomWts', splitchunks=10, nWorkers=40, stepSize="2048 MB",
               jetType='Central',forceProduction=False,onlyParquet=False, 
               minLeadPtThreshold=200., minSubLeadPtThreshold=200.,onlyControlHistos=False):
    
    nchunk=copy.deepcopy(splitchunks)
    cz=0
    onlyParquet=saveParquet

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
                                     onlyUnc='' if 'jes' in onlyUnc and not('HEM' in onlyUnc) else onlyUnc,era=era,
                                     verbose=False,splitCount=SC,jetType=jetType, 
                                     minLeadPt=minLeadPtThreshold, minSubLeadPt=minSubLeadPtThreshold,
                                     parquetExt=ext, onlyControlHistos=onlyControlHistos )
            if cz==0:
                #if verbose: print("Branches being read \n", my_processor._branchesToRead)
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

                #if verbose: 
                print(sample,year,f"Splitting input filelist of {sum(n_tosplit)} files, into sublists of {len(n_tosplit)} file chunks:{n_tosplit[0]}")
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
                    if 'pt' in sample.lower() and not ('flat' in sample.lower()): 
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]+sample.split("_Tune")[0].split("_")[2]}_UL{year}{ext}'
                    elif ('flat' in sample.lower()):
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]}_UL{year}{ext}'.replace('-','')
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
                    
                    
                    #my_processor=myProcessor(sampleName=sample, sampleDict=sampleDict_PFNano,
                    #                                     isMC=isMC, isSigMC=isSigMC,year=year,
                    #                                     saveParquet=saveParquet,onlyParquet=onlyParquet,
                    #                                     sysSource=sysSource,
                    #                                    wtUnc=wtUnc,sysUnc=sysUnc,
                    #                                     onlyUnc='' if 'jes' in onlyUnc else onlyUnc,era=era,
                    #                                     verbose=False,splitCount=SC,jetType=jetType, parquetExt=ext );
                
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
                        print(f"Time taken to load {len(events)} events from {n_tosplit[c]} files in {sample} = {elapsed1}, for {ext}, {c}")

                    if verbose: 
                        print(f'nEvents in this file chunk: {len(events)}')

                    if verbose: 
                        tstart2=time.time()

                    if verbose: print("#####Processing events#####")

                    processed_events=my_processor.process(events)
                    gc.collect()
                    if verbose:
                        elapsed2 = time.time()-tstart2
                        print(f"Time taken to build histos from {len(events)} events, for {fnstem,jetType} = {elapsed2}")

                    
                    #SC=chr(c)
                    if saveParquet: 
                        print(f"saving parquet for {c}")
                        my_processor.splitCount=f'{c+1}'
                        c=c+1
                        del(processed_events)
                        del(events)
                        gc.collect()
                        continue
                        
                    if c==0: 
                        hists=processed_events
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of reco phi,C nominal hist: {hists['recoJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Final integral of reco tau_1_1,C nominal hist: {hists['recoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of truereco tau_1_1,C nominal hist: {hists['truerecoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")


                                    print(f"Final integral of gen phi,C nominal hist: {hists['genJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    #if isMC:
                                    print(f"Final integral of gen tau_1_1,C nominal hist: {hists['genJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of accepgen tau_1_1,C nominal hist: {hists['accepgenJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")

                                    print(f"Final no flow integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=False)}, for {ext}, {c}")
                                    print(f"Final integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                            elif jetType=="Forward":
                                print(f"Current integral of reco phi,F nominal hist: {hists['recoJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Current integral of tau_1_1,F nominal hist: {hists['recoJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                    elif c>0 and writeAccumulatedORhadd==0: 
                        hists=accumulate([hists,processed_events])
                        if verbose: 
                            if jetType=='Central':
                                print(f"Current integral of reco phi,C nominal hist: {hists['recoJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Final integral of reco tau_1_1,C nominal hist: {hists['recoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of truereco tau_1_1,C nominal hist: {hists['truerecoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")


                                    print(f"Final integral of gen phi,C nominal hist: {hists['genJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    #if isMC:
                                    print(f"Final integral of gen tau_1_1,C nominal hist: {hists['genJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of accepgen tau_1_1,C nominal hist: {hists['accepgenJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")

                                    print(f"Final no flow integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=False)}, for {ext}, {c}")
                                    print(f"Final integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                            elif jetType=="Forward":
                                print(f"Current integral of reco phi,F nominal hist: {hists['recoJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Current integral of tau_1_1,F nominal hist: {hists['recoJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                    
                    if not sysUnc: 
                        if 'pt' in sample.lower() and not('flat' in sample.lower()): 
                            string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]+sample.split("_Tune")[0].split("_")[2]}_UL{year}{ext}'
                        elif ('flat' in sample.lower()):
                            string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]}_UL{year}{ext}'.replace('-','')
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
                        if verbose and c>0: 
                            if jetType=='Central':
                                print(f"Current integral of reco phi,C nominal hist: {hists['recoJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Final integral of reco tau_1_1,C nominal hist: {hists['recoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of truereco tau_1_1,C nominal hist: {hists['truerecoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")


                                    print(f"Final integral of gen phi,C nominal hist: {hists['genJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    #if isMC:
                                    print(f"Final integral of gen tau_1_1,C nominal hist: {hists['genJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of accepgen tau_1_1,C nominal hist: {hists['accepgenJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")

                                    print(f"Final no flow integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=False)}, for {ext}, {c}")
                                    print(f"Final integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                            elif jetType=="Forward":
                                print(f"Current integral of reco phi,F nominal hist: {hists['recoJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Current integral of tau_1_1,F nominal hist: {hists['recoJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    

                        #print (processed_events)
                        print ("Events proccesed, now writing output file(s)",fn,c,len(i),len(ifls))

                        with uproot.recreate(fn) as fout:#outputTest_{sample.split("_Tune")[0]}_{year}.root'
                            for key in processed_events.keys():
                                fout[key]=processed_events[key]
                            print (f'Done with creating output file:{fn}')
                            fout.close()
                    elif writeAccumulatedORhadd==0:# or len(n_tosplit)==1:
                        if verbose and c>0: 
                            if jetType=='Central':
                                print(f"Current integral of reco phi,C nominal hist: {hists['recoJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Final integral of reco tau_1_1,C nominal hist: {hists['recoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of truereco tau_1_1,C nominal hist: {hists['truerecoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")


                                    print(f"Final integral of gen phi,C nominal hist: {hists['genJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    #if isMC:
                                    print(f"Final integral of gen tau_1_1,C nominal hist: {hists['genJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    print(f"Final integral of accepgen tau_1_1,C nominal hist: {hists['accepgenJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")

                                    print(f"Final no flow integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=False)}, for {ext}, {c}")
                                    print(f"Final integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    
                            elif jetType=="Forward":
                                print(f"Current integral of gen phi,F nominal hist: {hists['genJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                if isMC:
                                    print(f"Current integral of gen tau_1_1,F nominal hist: {hists['genJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                    
                                    

                        #print (processed_events)
                        #print ("Events proccesed, now writing output file(s)",fn)   
                    my_processor.splitCount=f'{c+1}'
                    c=c+1
                    del(processed_events)
                    del(events)
                    gc.collect()
            if not(pathExistsFlag):
                if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):

                    if verbose: 
                        if jetType=='Central':
                            print(f"Final integral of reco phi,C nominal hist: {hists['recoJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                            if isMC:
                                print(f"Final integral of reco tau_1_1,C nominal hist: {hists['recoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                print(f"Final integral of truereco tau_1_1,C nominal hist: {hists['truerecoJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                
                                
                                print(f"Final integral of gen phi,C nominal hist: {hists['genJet_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                #if isMC:
                                print(f"Final integral of gen tau_1_1,C nominal hist: {hists['genJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                print(f"Final integral of accepgen tau_1_1,C nominal hist: {hists['accepgenJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                
                                print(f"Final no flow integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=False)}, for {ext}, {c}")
                                print(f"Final integral of respWithMissJet tau_1_1,C nominal hist: {hists['respWithMissJet_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                
                                

                        elif jetType=="Forward":
                            print(f"Final integral of reco phi,F nominal hist: {hists['recoJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                            if isMC:
                                print(f"Final integral of reco tau_1_1,F nominal hist: {hists['recoJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                print(f"Final integral of gen phi,F nominal hist: {hists['genJetF_phi_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                print(f"Final integral of gen tau_1_1,F nominal hist: {hists['genJetF_tau_1_1_nom_dijetSel'].sum(flow=True)}, for {ext}, {c}")
                                
                                

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




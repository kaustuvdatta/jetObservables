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
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
from collections import OrderedDict
import numpy as np
import awkward as ak
import uproot
import hist
from hist import Hist

from collections import defaultdict
import numba



class nSubBasis_unfoldingHistoProd_DijetsCentral(processor.ProcessorABC):
    
    def __init__(self, sampleName, sysSource=[],year='2017', era='', isMC=True, isSigMC=True, 
                 onlyUnc='', wtUnc=False, sampleDict=dictSamples,test=False, sysUnc=False, verbose=True):#isaltSigMC=False,
        self.test=test
        self.year = year
        self.isMC = isMC
        self.isSigMC = isSigMC
        self.era = era
        #self.isaltSigMC = isaltSigMC
        self.verbose=verbose
        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc
        
        
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        
        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with data' )
        
        
        ### Helpers
        self.listOfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        #self.outputDir = 'UL17and18_nano_tests/'
        
        self.nJet = [ 'Jet'] #1', 'Jet2' ]
        self.triggerTable = OrderedDict()
        self.triggerTable[ 'AK8PFJet80' ] = {# from the list below, only the first two numbers (trigger turn on) are used. The others were a test - Alejandro
                    '2016_preVFP' : [ 200,   242 ],
                    '2016' : [ 200,   245 ],
                    '2017' : [ 200,   257, 50510.97, 6727.43, 6854.43, 23381.99, 38371.17 ],
                    '2018' : [ 200,   267, 17046.53, 34000.07, 33988.86, 33998.78 ],
                    }
        self.triggerTable[ 'AK8PFJet140' ] = {
                    '2016_preVFP' : [ 242,   308 ],
                    '2016' : [ 245,   310 ],
                    '2017' : [ 257,   323, 672.78, 1644.76, 1255.72, 2027.79, 2315.59 ],
                    '2018' : [ 267,   332, 1601.95, 1207.43, 1220.84, 1184.09 ],
                    }
        self.triggerTable[ 'AK8PFJet200' ] = {
                    '2016_preVFP' : [ 308,   373 ],
                    '2016' : [ 310,   375 ],
                    '2017' : [ 323,   389, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 332,   398, 368.33, 281.73, 284.86, 276.30 ],
                    }
        self.triggerTable[ 'AK8PFJet260' ] = {
                    '2016_preVFP' : [ 373,   439 ],
                    '2016' : [ 375,   440 ],
                    '2017' : [ 389,   455, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 398,   464, 132.00, 128.00, 122.91, 127.42 ],
                    }
        self.triggerTable[ 'AK8PFJet320' ] = {
                    '2016_preVFP' : [ 439,   526 ],
                    '2016' : [ 440,   525 ],
                    '2017' : [ 455,   543, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 464,   551, 49.29, 48.00, 47.28, 47.92 ],
                    }
        self.triggerTable[ 'AK8PFJet400' ] = {
                    '2016_preVFP' : [ 526,   580 ],
                    '2016' : [ 525,   580 ],
                    '2017' : [ 543,   598, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 551,   606, 16.39, 16.00, 15.92, 15.99 ],
                    }
        self.triggerTable[ 'AK8PFJet450' ] = {
                    '2016_preVFP' : [ 580,   635 ],
                    '2016' : [ 580,   635 ],
                    '2017' : [ 598,   653, 54.41, 382.72, 263.79, 380.97, 410.70  ],
                    '2018' : [ 606,   661, 8.43, 8.00, 7.98, 8.00 ],
                    }
        self.triggerTable[ 'AK8PFJet500' ] = {
                    '2016_preVFP' : [ 635,   10000000. ],
                    '2016' : [ 635,   1000000. ],
                    '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                    '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
                    }
        #self.triggerTable[ 'AK8PFJet550' ] = {
        #            #'2016_preVFP' : [ 635,   10000000. ],
        #            #'2016' : [ 635,   1000000. ],
        #            '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
        #            '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
        #            }
        
        
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
                    #"_mass": np.array([(i/2.0) for i in np.arange(0.*2, 300*2.1)]),
                    #"_msoftdrop": np.array([(i/2.0) for i in np.arange(0.*2, 200*2.1)]),
                    #"_pt": np.array([(i/2.0) for i in np.arange(170.*2, 2500.*2.1)]),
                        }
        self.kinematic_labels=['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels=['_good_nPVs']

        self.dict_variables_kinematics = {

                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),
                            "_msoftdrop": np.array([i for i in np.arange(0., 300, 5.)]),
                            "_good_nPVs": np.array([i for i in np.arange(0., 100, 1.)]),

                              }
        
        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not( isys.endswith('nom') or self.sysUnc) ]
        self.sysWeightList = ( '_pu', '_pdf', '_isr', '_fsr', '_l1prefiring' ) #'_ps',
        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight'] if self.wtUnc else [] 
        if self.onlyUnc: self.sysSource = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] #################### Not using for right now
        
        #puWeights used only to change reco event weight (up/down) without application to gen weight, others applied to modify the genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)]
        self.selList = '_dijetSel' 
        
                
        
        self._branchesToRead = self.getBranchesToRead(
                                                     dirname=self.inputDir, year=self.year,
                                                     kinematic_labels=self.kinematic_labels,
                                                     reco_only_labels=self.reco_only_labels,
                                                     nSub_labels=self.dict_variables_toUnfold,
                                                    )# isMC=True, isSigMC=False, wtUnc=False, sysUnc=False,
        #self._accumulator = self.buildDictOfHistograms()#processor.dict_accumulator()
    
    #@property
    #def accumulator(self):
    #    return self._accumulator    
        
    def process(self, events):
        '''fill Hist histograms to accumulate event stats, convert to root or whatever else after returned by processor'''

        output = self.buildDictOfHistograms()
        branches = self._branchesToRead
        events = events[branches]
        
        #print (f'systematic sources being considered:{self.sysSource}')
        for sys in self.sysSource:
            if self.isMC and self.isSigMC:# or self.isaltSigMC:
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys
                if self.verbose: print(sys,s, len(events[events[f'totalRecoWeight{s}']!=0]))

                selRecoMask = events[f'totalRecoWeight{s}']!=0.
                
                trueRecoMask = (events[f'trueRecoJets{s}_pt']>0.) & (selRecoMask)# & (events[f'trueRecoJetsF{s}_pt']>0.)
                fakeRecoMask = (events[f'trueRecoJets{s}_pt']<0.) & (selRecoMask)# & (events[f'trueRecoJetsF{s}_pt']<0.)
                selGenMask = events[f'evtGenWeight_nom'] !=0.
                accepGenMask = (events[f'accepGenJets{s}_pt']>0.) & (selGenMask)# & (events[f'accepGenJetsF{s}_pt']>0.)
                missGenMask = (events[f'accepGenJets{s}_pt']<0.) & (selGenMask)# & (events[f'accepGenJetsF{s}_pt']<0.)
                
            elif self.isMC and (not self.isSigMC): #not so relevant for dijets but for background MC's in W/top
                
                if sys.endswith('nom'):
                    selRecoMask = events[f'totalRecoWeight{sys}']!=0.
                    selGenMask = events[f'evtGenWeight{sys}']!=0.
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                selRecoMasks = OrderedDict()
                #masking now trigger dependent to create different histos for data events passing a given prescale trigger exclusively or an unprescaled trigger
                print("building trigger masks")
                for itrigger, itrValues in self.triggerTable.items():    
                    triggerList=list(self.triggerTable.keys())
                    thistrigInd=triggerList.index(itrigger)
                    othertrigsInds=[i for i in range(0,len(triggerList)) if i!=thistrigInd]
                    
                    selRecoMasks[itrigger]=(events[f'totalRecoWeight{sys}']!=0.) & (events[f'passHLT_{itrigger}']==1 ) & ((events[f'selRecoJets{sys}_pt']>itrValues[self.year][0]) | (events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0])) & ((events[f'selRecoJets{sys}_pt']  < itrValues[self.year][1]) & (events[f'selRecoJetsF{sys}_pt']  < itrValues[self.year][1]))
                    
                    #if not(('500' in itrigger) and (self.year!='2017' or self.year!='2018')) and not(('450' in itrigger or '500' in itrigger) and (self.year.startswith('2016'))):
                        #avoiding doing the below for unprescaled triggers in each year (should I use the custom pass_HLT flags instead?check diff) to prevent duplicate events passing the prescaled triggers
                    #for t in othertrigsInds:
                    #    #if not(('500' in triggerList[t]) and (self.year!='2017' or self.year!='2018')) and not(('450' in triggerList[t] or '500' in itrigger) and (self.year.startswith('2016'))):
                    #    selRecoMasks[itrigger]=(selRecoMasks[itrigger]) & (events[f'passHLT_{triggerList[t]}']==0 )

            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            
            for isel in self.selList:
                
                
                listOfOutputHistosNames = [k for k,h in output.items() if ((sys in k) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                if self.verbose: print(listOfOutputHistosNames)
                for k in listOfOutputHistosNames:
                    key=k
                    
                    ############### Safety checks ##################
                    if not(sys in key) and not('resol' in key): 
                        if self.verbose: print(sys, key)
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

                    elif ('_tau' in key):# and 'nom' in key):
                        temp='WTA' if 'WTA' in key else 'tau' #hack, :(, to fix tau21_nonOPkT being histogrammed incorrectly
                        temp='exkT' if 'exkT' in key else 'tau'
                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                        varToFill = whichnSub[0]
                        

                    ################## Filling histos from accumulated event arrays ##################

                    # Decide what weights to use
                    if not (sys.startswith(self.sysWeightList)):
                        s=sys
                        #if self.verbose and sys.endswith('nom'): 
                        #    #print (s, events[f'totalRecoWeight_nom'],events[f'puWeightNom{s}'],events[f'l1prefiringWeightNom{s}'])
                        totalRecoWeight = events[f'evtGenWeight_nom']*events[f'puWeightNom{s}'] if self.isMC else events[f'totalRecoWeight_nom']
                        totalRecoWeight = totalRecoWeight*events[f'l1prefiringWeightNom{s}'] if (self.isMC and not('h7'in self.sampleName.lower())) else totalRecoWeight
                        totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None
                    else:
                        s='_nom'
                        if 'pu' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'l1prefiringWeightNom{s}']
                            
                        elif 'l1' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']
                            
                        elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                            totalGenWeight = events[f'evtGenWeight{s}']*events[f'{sys.split("_")[1]}{s}']
                            totalRecoWeight = totalGenWeight*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']

                    if (key.lower().startswith(('reco','good'))):
                        if 'nPV' in key:
                            output[key].fill(events[f'good_nPVs{s}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)
                        else: 
                            if self.verbose: print(s,sys, key,"filling reco", events[f'selRecoJets{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask])
                            output[key].fill(events[f'selRecoJets{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)

                    # Stuff below should have diff wts already available to them so no need to touch anything down here
                    if self.isMC and (not varToFill.startswith(tuple(self.reco_only_labels))):

                        if (key.lower().startswith('true')):
                            output[key].fill(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                        elif (key.lower().startswith('fake')):
                            output[key].fill(events[f'selRecoJets{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                             threads=8)

                        if (key.lower().startswith('genjet')):
                            output[key].fill(events[f'selGenJets_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)#=self._listofHistograms(histoName) #hist.Hist("Events", hist.Cat("branch", branch), hist.Bin("value", branch, 100, 0, 1000))    

                        elif (key.lower().startswith('accepgenjet')):# and not key.startswith(('accep','miss')):
                            output[key].fill(events[f'accepGenJets{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                             threads=8)

                        elif (key.lower().startswith('missgenjet')):
                            output[key].fill(events[f'selGenJets_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                             threads=8)
                        
                        elif ( 'resp' in key.lower() and not('miss' in key.lower())):
                            if self.verbose: print("filling resp")
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)     
                        
                        elif ('respwithmiss' in key.lower()):
                            if self.verbose: print("filling resp with miss")
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
                            

                        elif ('resol' in key.lower()) and sys.endswith('nom'):
                            zeroMask=(events[f'accepGenJets{s}{varToFill}']!=0.)&(accepGenMask)

                            response = events[f'trueRecoJets{s}{varToFill}'][zeroMask]/events[f'accepGenJets{s}{varToFill}'][zeroMask]
                            response = np.nan_to_num(response,nan=-999.)
                            
                            output[key].fill(response, weight=totalRecoWeight[zeroMask])
                            if 'noWt_' in key: output[key].fill(response)
        l=[]
        for x,y in output.items(): #y.SetDirectory(0)
        
            if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss','true','fake' )) and '_nom' in x:
                l.append(x)
        for k in l:
            del(output[k])
        
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
        print(self.selList)
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
                                dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())  
                                
                    if itype.startswith('truereco') and self.isMC and sysUnc.endswith('nom'):
                            dictOfHists['resol'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_pt'+isel, label='AK8 reco/gen jet pt', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD}/gen jet m_{SD}', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_mass'+isel, label='AK8 reco inv. m/gen jet inv. m', underflow=True,overflow=True).Weight())
                    #print('building unfolding histos')     
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
        
                
        return dictOfHists
    
    
    
    def getBranchesToRead(self, dirname, year, kinematic_labels, 
                          reco_only_labels, nSub_labels, 
                          jesSources=[], jerSources=[]):#,'pdfWeightAll'
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
            
            if 'pdfWeightAll' in self.wtSources: self.sysSources = self.sysSources+['pdfWeightAll'] 
        #else: self.sysSource = ['_nom']
            
        #print (sysSources)
        for sys in self.sysSource:
            #if not wtUnc and not sysUnc:
            reco_list.append('selRecoJetsF_nom_pt')
            #if self.isMC: reco_list.append('trueRecoJetsF'+sys+'_pt')
            #reco_list.append('accepGenJetsF'+sys+'_pt')
            for i in kinematic_labels+list(nSub_labels.keys()): 
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                    reco_list.append('selRecoJets'+sys+i)
                    
                if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    reco_list.append('trueRecoJets'+sys+i)
                    if not i in reco_only_labels: gen_list.append('accepGenJets'+sys+i)
                        
                if self.isMC and sys.endswith('nom'): 
                    if not i in reco_only_labels: gen_list.append('selGenJets'+sys+i)


            if self.isMC and (sys.endswith('nom') or self.sysUnc):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("totalRecoWeight"+sys)
                gen_list.append("evtGenWeight_nom")#+sys)
                reco_list.append("good_nPVs"+sys)
                
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("good_nPVs"+sys)
                for itrigger in self.triggerTable.keys():    
                    triggerBit_list.append(f'HLT_{itrigger}')
                    triggerBit_list.append(f'passHLT_{itrigger}')
                
                if '2018' in self.year or '2017' in self.year: 
                    triggerBit_list.append(f'passHLT_AK8PFJet550')
                    triggerBit_list.append(f'HLT_AK8PFJet550')
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or  self.sysUnc):
                if 'pu' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list
        elif (not self.isMC): branchesToRead=reco_list+triggerBit_list#gen_list+
        if self.verbose: print(branchesToRead)        
        return branchesToRead
    




class nSubBasis_unfoldingHistoProd_DijetsForward(processor.ProcessorABC):
    
    def __init__(self, sampleName, sysSource=[],year='2017', era='', isMC=True, isSigMC=True, 
                 onlyUnc='', wtUnc=False, sampleDict=dictSamples,test=False, sysUnc=False,verbose=False):#isaltSigMC=False,
        self.test=test
        self.year = year
        self.isMC = isMC
        self.isSigMC = isSigMC
        self.era = era
        #self.isaltSigMC = isaltSigMC
        self.verbose=verbose
        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc
        
        
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        
        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with data' )
        
        
        ### Helpers
        self.listOfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        #self.outputDir = 'UL17and18_nano_tests/'
        
        self.nJet = [ 'Jet'] #1', 'Jet2' ]
        self.triggerTable = OrderedDict()
        self.triggerTable[ 'AK8PFJet80' ] = {# from the list below, only the first two numbers (trigger turn on) are used. The others were a test - Alejandro
                    '2016_preVFP' : [ 200,   242 ],
                    '2016' : [ 200,   245 ],
                    '2017' : [ 200,   257, 50510.97, 6727.43, 6854.43, 23381.99, 38371.17 ],
                    '2018' : [ 200,   267, 17046.53, 34000.07, 33988.86, 33998.78 ],
                    }
        self.triggerTable[ 'AK8PFJet140' ] = {
                    '2016_preVFP' : [ 242,   308 ],
                    '2016' : [ 245,   310 ],
                    '2017' : [ 257,   323, 672.78, 1644.76, 1255.72, 2027.79, 2315.59 ],
                    '2018' : [ 267,   332, 1601.95, 1207.43, 1220.84, 1184.09 ],
                    }
        self.triggerTable[ 'AK8PFJet200' ] = {
                    '2016_preVFP' : [ 308,   373 ],
                    '2016' : [ 310,   375 ],
                    '2017' : [ 323,   389, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 332,   398, 368.33, 281.73, 284.86, 276.30 ],
                    }
        self.triggerTable[ 'AK8PFJet260' ] = {
                    '2016_preVFP' : [ 373,   439 ],
                    '2016' : [ 375,   440 ],
                    '2017' : [ 389,   455, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 398,   464, 132.00, 128.00, 122.91, 127.42 ],
                    }
        self.triggerTable[ 'AK8PFJet320' ] = {
                    '2016_preVFP' : [ 439,   526 ],
                    '2016' : [ 440,   525 ],
                    '2017' : [ 455,   543, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 464,   551, 49.29, 48.00, 47.28, 47.92 ],
                    }
        self.triggerTable[ 'AK8PFJet400' ] = {
                    '2016_preVFP' : [ 526,   580 ],
                    '2016' : [ 525,   580 ],
                    '2017' : [ 543,   598, 54.41, 382.72, 263.79, 380.97, 410.70 ],
                    '2018' : [ 551,   606, 16.39, 16.00, 15.92, 15.99 ],
                    }
        self.triggerTable[ 'AK8PFJet450' ] = {
                    '2016_preVFP' : [ 580,   635 ],
                    '2016' : [ 580,   635 ],
                    '2017' : [ 598,   653, 54.41, 382.72, 263.79, 380.97, 410.70  ],
                    '2018' : [ 606,   661, 8.43, 8.00, 7.98, 8.00 ],
                    }
        self.triggerTable[ 'AK8PFJet500' ] = {
                    '2016_preVFP' : [ 635,   10000000. ],
                    '2016' : [ 635,   1000000. ],
                    '2017' : [ 653.,  1000000., 1.00, 1.00, 1.00, 1.00, 1.00 ],
                    '2018' : [ 661,   1000000., 1, 1, 1, 1 ],
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
                    #"_mass": np.array([(i/2.0) for i in np.arange(0.*2, 300*2.1)]),
                    #"_msoftdrop": np.array([(i/2.0) for i in np.arange(0.*2, 200*2.1)]),
                    #"_pt": np.array([(i/2.0) for i in np.arange(170.*2, 2500.*2.1)]),
                        }
        self.kinematic_labels=['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels=['_good_nPVs']

        self.dict_variables_kinematics = {

                            "_pt": np.array([i for i in np.arange(170., 2570., 10.)]),
                            "_eta": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                            "_y": np.array([i for i in np.arange(-2.2, 2.2, 0.2)]),
                            "_phi": np.array([i for i in np.arange(-3.2, 3.2, 0.2)]),
                            "_mass": np.array([i for i in np.arange(0., 300., 5.)]),
                            "_msoftdrop": np.array([i for i in np.arange(0., 300, 5.)]),
                            "_good_nPVs": np.array([i for i in np.arange(0., 100, 1.)]),

                              }
        
        ### Uncertainties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not( isys.endswith('nom') or self.sysUnc) ]
        self.sysWeightList = ( '_pu', '_pdf', '_isr', '_fsr', '_l1prefiring' ) #'_ps',
        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight'] if self.wtUnc else [] 
        if self.onlyUnc: self.sysSource = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] #################### Not using for right now
        
        #puWeights used only to change reco event weight (up/down) without application to gen weight, others applied to modify the genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)]
        self.selList = '_dijetSel' 
        
                
        
        self._branchesToRead = self.getBranchesToRead(
                                                     dirname=self.inputDir, year=self.year,
                                                     kinematic_labels=self.kinematic_labels,
                                                     reco_only_labels=self.reco_only_labels,
                                                     nSub_labels=self.dict_variables_toUnfold,
                                                    )# isMC=True, isSigMC=False, wtUnc=False, sysUnc=False,
        #self._accumulator = self.buildDictOfHistograms()#processor.dict_accumulator()
    
    #@property
    #def accumulator(self):
    #    return self._accumulator    
        
    def process(self, events):
        '''fill Hist histograms to accumulate event stats, convert to root or whatever else after returned by processor'''

        output = self.buildDictOfHistograms()
        branches = self._branchesToRead
        events = events[branches]
        
        #print (f'systematic sources being considered:{self.sysSource}')
        for sys in self.sysSource:
            
            if self.isMC and self.isSigMC:# or self.isaltSigMC:
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys
                selRecoMask = events[f'totalRecoWeight{s}']!=0.
                trueRecoMask = (events[f'trueRecoJetsF{s}_pt']>0.) & (selRecoMask) #& (events[f'trueRecoJets{s}_pt']>0.) 
                fakeRecoMask = (events[f'trueRecoJetsF{s}_pt']<0.) & (selRecoMask) #& (events[f'trueRecoJets{s}_pt']<0.) 
                selGenMask = events[f'evtGenWeight_nom'] !=0.
                accepGenMask = (events[f'accepGenJetsF{s}_pt']>0.) & (selGenMask) #& (events[f'accepGenJets{s}_pt']>0.)
                missGenMask = (events[f'accepGenJetsF{s}_pt']<0.) & (selGenMask) #& (events[f'accepGenJets{s}_pt']<0.)
                
            elif self.isMC and (not self.isSigMC): #not so relevant for dijets but for background MC's in W/top
                
                if sys.endswith('nom'):
                    selRecoMask = events[f'totalRecoWeight{sys}']!=0.
                    selGenMask = events[f'evtGenWeight{sys}']!=0.
                    
            elif not(self.isMC) and sys.endswith('nom'): 
                selRecoMasks = OrderedDict()
                #masking now trigger dependent to create different histos for data events passing a given prescale trigger exclusively or an unprescaled trigger
                print("building trigger masks")
                for itrigger, itrValues in self.triggerTable.items():    
                    triggerList=list(self.triggerTable.keys())
                    thistrigInd=triggerList.index(itrigger)
                    othertrigsInds=[i for i in range(0,len(triggerList)) if i!=thistrigInd]
                    
                    selRecoMasks[itrigger]=(events[f'totalRecoWeight{sys}']!=0.) & (events[f'passHLT_{itrigger}']==1 ) & ((events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0]) | (events[f'selRecoJetsF{sys}_pt']>itrValues[self.year][0])) & ((events[f'selRecoJets{sys}_pt']  < itrValues[self.year][1]) & (events[f'selRecoJetsF{sys}_pt']  < itrValues[self.year][1]))
                    
                    #if not(('500' in itrigger) and (self.year!='2017' or self.year!='2018')) and not(('450' in itrigger or '500' in itrigger) and (self.year.startswith('2016'))):
                    #    #avoiding doing the below for unprescaled triggers in each year (should I use the custom pass_HLT flags instead?check diff) to prevent duplicate events passing the prescaled triggers
                    #for t in othertrigsInds:
                    #    selRecoMasks[itrigger]=(selRecoMasks[itrigger]) & (events[f'HLT_{triggerList[t]}']==0 )
                    
            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config') 
            
            
            for isel in self.selList:
                
                
                listOfOutputHistosNames = [k for k,h in output.items() if ((sys in k) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                for k in listOfOutputHistosNames:
                    key=k
                    
                    ############### Safety checks ##################
                    if not(sys in key) and not('resol' in key): 
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

                    elif ('_tau' in key):# and 'nom' in key):
                        temp='WTA' if 'WTA' in key else 'tau' #hack, :(, to fix tau21_nonOPkT being histogrammed incorrectly
                        temp='exkT' if 'exkT' in key else 'tau'
                        whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                        varToFill = whichnSub[0]
                        

                    ################## Filling histos from accumulated event arrays ##################

                    # Decide what weights to use
                    if not (sys.startswith(self.sysWeightList)):
                        s=sys
                        totalRecoWeight = events[f'evtGenWeight_nom']*events[f'puWeightNom{s}'] if self.isMC else events[f'totalRecoWeight_nom']
                        totalRecoWeight = totalRecoWeight*events[f'l1prefiringWeightNom{s}'] if (self.isMC and not('h7'in self.sampleName.lower())) else totalRecoWeight
                        totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None
                    else:
                        s='_nom'
                        if 'pu' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'l1prefiringWeightNom{s}']
                            
                        elif 'l1' in sys:
                            totalGenWeight = events[f'evtGenWeight{s}'] 
                            totalRecoWeight = totalGenWeight*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']
                            
                        elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                            totalGenWeight = events[f'evtGenWeight{s}']*events[f'{sys.split("_")[1]}{s}']
                            totalRecoWeight = totalGenWeight*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']

                    if (key.lower().startswith(('reco','good'))):
                        if 'nPV' in key:
                            output[key].fill(events[f'good_nPVs{s}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)
                        else: 
                            output[key].fill(events[f'selRecoJetsF{s}{varToFill}'][selRecoMasks[whichTrig] if not self.isMC else selRecoMask],weight=totalRecoWeight[selRecoMasks[whichTrig] if not self.isMC else selRecoMask],
                                             threads=8)

                    # Stuff below should have diff wts already available to them so no need to touch anything down here
                    if self.isMC and (not varToFill.startswith(tuple(self.reco_only_labels))):

                        if (key.lower().startswith('true')):
                            output[key].fill(events[f'trueRecoJetsF{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                        elif (key.lower().startswith('fake')):
                            output[key].fill(events[f'selRecoJetsF{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                             threads=8)

                        if (key.lower().startswith('genjet')):
                            output[key].fill(events[f'selGenJetsF_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                             threads=8)#=self._listofHistograms(histoName) #hist.Hist("Events", hist.Cat("branch", branch), hist.Bin("value", branch, 100, 0, 1000))    

                        elif (key.lower().startswith('accepgenjet')):# and not key.startswith(('accep','miss')):
                            output[key].fill(events[f'accepGenJetsF{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                             threads=8)

                        elif (key.lower().startswith('missgenjet')):
                            output[key].fill(events[f'selGenJetsF_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                             threads=8)
                        
                        elif ( 'resp' in key.lower() and not('miss' in key.lower())):
                            
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJetsF{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJetsF{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJetsF{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJetsF{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)     
                        
                        elif ('respwithmiss' in key.lower()):
                            
                            #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                            output[key].fill(gen=events[f'accepGenJetsF{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJetsF{s}{varToFill}'][trueRecoMask], weight=totalRecoWeight[trueRecoMask],
                                             threads=8)

                            #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                            output[key].fill(gen=events[f'accepGenJetsF{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJetsF{s}{varToFill}'][trueRecoMask])),
                                             weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                             threads=8)
                            
                            #fill missgen weight
                            output[key].fill(gen=events[f'selGenJetsF_nom{varToFill}'][missGenMask], reco=-1.*np.ones(len(events[f'selGenJetsF_nom{varToFill}'][missGenMask])),
                                             weight=totalGenWeight[missGenMask],
                                             threads=8)
                            

                        elif ('resol' in key.lower()) and sys.endswith('nom'):
                            zeroMask=(events[f'accepGenJetsF{s}{varToFill}']!=0.)&(accepGenMask)

                            response = events[f'trueRecoJetsF{s}{varToFill}'][zeroMask]/events[f'accepGenJetsF{s}{varToFill}'][zeroMask]
                            response = np.nan_to_num(response,nan=-999.)
                            
                            output[key].fill(response, weight=totalRecoWeight[zeroMask])
                            if 'noWt_' in key: output[key].fill(response)
        l=[]
        for x,y in output.items(): #y.SetDirectory(0)
        
            if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss','true','fake' )) and '_nom' in x:
                l.append(x)
        for k in l:
            del(output[k])
        
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
        print(self.selList)
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
                                dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())  
                                
                    if itype.startswith('truereco') and self.isMC and sysUnc.endswith('nom'):
                            dictOfHists['resol'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_pt'+isel, label='AK8 reco/gen jet pt', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD}/gen jet m_{SD}', underflow=True,overflow=True).Weight())
                            dictOfHists['resol'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(500, 0, 5, name='resol'+iJ+'_mass'+isel, label='AK8 reco inv. m/gen jet inv. m', underflow=True,overflow=True).Weight())
                    #print('building unfolding histos')     
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
        
                
        return dictOfHists
    
    
    
    def getBranchesToRead(self, dirname, year, kinematic_labels, 
                          reco_only_labels, nSub_labels, 
                          jesSources=[], jerSources=[]):#,'pdfWeightAll'
        reco_list = []
        gen_list = [] if self.isMC else None
        gen_reweights_list = [] if self.isSigMC else None
        reco_reweights_list = [] if self.isSigMC else None
        triggerBit_list = [] if not(self.isMC) else None

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
            reco_list.append('selRecoJets_nom_pt')
            #eco_list.append('trueRecoJets'+sys+'_pt')
            #reco_list.append('accepGenJets'+sys+'_pt')
            for i in kinematic_labels+list(nSub_labels.keys()): 
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                    reco_list.append('selRecoJetsF'+sys+i)
                    
                if (self.isSigMC and (sys.endswith('nom') or self.sysUnc)): 
                    reco_list.append('trueRecoJetsF'+sys+i)
                    if not i in reco_only_labels: gen_list.append('accepGenJetsF'+sys+i)
                        
                if self.isMC and sys.endswith('nom'): 
                    if not i in reco_only_labels: gen_list.append('selGenJetsF'+sys+i)


            if self.isMC and (sys.endswith('nom') or self.sysUnc):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("totalRecoWeight"+sys)
                gen_list.append("evtGenWeight_nom")#+sys)
                reco_list.append("good_nPVs"+sys)
                
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("good_nPVs"+sys)
                for itrigger in self.triggerTable.keys():    
                    triggerBit_list.append(f'HLT_{itrigger}')
                    triggerBit_list.append(f'passHLT_{itrigger}')
                
                if '2018' in self.year: triggerBit_list.append(f'passHLT_{itrigger}')
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or  self.sysUnc):
                if 'pu' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list
        elif self.isMC and (not self.isSigMC): branchesToRead=gen_list+reco_list
        elif (not self.isMC): branchesToRead=reco_list+triggerBit_list#gen_list+
        if self.verbose: print(branchesToRead)
        return branchesToRead
    
    
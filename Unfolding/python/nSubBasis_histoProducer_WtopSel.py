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
import correctionlib

from collections import defaultdict
import numba
import gc

"""
  #V3  vs V2 changes: 
  (a) added in on-the-fly btag SF and SF variation calculations (correlated and uncorrelated (b/w years) up&down for now)
  
  
"""

class nSubBasis_unfoldingHistoProd_WtopSel(processor.ProcessorABC):

    #####################################################################################################################################
    #####################################################################################################################################     #                                                                                                                                   # 
    #Relevant info wrt. what mode to run the W/top sample processing on, pragmatic choice: noDeltaR, 'purest'=massAndDeltaR (a la AL,SM)#
    # _________________________________________________________________________________________________________________________________ #
    #|      selection mode      |                                      details                                                         |#
    #|---------------------------------------------------------------------------------------------------------------------------------|#
    #|       massAndDeltaR      |                  AK8jet softdropped mass cut with dR(AK8,had.hem. b-jet)                             |#
    #|          noDeltaR        |           mass cut, and requiring exactly 1 b-tagged had. hem. AK4 and overall >1 btags              |#
    #|         onlyDeltaR       |   no mass cut, req. dR(AK8,had.hem. b-jet) and exactly 1 b-tagged had. hem. AK4 and overall >1 btags |#
    #|---------------------------------------------------------------------------------------------------------------------------------|#
    #|_________________________________________________________________________________________________________________________________|#
    #####################################################################################################################################
    #####################################################################################################################################
    
    def __init__(self, sampleName, selection='_topSel', withLeptHemBtag=False, sysSources=[],year='2017', era='', 
                 isMC=True, isSigMC=True, onlyUnc='', wtUnc=False, verbose=False, saveParquet=False, onlyParquet=False,
                 splitCount='0', sampleDict=dictSamples, test=False, sysUnc=False, selectionmode='massAndDeltaR',
                 btvJSON=None, effMap_b=None, effMap_c=None, effMap_udsg=None, btagWP="M"):
    
        self.test=test
        self.year = year
        self.isMC = isMC or isSigMC
        self.isSigMC = isSigMC
        self.era = era
        self.withLeptHemBtag = withLeptHemBtag
        self.splitCount=splitCount
        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc       
        self.dictSamples = sampleDict
        self.sampleName = sampleName
        self.verbose=verbose
        self.saveParquet=saveParquet
        self.onlyParquet=onlyParquet
        self.mode=selectionmode
        
        self.btvJSON = btvJSON
        self.effMap_b = effMap_b
        self.effMap_c = effMap_c
        self.effMap_udsg = effMap_udsg
        self.btagWP = btagWP
        
        self.WPBDisc_dict = OrderedDict()
        self.WPBDisc_dict['2016_preVFP'] = {
                                            "L": 0.0508, 
                                            "M": 0.2598, 
                                            "T": 0.6502,
                                           }
        self.WPBDisc_dict['2016'] = {
                                     "L": 0.0480, 
                                     "M": 0.2489, 
                                     "T": 0.6377,
                                    }
        self.WPBDisc_dict['2017'] = {
                                     "L": 0.0532, 
                                     "M": 0.3040, 
                                     "T": 0.7476,
                                    }
        self.WPBDisc_dict['2018'] = {
                                     "L": 0.0490, 
                                     "M": 0.2783, 
                                     "T": 0.7100,
                                    }
        self.WPBDisc = self.WPBDisc_dict[self.year][btagWP] 
        
        if self.isMC and not(self.onlyParquet):
            if self.btvJSON:
                self.btvjson = correctionlib.CorrectionSet.from_file(self.btvJSON)
            else:
                raise ValueError("btvJSON is not provided!")

            if self.effMap_b:
                effFile_b = uproot.open(self.effMap_b)
                self.effHist_b = effFile_b[f'efficiency_b_b_{self.year}{selection}']
            else:
                raise ValueError("effMap_b is not provided!")

            if self.effMap_c:
                effFile_c = uproot.open(self.effMap_c)
                self.effHist_c = effFile_c[f'efficiency_b_c_{self.year}{selection}']
            else:
                raise ValueError("effMap_c is not provided!")

            if self.effMap_udsg:
                effFile_udsg = uproot.open(self.effMap_udsg)
                self.effHist_udsg = effFile_udsg[f'efficiency_b_udsg_{self.year}{selection}']
            else:
                raise ValueError("effMap_udsg is not provided!")

            
        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with while handling data' )
        
        
        ### Helpers
        self.listOfUnfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 

        self.listOfControlHistTypes =  [ 'gen',  'reco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        
        self.nJet = [ 'Jet']
        self.selObjects = ['Mu','LeptW','MET','AK4bjetHadHem']

        self.kinematic_labels=['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop']
        self.reco_only_labels=['good_nPVs','nRecoBtags','nRecoHadBtags','nRecoLepBtags', 'selRecoHadHemDeltaR', 'selRecoHadHemDeltaRap', 'selRecoHadHemDeltaPhi']
        #self.gen_only_labels=['nGenBtags','nGenHadBtags','nGenLepBtags','selGenHadHemDeltaR','selGenHadHemDeltaRap','selGenHadHemDeltaPhi']

        self.dict_variables_kinematics_AK8 = {

                                            "_pt": np.array([i for i in np.arange(170., 3010., 10.)]) if 'WSel' in selection else  np.array([i for i in np.arange(300., 3010., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 165., 5.)]) if 'WSel' in selection else  np.array([i for i in np.arange(80., 325., 5.)]),
                                            "_msoftdrop": np.array([i for i in np.arange(0., 145., 5.)]) if 'WSel' in selection else  np.array([i for i in np.arange(100., 305., 5.)]),

                                        }  

        self.dict_variables_kinematics_Muon = {

                                            "_pt": np.array([i for i in np.arange(40., 1810., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 5., 0.05)]),

                                        } 

        self.dict_variables_kinematics_LeptW = {

                                            "_pt": np.array([i for i in np.arange(100., 2510., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 505., 5.)]),

                                        } 

        self.dict_variables_kinematics_AK4Had = {

                                            "_pt": np.array([i for i in np.arange(20., 2030., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 205., 5.)]),

                                        }

        self.dict_variables_kinematics_MET = {

                                            "_pt": np.array([i for i in np.arange(30., 1500., 10.)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),

                                        } 

        self.dict_variables_reco = { 
                                        "good_nPVs": np.array([i for i in np.arange(0., 101, 1.)]),
                                        "nRecoBtags": np.array([i for i in np.arange(0., 6, 1.)]),
                                        "nRecoHadBtags": np.array([i for i in np.arange(0., 6, 1.)]),
                                        "nRecoLepBtags": np.array([i for i in np.arange(0., 6, 1.)]),
                                        "selRecoHadHemDeltaR": np.array([i for i in np.arange(0., 1.85, 0.05)]),
                                        "selRecoHadHemDeltaRap": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                        "selRecoHadHemDeltaPhi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),

                                   }   
        
        self.dict_variables_gen = { 
                                        "nGenBtags": np.array([i for i in np.arange(0., 6., 1.)]),
                                        "nGenHadBtags": np.array([i for i in np.arange(0., 6., 1.)]),
                                        "nGenLepBtags": np.array([i for i in np.arange(0., 6., 1.)]),
                                        "selGenHadHemDeltaR": np.array([i for i in np.arange(0., 1.85, 0.05)]),
                                        "selGenHadHemDeltaRap": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                        "selGenHadHemDeltaPhi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),

                                   } 
        '''
        if self.isSigMC and not(self.wtUnc) and not(self.sysUnc) and not('var' in self.sampleName): 
            self.dict_variables_toUnfold = {

                            "_tau_0p25_1": np.array([(i/500) for i in np.arange(0.1*500, 1.*501)]),
                            "_tau_0p25_2": np.array([(i/1000) for i in np.arange(0.1*1000, 0.9*1001)]),
                            "_tau_0p25_3": np.array([(i/2000) for i in np.arange(0.1*2000, 0.8*2001)]),
                            "_tau_0p25_4": np.array([(i/5000) for i in np.arange(0.1*5000, 0.8*5001)]),
                            "_tau_0p25_5": np.array([(i/5000) for i in np.arange(0.1*5000, 0.8*5001)]),

                            "_tau_0p5_1": np.array([(i/500) for i in np.arange(0.*500, 0.9*501)]),
                            "_tau_0p5_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.8*1001)]),
                            "_tau_0p5_3": np.array([(i/5000) for i in np.arange(0.*5000, 0.6*5001)]),
                            "_tau_0p5_4": np.array([(i/5000) for i in np.arange(0.*5000, 0.6*5001)]),
                            "_tau_0p5_5": np.array([(i/5000) for i in np.arange(0.*5000, 0.54*5001)]),

                            "_tau_1_1": np.array([(i/500) for i in np.arange(0.*500, 0.7*501)]),
                            "_tau_1_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.54*1001)]),
                            "_tau_1_3": np.array([(i/5000) for i in np.arange(0.*5000, 0.4*5001)]),
                            "_tau_1_4": np.array([(i/5000) for i in np.arange(0.*5000, 0.3*5001)]),
                            "_tau_1_5": np.array([(i/5000) for i in np.arange(0.*5000, 0.28*5001)]),

                            "_tau_1p5_1": np.array([(i/500) for i in np.arange(0.*500, 0.7*501)]),
                            "_tau_1p5_2": np.array([(i/2000) for i in np.arange(0.*2000, 0.54*2001)]),
                            "_tau_1p5_3": np.array([(i/5000) for i in np.arange(0.*5000, 0.4*5001)]),
                            "_tau_1p5_4": np.array([(i/5000) for i in np.arange(0.*5000, 0.3*5001)]),
                            "_tau_1p5_5": np.array([(i/5000) for i in np.arange(0.*5000, 0.28*5001)]),

                            "_tau_2_1": np.array([(i/1000) for i in np.arange(0.*1000, 0.5*1001)]),
                            "_tau_2_2": np.array([(i/2000) for i in np.arange(0.*2000, 0.26*2001)]),
                            "_tau_2_3": np.array([(i/5000) for i in np.arange(0.*5000, 0.18*5001)]),
                            "_tau_2_4": np.array([(i/5000) for i in np.arange(0.*5000, 0.14*5001)]),
                            "_tau_2_5": np.array([(i/5000) for i in np.arange(0.*5000, 0.12*5001)]),

                            "_tau21": np.array([(i/1000) for i in np.arange(0.*1000, 1.2*1001)]),#for one-pass kT minimization as per CMS
                            "_tau32": np.array([(i/1000) for i in np.arange(0.*1000, 1.2*1001)]),#for one-pass kT minimization as per CMS

                            "_tau21_WTA": np.array([(i/1000) for i in np.arange(0.*1000, 1.2*1001)]),#for WTA-kT for comparison
                            "_tau32_WTA": np.array([(i/1000) for i in np.arange(0.*1000, 1.2*1001)]),#for WTA-kT for comparison

                            "_tau21_exkT": np.array([(i/2000) for i in np.arange(0.*2000, 1.4*2001)]),#for excl.-kT and E-scheme as per basis
                            "_tau32_exkT": np.array([(i/2000) for i in np.arange(0.*2000, 1.4*2001)]),#for excl.-kT and E-scheme as per basis

                            "_mass": np.array([(i/2) for i in np.arange(0.*2, 300*2.1)]),
                            "_msoftdrop": np.array([(i/2) for i in np.arange(0.*2, 300*2.1)]),
                            "_pt": np.array([(i/2) for i in np.arange(170.*2, 2500.*2.1)]),
                            }
        '''
        #else:
            

        if 'WSel' in selection:
            self.dict_variables_toUnfold = {

                            "_tau_0p25_1": np.array([(i/200) for i in np.arange(0.1*200, 1.05*200)]),
                            "_tau_0p25_2": np.array([(i/200) for i in np.arange(0.1*200, 0.915*200)]),
                            "_tau_0p25_3": np.array([(i/200) for i in np.arange(0.1*500, 0.862*500)]),
                            "_tau_0p25_4": np.array([(i/500) for i in np.arange(0.1*500, 0.812*500)]),
                            "_tau_0p25_5": np.array([(i/500) for i in np.arange(0.1*500, 0.802*500)]),

                            "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.1*200, 0.915*200)]),
                            "_tau_0p5_2": np.array([(i/200) for i in np.arange(0.*200, 0.805*200)]),
                            "_tau_0p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.701*1000)]),
                            "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.601*1000)]),
                            "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.511*1000)]),

                            "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.715*200)]),
                            "_tau_1_2": np.array([(i/200) for i in np.arange(0.*200, 0.515*200)]),
                            "_tau_1_3": np.array([(i/500) for i in np.arange(0.*500, 0.422*500)]),
                            "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.321*1000)]),
                            "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.301*1000)]),

                            "_tau_1p5_1": np.array([(i/200) for i in np.arange(0.*200, 0.615*200)]),
                            "_tau_1p5_2": np.array([(i/500) for i in np.arange(0.*500, 0.412*500)]),
                            "_tau_1p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.301*1000)]),
                            "_tau_1p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.251*1000)]),
                            "_tau_1p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.151*1000)]),

                            "_tau_2_1": np.array([(i/200) for i in np.arange(0.*200, 0.515*200)]),
                            "_tau_2_2": np.array([(i/500) for i in np.arange(0.*500, 0.262*500)]),
                            "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.161*1000)]),
                            "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.121*1000)]),
                            "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.101*1000)]),

                            "_tau21": np.array([(i/100) for i in np.arange(0.*100, 1.21*100)]),#for one-pass kT minimization as per CMS
                            "_tau32": np.array([(i/200) for i in np.arange(0.*200, 1.21*200)]),#for one-pass kT minimization as per CMS
                            "_tau21_WTA": np.array([(i/100) for i in np.arange(0.*100, 1.11*100)]),#for WTA-kT for comparison
                            "_tau32_WTA": np.array([(i/100) for i in np.arange(0.*100, 1.11*100)]),#for WTA-kT for comparison
                            "_tau21_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.41*200)]),#for excl.-kT and E-scheme as per basis
                            "_tau32_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.41*200)]),#for excl.-kT and E-scheme as per basis
                            #"_mass": np.array([(i*2) for i in np.arange(25., 81.)]),
                            #"_msoftdrop": np.array([(i*2) for i in np.arange(0., 81.)]),
                            #"_pt": np.array([(i*10) for i in np.arange(17., 251.)]),
                            }

        elif 'topSel' in selection:
            self.dict_variables_toUnfold = {

                            "_tau_0p25_1": np.array([(i/500) for i in np.arange(0.3*500, 0.922*500)]),
                            "_tau_0p25_2": np.array([(i/200) for i in np.arange(0.15*200, 0.875*200)]),
                            "_tau_0p25_3": np.array([(i/500) for i in np.arange(0.15*500, 0.802*500)]),
                            "_tau_0p25_4": np.array([(i/500) for i in np.arange(0.15*500, 0.722*500)]),
                            "_tau_0p25_5": np.array([(i/500) for i in np.arange(0.15*500, 0.712*500)]),

                            "_tau_0p5_1": np.array([(i/500) for i in np.arange(0.1*500, 0.882*500)]),
                            "_tau_0p5_2": np.array([(i/200) for i in np.arange(0.*200, 0.705*200)]),
                            "_tau_0p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.602*500)]),
                            "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.501*1000)]),
                            "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.451*1000)]),

                            "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.715*200)]),
                            "_tau_1_2": np.array([(i/200) for i in np.arange(0.*200, 0.515*200)]),
                            "_tau_1_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.321*1000)]),
                            "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.261*1000)]),
                            "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.221*1000)]),

                            "_tau_1p5_1": np.array([(i/200) for i in np.arange(0.*200, 0.525*200)]),
                            "_tau_1p5_2": np.array([(i/500) for i in np.arange(0.*500, 0.352*500)]),
                            "_tau_1p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.251*1000)]),
                            "_tau_1p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.151*1000)]),
                            "_tau_1p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.101*1000)]),

                            "_tau_2_1": np.array([(i/500) for i in np.arange(0.*500, 0.452*500)]),
                            "_tau_2_2": np.array([(i/500) for i in np.arange(0.*500, 0.252*500)]),
                            "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.121*1000)]),
                            "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.101*1000)]),
                            "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.081*1000)]),

                            "_tau21": np.array([(i/100) for i in np.arange(0.*100, 1.21*100)]),#for one-pass kT minimization as per CMS
                            "_tau32": np.array([(i/200) for i in np.arange(0.*200, 1.21*200)]),#for one-pass kT minimization as per CMS
                            "_tau21_WTA": np.array([(i/100) for i in np.arange(0.*100, 1.01*100)]),#for WTA-kT for comparison
                            "_tau32_WTA": np.array([(i/100) for i in np.arange(0.*100, 1.11*100)]),#for WTA-kT for comparison
                            "_tau21_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.21*200)]),#for excl.-kT and E-scheme as per basis
                            "_tau32_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.41*200)]),#for excl.-kT and E-scheme as per basis
                            #"_mass": np.array([(i*2) for i in np.arange(60., 151.)]),
                            #"_msoftdrop": np.array([(i*2) for i in np.arange(0., 151.)]),
                            #"_pt": np.array([(i*10) for i in np.arange(35., 251.)]),
                            }

        
        ### Uncertainties
        self.btaggingSystLabels = ["central", "up", "down"]
        
        self.sysSources = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSources if not( isys.endswith('nom'))]# or self.sysUnc) ]
        self.sysWeightList = ( '_pu', '_isr', '_fsr', '_pdf', '_l1prefiring', '_lepton', '_btag' ) #'_ps',
        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight','_btagWeightCorrelated', '_btagWeightUncorrelated', '_leptonWeightISO', '_leptonWeightID', '_leptonWeightTrig', '_leptonWeightRecoEff'] if self.wtUnc else [] #, '_leptonWeightAll'
        if self.onlyUnc: self.sysSources = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] #################### Not using for right now
        
        if self.wtUnc: #for the jes/jer cases sys self.sysSources is udpated  above
            self.sysSources = ['_nom'] + [ iwt+i for i in [ 'Up', 'Down' ] for iwt in self.wtSources if not iwt.endswith(('nom','pdfWeightAll')) ]
            
            if 'pdfWeightAll' in self.wtSources: self.sysSources = self.sysSources+['pdfWeightAll'] 
        
        #puWeights used only to change reco event weight (up/down) without application to gen weight, others applied to modify the genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)]
        
        #stupid hacks because this producer was initially planed to be run over multiple (W&top) selections, but is now only use one selection at a time... will change later if I have time
        self.selList = [ selection ] 
        self.sel = selection 
                    
        
        self._branchesToRead = self.getBranchesToRead(
                                                     dirname=self.inputDir, year=self.year,
                                                     kinematic_labels=self.kinematic_labels,
                                                     reco_only_labels=self.reco_only_labels,
                                                     nSub_labels=self.dict_variables_toUnfold,
                                                    )+['Jet_pt', 'Jet_eta', 'Jet_jetId', 'Jet_hadronFlavour', 'Jet_btagDeepFlavB']
        #load up Jet_* branches to calculate btaggingweights (sort of) 'on-the-fly'
        
        
        self.top_ptmin = 400.
        self.top_mSDmin = 140.
        self.top_mSDmax = 250.
        self.top_mINVmin = 140.
        self.top_mINVmax = 250.

        self.W_ptmin = 200.
        self.W_mSDmin = 65.
        self.W_mSDmax = 125.
        self.W_mINVmin = 65.
        self.W_mINVmax = 125.
        
        self.LeptW_ptmin = 150.
        self.Mu_ptmin = 55.
        self.MET_ptmin = 50.
        self.AK4_ptmin = 30.
        self.AK4_etamax = 2.4 #(abs)
        self.AK4_jetID = 2 #(cut on >=2; https://twiki.cern.ch/twiki/bin/view/CMS/JetID#nanoAOD%20Flags)
        self.genWeights = []
    

    def process(self, events):
        '''fill Hist histograms to accumulate event stats, convert to root or whatever else after returned by processor'''
        
        
        if not(self.onlyParquet): output = self.buildDictOfHistograms()
        
        # Get branches to be read 
        branches = self._branchesToRead
        
        # Currently the branches to be read are being defined somewhat circularly, the branches that are to be read are extracted from the class-level variable and only those columns are read in from the nanoskims, but keeping this functionality in for easy generalisability if needed; ie, if one wants to make this class self sufficient and add  in what I do in my histo production noteeboks as a postprocessor in this class. 
        events = events[branches]

        
        ############################### btagging weights, nom+variations ###############################
        
        if self.isMC and not(self.onlyParquet): 
            tstart=time.time()
            
            events_new_fields = self.compute_new_SF_fields(events)#, btvJSON, effHist_b, effHist_c, effHist_udsg)
        
            #if self.verbose: 

            print(f"Adding in new branches {events_new_fields.keys()} to access event weights from b-tagging")

            for key, values in events_new_fields.items():
                events[key] = values
            #if self.verbose: 
            elapsed = time.time() - tstart

            print (f'Finished adding in new b-tagging reweighting branches. Time taken:{elapsed}') 
        
        
        
        if not(self.onlyParquet): 
            print (f'Now, processing histograms; systematic sources being considered:{self.sysSources} for nevents={len(events)} (before masking)')
            if self.verbose: 
                print (output.keys())
        elif self.verbose and (self.onlyParquet):
            print (f'Now, producing .parquet files; for nevents={len(events)} (before masking)')
            
            
        for sys in self.sysSources:
            
            
            s='_nom' if (sys.startswith(self.sysWeightList)) else sys
            
            ############ building event masks ############
            
            if self.mode=='massAndDeltaR': 

                topRecoMask = (events[f'selRecoJets{s}_pt']>self.top_ptmin) & (events[f'selRecoJets{s}_msoftdrop']>self.top_mINVmin) & (events[f'selRecoJets{s}_msoftdrop']<self.top_mINVmax) & (events[f'selRecoHadHemDeltaR{s}']<0.8) & (events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & (events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'recoSelectedEventNumber_nom']>=1) #& (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 

                WRecoMask = (events[f'selRecoJets{s}_pt']>self.W_ptmin) & (events[f'selRecoJets{s}_msoftdrop']>self.W_mINVmin) & (events[f'selRecoJets{s}_msoftdrop']<self.W_mINVmax) & (events[f'selRecoHadHemDeltaR{s}']>0.8) & (events[f'selRecoHadHemDeltaR{s}']<1.6) & (events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & (events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'recoSelectedEventNumber_nom']>=1) #& (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 
                
                if not(self.isMC): 
                    topRecoMask = (topRecoMask) & (events[f'passHLT_Mu50{s}']==1)
                    WRecoMask = (WRecoMask) & (events[f'passHLT_Mu50{s}']==1)
            
                if self.isMC:# or self.isSigMC:    
                    
                    topGenMask = (events[f'selGenJets_nom_pt']>self.top_ptmin) & (events[f'selGenJets_nom_msoftdrop']>self.top_mINVmin) & (events[f'selGenJets_nom_msoftdrop']<self.top_mINVmax) & (events[f'selGenHadHemDeltaR_nom']<0.8) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) #& (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

                    WGenMask = (events[f'selGenJets_nom_pt']>self.W_ptmin) & (events[f'selGenJets_nom_msoftdrop']>self.W_mINVmin) & (events[f'selGenJets_nom_msoftdrop']<self.W_mINVmax) & (events[f'selGenHadHemDeltaR_nom']>0.8)  & (events[f'selGenHadHemDeltaR_nom']<1.6) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) #& (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

            elif self.mode=='onlyDeltaR':

                topRecoMask = (events[f'selRecoJets{s}_pt']>self.top_ptmin) & (events[f'selRecoHadHemDeltaR{s}']<0.8) & (events[f'selRecoLeptW{s}_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu{s}_pt']>self.Mu_ptmin ) & (events[f'selRecoMET{s}_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 

                WRecoMask = (events[f'selRecoJets{s}_pt']>self.W_ptmin) & (events[f'selRecoHadHemDeltaR{s}']>0.8) & (events[f'selRecoHadHemDeltaR{s}']<1.6) & (events[f'selRecoLeptW{s}_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu{s}_pt']>self.Mu_ptmin ) & (events[f'selRecoMET{s}_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 
            
                if self.isMC:# or self.isSigMC:    
                    
                    topGenMask = (events[f'selGenJets_nom_pt']>self.top_ptmin) & (events[f'selGenHadHemDeltaR_nom']<0.8) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) & (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

                    WGenMask = (events[f'selGenJets_nom_pt']>self.W_ptmin) & (events[f'selGenHadHemDeltaR_nom']>0.8)  & (events[f'selGenHadHemDeltaR_nom']<1.6) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) & (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

            elif self.mode=='noDeltaR':

                topRecoMask = (events[f'selRecoJets{s}_pt']>self.top_ptmin) & (events[f'selRecoJets{s}_msoftdrop']>self.top_mINVmin) & (events[f'selRecoJets{s}_msoftdrop']<self.top_mINVmax) & (events[f'selRecoLeptW{s}_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu{s}_pt']>self.Mu_ptmin ) & (events[f'selRecoMET{s}_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 

                WRecoMask = (events[f'selRecoJets{s}_pt']>self.W_ptmin) & (events[f'selRecoJets{s}_msoftdrop']>self.W_mINVmin) & (events[f'selRecoJets{s}_msoftdrop']<self.W_mINVmax) & (events[f'selRecoLeptW{s}_pt']>self.LeptW_ptmin ) & (events[f'selRecoMu{s}_pt']>self.Mu_ptmin ) & (events[f'selRecoMET{s}_pt']>self.MET_ptmin ) & (events[f'passRecoSel{s}']==1) & (events[f'totalRecoWeight{s}']!=0.) & (events[f'nRecoBtags_nom']>1) & (events[f'nRecoHadBtags_nom']==1) 
            
                if self.isMC:# or self.isSigMC:    
                    
                    topGenMask = (events[f'selGenJets_nom_pt']>self.top_ptmin) & (events[f'selGenJets_nom_msoftdrop']>self.top_mINVmin) & (events[f'selGenJets_nom_msoftdrop']<self.top_mINVmax) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) & (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

                    WGenMask = (events[f'selGenJets_nom_pt']>self.W_ptmin) & (events[f'selGenJets_nom_msoftdrop']>self.W_mINVmin) & (events[f'selGenJets_nom_msoftdrop']<self.W_mINVmax) & (events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & (events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & (events[f'selGenMET_nom_pt']>self.MET_ptmin ) & (events[f'passGenSel{s}']==1) & (events[f'evtGenWeight_nom']!=0.) & (events[f'nGenBtags_nom']>1) & (events[f'nGenHadBtags_nom']==1) 

            else: 
                print(f"######################## _WARNING_ ##############################")
                print(f"provided histoproducer mode, #{self.mode}#, is undefined; exiting")
                print(f"#################################################################")
                return -1

            if self.isMC:# or self.isSigMC:    
                if sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                    extraWeights = events[f'{sys.split("_")[1]}{s}']
                
                    self.genWeights = events['evtGenWeight_nom']*extraWeights
                    #genWeightsW = events['genWeight'][WRecoMask]*extraWeights[WRecoMask]
                else:
                    self.genWeights = events['evtGenWeight_nom']
                    #genWeightsW = events['genWeight'][WRecoMask]
           
            
            
            if self.verbose:
                print(s, sys, len(events[topRecoMask]), len(events[topGenMask]), len(events[WRecoMask]), len(events[WGenMask]), len(self.genWeights))

            if self.isSigMC:# and not(sys.startswith(self.sysWeightList)):
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys

                if self.verbose: 
                    print('making individual W and topSel mask',sys,s)

                #if self.withLeptHemBtag==False:

                
                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) #& (~(topRecoMask|topGenMask))
                    WtopGenMask = (WGenMask) #& (~(topRecoMask|topGenMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) #& (~(WRecoMask|WGenMask))
                    WtopGenMask = (topGenMask) #& (~(WGenMask|WRecoMask))    

            
                        
                selRecoMask = (WtopRecoMask) 
                selGenMask = (WtopGenMask) 

                trueRecoMask = (events[f'trueRecoJets{s}_pt']>0.) & ((selRecoMask) & (selGenMask)) 
                accepGenMask = (events[f'accepGenJets{s}_pt']>0.) & ((selGenMask) & (selRecoMask))

                fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) #((selRecoMask) ^ (trueRecoMask))#

                missGenMask =  ((selGenMask) & (~accepGenMask))#((selGenMask) ^ (accepGenMask))  #
                
                #selRecoMask = (trueRecoMask) | (fakeRecoMask) #& (event_mask)
                #selGenMask = (accepGenMask) | (missGenMask) #& (event_mask)
                
                if self.verbose:
                    print(s, sys, len(events[selRecoMask]), len(events[trueRecoMask]), len(events[fakeRecoMask]), len(events[selGenMask]), len(events[accepGenMask]), len(events[missGenMask]),len(self.genWeights))
                
                if self.saveParquet:
                    
                    #events['selRecoMask']=selRecoMask
                    #events['selGenMask']=selGenMask
                    #events['trueRecoMask']=trueRecoMask
                    #events['fakeRecoMask']=fakeRecoMask
                    #events['accepGenMask']=accepGenMask
                    #events['missGenMask']=missGenMask
    
                    if (sys.endswith('nom')): 
                        
                        
                        sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                        
                        print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_{self.splitCount}.parquet'}")

                        ak.to_parquet(events[selRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_reco_{self.splitCount}.parquet')
                        ak.to_parquet(events[selGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_gen_{self.splitCount}.parquet')
                        #ak.to_parquet(events[trueRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_truereco_{self.splitCount}.parquet')
                        #ak.to_parquet(events[accepGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_accepgen_{self.splitCount}.parquet')
                        #ak.to_parquet(events[fakeRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_fakereco_{self.splitCount}.parquet')
                        #ak.to_parquet(events[missGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_missgen_{self.splitCount}.parquet')
                        return 1

                    #elif self.sysUnc and self.onlyUnc: 
                    #    print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_{s}Unc{sel}_{self.splitCount}.parquet'}")

                    #    ak.to_parquet(events, self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_{s}Unc{sel}_{self.splitCount}.parquet')
                    #    return 1
                    

            elif self.isMC and not(self.isSigMC) and sys.endswith('nom'): #not so relevant for dijets but for background MC's in W/top
                

                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) #& (~(topRecoMask|topGenMask))
                    WtopGenMask = (WGenMask) #& (~(topRecoMask|topGenMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) #& (~(WRecoMask|WGenMask))
                    WtopGenMask = (topGenMask) #& (~(WGenMask|WRecoMask))    

            
                        
                selRecoMask = (WtopRecoMask) 
                selGenMask = (WtopGenMask) 

                trueRecoMask = (events[f'trueRecoJets{s}_pt']>0.) & ((selRecoMask) & (selGenMask)) 
                accepGenMask = (events[f'accepGenJets{s}_pt']>0.) & ((selGenMask) & (selRecoMask))

                fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) #((selRecoMask) ^ (trueRecoMask))#

                missGenMask =  ((selGenMask) & (~accepGenMask))#((selGenMask) ^ (accepGenMask))  #
                
                #selRecoMask = (trueRecoMask) | (fakeRecoMask) #& (event_mask)
                #selGenMask = (accepGenMask) | (missGenMask) #& (event_mask)
                
                
                sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                
               

                if self.saveParquet:
                    #events['selRecoMask']=selRecoMask
                    #events['selGenMask']=selGenMask
                    #events['trueRecoMask']=trueRecoMask
                    #events['fakeRecoMask']=fakeRecoMask
                    #events['accepGenMask']=accepGenMask
                    #events['missGenMask']=missGenMask
                        
                    print(f"Saving the following .parquet file with file-stem: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_reco/gen_{self.splitCount}.parquet'}")

                    ak.to_parquet(events[selRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_reco_{self.splitCount}.parquet')
                    ak.to_parquet(events[selGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_gen_{self.splitCount}.parquet')
                    return 1
            
            elif not(self.isMC) and sys.endswith('nom'): 
                
                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) & (~(topRecoMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) & (~(WRecoMask))

                       
                selRecoMask = (WtopRecoMask) 
                
                sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                
                if self.saveParquet:
                    events['selRecoMask']=selRecoMask
                    print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'{sel}_{self.splitCount}.parquet'}")

                    ak.to_parquet(events,self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'{sel}_{self.splitCount}.parquet')
                    return 1
                    
            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config', s,sys,self.isMC,self.isSigMC) 
            
            if self.verbose and self.isMC and not(self.onlyParquet):
                #    print (selRecoMask[0:20],selRecoMask[-20:-1])
                #    print (trueRecoMask[0:20],trueRecoMask[-20:-1])
                #    print (selGenMask[0:20],selGenMask[-20:-1])
                #    print (accepGenMask[0:20],accepGenMask[-20:-1])
                #    
                print(output.keys())
                print("Going to fill histos")
            
            if not(self.onlyParquet): 
                for isel in self.selList:

                    if self.verbose: print(sys,s)

                    if sys.endswith('nom'):
                        listOfMiscOutputHistosNames = [k for k,h in output.items() if ((not ('Jet' in k) and ('nom' in k)) and ('Mu' in k or 'LeptW' in k or 'MET' in k or 'AK4bjetHadHem' in k or 'nPV' in k or 'DeltaR' or 'DeltaPhi' in k or 'Btag' in k))]
                        #btagWeightNom{s}
                        totalRecoWeight = events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']*events[f'leptonWeightNom{s}']*self.genWeights if self.isMC else events[f'totalRecoWeight_nom']#events[f'evtGenWeight_nom']*
                        totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None

                        if self.verbose: print(listOfMiscOutputHistosNames)
                        
                        for histos in listOfMiscOutputHistosNames:
                            key=histos

                            if 'Mu' in key: 
                                dicto=self.dict_variables_kinematics_Muon
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(events[f'selRecoMu{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(events[f'selGenMu{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'MET' in key: 
                                dicto=self.dict_variables_kinematics_MET
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(events[f'selRecoMET{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(events[f'selGenMET{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'LeptW' in key: 
                                dicto=self.dict_variables_kinematics_LeptW
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(events[f'selRecoLeptW{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(events[f'selGenLeptW{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'AK4bjetHadHem' in key and not ('DeltaR' in key or 'DeltaPhi' in key): 
                                dicto=self.dict_variables_kinematics_AK4Had
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(events[f'selRecoAK4bjetHadHem{s}{i}'][selRecoMask], weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower() and self.isMC:
                                            output[key].fill(events[f'selGenAK4bjetHadHem{s}{i}'][selGenMask], weight=totalGenWeight[selGenMask], threads=8)

                            elif 'DeltaRap' in key:#
                                if 'reco' in key.lower():
                                    output[key].fill(events[f'selRecoHadHemDeltaRap{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(events[f'selGenHadHemDeltaRap{s}'][selGenMask], weight=totalGenWeight[selGenMask])

                            elif 'DeltaR' in key and not('DeltaRap' in key):
                                if 'reco' in key.lower():
                                    output[key].fill(events[f'selRecoHadHemDeltaR{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(events[f'selGenHadHemDeltaR{s}'][selGenMask], weight=totalGenWeight[selGenMask])

                            elif 'DeltaPhi' in key:
                                if 'reco' in key.lower() and self.isMC:
                                    output[key].fill(events[f'selRecoHadHemDeltaPhi{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(events[f'selGenHadHemDeltaPhi{s}'][selGenMask], weight=totalGenWeight[selGenMask])


                            elif 'good_nPVs' in key:
                                output[key].fill(events[f'good_nPVs{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nRecoBtags' in key:
                                output[key].fill(events[f'nRecoBtags{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenBtags' in key and self.isMC:
                                output[key].fill(events[f'nGenBtags{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                            elif 'nRecoHadBtags' in key:
                                output[key].fill(events[f'nRecoHadBtags{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenHadBtags' in key and self.isMC:
                                output[key].fill(events[f'nGenHadBtags{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                            elif 'nRecoLepBtags' in key:
                                output[key].fill(events[f'nRecoLepBtags{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenLepBtags' in key and self.isMC:
                                output[key].fill(events[f'nGenLepBtags{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                    listOfJetOutputHistosNames = [k for k,h in output.items() if ((('Jet' in k) and (sys in k)) or ('residual' in k and not(sys in k)) or ('resol' in k and not(sys in k)))] #prevent needless iterations
                    if self.verbose: print (listOfJetOutputHistosNames[0:10],sys,s)

                    for k in listOfJetOutputHistosNames:
                        key=k

                        ############### Safety checks ##################
                        if not(sys in key) and not('residual' in key) and not('resol' in key): 
                            if self.verbose: 
                                print(key,sys,s)
                            continue

                        if not('_tau' in key):
                            whichKinVar = [var for var in self.dict_variables_kinematics_AK8.keys() if var[1:] in key]#[0]
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
                            if self.sampleName.startswith('TTToSemi') and self.year=='2016_preVFP' and self.verbose:
                                print(whichnSub,key,keep,temp)
                            varToFill = keep
                        
                        elif (('_tau' in key and not('21' in key or '32' in key)) or '__' in key):# and 'nom' in key):
                            temp='tau'
                            #temp='exkT' if 'exkT' in key else 'tau'
                            whichnSub = [var for var in self.dict_variables_toUnfold.keys() if (var in key and temp in var)]
                            varToFill = whichnSub[0]


                        ################## Filling histos from accumulated event arrays ##################
                        if self.verbose: print(key,sys,s,varToFill)
                        # Decide what weights to use
                        if not (sys.startswith(self.sysWeightList)):
                            s=sys
                            totalRecoWeight = self.genWeights*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']*events[f'leptonWeightNom{s}'] if self.isMC else events[f'totalRecoWeight_nom']
                            totalGenWeight = events[f'evtGenWeight_nom'] if self.isMC else None
                        else:
                            if self.isSigMC:
                                s='_nom'
                                #print(f'{sys.split("_")[1]}{s}')
                                if 'pu' in sys:
                                    totalGenWeight = events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*events[f'{sys.split("_")[1]}{s}']*events[f'l1prefiringWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']*events[f'leptonWeightNom{s}']

                                elif 'l1' in sys:
                                    totalGenWeight = events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']*events[f'leptonWeightNom{s}']

                                elif '_btag' in sys:
                                    totalGenWeight = events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*events[f'new{self.sel}_{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'leptonWeightNom{s}']

                                elif 'lepton' in sys:
                                    totalGenWeight = events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*events[f'{sys.split("_")[1]}{s}']*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']

                                elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                                    totalGenWeight = events[f'evtGenWeight{s}']*events[f'{sys.split("_")[1]}{s}']
                                    totalRecoWeight = self.genWeights*events[f'puWeightNom{s}']*events[f'l1prefiringWeightNom{s}']*events[f'new{self.sel}_btagWeightNom{s}']*events[f'leptonWeightNom{s}']

                        if (key.lower().startswith('reco')):
                            if self.verbose: 
                                print('all reco', key,s,sys, len(events[f'selRecoJets{s}{varToFill}'][selRecoMask]))
                            output[key].fill(events[f'selRecoJets{s}{varToFill}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)

                        # Stuff below should have diff wts already available to them so no need to touch anything down here
                        if self.isMC:# and (not varToFill.startswith(tuple(self.reco_only_labels))):

                            if (key.lower().startswith('true')) and self.isSigMC :
                                if self.verbose: print('true', key,s,sys, len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask]))
                                output[key].fill(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                            elif (key.lower().startswith('fake')) and self.isSigMC:
                                if self.verbose: print('fake', key,s,sys,len(events[f'selRecoJets{s}{varToFill}'][fakeRecoMask]))
                                output[key].fill(events[f'selRecoJets{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                                 threads=8)

                            if (key.lower().startswith('genjet')):
                                if self.verbose: 
                                    print('all gen', key,s,sys, len(events[f'selGenJets_nom{varToFill}'][selGenMask]))

                                output[key].fill(events[f'selGenJets_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                                 threads=8)

                            elif (key.lower().startswith('accepgen')):
                                if self.verbose: 
                                    print('accepgen', key,s,sys, len(events[f'accepGenJets{s}{varToFill}'][accepGenMask]))
                                output[key].fill(events[f'accepGenJets{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                                 threads=8)

                            elif (key.lower().startswith('missgenjet')):
                                if self.verbose: 
                                    print('missgen', key,s,sys, len(events[f'selGenJets_nom{varToFill}'][missGenMask]))
                                output[key].fill(events[f'selGenJets_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                                 threads=8)

                            elif ( 'resp' in key.lower() and not('withmiss' in key.lower())) and self.isSigMC:
                                #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                                output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                                output[key].fill(gen=events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                                 weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                                 threads=8)     

                            elif ('respwithmiss' in key.lower()) and self.isSigMC:
                                if self.verbose: print("filling resp w. misses")
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


                            elif ('residual' in key.lower() or 'resol' in key.lower()) and sys.endswith('nom') and self.isSigMC:
                                zeroMask=(events[f'accepGenJets{s}{varToFill}']!=0.)&(accepGenMask)
                                
                                response = events[f'trueRecoJets{s}{varToFill}'][zeroMask]/events[f'accepGenJets{s}{varToFill}'][zeroMask]
                                response = np.nan_to_num(response,nan=-999.)
                                residual = events[f'trueRecoJets{s}{varToFill}'][zeroMask]-events[f'accepGenJets{s}{varToFill}'][zeroMask]
                                #residual = np.nan_to_num(residual, nan=-929.)
                                
                                if 'noWt_' in key: output[key].fill(response)#, weight=totalRecoWeight[zeroMask])
                                elif 'residual' in key: output[key].fill(residual)
                            

                print("Histos filled!",self.selList,sys)
                l=[]
                for x,y in output.items(): #y.SetDirectory(0)

                    if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss')) and '_nom' in x:
                        l.append(x)
                    #elif self.isMC and not(self.isSigMC) and x.startswith(('accepgen','miss','true','fake','residual' )) and '_nom' in x:
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

                for sysUnc in self.sysSources:
                    
                    for x, y in self.dict_variables_kinematics_AK8.items():
                        binning = y
                        if sysUnc.endswith('nom') or self.onlyUnc:
                            if not x in tuple(self.reco_only_labels): 
                                dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())        
                            else: 
                                dictOfHists[x[1:]+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=x[1:]+sysUnc+isel, label=f'AK8 {itype} jet {x}', underflow=True,overflow=True).Weight())  
                    if itype.startswith('truereco')  and self.isSigMC and sysUnc.endswith('nom'):
                        dictOfHists['residual'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_pt'+isel, label='AK8 reco - gen jet pt').Weight())
                        dictOfHists['residual'+iJ+'_msoftdrop'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_msoftdrop'+isel, label='AK8 reco m_{SD} - gen jet m_{SD}').Weight())
                        dictOfHists['residual'+iJ+'_mass'+isel] = (hist.Hist.new.Regular(2000, -10, 10, name='residual'+iJ+'_mass'+isel, label='AK8 reco inv. m - gen jet inv. m').Weight())
                    
                    #print('building unfolding histos')     
                    for x, y in self.dict_variables_toUnfold.items():
                        binning = y
                        
                        #binning_coarse=np.array([binning[i] for i in range(len(binning)) if i%5==0])
                        
                        dictOfHists[itype+iJ+x+sysUnc+isel] = (hist.Hist.new.Variable(binning,name=itype+iJ+x+sysUnc+isel, label='AK8 '+itype+' jet #tau', underflow=True,overflow=True).Weight())
                        
                        if itype.startswith('truereco') and self.isMC: #and self.isSigMC:

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
                                dictOfHists['residual'+iJ+x+isel] = (hist.Hist.new.Regular(6000, -1.5, 1.5, name='residual'+iJ+x+isel, label='(Ak8 reco - gen jet) ', underflow=True,overflow=True).Weight())
                                dictOfHists['noWt_resol'+iJ+x+isel] = (hist.Hist.new.Regular(1000, 0, 8, name='noWt_resol'+iJ+x+isel, label='AK8 reco/gen jet ', underflow=True,overflow=True).Weight())
        
        
        #build from list of histos for control plots
        if not(self.onlyUnc or self.sysUnc):
            for isel in self.selList:
                sysUnc= '_nom'

                #if self.verbose:    print (self.listOfControlHistTypes)

                for itype in self.listOfControlHistTypes:
                    for iO in self.selObjects:

                        #making control histos only for nominal case for now
                            
                        if 'mu' in iO.lower():
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

                #if self.verbose: print (self.dict_variables_reco.keys())

                for itype in self.dict_variables_reco:
                    bins=self.dict_variables_reco[itype]
                    
                    if 'good_nPVs' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'Good nPVs').Weight())  
                    elif 'nRecoBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco').Weight())  
                    elif 'nRecoHadBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco (hadr. hem.)').Weight())  
                    elif 'nRecoLepBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco (lept. hem.)').Weight())  
                    elif 'selRecoHadHemDeltaRap' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dy(AK8, AK4)').Weight())  
                    elif 'selRecoHadHemDeltaR' in itype and not('DeltaRap' in itype): dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, AK4)').Weight())  
                    elif 'selRecoHadHemDeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(AK8, AK4)').Weight())  

                if self.isMC or self.isSigMC:
                    #if self.verbose: print (self.dict_variables_gen.keys())

                    for itype in self.dict_variables_gen:
                        bins=self.dict_variables_gen[itype]
                        
                        if 'nGenBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen').Weight())  
                        elif 'nGenHadBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen (hadr. hem.)').Weight())  
                        elif 'nGenLepBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen (lept. hem.)').Weight())  
                        elif 'selGenHadHemDeltaRap' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dy(AK8, AK4)').Weight())  
                        elif 'selGenHadHemDeltaR' in itype and not('DeltaRap' in itype): dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, AK4)').Weight())  
                        elif 'selGenHadHemDeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(AK8, AK4)').Weight())  
        
        #if self.verbose: print (dictOfHists)
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
            self.wtUnc=False
                    
        print (self.sysSources)

        for sys in self.sysSources:
            if 'btag' in sys: continue #this is since the btagWeight branches are now generated on the fly vs. being read in
            #if not wtUnc and not sysUnc:
            for i in kinematic_labels+list(nSub_labels.keys()): 
                
                if sys.endswith('nom') or (self.isSigMC and self.sysUnc): 
                    reco_list.append('selRecoJets'+sys+i)
                
                if not ('softdrop' in i) and not('tau' in i):
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                        if i in self.dict_variables_kinematics_Muon and sys.endswith('nom'):
                            reco_list.append('selRecoMu'+sys+i)
                        if i in self.dict_variables_kinematics_MET and sys.endswith('nom'):
                            reco_list.append('selRecoMET'+sys+i)
                        if i in self.dict_variables_kinematics_LeptW and sys.endswith('nom'):
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
                                if i in self.dict_variables_kinematics_Muon and sys.endswith('nom'):
                                    gen_list.append('selGenMu'+sys+i)
                                if i in self.dict_variables_kinematics_MET and sys.endswith('nom'):
                                    gen_list.append('selGenMET'+sys+i)
                                if i in self.dict_variables_kinematics_LeptW and sys.endswith('nom'):
                                    gen_list.append('selGenLeptW'+sys+i)
                                if i in self.dict_variables_kinematics_AK4Had:
                                    gen_list.append('selGenAK4bjetHadHem'+sys+i)

            if self.isMC and (sys.endswith('nom') or (self.isSigMC and self.sysUnc)):
                reco_list.append("puWeightNom"+sys)
                reco_list.append("btagWeightNom"+sys)
                
                reco_list.append("leptonWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("pdfWeightNom"+sys)

                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("passRecoSel"+sys)
                
                ###
                if not(self.onlyParquet) and (self.isSigMC and self.wtUnc) and not(self.sysUnc): 
                    for ud in ['Up','Down']:
                        reco_list.append("leptonWeightISO"+ud+sys)
                        reco_list.append("leptonWeightID"+ud+sys)
                        reco_list.append("leptonWeightTrig"+ud+sys)
                        reco_list.append("leptonWeightRecoEff"+ud+sys)
                ###
                
                if self.isMC: 
                    reco_list.append(f"selRecoJets{sys}_msoftdrop_corr_PUPPI")#+sys)
                    
                reco_list.append("FlagRecoLeptHemBjet"+sys)
                reco_list.append("selRecoHadHemDeltaR"+sys)
                reco_list.append("selRecoHadHemDeltaPhi"+sys)
                reco_list.append("selRecoHadHemDeltaRap"+sys)

                gen_list.append("evtGenWeight"+sys)
                gen_list.append("passGenSel"+sys)
                gen_list.append("FlagGenLeptHemBjet"+'_nom')
                gen_list.append("selGenHadHemDeltaR"+sys)
                gen_list.append("selGenHadHemDeltaPhi"+sys)
                gen_list.append("selGenHadHemDeltaRap"+sys)
                
                
                for i in self.dict_variables_reco:
                    reco_list.append(i+sys)
                for i in self.dict_variables_gen:
                    gen_list.append(i+(sys if not('nGen' in i) else '_nom') )
                    
            elif (not self.isMC) and sys.endswith('nom'):
                reco_list.append("totalRecoWeight"+sys)               
                reco_list.append("passRecoSel"+sys)               
                reco_list.append("passHLT_Mu50"+sys)               
                reco_list.append("selRecoHadHemDeltaR"+sys)
                reco_list.append("selRecoHadHemDeltaPhi"+sys)
                reco_list.append("selRecoHadHemDeltaRap"+sys) 
                for i in self.dict_variables_reco:
                    reco_list.append(i+sys)
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or  self.sysUnc):
                if 'pu' in sys or 'l1' in sys or 'lepton' in sys or 'btag' in sys: reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: gen_reweights_list.append(sys.split('_')[1]+'_nom')

        if self.isSigMC: 
            branchesToRead=gen_list+reco_list+gen_reweights_list+reco_reweights_list+['genWeight','recoSelectedEventNumber_nom','genSelectedEventNumber_nom']#,'recoEventCategory_nom','genEventCategory_nom']
        elif self.isMC and (not self.isSigMC): 
            branchesToRead=gen_list+reco_list+['genWeight','recoSelectedEventNumber_nom','genSelectedEventNumber_nom']#,'recoEventCategory_nom','genEventCategory_nom']
        elif (not self.isMC): 
            branchesToRead=reco_list+['recoSelectedEventNumber_nom']#,'recoEventCategory_nom']

        
        if self.verbose: 
            #print(branchesToRead)
            print(reco_reweights_list)
            print(gen_reweights_list)

        return branchesToRead
    
    def compute_new_SF_fields(self, events):#, btvjson, effHist_b, effHist_c, effHist_udsg):
        
        # the setting of the keys for the systematic suffix map can be raised to the global/class level but keeping things hard-coded here since my needs are fixed and I use custom nomenclature for self-created branches in ntuples, and it meshes better with how I do my histo productions
        
        #for non-signal MC, signal MC sys variations on jes/jer, or signal MC when run without exp/theory weight unc. variations
        if (self.isMC and not(self.isSigMC)) or (self.isSigMC and not(self.wtUnc)): 
            
            systematic_suffix_map = {
                                     'central': 'Nom',
                                    }
        
        elif (self.isSigMC and self.wtUnc):
            
            systematic_suffix_map = {
                                     'central': 'Nom',
                                     'up_correlated': 'CorrelatedUp',
                                     'down_correlated': 'CorrelatedDown',
                                     'up_uncorrelated': 'UncorrelatedUp',
                                     'down_uncorrelated': 'UncorrelatedDown'
                                    }
        
        
        new_fields = OrderedDict()
        
        
        for s in systematic_suffix_map.keys(): 
            if not self.sysUnc: 
                new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[s]}_nom'] = []
            else:
                if s!='central': continue
                    
                for sys in self.sysSources:
                    new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[s]}{sys}'] = []
            
        #new_fields = {
        #                'new{self.sel}_btagWeightNom_nom': [],
        #                'new{self.sel}_btagWeightUp_nom': [],
        #                'new{self.sel}_btagWeightDown_nom': []
        #             }

        if not(self.sysUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        print(f"Computing btagging weights ({systematic_suffix_map.keys()}) for following sys sources: {sysList}")
        
        
        for i in range(0, len(events)):
                
            for sys in sysList:     
                
                if events[f'passRecoSel{sys}'][i]!=1: 
                    for systematic in systematic_suffix_map.keys():
                        new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(1.)
                else:

                    jet_pt = events['Jet_pt'][i]
                    jet_eta = events['Jet_eta'][i]
                    jet_ID = events['Jet_jetId'][i]
                    jet_flav = events['Jet_hadronFlavour'][i]
                    jet_discr = events['Jet_btagDeepFlavB'][i]

                    # Define the masks
                    ID_mask = jet_ID>=self.AK4_jetID
                    discr_mask = (jet_discr >= self.WPBDisc)
                    ptMask = (jet_pt > self.AK4_ptmin)
                    etaMask =((jet_eta<=self.AK4_etamax) & (jet_eta>=-1.*self.AK4_etamax))# (np.abs(jet_eta) <= self.AK4_etamax)
                    jet_eta = np.abs(jet_eta)
                    
                    # Apply the masks
                    pt, eta, flav = mask(ptMask, etaMask, discr_mask, ID_mask, jet_pt, jet_eta, jet_flav)

                    # Identify b, c, and light jets
                    b_jets = np.where(flav == 5)
                    c_jets = np.where(flav == 4)
                    light_jets = np.where(flav == 0)

                    # Get efficiencies
                    b_efficiencies = np.array([get_efficiency(self.effHist_b, p, e) for p, e in zip(pt[b_jets], eta[b_jets])])
                    c_efficiencies = np.array([get_efficiency(self.effHist_c, p, e) for p, e in zip(pt[c_jets], eta[c_jets])])
                    light_efficiencies = np.array([get_efficiency(self.effHist_udsg, p, e) for p, e in zip(pt[light_jets], eta[light_jets])])

                    # Compute SFs for each systematic
                    for systematic in systematic_suffix_map.keys():#self.btaggingSystLabels:
                        _, _, _, SFs = compute_btagSFAndWeight(self.btvjson, systematic, self.btagWP, 
                                                               flav, eta, pt, 
                                                               b_jets, c_jets, light_jets,
                                                               b_efficiencies, c_efficiencies, light_efficiencies)

                        new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(np.prod(SFs))
    

        # Convert lists to arrays
        for key in new_fields:
            new_fields[key] = np.array(new_fields[key])

        return new_fields
        
        
#####################################################
###################### Helpers ######################
#####################################################
def compute_btagSFAndWeight(btvjson, systematic, WP,
                            flav, eta, pt, 
                            b_jets, c_jets, light_jets,
                            b_efficiencies, c_efficiencies, light_efficiencies
                           ):
    bJet_SFs = btvjson["deepJet_comb"].evaluate(systematic, WP, 
                                               np.array(flav[b_jets]), 
                                               np.array(eta[b_jets]), 
                                               np.array(pt[b_jets]))
    cJet_SFs = btvjson["deepJet_comb"].evaluate(systematic, WP, 
                                               np.array(flav[c_jets]), 
                                               np.array(eta[c_jets]), 
                                               np.array(pt[c_jets]))
    lightJet_SFs = btvjson["deepJet_incl"].evaluate(systematic, WP, 
                                                    np.array(flav[light_jets]), 
                                                    np.array(eta[light_jets]), 
                                                    np.array(pt[light_jets]))

    bJet_SFs *= b_efficiencies 
    cJet_SFs *= c_efficiencies
    lightJet_SFs *= light_efficiencies
    weight = np.concatenate((bJet_SFs, cJet_SFs, lightJet_SFs))

    return bJet_SFs, cJet_SFs, lightJet_SFs, weight

def mask(ptMask, etaMask, discrMask, idMask, pt, eta, flav):
    
    #print(ptMask, etaMask, discrMask, pt, eta, flav)
    
    try:
        masking = ptMask & etaMask & discrMask & idMask
        return pt[masking],eta[masking],flav[masking]
    except TypeError:
        if type(pt)==float and (type(eta)==float or type(eta)==np.float64) and type(flav)==int:
            #masking = ptMask & etaMask & discrMask & idMask
            #print(pt, eta, flav)
            ptMask, etaMask, discrMask, idMask, pt, eta, flav = np.array([ptMask]), np.array([etaMask]), np.array([discrMask]), np.array([idMask]), np.array([pt]), np.array([eta]), np.array([flav])
            
            masking = ptMask & etaMask & discrMask & idMask

            return pt[masking],eta[masking],flav[masking]
        else:
            
            print("Unhandled error occured for these inputs", ptMask, etaMask, discrMask, idMask, pt,eta,flav)
            print(type(ptMask), type(etaMask), type(discrMask), type(idMask), type(pt),type(eta),type(flav))
            return [],[],[]
    
    #masking = ptMask & etaMask & discrMask & idMask
    
    #if len(pt)>=1 and len(eta)>=1 and len(flav)>=1:
    #    masking = ptMask & etaMask & discrMask & idMask
    #else: 
    #    if type(pt)==float and type(eta)==float and type(flav)==float:
    #        #masking = ptMask & etaMask & discrMask & idMask
    #        pt, eta, flav = list(pt),list(eta),list(flav)
    #try:
    
    #return pt[masking],eta[masking],flav[masking]
    
    #except TypeError:
    #    try: 
    #        print("Shapes:", pt.shape, eta.shape, flav.shape)
    #    except AttributeError:
    #        print("Values:", pt, eta, flav)
    #    return [pt], [eta], [flav]
def get_efficiency(hist, pt, eta):
    # Extract the bin edges for eta and pt
    eta_edges = hist.axis('x').edges()
    pt_edges = hist.axis('y').edges()

    # Find the bin index for eta and pt
    binx = np.digitize(eta, eta_edges) - 1
    biny = np.digitize(pt, pt_edges) - 1

    # Handle overflow bins <=need to check this again maybe?
    binx = np.clip(binx, 0, len(eta_edges) - 2)
    biny = np.clip(biny, 0, len(pt_edges) - 2)

    # Extract relevant the bin contents
    bin_contents_eff = hist.values()

    # Return the said bin content
    return bin_contents_eff[binx, biny]


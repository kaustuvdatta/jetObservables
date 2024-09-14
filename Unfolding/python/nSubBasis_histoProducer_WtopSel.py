import time,re,uproot,hist,correctionlib,numba,gc,copy,os,glob
from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict
from collections import OrderedDict
import numpy as np
import awkward as ak
from hist import Hist
import ROOT 
from collections import defaultdict


class nSubBasis_unfoldingHistoProd_WtopSel():#processor.ProcessorABC):

    #####################################################################################################################################
    #####################################################################################################################################     #                                                                                                                                   # 
    #Relevant info wrt. what mode to run the W/top sample processing on, pragmatic choice: noDeltaR, 'purest'=massAndDeltaR (a la AL,SM)#
    # _________________________________________________________________________________________________________________________________ #
    #|      selection mode      |                                      details                                                         |#
    #|---------------------------------------------------------------------------------------------------------------------------------|#
    #|       massAndDeltaR      |              AK8jet pt and softdropped mass and dR(muon,lept. hem. b-jet) cut                        |#
    #|          noDeltaR        |        AK8jet pt and softdropped mass cut [& dR(muon,lept. hem. b-jet)<=1.6 by default]              |#
    #|---------------------------------------------------------------------------------------------------------------------------------|#
    #|_________________________________________________________________________________________________________________________________|#
    #####################################################################################################################################
    #####################################################################################################################################
    
    def __init__(
                 self, sampleName, selection='_topSel', sysSources=[], year='2017', era='', 
                 isMC=True, isSigMC=True, onlyUnc='', wtUnc=False, 
                 verbose=False, saveParquet=False, onlyParquet=False, test=False,
                 splitCount='0', sampleDict=dictSamples, sysUnc=False, 
                 selectionmode='massAndDeltaR', massToCutOn = 'msoftdrop',
                 btvJSON_fn=None, effMap=None,btagWP="M", withBTaggingWeights=True, #_b=None, effMap_c=None, effMap_udsg=None, 
                 withLeptonWeights=True,
                 muIDISO_JSON_fn=None,muRecoEff_JSON_fn=None,muTrig_JSON_fn=None,profile_time=True,
                 top_ptmin = 400.,top_mSDmin = 140.,top_mSDmax = 300.,top_mINVmin = 140.,top_mINVmax = 300.,
                 W_ptmin = 200.,W_mSDmin = 65.,W_mSDmax = 125.,W_mINVmin = 65.,W_mINVmax = 125.,
                 LeptW_ptmin = 100.,Mu_ptmin = 55., Mu_tkRelIso=0.05, MET_ptmin = 30.,
                 AK4_ptmin = 30.,AK4_etamax = 2.4 ,AK4_jetID = 2 ,dR_max = 1.6,dR_min = 0.2,                 
                ):
    
        self.test=test
        self.events = None

        self.profile_time=profile_time
        self.verbose=verbose
        self.saveParquet=saveParquet
        self.onlyParquet=onlyParquet
        
        self.year = year
        self.isMC = isMC or isSigMC
        self.isSigMC = isSigMC
        self.era = era
        self.splitCount=splitCount

        self.onlyUnc = onlyUnc
        self.wtUnc = wtUnc
        self.sysUnc = sysUnc  #somewhat confusing maybe, but this refers to the mode of the producer for jer/jes unc. variation histos
        
        self.massToCutOn = massToCutOn #choices: 'msoftdrop', 'mass' (i.e., PUPPI corrected softdrop mass or the invariant mass of the AK8)
        self.mode=selectionmode

        self.dictSamples = sampleDict
        self.sampleName = sampleName
        
        
        self.muIDISO_JSON_fn = muIDISO_JSON_fn
        self.muRecoEff_JSON_fn = muRecoEff_JSON_fn
        self.muTrig_JSON_fn = muTrig_JSON_fn


        self.btvJSON_fn = btvJSON_fn
        self.effMap = effMap
        #self.effFile = correctionlib.CorrectionSet.from_file(self.effMap)
        #self.effMap_c = effMap_c
        #self.effMap_udsg = effMap_udsg
        self.btagWP = btagWP
        self.withBTaggingWeights = withBTaggingWeights
        self.withLeptonWeights = withLeptonWeights
        
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
        if self.isMC: # and not(self.onlyParquet)
            
            if self.withLeptonWeights:
                if self.muIDISO_JSON_fn:
                    self.muIDISOjson = correctionlib.CorrectionSet.from_file(self.muIDISO_JSON_fn)
                else:
                    raise ValueError("muIDISO_JSON_fn is not provided!")

                if self.muRecoEff_JSON_fn:
                    self.muRecoEffjson = correctionlib.CorrectionSet.from_file(self.muRecoEff_JSON_fn)
                else:
                    raise ValueError("muRecoEff_JSON_fn is not provided!")

                if self.muTrig_JSON_fn:
                    self.muTrigjson = correctionlib.CorrectionSet.from_file(self.muTrig_JSON_fn)
                else:
                    raise ValueError("muTrig_JSON_fn is not provided!")
            
            if self.withBTaggingWeights: #
                if self.btvJSON_fn:
                    self.btvjson = correctionlib.CorrectionSet.from_file(self.btvJSON_fn)
                else:
                    raise ValueError("btvJSON_fn is not provided!")

                if self.effMap:#_b:
                    self.effCorr_b = correctionlib.CorrectionSet.from_file(self.effMap)#.Get(f'efficiency_b_b_{self.year}{selection}') #effFile_b[f'efficiency_b_b_{self.year}{selection}']
                    self.effCorr_c = correctionlib.CorrectionSet.from_file(self.effMap)#.Get(f'efficiency_b_c_{self.year}{selection}') #effFile_c[f'efficiency_b_c_{self.year}{selection}']
                    self.effCorr_udsg = correctionlib.CorrectionSet.from_file(self.effMap)#.Get(f'efficiency_b_udsg_{self.year}{selection}') #effFile_udsg[f'efficiency_b_udsg_{self.year}{selection}']
                else:
                    raise ValueError("effMap JSON is not provided!")

                #if self.effMap_c:
                #    self.effFile_c = self.effMap#uproot.open(self.effMap_c)
                    
                #else:
                #    raise ValueError("effMap_c JSON is not provided!")

                #if self.effMap_udsg:
                #    self.effFile_udsg = self.effMap#uproot.open(self.effMap_udsg)
                    
                #else:
                #    raise ValueError("effMap_udsg JSON is not provided!")


        
        if (not self.isMC) and self.era=='': print (f'Data-loading error: You need to specify what era if you want to work with while handling data' )
        
        
        ### Helpers
        self.listOfUnfHistTypes =  [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] 

        self.listOfControlHistTypes =  [ 'gen',  'reco' ] if self.isMC else [ 'reco' ] 
        
        self.inputDir = checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'] if self.isMC else checkDict(sampleName,self.dictSamples)[self.year]['t3_dirs'][self.era]
        
        self.nJet = [ 'Jet']
        self.selObjects = ['Mu','LeptW', 'LeptTop','MET','AK4bjetLeptHem']

        self.kinematic_labels=['_pt','_eta', '_y', '_phi', '_mass', '_msoftdrop', '_mt', '_mt1', '_p', '_ptRel', '_tkRelIso', '_jetId', '_hadronFlavour']
        self.reco_only_labels=['good_nPVs','nRecoBtags','nRecoLeptBtags', 
                               'selRecoLeptHemDeltaR', 'selRecoLeptHemDeltaRap', 'selRecoLeptHemDeltaPhi',
                               'selRecoLeptHemAK8DeltaPhi','selRecoLeptHemAK8DeltaR']
        #self.gen_only_labels=['nGenBtags','nGenHadBtags','nGenLepBtags','selGenHadHemDeltaR','selGenHadHemDeltaRap','selGenHadHemDeltaPhi']

        self.dict_variables_kinematics_AK8 = {

                                            "_pt": np.array([i for i in np.arange(170., 3010., 10.)]) if 'WSel' in selection else  np.array([i for i in np.arange(300., 3010., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 205., 5.)]) if 'WSel' in selection else  np.array([i for i in np.arange(80., 405., 5.)]),
                                            "_msoftdrop": np.array([i for i in np.arange(0., 145., 5.)]) if 'WSel' in selection else  np.array([i for i in np.arange(100., 305., 5.)]),

                                        }  

        self.dict_variables_kinematics_Muon = {

                                            "_pt": np.array([i for i in np.arange(50., 1210., 5.)]),
                                            "_p": np.array([i for i in np.arange(0., 1210., 5.)]),
                                            "_ptRel": np.array([i for i in np.arange(0., 410., 5.)]),
                                            "_tkRelIso": np.array([i for i in np.arange(0., 0.1, 0.002)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 5., 0.05)]),

                                        } 

        self.dict_variables_kinematics_LeptW = {

                                            "_pt": np.array([i for i in np.arange(80., 2510., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 705., 5.)]),
                                            "_mt": np.array([i for i in np.arange(0., 1810., 5.)]),
                                            "_mt1": np.array([i for i in np.arange(0., 1810., 5.)]),

                                        } 

        self.dict_variables_kinematics_LeptTop = {

                                            "_pt": np.array([i for i in np.arange(100., 2510., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 805., 5.)]),
                                            #"_mt": np.array([i for i in np.arange(0., 1810., 5.)]),

                                        } 

        self.dict_variables_kinematics_AK4Lept = {

                                            "_pt": np.array([i for i in np.arange(20., 2030., 10.)]),
                                            "_eta": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_y": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                            "_mass": np.array([i for i in np.arange(0., 52.5, 2.5)]),
                                            "_jetId": np.array([i for i in np.arange(0, 10, 1)]),
                                            "_hadronFlavour": np.array([i for i in np.arange(-7,8,1)]),


                                        }

        self.dict_variables_kinematics_MET = {

                                            "_pt": np.array([i for i in np.arange(10., 1500., 10.)]),
                                            "_phi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),

                                        } 

        self.dict_variables_reco = { 
                                        "good_nPVs": np.array([i for i in np.arange(0., 101, 1.)]),
                                        "nRecoBtags": np.array([i for i in np.arange(0., 6, 1.)]),
                                        "nRecoAK4s": np.array([i for i in np.arange(0, 14, 1)]),
                                        "nRecoLeptBtags": np.array([i for i in np.arange(0., 6, 1.)]),
                                        "nselRecoJets": np.array([i for i in np.arange(0., 10, 1.)]),
                                        "selRecoLeptHemDeltaR": np.array([i for i in np.arange(0., 1.85, 0.05)]),
                                        "selRecoLeptHemDeltaRap": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                        "selRecoLeptHemDeltaPhi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                        "selRecoLeptHemAK8DeltaPhi": np.array([i for i in np.arange(-3.3, 3.4, 0.1)]),
                                        "selRecoLeptHemAK8DeltaR": np.array([i for i in np.arange(0., 4.1, 0.1)])

                                   }   
        
        self.dict_variables_gen = { 
                                        "nGenBtags": np.array([i for i in np.arange(0., 6., 1.)]),
                                        "nGenAK4s": np.array([i for i in np.arange(0, 14, 1)]),
                                        "nGenLeptBtags": np.array([i for i in np.arange(0., 6., 1.)]),
                                        "nselGenJets": np.array([i for i in np.arange(0., 10., 1.)]),
                                        "selGenLeptHemDeltaR": np.array([i for i in np.arange(0., 1.85, 0.05)]),
                                        "selGenLeptHemDeltaRap": np.array([i for i in np.arange(-2.4, 2.6, 0.2)]),
                                        "selGenLeptHemDeltaPhi": np.array([i for i in np.arange(-3.4, 3.6, 0.2)]),
                                        "selGenLeptHemAK8DeltaPhi": np.array([i for i in np.arange(-3.3, 3.4, 0.1)]),
                                        "selGenLeptHemAK8DeltaR": np.array([i for i in np.arange(0., 4.1, 0.1)])
                                   } 
        
        self.recoMask = None
        if self.isMC:# and self.onlyParquet:
            self.genMask = None
            self.acceptedGenMask = None
            self.trueRecoMask = None
            self.recoWeight = None
            self.genWeight = None
            
            

        if 'WSel' in selection:
            self.dict_variables_toUnfold = {

                                            "_tau_0p25_1": np.array([(i/200) for i in np.arange(0.*200, 1.*200)]),
                                            "_tau_0p25_2": np.array([(i/500) for i in np.arange(0.*500, 1.*500)]),
                                            "_tau_0p25_3": np.array([(i/500) for i in np.arange(0.*500, 0.98*500)]),
                                            "_tau_0p25_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.96*1000)]),
                                            "_tau_0p25_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.94*1000)]),

                                            "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.*200, 0.96*200)]),
                                            "_tau_0p5_2": np.array([(i/500) for i in np.arange(0.*500, 0.92*500)]),
                                            "_tau_0p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.82*500)]),
                                            "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.74*1000)]),
                                            "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.64*1000)]),

                                            "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.9*200)]),
                                            "_tau_1_2": np.array([(i/500) for i in np.arange(0.*500, 0.8*500)]),
                                            "_tau_1_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.6*1000)]),
                                            "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.54*1000)]),
                                            "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.5*1000)]),

                                            "_tau_1p5_1": np.array([(i/500) for i in np.arange(0.*500, 0.74*500)]),
                                            "_tau_1p5_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.58*1000)]),
                                            "_tau_1p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.46*1000)]),
                                            "_tau_1p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.42*1000)]),
                                            "_tau_1p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.30*1000)]),

                                            "_tau_2_1": np.array([(i/500) for i in np.arange(0.*500, 0.64*500)]),
                                            "_tau_2_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.46*1000)]),
                                            "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1000)]),
                                            "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.28*1000)]),
                                            "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.26*1000)]),

                                            "_tau21": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for one-pass kT minimization as per CMS
                                            "_tau32": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for one-pass kT minimization as per CMS
                                            "_tau21_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for WTA-kT for comparison
                                            "_tau32_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for WTA-kT for comparison
                                            "_tau21_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.4*200)]),#for excl.-kT and E-scheme as per basis
                                            "_tau32_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.6*200)]),#for excl.-kT and E-scheme as per basis
                                            "_mass": np.array([(i*2) for i in np.arange(0., 80.)]),
                                            #"_msoftdrop": np.array([(i*2) for i in np.arange(0., 80.)]),
                                            #"_pt": np.array([(i*20) for i in np.arange(9., 121.)]),
                                            }

        elif 'topSel' in selection:
            self.dict_variables_toUnfold = {

                                            "_tau_0p25_1": np.array([(i/200) for i in np.arange(0.2*200, 1.*200)]),
                                            "_tau_0p25_2": np.array([(i/500) for i in np.arange(0.*500, 1.*500)]),
                                            "_tau_0p25_3": np.array([(i/500) for i in np.arange(0.*500, 0.96*500)]),
                                            "_tau_0p25_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.94*1000)]),
                                            "_tau_0p25_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.9*1000)]),

                                            "_tau_0p5_1": np.array([(i/200) for i in np.arange(0.1*200, 0.96*200)]),
                                            "_tau_0p5_2": np.array([(i/500) for i in np.arange(0.*500, 0.9*500)]),
                                            "_tau_0p5_3": np.array([(i/500) for i in np.arange(0.*500, 0.8*500)]),
                                            "_tau_0p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.7*1000)]),
                                            "_tau_0p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.64*1000)]),

                                            "_tau_1_1": np.array([(i/200) for i in np.arange(0.*200, 0.84*200)]),
                                            "_tau_1_2": np.array([(i/500) for i in np.arange(0.*500, 0.72*500)]),
                                            "_tau_1_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.52*1000)]),
                                            "_tau_1_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.42*1000)]),
                                            "_tau_1_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.38*1000)]),

                                            "_tau_1p5_1": np.array([(i/500) for i in np.arange(0.*500, 0.7*500)]),
                                            "_tau_1p5_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.52*1000)]),
                                            "_tau_1p5_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.42*1000)]),
                                            "_tau_1p5_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1000)]),
                                            "_tau_1p5_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.28*1000)]),

                                            "_tau_2_1": np.array([(i/500) for i in np.arange(0.*500, 0.62*500)]),
                                            "_tau_2_2": np.array([(i/1000) for i in np.arange(0.*1000, 0.42*1000)]),
                                            "_tau_2_3": np.array([(i/1000) for i in np.arange(0.*1000, 0.3*1000)]),
                                            "_tau_2_4": np.array([(i/1000) for i in np.arange(0.*1000, 0.26*1000)]),
                                            "_tau_2_5": np.array([(i/1000) for i in np.arange(0.*1000, 0.24*1000)]),

                                            "_tau21": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for one-pass kT minimization as per CMS
                                            "_tau32": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for one-pass kT minimization as per CMS
                                            "_tau21_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for WTA-kT for comparison
                                            "_tau32_WTA": np.array([(i/200) for i in np.arange(0.*200, 1.2*200)]),#for WTA-kT for comparison
                                            "_tau21_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.4*200)]),#for excl.-kT and E-scheme as per basis
                                            "_tau32_exkT": np.array([(i/200) for i in np.arange(0.*200, 1.6*200)]),#for excl.-kT and E-scheme as per basis
                                            "_mass": np.array([(i*2) for i in np.arange(0., 180.)]),
                                            #"_msoftdrop": np.array([(i*2) for i in np.arange(0., 160.)]),
                                            #"_pt": np.array([(i*20) for i in np.arange(19., 141.)]),
                                            }

        
        ### Uncertainties
        
        self.sysSources = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSources if not( isys.endswith('nom'))]# or self.sysUnc) ]
        self.sysWeightList = ( '_pu', '_isr', '_fsr', '_pdf', '_l1prefiring', '_lepton', '_btag' ) #'_ps',
        self.constJESList = ( '_constituentJES_neutral', '_constituentJES_charged', '_constituentJES_photon' ) #'_ps',
        self.wtSources=['_puWeight','_isrWeight','_fsrWeight','_pdfWeight', '_l1prefiringWeight', 
                        '_leptonWeightISOStat', '_leptonWeightIDStat', '_leptonWeightTrigStat', '_leptonWeightRecoEffStat', 
                        '_leptonWeightISOSyst', '_leptonWeightIDSyst', '_leptonWeightTrigSyst', '_leptonWeightRecoEffSyst', 
                        '_btagWeightEfficiency', '_btagWeightCorrelated', '_btagWeightUncorrelated',] if self.wtUnc else [] #, '_leptonWeightAll', 
        if self.onlyUnc: self.sysSources = ['_nom'] + [ onlyUnc+i for i in [ 'Up', 'Down' ] ] #################### Not using for right now
        
        if self.wtUnc: #for the jes/jer cases sys self.sysSources is udpated  above
            self.sysSources = ['_nom'] + [ iwt+i for i in [ 'Up', 'Down' ] for iwt in self.wtSources if not iwt.endswith(('nom','pdfWeightAll')) ]
            
            if 'pdfWeightAll' in self.wtSources: self.sysSources = self.sysSources+['pdfWeightAll'] 
        
        #puWeights used only to change reco event weight (up/down) without application to gen weight, others applied to modify the genweight [and thereby also the recoWeight(=self.genWeight*self.puWeights)]
        
        
        self.selList = [ selection ] 
        self.sel = selection 
                    
        
        self._branchesToRead = self.getBranchesToRead(
                                                     dirname=self.inputDir, year=self.year,
                                                     kinematic_labels=self.kinematic_labels,
                                                     reco_only_labels=self.reco_only_labels,
                                                     nSub_labels=self.dict_variables_toUnfold,
                                                    )
        if self.isMC and self.withBTaggingWeights: 
            self._branchesToRead+=['selRecoAK4bjetLeptHem_nom_pt', 'selRecoAK4bjetLeptHem_nom_eta', 'selRecoAK4bjetLeptHem_nom_jetId', 'selRecoAK4bjetLeptHem_nom_hadronFlavour'] 
        if self.isMC and self.withLeptonWeights:     
            self._branchesToRead+=[ 'selRecoAK4bjetLeptHem_nom_btagDeepFlavB', 'selRecoMu_nom_tkRelIso']
        
        
        self.genWeights = []
        #self.initialise_cut_values(
        #                             top_ptmin = 400.,top_mSDmin = 140.,top_mSDmax = 300.,top_mINVmin = 140.,top_mINVmax = 300.,
        #                             W_ptmin = 200.,W_mSDmin = 65.,W_mSDmax = 120.,W_mINVmin = 65.,W_mINVmax = 120.,
        #                             LeptW_ptmin = 120.,Mu_ptmin = 55.,MET_ptmin = 40.,
        #                             AK4_ptmin = 30.,AK4_etamax = 2.4 ,AK4_jetID = 2 ,dR_max = 1.6,dR_min = 0.4,         )
        
        #def initialise_cut_values(self, 
        #                      top_ptmin = 400.,top_mSDmin = 140.,top_mSDmax = 300.,top_mINVmin = 140.,top_mINVmax = 300.,
        #                      W_ptmin = 200.,W_mSDmin = 65.,W_mSDmax = 120.,W_mINVmin = 65.,W_mINVmax = 120.,
        #                      LeptW_ptmin = 120.,Mu_ptmin = 55.,MET_ptmin = 40.,
        #                      AK4_ptmin = 30.,AK4_etamax = 2.4 ,AK4_jetID = 2 ,dR_max = 1.6,dR_min = 0.4,         ):
        
        self.top_ptmin = top_ptmin
        self.top_mSDmin = top_mSDmin
        self.top_mSDmax = top_mSDmax
        self.top_mINVmin = top_mINVmin
        self.top_mINVmax = top_mINVmax

        self.W_ptmin = W_ptmin
        self.W_mSDmin = W_mSDmin
        self.W_mSDmax = W_mSDmax
        self.W_mINVmin = W_mINVmin
        self.W_mINVmax = W_mINVmax

        self.LeptW_ptmin = LeptW_ptmin
        self.Mu_ptmin = Mu_ptmin
        self.MET_ptmin = MET_ptmin
        self.AK4_ptmin = AK4_ptmin
        self.AK4_etamax = AK4_etamax
        self.AK4_jetID = AK4_jetID
        self.dR_max = dR_max
        self.dR_min = dR_min
        self.Mu_tkRelIso = Mu_tkRelIso
        
    def add_btagAndLeptonWeights(self):
        if self.withLeptonWeights:
             
            tstart=time.time()

            events_leptonWt_fields = self.include_leptonWeights()#events)

            if self.verbose: print(f"Adding in new branches {events_leptonWt_fields.keys()} to access event weights+variations from lepton SFs")

            for key, values in events_leptonWt_fields.items():
                self.events[key] = values
            elapsed = time.time() - tstart
            if self.verbose: print (f'Finished adding in new lepton weight branches. Time taken:{elapsed}') 

        else:
            #for testing purposes
            tstart=time.time()

            events_leptonWt_fields = self.include_dummy_leptonWeights()#simply returns 1'/0's for pass/fail in reco-sel events

            if self.verbose: print(f"Adding in dummy weights {events_leptonWt_fields.keys()} to access event weights+variations from lepton SFs")

            for key, values in events_leptonWt_fields.items():
                self.events[key] = values
            elapsed = time.time() - tstart
            if self.verbose: print (f'Finished adding in dummy lepton weight branches. Time taken:{elapsed}') 


        if self.withBTaggingWeights:

            tstart=time.time()

            events_btagWt_fields = self.include_btagWeights()
            if self.verbose: print(f"Adding in new branches {events_btagWt_fields.keys()} to access event weights from b-tagging")

            for key, values in events_btagWt_fields.items():
                self.events[key] = values
            elapsed = time.time() - tstart
            if self.verbose: print (f'Finished adding in new b-tagging reweighting branches. Time taken:{elapsed}') 
            
        else:
            tstart=time.time()

            if self.verbose: print(f"Adding in dummy weights for b-tagging")
            events_btagWt_fields = self.include_dummy_btagWeights()

            for key, values in events_btagWt_fields.items():
                self.events[key] = values
            elapsed = time.time() - tstart
            if self.verbose: print (f'Finished adding in dummy b-tagging weight branches. Time taken:{elapsed}') 
        
        if self.test: #for testing selections, without running the histo processor
            s='_nom'
            self.recoWeight = self.events['evtGenWeight_nom']*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}']))) if self.isMC else self.events[f'totalRecoWeight_nom']
            
    def process(self, events):
        '''fill Hist histograms to accumulate event stats, convert to root or whatever else after returned by processor'''
        
        
        if not(self.onlyParquet): output = self.buildDictOfHistograms()
        
        # Get branches to be read 
        branches = self._branchesToRead
        
        # Currently the branches to be read are being defined somewhat circularly, the branches that are to be read are extracted from the class-level variable and only those columns are read in from the nanoskims, but keeping this functionality in for easy generalisability if needed; ie, if one wants to make this class self sufficient and add  in what I do in my histo production noteeboks as a postprocessor in this class. 
        #events = events[branches]
        self.events=events                
        
        ############################### btagging and lepton weights, nom+variations ###############################
        
        if self.isMC and self.test: self.genWeight = self.events['evtGenWeight_nom']
        
        if self.isMC and 'passRecoSel_nom' in self.events.fields:# and not(self.onlyParquet): 

            self.add_btagAndLeptonWeights()
        
        
        if not(self.onlyParquet): 
            if self.verbose: print (f'Now, processing histograms; systematic sources being considered:{self.sysSources} for nevents={len(events)} (before masking)')
            if self.verbose: 
                print (output.keys())
        elif self.verbose and (self.onlyParquet):
            print (f'Now, producing .parquet files; for nevents={len(events)} (before masking)')
            
            
        for sys in self.sysSources:
            
            
            s='_nom' if (sys.startswith(self.sysWeightList) or 'const' in sys) else sys
            
            ############ building event masks ############
            
            if self.mode=='noDeltaR': 

                topRecoMask = ( (self.events[f'selRecoJets{s}_pt']>self.top_ptmin) & 
                                (self.events[f'selRecoJets{s}_{self.massToCutOn}']>self.top_mSDmin) & 
                                (self.events[f'selRecoJets{s}_{self.massToCutOn}']<self.top_mSDmax) & 
                                (self.events[f'selRecoLeptHemDeltaR_nom']<self.dR_max) & 
                                (self.events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                (self.events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & 
                                (self.events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & 
                                (self.events[f'passRecoSel{s}']==1) & 
                                (self.events[f'totalRecoWeight{s}']!=0.) & 
                                (self.events[f'recoSelectedEventNumber_nom']>=1) ) #& (self.events[f'nRecoBtags_nom']>1) & (self.events[f'nRecoHadBtags_nom']==1) #0.04+10/self.events[f'selRecoMu_nom_pt'])

                WRecoMask = ( (self.events[f'selRecoJets{s}_pt']>self.W_ptmin) & 
                              (self.events[f'selRecoJets{s}_{self.massToCutOn}']>self.W_mSDmin) & 
                              (self.events[f'selRecoJets{s}_{self.massToCutOn}']<self.W_mSDmax) & 
                              (self.events[f'selRecoLeptHemDeltaR_nom']<self.dR_max) & 
                              (self.events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & 
                              (self.events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & 
                              (self.events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & 
                              (self.events[f'passRecoSel{s}']==1) & 
                              (self.events[f'totalRecoWeight{s}']!=0.) & 
                              (self.events[f'recoSelectedEventNumber_nom']>=1) )#& (self.events[f'nRecoBtags_nom']>1) & (self.events[f'nRecoHadBtags_nom']==1) 

                if not(self.isMC): 
                    topRecoMask = (topRecoMask) & (self.events[f'passHLT_Mu50{s}']==1)
                    WRecoMask = (WRecoMask) & (self.events[f'passHLT_Mu50{s}']==1)

                if self.isMC:# or self.isSigMC:    

                    topGenMask = ( (self.events[f'selGenJets_nom_pt']>self.top_ptmin) & 
                                   (self.events[f'selGenJets_nom_{self.massToCutOn}']>self.top_mSDmin) & 
                                   (self.events[f'selGenJets_nom_{self.massToCutOn}']<self.top_mSDmax) & 
                                   (self.events[f'selGenLeptHemDeltaR_nom']<self.dR_max) & 
                                   (self.events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                   (self.events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & 
                                   (self.events[f'selGenMET_nom_pt']>self.MET_ptmin ) & 
                                   (self.events[f'passGenSel{s}']==1) & 
                                   (self.events[f'evtGenWeight_nom']!=0.) )#& selRecoLeptHemDeltaR  (self.events[f'nGenBtags_nom']>1) & (self.events[f'nGenHadBtags_nom']==1) 

                    WGenMask = ( (self.events[f'selGenJets_nom_pt']>self.W_ptmin) & 
                                 (self.events[f'selGenJets_nom_{self.massToCutOn}']>self.W_mSDmin) & 
                                 (self.events[f'selGenJets_nom_{self.massToCutOn}']<self.W_mSDmax) & 
                                 (self.events[f'selGenLeptHemDeltaR_nom']<self.dR_max) & 
                                 (self.events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                 (self.events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & 
                                 (self.events[f'selGenMET_nom_pt']>self.MET_ptmin ) & 
                                 (self.events[f'passGenSel{s}']==1) & 
                                 (self.events[f'evtGenWeight_nom']!=0.) )#& (self.events[f'nGenBtags_nom']>1) & (self.events[f'nGenHadBtags_nom']==1) 

            
                    
            elif self.mode=='massAndDeltaR':# and self.massToCutOn=='msoftdrop': #new default

                topRecoMask = ( (self.events[f'selRecoJets{s}_pt']>self.top_ptmin) & 
                                (self.events[f'selRecoJets{s}_{self.massToCutOn}']>=self.top_mSDmin) & 
                                (self.events[f'selRecoJets{s}_{self.massToCutOn}']<self.top_mSDmax) & 
                                (self.events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                (self.events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & 
                                (self.events[f'selRecoMu_nom_tkRelIso']<=self.Mu_tkRelIso ) & 
                                (self.events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & 
                                (self.events[f'passRecoSel{s}']==1) & 
                                #(self.events[f'nRecoAK4s{s}']>=3) & 
                                (self.events[f'totalRecoWeight{s}']!=0.) & 
                                (self.events[f'selRecoLeptHemDeltaR_nom']<=self.dR_max) &
                                (self.events[f'selRecoAK4bjetLeptHem_nom_jetId']>=self.AK4_jetID) & 
                                (self.events[f'selRecoLeptHemDeltaR_nom']>self.dR_min) ) # & (self.events[f'nRecoBtags_nom']>1)# & (self.events[f'nRecoHadBtags_nom']==1) 

                WRecoMask = ( (self.events[f'selRecoJets{s}_pt']>self.W_ptmin) & 
                              (self.events[f'selRecoJets{s}_{self.massToCutOn}']>=self.W_mSDmin) & 
                              (self.events[f'selRecoJets{s}_{self.massToCutOn}']<self.W_mSDmax) & 
                              (self.events[f'selRecoLeptW_nom_pt']>self.LeptW_ptmin ) & 
                              (self.events[f'selRecoMu_nom_pt']>self.Mu_ptmin ) & 
                              (self.events[f'selRecoMu_nom_tkRelIso']<=self.Mu_tkRelIso ) & 
                              (self.events[f'selRecoMET_nom_pt']>self.MET_ptmin ) & 
                              (self.events[f'passRecoSel{s}']==1) & 
                              #(self.events[f'nRecoAK4s{s}']>=3) & 
                              (self.events[f'totalRecoWeight{s}']!=0.) & 
                              (self.events[f'selRecoLeptHemDeltaR_nom']<=self.dR_max) &
                              (self.events[f'selRecoAK4bjetLeptHem_nom_jetId']>=self.AK4_jetID) & 
                              (self.events[f'selRecoLeptHemDeltaR_nom']>self.dR_min) )# & (self.events[f'nRecoBtags_nom']>1)# & (self.events[f'nRecoHadBtags_nom']==1) 

                if self.isMC:# or self.isSigMC:    

                    topGenMask = ( (self.events[f'selGenJets_nom_pt']>self.top_ptmin) & 
                                   (self.events[f'selGenJets_nom_{self.massToCutOn}']>=self.top_mSDmin) & 
                                   (self.events[f'selGenJets_nom_{self.massToCutOn}']<self.top_mSDmax) & 
                                   (self.events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                   (self.events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & 
                                   (self.events[f'selGenMET_nom_pt']>self.MET_ptmin ) & 
                                   (self.events[f'passGenSel{s}']==1) & 
                                   #(self.events[f'nGenAK4s_nom']>=3) & 
                                   (self.events[f'evtGenWeight_nom']!=0.) & 
                                   (self.events[f'selGenLeptHemDeltaR_nom']<=self.dR_max) & 
                                   (self.events[f'selGenLeptHemDeltaR_nom']>self.dR_min)) #& (self.events[f'nGenBtags_nom']>1) & (self.events[f'nGenHadBtags_nom']==1) 

                    WGenMask = ( (self.events[f'selGenJets_nom_pt']>self.W_ptmin) & 
                                 (self.events[f'selGenJets_nom_{self.massToCutOn}']>=self.W_mSDmin) & 
                                 (self.events[f'selGenJets_nom_{self.massToCutOn}']<self.W_mSDmax) & 
                                 (self.events[f'selGenLeptW_nom_pt']>self.LeptW_ptmin ) & 
                                 (self.events[f'selGenMu_nom_pt']>self.Mu_ptmin ) & 
                                 (self.events[f'selGenMET_nom_pt']>self.MET_ptmin ) & 
                                 (self.events[f'passGenSel{s}']==1) & 
                                 #(self.events[f'nGenAK4s_nom']>=3) & 
                                 (self.events[f'evtGenWeight_nom']!=0.) & 
                                 (self.events[f'selGenLeptHemDeltaR_nom']<=self.dR_max) & 
                                 (self.events[f'selGenLeptHemDeltaR_nom']>self.dR_min)) #& (self.events[f'nGenBtags_nom']>1) & (self.events[f'nGenHadBtags_nom']==1) 
                               
                    

            else: 
                print(f"############################################ _WARNING_ ##################################################")
                print(f"provided histoproducer mode, #{self.mode}#, is undefined for at least the {self.massToCutOn} cut; exiting")
                print(f"#########################################################################################################")
                return -1

            if self.isMC:# or self.isSigMC:    
                if sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                    extraWeights = self.events[f'{sys.split("_")[1]}{s}']
                
                    self.genWeights = self.events['evtGenWeight_nom']*extraWeights
                    #genWeightsW = self.events['genWeight'][WRecoMask]*extraWeights[WRecoMask]
                else:
                    self.genWeights = self.events['evtGenWeight_nom']
                    #genWeightsW = self.events['genWeight'][WRecoMask]
           
            
            
            if self.verbose:
                print(s, sys, len(self.events[topRecoMask]), len(self.events[topGenMask]), len(self.events[WRecoMask]), len(self.events[WGenMask]), len(self.genWeights))

            if self.isSigMC:# and not(sys.startswith(self.sysWeightList)):
                
                s='_nom' if (sys.startswith(self.sysWeightList)) else sys

                if self.verbose: 
                    print('making individual W and topSel mask',sys,s)


                
                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) #& (~(topRecoMask|topGenMask))
                    WtopGenMask = (WGenMask) #& (~(topRecoMask|topGenMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) #& (~(WRecoMask|WGenMask))
                    WtopGenMask = (topGenMask) #& (~(WGenMask|WRecoMask))    
                    

                self.recoMask = WtopRecoMask
                self.genMask = WtopGenMask
                        
                selRecoMask = (WtopRecoMask) 
                selGenMask = (WtopGenMask) 

                trueRecoMask = (self.events[f'trueRecoJets{s}_pt']>0.) & ((selRecoMask) & (selGenMask)) 
                accepGenMask = (self.events[f'accepGenJets{s}_pt']>0.) & ((selGenMask) & (selRecoMask))
                
                self.accepGenMask = accepGenMask
                self.trueRecoMask = trueRecoMask

                fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) #((selRecoMask) ^ (trueRecoMask))#

                missGenMask =  ((selGenMask) & (~accepGenMask))#((selGenMask) ^ (accepGenMask))  #
                
                
                
                if self.verbose:
                    print(s, sys, len(self.events[selRecoMask]), len(self.events[trueRecoMask]), len(self.events[fakeRecoMask]), len(self.events[selGenMask]), len(self.events[accepGenMask]), len(self.events[missGenMask]),len(self.genWeights))
                
                if self.saveParquet:
                    
    
                    if (sys.endswith('nom')): 
                        
                        
                        sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                        
                        print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_{self.splitCount}.parquet'}")

                        ak.to_parquet(self.events[selRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_reco_{self.splitCount}.parquet')
                        ak.to_parquet(self.events[selGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_OC_gen_{self.splitCount}.parquet')
                        
                        return 1

                    
            elif self.isMC and not(self.isSigMC) and sys.endswith('nom'): 
                

                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) #& (~(topRecoMask|topGenMask))
                    WtopGenMask = (WGenMask) #& (~(topRecoMask|topGenMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) #& (~(WRecoMask|WGenMask))
                    WtopGenMask = (topGenMask) #& (~(WGenMask|WRecoMask))    

                
                self.recoMask = WtopRecoMask
                self.genMask = WtopGenMask
                
                selRecoMask = (WtopRecoMask) 
                selGenMask = (WtopGenMask) 

                trueRecoMask = (self.events[f'trueRecoJets{s}_pt']>0.) & ((selRecoMask) & (selGenMask)) 
                accepGenMask = (self.events[f'accepGenJets{s}_pt']>0.) & ((selGenMask) & (selRecoMask))
                
                self.accepGenMask = accepGenMask
                self.trueRecoMask = trueRecoMask
                
                fakeRecoMask = ((selRecoMask) & (~trueRecoMask)) #((selRecoMask) ^ (trueRecoMask))#

                missGenMask =  ((selGenMask) & (~accepGenMask))#((selGenMask) ^ (accepGenMask))  #
                
                #selRecoMask = (trueRecoMask) | (fakeRecoMask) #& (event_mask)
                #selGenMask = (accepGenMask) | (missGenMask) #& (event_mask)
                
                
                sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                
               

                if self.saveParquet:
                    
                        
                    print(f"Saving the following .parquet file with file-stem: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_reco/gen_{self.splitCount}.parquet'}")

                    ak.to_parquet(self.events[selRecoMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/recoMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_reco_{self.splitCount}.parquet')
                    ak.to_parquet(self.events[selGenMask],self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/genMasked/'+self.year+'/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+f'_nomWts{sel}_gen_{self.splitCount}.parquet')
                    return 1
            
            elif not(self.isMC) and sys.endswith('nom'): 
                
                if '_WSel' in self.selList:
                    WtopRecoMask = (WRecoMask) & (~(topRecoMask))

                elif '_topSel' in self.selList: 
                    WtopRecoMask = (topRecoMask) & (~(WRecoMask))

                self.recoMask = WtopRecoMask
                
                selRecoMask = (WtopRecoMask) 
                
                sel = '_WSel' if '_WSel' in self.selList else '_topSel'
                
                if self.saveParquet:
                    self.events['selRecoMask']=selRecoMask
                    print(f"Saving the following .parquet file: {self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'{sel}_{self.splitCount}.parquet'}")

                    ak.to_parquet(events,self.inputDir[0].split('jetObservables/')[0]+'jetObservables/Wtop_rootToParquet/'+self.inputDir[0].split('kadatta/jetObservables/')[1].split('/')[0]+'_UL'+self.year+self.era+f'{sel}_{self.splitCount}.parquet')
                    return 1
                    
            else: 
                print (f'something fishy in what type of sample you want me to load, recheck input config', s,sys,self.isMC,self.isSigMC) 
            
            if self.verbose and self.isMC and not(self.onlyParquet):
                
                print(output.keys())
                print("Going to fill histos")
            
            if not(self.onlyParquet): 
                for isel in self.selList:

                    if self.verbose: print(sys,s)

                    if sys.endswith('nom'):
                        listOfMiscOutputHistosNames = [k for k,h in output.items() if ((not ('Jet' in k and ('nsel' in k or 'ntrue' in k or 'naccep' in k)) and ('nom' in k)) and ('Mu' in k or 'LeptW' in k or 'LeptTop' in k or 'MET' in k or 'AK4bjetHadHem' in k or 'nPV' in k or 'DeltaR' or 'DeltaPhi' in k or 'Btag' in k))]
                        #btagWeightNom{s}
                        totalRecoWeight = self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*self.genWeights*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}']))) if self.isMC else self.events[f'totalRecoWeight_nom']#self.events[f'evtGenWeight_nom']*#
                        
                        self.recoWeight = totalRecoWeight
                        totalGenWeight = self.genWeights if self.isMC else None

                        if self.verbose: print(listOfMiscOutputHistosNames)
                        
                        for histos in listOfMiscOutputHistosNames:
                            key=histos

                            if 'Mu' in key: 
                                dicto=self.dict_variables_kinematics_Muon
                                for i in dicto:
                                    if 'gen' in key.lower() and i.endswith(('_p','_ptRel','_tkRelIso')): continue

                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(self.events[f'selRecoMu{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(self.events[f'selGenMu{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'MET' in key: 
                                dicto=self.dict_variables_kinematics_MET
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(self.events[f'selRecoMET{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(self.events[f'selGenMET{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'LeptW' in key: 
                                dicto=self.dict_variables_kinematics_LeptW
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(self.events[f'selRecoLeptW{s}{i}'][selRecoMask], weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(self.events[f'selGenLeptW{s}{i}'][selGenMask], weight=totalGenWeight[selGenMask], threads=8)
                            elif 'LeptTop' in key: 
                                dicto=self.dict_variables_kinematics_LeptTop
                                for i in dicto:
                                    if i in key:
                                        if 'reco' in key.lower(): 
                                            output[key].fill(self.events[f'selRecoLeptTop{s}{i}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower():
                                            output[key].fill(self.events[f'selGenLeptTop{s}{i}'][selGenMask],weight=totalGenWeight[selGenMask], threads=8)

                            elif 'AK4bjetLeptHem' in key and not ('DeltaR' in key or 'DeltaPhi' in key): 
                                dicto=self.dict_variables_kinematics_AK4Lept
                                for i in dicto:
                                    if 'gen' in key.lower() and i.endswith(('_jetId','_hadronFlavour')): continue
                                    if 'reco' in key.lower() and not(self.isMC) and i.endswith('_hadronFlavour'): continue
                                
                                    if i in key:
                                        
                                        if 'reco' in key.lower(): 
                                            output[key].fill(self.events[f'selRecoAK4bjetLeptHem{s}{i}'][selRecoMask], weight=totalRecoWeight[selRecoMask], threads=8)
                                        elif 'gen' in key.lower() and self.isMC:
                                            output[key].fill(self.events[f'selGenAK4bjetLeptHem{s}{i}'][selGenMask], weight=totalGenWeight[selGenMask], threads=8)

                            elif 'DeltaRap' in key:#
                                if 'reco' in key.lower():
                                    output[key].fill(self.events[f'selRecoLeptHemDeltaRap{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selGenLeptHemDeltaRap{s}'][selGenMask], weight=totalGenWeight[selGenMask])

                            elif 'DeltaR' in key and not('DeltaRap' in key) and not ('AK8' in key):
                                if 'reco' in key.lower():
                                    output[key].fill(self.events[f'selRecoLeptHemDeltaR{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selGenLeptHemDeltaR{s}'][selGenMask], weight=totalGenWeight[selGenMask])

                            elif 'DeltaPhi' in key and not ('AK8' in key):
                                if 'reco' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selRecoLeptHemDeltaPhi{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selGenLeptHemDeltaPhi{s}'][selGenMask], weight=totalGenWeight[selGenMask])
                            
                            elif 'DeltaPhi' in key and ('AK8' in key):
                                print(key)
                                if 'reco' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selRecoLeptHemAK8DeltaPhi{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selGenLeptHemAK8DeltaPhi{s}'][selGenMask], weight=totalGenWeight[selGenMask])
                            
                            elif 'DeltaR' in key and not('DeltaRap' in key) and ('AK8' in key):
                                print(key)
                                if 'reco' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selRecoLeptHemAK8DeltaR{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )
                                elif 'gen' in key.lower() and self.isMC:
                                    output[key].fill(self.events[f'selGenLeptHemAK8DeltaR{s}'][selGenMask], weight=totalGenWeight[selGenMask])


                            elif 'good_nPVs' in key:
                                output[key].fill(self.events[f'good_nPVs{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nRecoBtags' in key:
                                output[key].fill(self.events[f'nRecoBtags{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenBtags' in key and self.isMC:
                                output[key].fill(self.events[f'nGenBtags{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                            elif 'nselRecoJets' in key:
                                output[key].fill(self.events[f'nselRecoJets{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nselGenJets' in key and self.isMC:
                                output[key].fill(self.events[f'nselGenJets{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )
                                
                            elif 'nRecoAK4s' in key:
                                output[key].fill(self.events[f'nRecoAK4s{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenAK4s' in key and self.isMC:
                                output[key].fill(self.events[f'nGenAK4s{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                            elif 'nRecoLeptBtags' in key:
                                output[key].fill(self.events[f'nRecoLeptBtags{s}'][selRecoMask], weight=totalRecoWeight[selRecoMask] )

                            elif 'nGenLeptBtags' in key and self.isMC:
                                output[key].fill(self.events[f'nGenLeptBtags{s}'][selGenMask], weight=totalRecoWeight[selGenMask] )

                    listOfJetOutputHistosNames = [k for k,h in output.items() if ((('Jet' in k and not('nsel' in k or 'ntrue' in k or 'naccep' in k)) and (sys in k)) or ('residual' in k and not(sys in k)) or ('resol' in k and not(sys in k)))] 
                    if self.verbose: print (listOfJetOutputHistosNames[0:10],sys,s)

                    for k in listOfJetOutputHistosNames:
                        key=k
                        #print(key)
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
                        if not (sys.startswith(self.sysWeightList) or sys.startswith(self.constJESList)):
                            s=sys
                            totalRecoWeight = self.genWeights*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}']))) if self.isMC else self.events[f'totalRecoWeight_nom']
                            totalGenWeight = self.events[f'evtGenWeight_nom'] if self.isMC else None
                        else:
                            if self.isSigMC:
                                s='_nom'
                                #print(f'{sys.split("_")[1]}{s}')
                                if 'pu' in sys:
                                    totalGenWeight = self.events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*self.events[f'{sys.split("_")[1]}{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}'])))

                                elif 'l1' in sys:
                                    totalGenWeight = self.events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*self.events[f'{sys.split("_")[1]}{s}']*self.events[f'puWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}'])))

                                elif '_btag' in sys:
                                    totalGenWeight = self.events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*self.events[f'new{self.sel}_{sys.split("_")[1]}{s}']*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))

                                elif 'lepton' in sys:
                                    #print(sys,f'new{self.sel}_{sys.split("_")[1]}{s}')
                                    totalGenWeight = self.events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*self.events[f'new{self.sel}_{sys.split("_")[1]}{s}']*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}'])))

                                elif sys.startswith(('_isr','_fsr','_pdf')):#,'_pdf''):
                                    totalGenWeight = self.genWeights#self.events[f'evtGenWeight{s}']*self.events[f'{sys.split("_")[1]}{s}']
                                    totalRecoWeight = self.genWeights*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_leptonWeightNom{s}'] if self.withLeptonWeights else np.ones(len(self.events[f'puWeightNom{s}'])))*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}'])))
                                
                                elif 'const' in sys:
                                    totalGenWeight = self.events[f'evtGenWeight{s}'] 
                                    totalRecoWeight = self.genWeights*self.events[f'new{self.sel}_leptonWeightNom{s}']*self.events[f'puWeightNom{s}']*self.events[f'l1prefiringWeightNom{s}']*(self.events[f'new{self.sel}_btagWeightNom{s}'] if self.withBTaggingWeights else np.ones(len(self.events[f'puWeightNom{s}'])))
                                
                        if self.isSigMC and 'const' in sys:
                            s=sys
                        
                        if (key.lower().startswith('reco')):
                            if self.verbose: 
                                print('all reco', key,s,sys, len(self.events[f'selRecoJets{s}{varToFill}'][selRecoMask]))
                            output[key].fill(self.events[f'selRecoJets{s}{varToFill}'][selRecoMask],weight=totalRecoWeight[selRecoMask], threads=8)

                        # Stuff below should have diff wts already available to them so no need to touch anything down here
                        if self.isMC:# and (not varToFill.startswith(tuple(self.reco_only_labels))):

                            if (key.lower().startswith('true')) and self.isSigMC :
                                if self.verbose: print('true', key,s,sys, len(self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask]))
                                output[key].fill(self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                            elif (key.lower().startswith('fake')) and self.isSigMC:
                                if self.verbose: print('fake', key,s,sys,len(self.events[f'selRecoJets{s}{varToFill}'][fakeRecoMask]))
                                output[key].fill(self.events[f'selRecoJets{s}{varToFill}'][fakeRecoMask],weight=totalRecoWeight[fakeRecoMask],
                                                 threads=8)

                            if (key.lower().startswith('genjet')):
                                if self.verbose: 
                                    print('all gen', key,s,sys, len(self.events[f'selGenJets_nom{varToFill}'][selGenMask]))

                                output[key].fill(self.events[f'selGenJets_nom{varToFill}'][selGenMask],weight=totalGenWeight[selGenMask],
                                                 threads=8)

                            elif (key.lower().startswith('accepgen')):
                                if self.verbose: 
                                    print('accepgen', key,s,sys, len(self.events[f'accepGenJets{s}{varToFill}'][accepGenMask]))
                                output[key].fill(self.events[f'accepGenJets{s}{varToFill}'][accepGenMask],weight=totalGenWeight[accepGenMask],
                                                 threads=8)

                            elif (key.lower().startswith('missgenjet')):
                                if self.verbose: 
                                    print('missgen', key,s,sys, len(self.events[f'selGenJets_nom{varToFill}'][missGenMask]))
                                output[key].fill(self.events[f'selGenJets_nom{varToFill}'][missGenMask],weight=totalGenWeight[missGenMask],
                                                 threads=8)

                            elif ( 'resp' in key.lower() and not('withmiss' in key.lower())) and self.isSigMC:
                                #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                                output[key].fill(gen=self.events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask],weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                                output[key].fill(gen=self.events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                                 weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                                 threads=8)     

                            elif ('respwithmiss' in key.lower()) and self.isSigMC:
                                if self.verbose: print("filling resp w. misses")
                                #fill matched entries with weight wgen*wxreco='totalrecoweight' (x=>excl.to reco)
                                output[key].fill(gen=self.events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask], weight=totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill counter weight wgen(1-wxrec)= wgen-'totalrecoWeight'
                                output[key].fill(gen=self.events[f'accepGenJets{s}{varToFill}'][accepGenMask], reco=-1.*np.ones(len(self.events[f'trueRecoJets{s}{varToFill}'][trueRecoMask])),
                                                 weight=totalGenWeight[accepGenMask]-totalRecoWeight[trueRecoMask],
                                                 threads=8)

                                #fill missgen weight
                                output[key].fill(gen=self.events[f'selGenJets_nom{varToFill}'][missGenMask], reco=-1.*np.ones(len(self.events[f'selGenJets_nom{varToFill}'][missGenMask])),
                                                 weight=totalGenWeight[missGenMask],
                                                 threads=8)


                            elif ('residual' in key.lower() or 'resol' in key.lower()) and sys.endswith('nom') and self.isSigMC:
                                zeroMask=(self.events[f'accepGenJets{s}{varToFill}']!=0.)&(accepGenMask)
                                
                                response = self.events[f'trueRecoJets{s}{varToFill}'][zeroMask]/self.events[f'accepGenJets{s}{varToFill}'][zeroMask]
                                response = np.nan_to_num(response,nan=-999.)
                                residual = self.events[f'trueRecoJets{s}{varToFill}'][zeroMask]-self.events[f'accepGenJets{s}{varToFill}'][zeroMask]
                                #residual = np.nan_to_num(residual, nan=-929.)
                                
                                if 'noWt_' in key: output[key].fill(response)#, weight=totalRecoWeight[zeroMask])
                                elif 'residual' in key: output[key].fill(residual)
                            

                if self.verbose: print("Histos filled!",self.selList,sys)
                l=[]
                for x,y in output.items(): #y.SetDirectory(0)

                    if self.sysUnc and self.onlyUnc and x.startswith(('accepgen','miss')) and '_nom' in x:
                        l.append(x)
                    #elif self.isMC and not(self.isSigMC) and x.startswith(('accepgen','miss','true','fake','residual' )) and '_nom' in x:
                    #    l.append(x)
                for k in l:
                    del(output[k])
        print("Histos filled!", self.selList, self.year)

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
                        dictOfHists['residual'+iJ+'_pt'+isel] = (hist.Hist.new.Regular(1000, -2, 2, name='residual'+iJ+'_pt'+isel, label='AK8 reco - gen jet pt').Weight())
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
                                if 'gen' in itype and x.endswith(('_p','_ptRel','_tkRelIso')): continue
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} muon {x}').Weight())        
                                
                        elif 'met' in iO.lower():
                            for x, y in self.dict_variables_kinematics_MET.items():
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} MET {x}').Weight())        
                                
                        elif 'leptw' in iO.lower():
                            for x, y in self.dict_variables_kinematics_LeptW.items():
                                #print(itype+iO+x+sysUnc+isel)
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} leptonic W {x}').Weight())    
                                
                        elif 'lepttop' in iO.lower():
                            for x, y in self.dict_variables_kinematics_LeptTop.items():
                                #print(itype+iO+x+sysUnc+isel)
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} leptonic top {x}').Weight())    
                                
                        elif 'ak4' in iO.lower():
                            for x, y in self.dict_variables_kinematics_AK4Lept.items():
                                if 'gen' in itype and x.endswith(('_jetId','_hadronFlavour')): continue
                                if 'reco' in itype and not(self.isMC) and x.endswith('_hadronFlavour'): continue
                                dictOfHists[itype+iO+x+sysUnc+isel] = (hist.Hist.new.Variable(y,name=itype+iO+x+sysUnc+isel, label=f'{itype} AK4 {x}').Weight())        

                #if self.verbose: print (self.dict_variables_reco.keys())

                for itype in self.dict_variables_reco:
                    bins=self.dict_variables_reco[itype]
                    
                    if 'good_nPVs' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'Good nPVs').Weight())  
                    elif 'nRecoBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco').Weight())  
                    elif 'nRecoAK4s' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nAK4 (reco)').Weight())  
                    elif 'nselRecoJets' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nAK8recojets').Weight())  
                    elif 'nRecoLeptBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBtags_reco (lept. hem.)').Weight())  
                    elif 'selRecoLeptHemDeltaRap' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dy(#mu, b-cand)').Weight())  
                    elif 'selRecoLeptHemDeltaR' in itype and not('DeltaRap' in itype): dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(#mu, b-cand)').Weight())  
                    elif 'selRecoLeptHemDeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(#mu, b-cand)').Weight())  
                    elif 'selRecoLeptHemAK8DeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(AK8, lept. b-cand)').Weight())  
                    elif 'selRecoLeptHemAK8DeltaR' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, lept. b-cand)').Weight())  

                if self.isMC or self.isSigMC:
                    #if self.verbose: print (self.dict_variables_gen.keys())

                    for itype in self.dict_variables_gen:
                        bins=self.dict_variables_gen[itype]
                        
                        if 'nGenBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen').Weight())  
                        elif 'nselGenJets' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nAK8genjets').Weight())  
                        elif 'nGenLeptBtags' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nBmatched_Gen (lept. hem.)').Weight())  
                        elif 'nGenAK4s' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'nAK4 (gen)').Weight())  
                        elif 'selGenLeptHemDeltaRap' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dy(#mu, b-cand)').Weight())  
                        elif 'selGenLeptHemDeltaR' in itype and not('DeltaRap' in itype): dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(#mu, b-cand)').Weight())  
                        elif 'selGenLeptHemDeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(#mu, b-cand)').Weight())  
                        elif 'selGenLeptHemAK8DeltaPhi' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dPhi(AK8, lept. b-cand)').Weight())  
                        elif 'selGenLeptHemAK8DeltaR' in itype: dictOfHists[itype+sysUnc+isel] = (hist.Hist.new.Variable(bins,name=itype+sysUnc+isel, label=f'dR(AK8, lept. b-cand)').Weight())  
        
        #if self.verbose: print ("Dictionary of hists to be filled:",dictOfHists)
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
                    
        if self.verbose: print (self.sysSources)

        for sys in self.sysSources:
            if ('btag' in sys) or ('lepton' in sys): 
                continue #this is since the btag/leptonWeight branches are now generated on the fly vs. being read in
            #if not wtUnc and not sysUnc:
            
            
            for i in kinematic_labels+list(nSub_labels.keys()): 
                
                if (sys.endswith('nom') or (self.isSigMC and self.sysUnc) ) and (i in self.dict_variables_kinematics_AK8 or i in nSub_labels.keys()):#not('mt' in i or '_p' in i or 'jetId' in i or 'jetId' in i): 
                    reco_list.append('selRecoJets'+sys+i)
                
                if not ('softdrop' in i) and not('tau' in i):
                    if sys.endswith('nom') or (self.isSigMC and self.sysUnc and not('const' in sys)):
                        #if not('mt' in i): 
                        
                        if i in self.dict_variables_kinematics_AK4Lept:
                            if not(self.isMC) and i.endswith('_hadronFlavour'): continue
                            reco_list.append('selRecoAK4bjetLeptHem'+sys+i)

                        if i in self.dict_variables_kinematics_Muon and sys.endswith('nom'):
                            reco_list.append('selRecoMu'+sys+i)

                        if i in self.dict_variables_kinematics_MET and sys.endswith('nom'):
                            reco_list.append('selRecoMET'+sys+i)

                        if i in self.dict_variables_kinematics_LeptW and sys.endswith('nom'):
                            reco_list.append('selRecoLeptW'+sys+i)

                        if i in self.dict_variables_kinematics_LeptTop and sys.endswith('nom'):
                            reco_list.append('selRecoLeptTop'+sys+i)
                        

                if (self.isMC and (sys.endswith('nom') or self.sysUnc)) and (i in self.dict_variables_kinematics_AK8 or i in nSub_labels.keys()):#not('mt' in i): 
                    reco_list.append('trueRecoJets'+sys+i)
                    if not i in reco_only_labels: gen_list.append('accepGenJets'+sys+i)
                        
                if self.isMC and sys.endswith('nom') : 
                    if not i in reco_only_labels and (i in self.dict_variables_kinematics_AK8 or i in nSub_labels.keys()): 
                        #if not('mt' in i): 
                        gen_list.append('selGenJets'+sys+i)
                            
                    if not 'softdrop' in i and not('tau' in i):
                        if sys.endswith('nom') or (self.isSigMC and self.sysUnc):
                            if i in self.dict_variables_kinematics_Muon and sys.endswith('nom'):
                                if (i.endswith(('_p','_ptRel','_tkRelIso'))): continue
                                gen_list.append('selGenMu'+sys+i)

                            if i in self.dict_variables_kinematics_MET and sys.endswith('nom'):
                                gen_list.append('selGenMET'+sys+i)

                            if i in self.dict_variables_kinematics_AK4Lept:
                                if i.endswith(('_jetId','_hadronFlavour')): continue
                                gen_list.append('selGenAK4bjetLeptHem'+sys+i)

                            if i in self.dict_variables_kinematics_LeptW and sys.endswith('nom'):
                                gen_list.append('selGenLeptW'+sys+i)

                            if i in self.dict_variables_kinematics_LeptTop and sys.endswith('nom'):
                                gen_list.append('selGenLeptTop'+sys+i)


                

            if self.isMC and (sys.endswith('nom') or (self.isSigMC and self.sysUnc and not('const' in sys))):
                
                
                reco_list.append('nselRecoJets'+sys)
                reco_list.append('ntrueRecoJets'+sys)
                reco_list.append('naccepGenJets'+sys)
                

                if sys.endswith('nom'): 
                    #print(sys,'nselGenJets'+sys)
                    gen_list.append('nselGenJets'+sys)
                    gen_list.append("nGenAK4s"+sys)
                    
                reco_list.append("puWeightNom"+sys)
                reco_list.append("leptonWeightNom"+sys)
                reco_list.append("l1prefiringWeightNom"+sys)
                reco_list.append("pdfWeightNom"+sys)

                reco_list.append("totalRecoWeight"+sys)
                reco_list.append("passRecoSel"+sys)
                reco_list.append("nRecoAK4s"+sys)               
                
                reco_list.append(f"selRecoJets{sys}_msoftdrop_corr_PUPPI")#+sys)
                    
                
                gen_list.append("evtGenWeight"+sys)
                gen_list.append("passGenSel"+sys)
                
                for i in self.dict_variables_reco:
                    reco_list.append(i+(sys if not ('LeptHemDelta' in i) else '_nom'))
                for i in self.dict_variables_gen:
                    gen_list.append(i+(sys if not('nGen' in i or 'nselGen' in i) and not('LeptHem' in i) else '_nom') )
                    
            elif (not self.isMC) and sys.endswith('nom'):

                reco_list.append('nselRecoJets'+sys)

                reco_list.append("totalRecoWeight"+sys)               
                reco_list.append("passRecoSel"+sys)               
                reco_list.append("nRecoAK4s"+sys)               
                
                reco_list.append("passHLT_Mu50"+sys)               
                reco_list.append(f"selRecoJets{sys}_msoftdrop_corr_PUPPI")#+sys)

                for i in self.dict_variables_reco:
                    reco_list.append(i+sys)
                    
            elif (self.isMC and self.isSigMC) and self.wtUnc and not(sys.endswith('_nom') or  self.sysUnc):

                if 'pu' in sys or 'l1' in sys or 'lepton' in sys or 'btag' in sys: 
                    reco_reweights_list.append(sys.split('_')[1]+'_nom')
                else: 
                    gen_reweights_list.append(sys.split('_')[1]+'_nom')

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
    
    
    def include_dummy_leptonWeights(self):
        
        systematic_suffix_map = {
                                 'nominal': 'Nom',
                                }
        
        new_fields = OrderedDict()
        #new_fields[f'new{self.sel}_leptonWeightNom_nom'] = []
        
        if not(self.sysUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        if self.verbose: print(f"Initialising dummy lepton weights (0=fail,1=pass)")
        
    
        for sys in sysList:
            passMask = (self.events[f'passRecoSel{sys}']==1)
            failMask = (self.events[f'passRecoSel{sys}']!=1)
            
            dummyOnes = np.ones(len(self.events))
            ID_SFs, ISO_SFs, RecoEff_SFs, Trig_SFs = dummyOnes,dummyOnes,dummyOnes,dummyOnes
            
            ID_SFs[failMask], ISO_SFs[failMask], RecoEff_SFs[failMask], Trig_SFs[failMask] = np.zeros(len(ID_SFs[failMask])), np.zeros(len(ISO_SFs[failMask])), np.zeros(len(RecoEff_SFs[failMask])), np.zeros(len(Trig_SFs[failMask]))
            
            new_fields[f'new{self.sel}_leptonWeightNom{sys}'] = ID_SFs*ISO_SFs*RecoEff_SFs*Trig_SFs
         
        # Convert lists to arrays
        for key in new_fields:
            #print(key)
            if not(type(new_fields[key])==np.ndarray): new_fields[key] = np.array(new_fields[key])

        return new_fields
    
    
    def include_dummy_btagWeights(self):
        # include dummy b-tagging weights when running the histoProducer for tests/whatever other reason
        
        
        new_fields = OrderedDict()
       # new_fields[f'new{self.sel}_btagWeightNom_{sys}'] = np.ones(len(self.events))
        
        if not(self.sysUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        if self.verbose: print(f"Initialising dummy btagging weights (0=fail,1=pass)")
        
        for sys in sysList:
            passMask = (self.events[f'passRecoSel{sys}']==1)
            failMask = (self.events[f'passRecoSel{sys}']!=1)
            
            new_fields[f'new{self.sel}_btagWeightNom{sys}'] = np.ones(len(self.events))
            new_fields[f'new{self.sel}_btagWeightNom{sys}'][failMask] = np.zeros(len(new_fields[f'new{self.sel}_btagWeightNom{sys}'][failMask])) 
                            

        # Convert lists to arrays
        #for key in new_fields:
        #    if not(type(new_fields[key])==np.ma.core.MaskedArray): new_fields[key] =np.ma.core.MaskedArray(new_fields[key],dtype=np.float32)

        return new_fields
    
    def include_btagWeights(self):
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
                                     'down_uncorrelated': 'UncorrelatedDown',
                                     'up_efficiency': 'EfficiencyUp',
                                     'down_efficiency': 'EfficiencyDown',
                                
                                    }
        
        
        new_fields = OrderedDict()
        
        
        for s in systematic_suffix_map.keys(): 
            if not (self.sysUnc) or 'const' in self.onlyUnc: 
                new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[s]}_nom'] = []
            else:
                if s!='central': continue
                    
                for sys in self.sysSources:
                    new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[s]}{sys}'] = []
            
       
        if not(self.sysUnc) or ('const' in self.onlyUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        if self.verbose: print(f"Computing btagging weights ({systematic_suffix_map.keys()}) for following sys sources: {sysList}")
        
                        
        for sys in sysList:     
            
            

            jet_pt = self.events['selRecoAK4bjetLeptHem_nom_pt']
            jet_eta = self.events['selRecoAK4bjetLeptHem_nom_eta']
            jet_ID = self.events['selRecoAK4bjetLeptHem_nom_jetId']
            jet_flav = self.events['selRecoAK4bjetLeptHem_nom_hadronFlavour']
            jet_discr = self.events['selRecoAK4bjetLeptHem_nom_btagDeepFlavB']
            
            
            # Define the masks to be used on the selected lept. hem. AK4
            ID_mask = jet_ID>=self.AK4_jetID
            discr_mask = (jet_discr >= self.WPBDisc)
            ptMask = (jet_pt > self.AK4_ptmin)
            etaMask =((jet_eta<=self.AK4_etamax) & (jet_eta>=-1.*self.AK4_etamax))# (np.abs(jet_eta) <= self.AK4_etamax)
            jet_eta = np.abs(jet_eta)
            
            
            # Identify b, c, and light jets
            b_jetsMask = (jet_flav == 5)
            c_jetsMask = (jet_flav == 4)
            udsg_jetsMask = (jet_flav == 0)


            passMask = (self.events[f'passRecoSel{sys}']==1) & (ID_mask) & (discr_mask) & (ptMask) & (etaMask)
            failMask = (self.events[f'passRecoSel{sys}']!=1)  & ( (~ID_mask) | (~discr_mask) | (~ptMask) | (~etaMask))

            
            dummyOnes = np.ones(len(self.events))
            dummyZeros = np.zeros(len(self.events['selRecoAK4bjetLeptHem_nom_pt'][failMask]))

            
            # Compute wt=eff*SFs for each systematic
            for systematic in systematic_suffix_map.keys():#self.btaggingSystLabels:
                
                efficiency = np.ones_like(dummyOnes)
                SF = np.ones_like(dummyOnes)

                if not('efficiency' in systematic):


                    efficiency[(b_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_b_b_{self.year}{self.sel}'].evaluate(jet_pt[(b_jetsMask) & (passMask)],
                                                                                                                               np.abs(jet_eta[(b_jetsMask) & (passMask)]))

                    efficiency[(c_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_b_c_{self.year}{self.sel}'].evaluate(jet_pt[(c_jetsMask) & (passMask)],
                                                                                                                               np.abs(jet_eta[(c_jetsMask) & (passMask)]))

                    efficiency[(udsg_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_b_udsg_{self.year}{self.sel}'].evaluate(jet_pt[(udsg_jetsMask) & (passMask)],
                                                                                                                               np.abs(jet_eta[(udsg_jetsMask) & (passMask)]))
                    efficiency[failMask] = np.zeros_like(dummyZeros)

                    SF[(b_jetsMask) & (passMask) ] = self.btvjson['deepJet_comb'].evaluate(systematic, self.btagWP, 
                                                                                           jet_flav[(b_jetsMask) & (passMask)],
                                                                                           np.abs(jet_eta[(b_jetsMask) & (passMask)]),
                                                                                           jet_pt[(b_jetsMask) & (passMask)],
                                                                                           
                                                                                           )

                    SF[(c_jetsMask) & (passMask) ] = self.btvjson['deepJet_comb'].evaluate(systematic, self.btagWP, 
                                                                                           jet_flav[(c_jetsMask) & (passMask)],
                                                                                           np.abs(jet_eta[(c_jetsMask) & (passMask)]),
                                                                                           jet_pt[(c_jetsMask) & (passMask)],

                                                                                           )

                    SF[(udsg_jetsMask) & (passMask) ] = self.btvjson['deepJet_incl'].evaluate(systematic, self.btagWP, 
                                                                                              jet_flav[(udsg_jetsMask) & (passMask)],
                                                                                              np.abs(jet_eta[(udsg_jetsMask) & (passMask)]),
                                                                                              jet_pt[(udsg_jetsMask) & (passMask)],

                                                                                             )
                    SF[failMask] = np.zeros_like(dummyZeros)


                    new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'] = SF*efficiency

                else:
                    
                    SF[(b_jetsMask) & (passMask) ] = self.btvjson['deepJet_comb'].evaluate('central', self.btagWP, 
                                                                                           jet_flav[(b_jetsMask) & (passMask)],
                                                                                           np.abs(jet_eta[(b_jetsMask) & (passMask)]),
                                                                                           jet_pt[(b_jetsMask) & (passMask)],
                                                                                           
                                                                                           )

                    SF[(c_jetsMask) & (passMask) ] = self.btvjson['deepJet_comb'].evaluate('central', self.btagWP, 
                                                                                           jet_flav[(c_jetsMask) & (passMask)],
                                                                                           np.abs(jet_eta[(c_jetsMask) & (passMask)]),
                                                                                           jet_pt[(c_jetsMask) & (passMask)],

                                                                                           )

                    SF[(udsg_jetsMask) & (passMask) ] = self.btvjson['deepJet_incl'].evaluate('central', self.btagWP, 
                                                                                              jet_flav[(udsg_jetsMask) & (passMask)],
                                                                                              np.abs(jet_eta[(udsg_jetsMask) & (passMask)]),
                                                                                              jet_pt[(udsg_jetsMask) & (passMask)],

                                                                                             )
                    SF[failMask] = np.zeros_like(dummyZeros)
                    
                    if 'up' in systematic.lower():
                        
                        efficiency[(b_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_up_b_b_{self.year}{self.sel}'].evaluate(jet_pt[(b_jetsMask) & (passMask)],
                                                                                                                                      np.abs(jet_eta[(b_jetsMask) & (passMask)]))

                        efficiency[(c_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_up_b_c_{self.year}{self.sel}'].evaluate(jet_pt[(c_jetsMask) & (passMask)],
                                                                                                                                      np.abs(jet_eta[(c_jetsMask) & (passMask)]))

                        efficiency[(udsg_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_up_b_udsg_{self.year}{self.sel}'].evaluate(jet_pt[(udsg_jetsMask) & (passMask)],
                                                                                                                                            np.abs(jet_eta[(udsg_jetsMask) & (passMask)]))
                        efficiency[failMask] = np.zeros_like(dummyZeros)

                    elif 'down' in systematic.lower():
                        efficiency[(b_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_down_b_b_{self.year}{self.sel}'].evaluate(jet_pt[(b_jetsMask) & (passMask)],
                                                                                                                                   np.abs(jet_eta[(b_jetsMask) & (passMask)]))

                        efficiency[(c_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_down_b_c_{self.year}{self.sel}'].evaluate(jet_pt[(c_jetsMask) & (passMask)],
                                                                                                                                   np.abs(jet_eta[(c_jetsMask) & (passMask)]))

                        efficiency[(udsg_jetsMask) & (passMask) ] = self.effCorr_b[f'efficiency_down_b_udsg_{self.year}{self.sel}'].evaluate(jet_pt[(udsg_jetsMask) & (passMask)],
                                                                                                                                         np.abs(jet_eta[(udsg_jetsMask) & (passMask)]))
                        efficiency[failMask] = np.zeros_like(dummyZeros)
                        
                    

                    new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}']=SF*efficiency
                    

        
        return new_fields
    
    def include_leptonWeights(self):#, btvjson, effHist_b, effHist_c, effHist_udsg):
        
        # the setting of the keys for the systematic suffix map can be raised to the global/class level but keeping things hard-coded here since my needs are fixed and I use custom nomenclature for self-created branches in ntuples, and it meshes better with how I do my histo productions
        
        #for non-signal MC, signal MC sys variations on jes/jer, or signal MC when run without exp/theory weight unc. variations

        if self.isMC:
            muon_SF_labels = [  
                                'ID',
                                'ISO',
                                'Trig',
                                'RecoEff',
                             ]

        if (self.isMC and not(self.isSigMC)) or (self.isSigMC and not(self.wtUnc)): 

            systematic_suffix_map = {
                                     'nominal': 'Nom',
                                    }
        
        elif (self.isSigMC and self.wtUnc):
            
            systematic_suffix_map = {
                                     'nominal': 'Nom',
                                     'up_stat': 'StatUp',
                                     'down_stat': 'StatDown',
                                     'up_syst': 'SystUp',
                                     'down_syst': 'SystDown',
                                     
                                    }
        
        
        new_fields = OrderedDict()
        new_fields[f'new{self.sel}_leptonWeightNom_nom'] = []

        
        if not(self.sysUnc) or ('const' in self.onlyUnc): sysList = ['_nom']
        else: sysList=self.sysSources

        if self.verbose: print(f"Computing lepton weights ({systematic_suffix_map.keys()}) for following sys sources: {sysList}, using {muon_SF_labels} SFs")
        
    
        for sys in sysList:
            
            s='_nom' if sys.startswith(('_jes','_jer','_const')) else sys
            passMask = (self.events[f'passRecoSel{s}']==1)
            failMask = (self.events[f'passRecoSel{s}']!=1)

            ############### nominal ################
            dummyOnes = np.ones(len(self.events))
            dummyZeroes=np.zeros(len(self.events[f'selRecoMu_nom_pt'][failMask]))

            ID_SFs, ISO_SFs, RecoEff_SFs, Trig_SFs = np.ones_like(dummyOnes),np.ones_like(dummyOnes),np.ones_like(dummyOnes),np.ones_like(dummyOnes)
            
            ID_SFs[passMask], ISO_SFs[passMask], RecoEff_SFs[passMask], Trig_SFs[passMask] = compute_MuonSFs(self.muIDISOjson, 
                                                                                                             self.muRecoEffjson, self.muTrigjson,
                                                                                                             'nominal', 
                                                                                                             np.array(self.events[f'selRecoMu_nom_eta'][passMask]), 
                                                                                                             np.array(self.events[f'selRecoMu_nom_p'][passMask]), 
                                                                                                             np.array(self.events[f'selRecoMu_nom_pt'][passMask]))
            #ID_SFs_pass, ISO_SFs_pass, RecoEff_SFs_pass, Trig_SFs_pass
            ID_SFs[failMask], ISO_SFs[failMask], RecoEff_SFs[failMask], Trig_SFs[failMask] = np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes)#ID_SFs_fail, ISO_SFs_fail, RecoEff_SFs_fail, Trig_SFs_fail
            
            new_fields[f'new{self.sel}_leptonWeightNom{sys}'] = ID_SFs*ISO_SFs*RecoEff_SFs*Trig_SFs
            #print(new_fields[f'new{self.sel}_leptonWeightNom{sys}'][0:5], ID_SFs[0:5], ISO_SFs[0:5], RecoEff_SFs[0:5], Trig_SFs[0:5])    
            #print( sys, ID_SFs)
            #print( sys, ISO_SFs)
            #print( sys, RecoEff_SFs)
            #print( sys, Trig_SFs)
            #print( sys, new_fields[f'new{self.sel}_leptonWeightNom{sys}'])
            
            if not(self.sysUnc) and (self.wtUnc and self.isSigMC): 
                ############## Syst./stat. variations #############
                ID_SF_statUnc,ISO_SF_statUnc,RecoEff_SF_statUnc,Trig_SF_statUnc = np.ones_like(dummyOnes), np.ones_like(dummyOnes), np.ones_like(dummyOnes), np.ones_like(dummyOnes)
                ID_SF_systUnc, ISO_SF_systUnc, RecoEff_SF_systUnc, Trig_SF_systUnc = np.ones_like(dummyOnes), np.ones_like(dummyOnes), np.ones_like(dummyOnes), np.ones_like(dummyOnes)
 

                

                ID_SF_statUnc[passMask], ISO_SF_statUnc[passMask], RecoEff_SF_statUnc[passMask], Trig_SF_statUnc[passMask] = compute_MuonSFs(self.muIDISOjson, self.muRecoEffjson, self.muTrigjson, 'stat', (self.events[f'selRecoMu_nom_eta'][passMask]), (self.events[f'selRecoMu_nom_p'][passMask]), (self.events[f'selRecoMu_nom_pt'][passMask]))

                ID_SF_statUnc[failMask], ISO_SF_statUnc[failMask], RecoEff_SF_statUnc[failMask], Trig_SF_statUnc[failMask] = np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes)
            
            
                ID_SF_systUnc[passMask], ISO_SF_systUnc[passMask], RecoEff_SF_systUnc[passMask], Trig_SF_systUnc[passMask] = compute_MuonSFs(self.muIDISOjson, self.muRecoEffjson, self.muTrigjson, 'syst', (self.events[f'selRecoMu_nom_eta'][passMask]), (self.events[f'selRecoMu_nom_p'][passMask]), (self.events[f'selRecoMu_nom_pt'][passMask]))

                ID_SF_systUnc[failMask], ISO_SF_systUnc[failMask], RecoEff_SF_systUnc[failMask], Trig_SF_systUnc[failMask] = np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes),np.zeros_like(dummyZeroes)
                
                #print( sys, ID_SFs,ID_SF_statUnc,ID_SF_systUnc)
                #print( sys, ISO_SFs,ISO_SF_statUnc,ISO_SF_systUnc)
                #print( sys, RecoEff_SFs,RecoEff_SF_statUnc,RecoEff_SF_systUnc)
                #print( sys, Trig_SFs,Trig_SF_statUnc,Trig_SF_systUnc)
                for corr in muon_SF_labels:
                    for s in systematic_suffix_map.keys():
                        if 'nominal' in s: continue


                        if 'up_' in s: 

                            if 'ID' in corr: 
                                #print( sys, corr, s , ID_SFs,ID_SF_statUnc,ID_SF_systUnc, )
                                #print( sys, corr, s , ID_SFs+ID_SF_statUnc,ID_SFs+ID_SF_systUnc )

                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (ID_SFs)+(ID_SF_statUnc if 'stat' in s else ID_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*ISO_SFs*RecoEff_SFs*Trig_SFs

                            elif 'ISO' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (ISO_SFs)+(ISO_SF_statUnc if 'stat' in s else ISO_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*RecoEff_SFs*Trig_SFs

                            elif 'Trig' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (Trig_SFs)+(Trig_SF_statUnc if 'stat' in s else Trig_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*ISO_SFs*RecoEff_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']#*Trig_SFs

                            elif 'RecoEff' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (RecoEff_SFs)+(RecoEff_SF_statUnc if 'stat' in s else RecoEff_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*ISO_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*Trig_SFs

                        elif 'down_' in s: 

                            if 'ID' in corr: 
                                #print( sys, corr, s , ID_SFs,ID_SF_statUnc,ID_SF_systUnc, )
                                #print( sys, corr, s , ID_SFs-ID_SF_statUnc,ID_SFs-ID_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (ID_SFs)-(ID_SF_statUnc if 'stat' in s else ID_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*ISO_SFs*RecoEff_SFs*Trig_SFs

                            elif 'ISO' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (ISO_SFs)-(ISO_SF_statUnc if 'stat' in s else ISO_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*RecoEff_SFs*Trig_SFs

                            elif 'Trig' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (Trig_SFs)-(Trig_SF_statUnc if 'stat' in s else Trig_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*ISO_SFs*RecoEff_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']

                            elif 'RecoEff' in corr: 
                                new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}'] = (RecoEff_SFs)-(RecoEff_SF_statUnc if 'stat' in s else RecoEff_SF_systUnc)
                                new_fields[f'new{self.sel}_leptonWeight{corr}{systematic_suffix_map[s]}{sys}'] = ID_SFs*ISO_SFs*new_fields[f'new{self.sel}_leptonSF{corr}{systematic_suffix_map[s]}{sys}']*Trig_SFs


           

        #print(f"New reweighting/SF branches for muons: {new_fields.keys()}")

        # Convert lists to arrays
        # for key in new_fields:
        #    print(key)
        #    #if not(type(new_fields[key])==np.ma.core.MaskedArray): new_fields[key] = np.ma.core.MaskedArray(new_fields[key],dtype=np.float32) #.ndarray
        #    print("Event 0-10 weight/SF, key", key, new_fields[key][0:10])

        return new_fields
    
        
#####################################################
###################### Helpers ######################
#####################################################

def compute_MuonSFs(mujsonIDISO, mujsonRecoEff, mujsonTrig, systematic, eta, momentum, pt):#, iEvt):
    
    ID_SFs = mujsonIDISO["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate( np.abs(eta), pt, systematic)
    ISO_SFs = mujsonIDISO["NUM_probe_TightRelTkIso_DEN_HighPtProbes"].evaluate( np.abs(eta), pt, systematic)
    RecoEff_SFs = mujsonRecoEff["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate( np.abs(eta), momentum, systematic)
    TriggerEff_SFs = mujsonTrig["NUM_HLT_DEN_HighPtTightRelIsoProbes"].evaluate( np.abs(eta), pt, systematic)

    #print("Event 0-10, syst, SFs ID,ISO,Reco eff., Trig. eff. =", systematic, ID_SFs[0:10], ISO_SFs[0:10], RecoEff_SFs[0:10], TriggerEff_SFs[0:10])

    return ID_SFs, ISO_SFs, RecoEff_SFs, TriggerEff_SFs

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
    
    


def make_simple_plot(plotDict,nBins = 20):
    canvas = ROOT.TCanvas('can_simple', 'can_simple', 800, 600)
    histDict = {}
    color=1
    for k in plotDict.keys():
        
        histDict[k] = ROOT.TH1F(k,k,nBins,min(plotDict[k]),max(plotDict[k]))
        histDict[k].SetDirectory(0)
        histDict[k].SetLineColor(color)
        if color==1:
            histDict[k].Draw('hist E')
        else:
            histDict[k].Draw('hist E same')
        color+=1
    return canvas
        
    
################# DEPRECATED FUNCTIONS #################
"""
def include_btagWeights_old(self):
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
                                     'down_uncorrelated': 'UncorrelatedDown',
                                     'up_efficiency': 'EfficiencyUp',
                                     'down_efficiency': 'EfficiencyDown',
                                
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
        
        
        for i in range(0, len(self.events)):
                
            for sys in sysList:     
                
                if self.events[f'passRecoSel{sys}'][i]!=1:# or self.events[f'totalRecoWeight{sys}'][i]==0. or (i>=1 and self.events['selRecoAK4bjetLeptHem_nom_pt'][i]==self.events['selRecoAK4bjetLeptHem_nom_pt'][i-1] and self.events[f'passRecoSel{sys}'][i]==self.events[f'passRecoSel{sys}'][i-1]):
                    for systematic in systematic_suffix_map.keys():
                        new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(0.)
                else:
                    
                    #not calculating weights based on more than one jet in consideration (pick leading lept hem AK4, check if b-tagged)
                    #so, btagging wt. is then given by eff*SF for just that jet, but keeping masking based approach since it's
                    #generalisable to a degree and can be modified easily for >1 b-jet
                    
                    
                    try: 
                        #if not(type(self.events['Jet_pt'][i]))==list:
                        #    #.index(self.events['selRecoAK4bjetLeptHem_nom_pt'][i])
                        j = list(self.events['Jet_pt'][i]).index(self.events['selRecoAK4bjetLeptHem_nom_pt'][i])
                        #else:
                        #    j = (self.events['Jet_pt'][i]).index(self.events['selRecoAK4bjetLeptHem_nom_pt'][i])

                    except TypeError: 
                        #print((self.events['Jet_pt'][i]))
                        print(self.events['selRecoAK4bjetLeptHem_nom_pt'][0], self.events['selRecoAK4bjetLeptHem_nom_pt'][i],self.events['selRecoAK4bjetLeptHem_nom_pt'][i-1 if i>0 else i])
                        print(self.events['Jet_pt'][0], self.events['Jet_pt'][i],self.events['Jet_pt'][i-1 if i>0 else i])

                        try: 
                            j = list([self.events['Jet_pt'][i]]).index(self.events['selRecoAK4bjetLeptHem_nom_pt'][i])
                        except ValueError:
                            for systematic in systematic_suffix_map.keys():
                                new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(1.)
                            continue
                    
                    #if type(self.events['Jet_pt'][i])==list:
                    jet_pt = np.array([self.events['Jet_pt'][i][j]])
                    jet_eta = np.array([self.events['Jet_eta'][i][j]])
                    jet_ID = np.array([self.events['Jet_jetId'][i][j]])
                    jet_flav = np.array([self.events['Jet_hadronFlavour'][i][j]])
                    jet_discr = np.array([self.events['Jet_btagDeepFlavB'][i][j]])
                    
                    #elif type(self.events['Jet_pt'][i])==float:
                        
                    #    jet_pt = np.array([self.events['Jet_pt'][i]])
                    #    jet_eta = np.array([self.events['Jet_eta'][i]])
                    #    jet_ID = np.array([self.events['Jet_jetId'][i]])
                    #    jet_flav = np.array([self.events['Jet_hadronFlavour'][i]])
                    #    jet_discr = np.array([self.events['Jet_btagDeepFlavB'][i]])

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
                    b_efficiencies = np.array([get_efficiency(self.effCorr_b, p, e)[0] for p, e in zip(pt[b_jets], eta[b_jets])])#np.array([1])#
                    c_efficiencies = np.array([get_efficiency(self.effCorr_c, p, e)[0] for p, e in zip(pt[c_jets], eta[c_jets])])#np.array([1])#
                    light_efficiencies = np.array([get_efficiency(self.effCorr_udsg, p, e)[0] for p, e in zip(pt[light_jets], eta[light_jets])])#np.array([1])#
                    
                    #Get efficiency errors
                    b_efficiencies_err = np.array([get_efficiency(self.effCorr_b, p, e)[1] for p, e in zip(pt[b_jets], eta[b_jets])])#np.array([1])#
                    c_efficiencies_err = np.array([get_efficiency(self.effCorr_c, p, e)[1] for p, e in zip(pt[c_jets], eta[c_jets])])#np.array([1])#
                    light_efficiencies_err = np.array([get_efficiency(self.effCorr_udsg, p, e)[1] for p, e in zip(pt[light_jets], eta[light_jets])])#np.array([1])#
                    
                    #if i%25000==0: 
                    #    print('\n',"Event #:", i)

                    # Compute SFs for each systematic
                    for systematic in systematic_suffix_map.keys():#self.btaggingSystLabels:
                        
                        if not('efficiency' in systematic):
                            _, _, _, SFs = compute_btagSFAndWeight(self.btvjson, systematic, self.btagWP, 
                                                                   flav, eta, pt, 
                                                                   b_jets, c_jets, light_jets,
                                                                   b_efficiencies, c_efficiencies, light_efficiencies,i)

                            #if i%10000==0: 
                            #    print("Event efficiencies b,c,udsg,and jet pt, jet eta =", b_efficiencies, c_efficiencies, light_efficiencies, pt, eta)

                            #    print(f"Event weight {systematic} = {SFs}, {np.prod(SFs)}")

                            #    #print("Event weight old =", i, self.events['btagWeightNom_nom'][i])

                            new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(np.prod(SFs))

                        else:
                            
                            
                            if 'up' in systematic.lower():
                                
                                _, _, _, SFs = compute_btagSFAndWeight(self.btvjson, 'central', self.btagWP, 
                                                                       flav, eta, pt, 
                                                                       b_jets, c_jets, light_jets,
                                                                       b_efficiencies + b_efficiencies_err, 
                                                                       c_efficiencies + c_efficiencies_err, 
                                                                       light_efficiencies + light_efficiencies_err,i)
                            elif 'down' in systematic.lower():
                                _, _, _, SFs = compute_btagSFAndWeight(self.btvjson, 'central', self.btagWP, 
                                                                       flav, eta, pt, 
                                                                       b_jets, c_jets, light_jets,
                                                                       b_efficiencies - b_efficiencies_err, 
                                                                       c_efficiencies - c_efficiencies_err, 
                                                                       light_efficiencies - light_efficiencies_err,i)
                                
                            #if i%10000==0: 
                            #    print("Event efficiencies b,c,udsg,and jet pt, jet eta =", b_efficiencies, c_efficiencies, light_efficiencies, pt, eta)

                            #    print(f"Event weight {systematic} = {SFs}, {np.prod(SFs)}")

                            #    #print("Event weight old =", i, self.events['btagWeightNom_nom'][i])

                            new_fields[f'new{self.sel}_btagWeight{systematic_suffix_map[systematic]}{sys}'].append(np.prod(SFs))
                            

        # Convert lists to arrays
        #for key in new_fields:
        #    if not(type(new_fields[key])==np.ma.core.MaskedArray): new_fields[key] = np.ma.core.MaskedArray(new_fields[key],dtype=np.float32)
        #    #if not(type(new_fields[key])==np.ndarray): new_fields[key] = np.array(new_fields[key])

        return new_fields
    
def compute_btagSFAndWeight(btvjson, systematic, WP,
                            flav, eta, pt, 
                            b_jets, c_jets, light_jets,
                            b_efficiencies, c_efficiencies, light_efficiencies, iEvt
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

    #if iEvt%25000==0:
    #    print(iEvt,"Event SFs b,c,udsg =", bJet_SFs, cJet_SFs, lightJet_SFs)

    bJet_SFs *= b_efficiencies 
    cJet_SFs *= c_efficiencies
    lightJet_SFs *= light_efficiencies
    all_SFs = np.concatenate((bJet_SFs, cJet_SFs, lightJet_SFs))

    return bJet_SFs, cJet_SFs, lightJet_SFs, all_SFs
    
def get_efficiency_old(hist, pt, eta):
    '''Reads efficiency maps for tagging a b/c/udsg as a b; implementation based on uproot and hist'''
    # deprecated since slower somehow a la %timeit
    
    eta_edges = hist.axis('y').edges()
    pt_edges = hist.axis('x').edges()

    biny = np.digitize(eta, eta_edges)-1
    binx = np.digitize(pt, pt_edges)-1

    biny = np.clip(biny, 0, len(eta_edges) - 2)
    binx = np.clip(binx, 0, len(pt_edges) - 2)

    bin_contents_eff = hist.values()

    return bin_contents_eff[binx, biny]

def get_efficiency(histogram, pt, eta):
    '''Reads efficiency maps for tagging a b/c/udsg as a b; implementation based on (py)ROOT'''
    pt_bin = histogram.GetXaxis().FindFixBin(pt)
    eta_bin = histogram.GetYaxis().FindFixBin(eta)
    
    # Handle overflows
    n_pt_bins = histogram.GetNbinsX()
    n_eta_bins = histogram.GetNbinsY()
    
    if pt_bin > n_pt_bins:
        pt_bin = n_pt_bins
    if eta_bin > n_eta_bins:
        eta_bin = n_eta_bins
    
    efficiency = histogram.GetBinContent(pt_bin, eta_bin)
    error = histogram.GetBinError(pt_bin, eta_bin)
    
    return [efficiency, error]
"""

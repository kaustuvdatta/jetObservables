from collections import OrderedDict
import copy
import pprint 
import ROOT
import numpy as np
import array
from array import array
import bisect
#from legend import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from root_numpy import array2hist, hist2array
#import histoHelpers
from histoHelpers import *
#import unfoldingPlottersAndHelpers
from unfoldingPlottersAndHelpers import *
import os
import glob
import sys
import math
import yoda
#sys.path.insert(0,'../test/')
#from DrawHistogram_dijetSel import *

sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
sys.path.insert(0,'../../')
from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict

#ROOT.gROOT.ForceStyle()
#tdrstyle.setTDRStyle()
import gc
#import pdb
lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)

##############################################################################################

#################################### MAIN Unfolding script ###################################

##############################################################################################



def runTUnfold(
                dataFile, sigFiles, bkgFiles, variables, sel, sysUncert, process, ext, lumi=1., sysSignalLabels=[], 
                year='2017',runMLU=False, sysSigFiles={}, varSigFiles={}, outputFolder='../Results',
                version='_Sept23',verbose=False, return_tunfolder_object = False, mainMC='MLM_HTbin', altMC='Ptbin',extraMC=False,
                noScale_dijets=False, include_FSR_in_unfolded_result=False, areaConstraint=None, dict_t3_samples=None
              ):
    """ 
    Response matrix for unfoldings have miss corrections applied already a la UF bin on y-axis;
    a) this means that the unfolded distribution's integral should match the overall gen-count in closure tests (not just accepgen counts), 
    b) and, the folded unfolded distribution from tunfold will match the true-reco in closure tests and data-bkg in data unfolding;
    """    
    #selection=sel
    
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    ROOT.TH1.StatOverflows(ROOT.kTRUE)
    ROOT.TH2.StatOverflows(ROOT.kTRUE)

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 41, 45, 48, 49, 37, 38, 33, 17]
    years_list = ['2016_preVFP','2016','2017','2018']
    print(variables.keys())
    import sys
    
    if not('dijet' in sel):
        sys.path.insert(0,'../../')
        from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict
        
        dictSamples = dict_t3_samples
        
        signalLabel = 'TTTo'#SemiLeptonic
        sigPlotLabel = 'Powheg+Pythia8'
        signalLabelBegin = 'TTTo'#SemiLeptonic
        varSignalLabelBegin = 'varTTTo'#Semileptonic
        sysSignalLabelBegin = 'sysTTTo'#Semileptonic
        fsrLabel = 'sysTTToSemiLeptonic_fsrWeight'#SemiLeptonic
        #altSignalLabelBegin = 'TTJets'
        #altSigPlotLabel = 'aMC@NLO-FXFX+Pythia8'
        #altSignalLabel = 'TTJets'
        if altMC.startswith('TT_'):
            altSignalLabelBegin = 'TT_TuneCH3'
            altSigPlotLabel = 'Powheg+Herwig7'
            altSignalLabel = 'TT_TuneCH3'
            alt1SignalLabelBegin = 'TTJets'
            alt1SigPlotLabel = 'aMC@NLO-FXFX+Pythia8'
            alt1SignalLabel = 'TTJets'
            alt1 = 'TTJets'
            alt1MC = 'TTJets'
            
        elif altMC.startswith('TTJets'):
            altSignalLabelBegin = 'TTJets'
            altSigPlotLabel = 'aMC@NLO-FXFX+Pythia8'
            altSignalLabel = 'TTJets'
            alt1SignalLabelBegin = 'TT_TuneCH3'
            alt1SigPlotLabel = 'Powheg+Herwig7'
            alt1SignalLabel = 'TT_TuneCH3'
            alt1 = 'TT_TuneCH3'
            alt1MC = 'TT_TuneCH3'
        
        alt2SignalLabelBegin = None
        alt2SigPlotLabel = None
        alt2SignalLabel = None
        
        
        varSignalLabels = ['varTTToSemileptonic_TuneCP5Up',
                           'varTTToSemileptonic_TuneCP5Down',
                           'varTTToSemileptonic_hdampUp_TuneCP5', 
                           'varTTToSemileptonic_hdampDown_TuneCP5', 
                           'varTTToSemileptonic_TuneCP5_erdON', 
                           'varTTToSemileptonic_TuneCP5CR1', 
                           'varTTToSemileptonic_TuneCP5CR2', 
                           'varTTToSemileptonic_mtop171p5_TuneCP5', 
                           'varTTToSemileptonic_mtop173p5_TuneCP5', 

                          ]

        bkgLabels = [
                     'TTToHadronic', 'TTTo2L2Nu', 
                     'WJetsToLNu', 
                     'ST_s-channel_4f_leptonDecays', 
                     'ST_t-channel_top_5f_InclusiveDecays', 'ST_t-channel_antitop_5f_InclusiveDecays', 
                     'ST_tW_top_5f_NoFullyHadronicDecays', 'ST_tW_antitop_5f_NoFullyHadronicDecays', 
                     'WW', 'ZZ', 'WZ', 
                     'QCD_Pt-1000_MuEnrichedPt5'
                    ]
    else:
        sys.path.insert(0,'../../')
        from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
        
        dictSamples = dict_t3_samples
        '''
        if '2016' in year:
            sysSignalLabels = [
                            #'sysMLMQCD_puWeight_HT2000toInf', 'sysMLMQCD_isrWeight_HT2000toInf',
                            'sysMLMQCD_fsrWeight_HT2000toInf',
                            'sysMLMQCD_jer_HT2000toInf', 
                            #'sysMLMQCD_pdfWeight_HT2000toInf', 'sysMLMQCD_l1prefiringWeight_HT2000toInf', 
                            #'sysMLMQCD_jesAbsolute_HT2000toInf', 'sysMLMQCD_jesFlavorQCD_HT2000toInf', 'sysMLMQCD_jesHF_HT2000toInf',
                            #'sysMLMQCD_jesRelativeBal_HT2000toInf', 'sysMLMQCD_jesBBEC1_HT2000toInf', 'sysMLMQCD_jesEC2_HT2000toInf',
                            #'sysMLMQCD_jesAbsolute_2016_HT2000toInf', 'sysMLMQCD_jesBBEC1_2016_HT2000toInf', 'sysMLMQCD_jesEC2_2016_HT2000toInf',
                            #'sysMLMQCD_jesHF_2016_HT2000toInf', 'sysMLMQCD_jesRelativeSample_2016_HT2000toInf'
                        ]
        elif '2017' in year:
            sysSignalLabels = [
                            #'sysMLMQCD_puWeight_HT2000toInf', 'sysMLMQCD_isrWeight_HT2000toInf',
                             'sysMLMQCD_fsrWeight_HT2000toInf',
                            'sysMLMQCD_jer_HT2000toInf', 
                            #'sysMLMQCD_pdfWeight_HT2000toInf', 'sysMLMQCD_l1prefiringWeight_HT2000toInf', 
                            #'sysMLMQCD_jesAbsolute_HT2000toInf', 'sysMLMQCD_jesFlavorQCD_HT2000toInf', 'sysMLMQCD_jesHF_HT2000toInf',
                            #'sysMLMQCD_jesRelativeBal_HT2000toInf', 'sysMLMQCD_jesBBEC1_HT2000toInf', 'sysMLMQCD_jesEC2_HT2000toInf',
                            #'sysMLMQCD_jesAbsolute_2017_HT2000toInf', 'sysMLMQCD_jesBBEC1_2017_HT2000toInf', 'sysMLMQCD_jesEC2_2017_HT2000toInf',
                            #'sysMLMQCD_jesHF_2017_HT2000toInf', 'sysMLMQCD_jesRelativeSample_2017_HT2000toInf'
                        ]
        elif '2018' in year:
            sysSignalLabels = [
                            #'sysMLMQCD_puWeight_HT2000toInf', 'sysMLMQCD_isrWeight_HT2000toInf',
                             'sysMLMQCD_fsrWeight_HT2000toInf',
                            'sysMLMQCD_jer_HT2000toInf', 
                            #'sysMLMQCD_pdfWeight_HT2000toInf', 'sysMLMQCD_l1prefiringWeight_HT2000toInf',
                            #'sysMLMQCD_jesAbsolute_HT2000toInf', 'sysMLMQCD_jesFlavorQCD_HT2000toInf', 'sysMLMQCD_jesHF_HT2000toInf',
                            #'sysMLMQCD_jesRelativeBal_HT2000toInf', 'sysMLMQCD_jesBBEC1_HT2000toInf', 'sysMLMQCD_jesEC2_HT2000toInf',
                            #'sysMLMQCD_jesAbsolute_2018_HT2000toInf', 'sysMLMQCD_jesBBEC1_2018_HT2000toInf', 'sysMLMQCD_jesEC2_2018_HT2000toInf',
                            #'sysMLMQCD_jesHF_2018_HT2000toInf', 'sysMLMQCD_jesRelativeSample_2018_HT2000toInf'
                        ]
        ''' 
        if extraMC:
            alt1MC='HTbin'
            alt2MC='Ptbin'
        else:
            alt1MC=None
            alt2MC=None
        
        if mainMC.startswith('MLM_HTbin'):
            signalLabelBegin = 'MLMQCD_HT'
            signalLabel = 'MLMQCD_HT2000toInf'

        elif mainMC.startswith('HTbin'):
            signalLabelBegin = 'QCD_HT'
            signalLabel = 'QCD_HT2000toInf'

        elif mainMC.startswith('H7MLM_HTbin'):
            signalLabelBegin = 'H7MLMQCD_HT'
            signalLabel = 'H7MLMQCD_HT2000toInf'

        if altMC.startswith('MLM_HTbin'):
            altSignalLabelBegin = 'MLMQCD_HT'
            altSignalLabel = 'MLMQCD_HT2000toInf'

        elif altMC.startswith('HTbin'):
            altSignalLabelBegin = 'QCD_HT'
            altSignalLabel = 'QCD_HT2000toInf'

        elif altMC.startswith('H7MLM_HTbin'):
            altSignalLabelBegin = 'H7MLMQCD_HT'
            altSignalLabel = 'H7MLMQCD_HT2000toInf'
            #altSignalLabel = 'QCD_Pt-15to7000'
            
        if extraMC:
            if alt1MC.startswith('HTbin'):
                alt1SignalLabelBegin = 'QCD_HT'
                alt1SignalLabel = 'QCD_HT2000toInf'

            elif alt1MC.startswith('Ptbin'):
                alt1SignalLabelBegin = 'QCD_Pt_'
                alt1SignalLabel = 'QCD_Pt_3200toInf'

            elif alt1MC.startswith('H7MLM_HTbin'):
                alt1SignalLabelBegin = 'H7MLMQCD_HT'
                alt1SignalLabel = 'H7MLMQCD_HT2000toInf'
                #altSignalLabel = 'QCD_Pt-15to7000'


            if alt2MC.startswith('Ptbin'):
                alt2SignalLabelBegin = 'QCD_Pt_'
                alt2SignalLabel = 'QCD_Pt_3200toInf'

            elif alt2MC.startswith('HTbin'):
                alt2SignalLabelBegin = 'QCD_HT'
                alt2SignalLabel = 'QCD_HT2000toInf'

            elif alt2MC.startswith('H7MLM_HTbin'):
                alt2SignalLabelBegin = 'H7MLMQCD_HT'
                alt2SignalLabel = 'H7MLMQCD_HT2000toInf'
                #altSignalLabel = 'QCD_Pt-15to7000'
            
            
        sysSignalLabelBegin = 'sysMLMQCD'  if mainMC.startswith('MLM') else 'sysQCD' 
        fsrLabel = 'sysMLMQCD_fsrWeight_HT2000toInf' if mainMC.startswith('MLM') else 'sysQCD_fsrWeight_HT2000toInf' 
        
    try: print ("Labels:", signalLabelBegin,altSignalLabelBegin,alt1SignalLabelBegin,alt2SignalLabelBegin,fsrLabel)
    except: pass
    
    
    
    for ivar in variables:
        gc.collect()
        if not('tau' in ivar) : continue
        print (ivar)
        outputDir=outputFolder+sel.split('_')[1]+'/'+year+'/Unfolding/'+ivar+'/'+process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        genBin = variables[ivar]['bins']
        recoBin = variables[ivar]['bins_reco']

        ###########################################################################################
        if year.startswith('all'):
            dataHistos = { }
            dataHistostrue = {}
            bkgHistos = { }
            sysUncs_added = []
            signalHistos = {
                    signalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_gen'+ivar+'_nom'+sel).Clone(),
                
                    signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin': dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin').Clone(),
                
                    signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin': dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin').Clone()

            }
            for yt in years_list[1:]:
                print(yt,signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), dataFile[ivar+f'_{yt}'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Integral(), 
                     (dataFile[ivar+f'_{yt}'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Integral()+signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Integral()))
                
                signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_gen'+ivar+'_nom'+sel).Clone() 
                )

                signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin').Clone() 
                )

                signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Add( 
                    
                    dataFile[ivar+f'_{yt}'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin').Clone() 
                )
            print(signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Integral())
            
            if not process.startswith('MC'):
                sysSignalHistos={}

                print ("Loading up all syst. variations from the following:", sysUncert)
                for isys,sys in enumerate(sysUncert):
                    #print (isys,sys)
                    if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                    for upDown in [ 'Up', 'Down' ]:
                        #print (isys,sys,upDown)
                        if sys+upDown in sysUncs_added: 
                            print("unc already added in, moving to next unc.",sys+upDown)
                            continue
                        s = [i for i in sysSignalLabels if sys in i]
                        #if verbose: 
                        #print ('s=',s)
                        s=[s[0]]    
                    
                        if ('2016' in sys or '2017' in sys or '2018' in sys):
                            if not ('2016' in s[0] or '2017' in s[0] or '2018' in s[0]):
                                print("Sys and s[0] mismatch, where each are (respectively):",sys,s[0])

                        #        if '2016' in s: s=[s[1]]
                        #        elif '2017' in s: s=[s[2]]
                        #        elif '2018' in s: s=[s[3]]
                        #    else:
                        #        s=[s[0]]
                        #    #print (s[0]+'_reco'+ivar+sys+upDown+sel)
                        #elif len(s)==1:
                        #    s=[s[0]]
                        #    #print (s[0]+'_reco'+ivar+sys+upDown+sel)
                        #else:                            
                        #    continue
                            
                        #print ('s=',s)
                        if (not ('2017' in sys) and not ( '2018' in sys ) and not ( '2016' in sys) ): 
                            #print(sys,s[0],s[0]+'_reco'+ivar+sys+upDown+sel,)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                            if 'fsr' in sys:
                                sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_gen'+ivar+sys+upDown+sel)
                                sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )
                                sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )
                                sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )

                            sysUncs_added.append(sys+upDown)
                        # dealing with uncorrelated jes unc sources below
                        elif '2016' in sys and not( '2017' in sys) and not('2018' in sys):# and s[0] in sys:#.endswith(sys): 
                            #print(sys+upDown,'s[0]=',s[0],s[0]+'_reco'+ivar+sys+upDown+sel)
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone() 
                            
                            #print(sys+upDown, s[0], sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Integral(), dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Integral())
                            
                            #add 2016 systematic to 2016_preVFP systematic
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone())  
                            
                            #print(sys+upDown, s[0], sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Integral(),                                 dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Integral())
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Clone())
                            #print(sys+upDown, s[0], sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Integral(),                                 dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Integral())
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Clone())
                            #print(sys+upDown, s[0], sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Integral(),                                 dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel).Integral())
                            
                            
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            
                            
                            sysUncs_added.append(sys+upDown)

                        elif '2017' in sys and not( '2016' in sys) and not('2018' in sys):# and s[0] in sys:#s[0].endswith(sys) : 
                            #print(sys,s[0],s[0]+'_reco'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysUncs_added.append(sys+upDown)
                            
                        elif '2018' in sys and not( '2016' in sys) and not('2017' in sys):# and s[0] in sys:#s[0].endswith(sys): 
                            #print(sys,s[0],s[0]+'_reco'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysUncs_added.append(sys+upDown)
                    #print(sysUncs_added)
                        ################ added recohistos & resp matrices from nominal_2018(/2017/2016) and jesUncorrUnc_2017(/2018/2016) ########################
                        
            if not('self' in process.lower()):
                
                altSignalHistos = {
                    altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'),
                    
                    altSignalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin'),
                    
                    altSignalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_missgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_missgen'+ivar+'_nom'+sel),
                    
                    altSignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel),
                    }
                print(altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Integral())
                for yt in years_list[1:]:
                    #print(yt,altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), 
                    #      dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel).Integral(), 
                    #      (dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel).Integral()+altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()))

                    altSignalHistos[ altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )

                    altSignalHistos[ altSignalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin') )

                    altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_missgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_missgen'+ivar+'_nom'+sel) )

                    altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel) )
                    print(yt,altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Integral())


                
            if not('mc' in process.lower()) and extraMC:
                alt1SignalHistos = {
                    alt1SignalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_respWithMiss'+ivar+'_nom'+sel),

                    alt1SignalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel),

                    alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'),

                    alt1SignalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_truereco'+ivar+'_nom'+sel),

                    alt1SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin'),

                    alt1SignalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_fakereco'+ivar+'_nom'+sel),

                    alt1SignalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_accepgen'+ivar+'_nom'+sel),
                    
                    alt1SignalLabel+'_missgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_missgen'+ivar+'_nom'+sel),

                    alt1SignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel),
                    }

                for yt in years_list[1:]:
                    #print(yt,alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), 
                    #      dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel).Integral(), 
                    #      (dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel).Integral()+alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()))

                    alt1SignalHistos[ alt1SignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )

                    alt1SignalHistos[ alt1SignalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_truereco'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin') )

                    alt1SignalHistos[ alt1SignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_fakereco'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_accepgen'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_missgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_missgen'+ivar+'_nom'+sel) )

                    alt1SignalHistos[ alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel) )
                
                if alt2SignalLabelBegin!=None:
                    alt2SignalHistos = {
                        alt2SignalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_respWithMiss'+ivar+'_nom'+sel),

                        alt2SignalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel),

                        alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'),

                        alt2SignalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_truereco'+ivar+'_nom'+sel),

                        alt2SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin'),

                        alt2SignalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_fakereco'+ivar+'_nom'+sel),

                        alt2SignalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_accepgen'+ivar+'_nom'+sel),
                        
                        alt2SignalLabel+'_missgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_missgen'+ivar+'_nom'+sel),

                        alt2SignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel),
                        }

                    for yt in years_list[1:]:
                        #print(yt,alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), 
                        #      dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel).Integral(), 
                        #      (dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel).Integral()+alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()))

                        alt2SignalHistos[ alt2SignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )

                        alt2SignalHistos[ alt2SignalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_truereco'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin') )

                        alt2SignalHistos[ alt2SignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_fakereco'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_accepgen'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_missgen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_missgen'+ivar+'_nom'+sel) )

                        alt2SignalHistos[ alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+f'_{yt}'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel) )
            
            if process.startswith('data') and sel.startswith(('_W','_top')):
                
                if verbose: print("Processing bkgs from amongst the following uncertainty sources:", bkgLabels)
                
                for ibkg in bkgLabels:
                    print(ibkg)
                    bkgHistos[ ibkg+'_reco'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(ibkg+'_reco'+ivar+'_nom'+sel)
                    bkgHistos[ ibkg+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(ibkg+'_reco'+ivar+'_nom'+sel) )
                    bkgHistos[ ibkg+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(ibkg+'_reco'+ivar+'_nom'+sel) )
                    bkgHistos[ ibkg+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(ibkg+'_reco'+ivar+'_nom'+sel) )

                    
                varSignalHistos={}
                s=[]

                
                if verbose: print ("Processing signal variations from amongst the following uncertainty sources: ", sysUncert)

                for sys in sysUncert:
                    print (sys)
                    s = [i for i in varSignalLabels if (sys.split('_')[1] in i)]
                    for j in s:
                        if 'Tune' in sys and not('TuneCP5Up' in j or 'TuneCP5Down' in j): continue
                        
                        #print (j+'_reco'+ivar+'_nom'+sel,s)
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_reco'+ivar+'_nom'+sel)
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_respWithMiss'+ivar+'_nom'+sel)
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )
            
            
            
            allHistos = {
                    'dataHisto' : dataFile[ivar+'_2016_preVFP'].Get( 'dataHisto' ),
                    'dataHistoGenBin' : dataFile[ivar+'_2016_preVFP'].Get( 'dataHistoGenBin' )
                    }
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2016'].Get( 'dataHisto' ) )
            allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2016'].Get( 'dataHistoGenBin' ) )
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2017'].Get( 'dataHisto' ) )
            allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2017'].Get( 'dataHistoGenBin' ) )
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2018'].Get( 'dataHisto' ) )
            allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2018'].Get( 'dataHistoGenBin' ) )

            dataHistos['data_reco'+ivar+'_nom'+sel] = allHistos[ 'dataHisto' ].Clone('data_reco'+ivar+'_nom'+sel)
            dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = allHistos[ 'dataHistoGenBin' ].Clone('data_reco'+ivar+'_nom'+sel+'_genBin')
            
            
            #dummy data copy: using truereco signal MC as data
            if 'MC' in process:
                dataHistostrue['data_reco'+ivar+'_nom'+sel] = signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Clone(signalLabel+'_data_truereco'+ivar+'_nom'+sel) 
                
                dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Clone( signalLabel+'_data_truereco'+ivar+'_nom'+sel+'_genBin' ) 
                
                dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin(len(genBin)-1, signalLabel+'_data_truereco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
            
            else:
                dataHistostrue['data_reco'+ivar+'_nom'+sel] = allHistos[ 'dataHisto' ].Clone()
                
                dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = allHistos[ 'dataHistoGenBin' ].Clone()
                
                    
            

        else:
            print('|-------> Running single year '+year)
            ### Getting input histos
            allHistos = {}
            
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
            
            tmp2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(altSignalLabelBegin)  }
            
            if not('self' in process.lower()): 
                altSignalHistos = loadHistograms( tmp2SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                  lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
                #for ih in altSignalHistos.keys(): print(altSignalHistos[ih].Integral(),altSignalHistos[ih].GetName())
                    
            if 'data' in process and extraMC:
                tmp3SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt1SignalLabelBegin)  }
                alt1SignalHistos = loadHistograms( tmp3SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                  lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
                if alt2SignalLabelBegin!=None:
                    tmp4SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt2SignalLabelBegin)  }
                    alt2SignalHistos = loadHistograms( tmp4SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                      lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )

            
            if process.startswith('data'): 
                sysSignalHistos = loadHistograms( sysSigFiles, ivar, sel, sysUnc=sysUncert, respOnly=False, isMC=True, lumi=lumi,
                                                  year=year, process=process, variables=variables,outputFolder=outputFolder)

                if sel.startswith(('_W','_top')): 
                    varSignalHistos = loadHistograms( varSigFiles, ivar, sel, sysUnc=[], respOnly=False, isMC=True, lumi=lumi,
                                                  year=year,process=process, variables=variables, outputFolder=outputFolder )
                    #print(varSignalHistos.keys())

            if sel.startswith(('_W','_top')): 

                bkgHistos = loadHistograms( bkgFiles, ivar, sel, sysUnc=[], respOnly=False, isMC=True, lumi=lumi, 
                                            year=year,process=process, variables=variables, outputFolder=outputFolder )
                
                #for ih in bkgHistos.keys():
                #    if 'recoJet' in ih and not('genBin' in ih) and not('true' in ih) and not('fake' in ih):
                #        print(ih, bkgHistos[ih].Integral())
            else: 
                
                bkgHistos = {}

            #if verbose: print ("All signal histos:", signalHistos)

              
                
            if process.startswith("MC"):
                dataHistostrue = { 'data_reco'+k.split(('_truereco'))[1] : v.Clone() for (k,v) in signalHistos.items() if ('_truereco' in k)}# and not ('genBin' in k ))}
                
               

                dataHistos = { 'data_reco'+k.split(('_reco'))[1] : v.Clone() for (k,v) in signalHistos.items() if ('_reco' in k)}# and not ('genBin' in k in k ))}
                
                if 'dijet' in sel:
                    dataHistosForRescaling = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process='data', variables=variables,outputFolder=outputFolder)
                    dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Clone('data_reco'+ivar+'_nom'+sel+'_genBin')
                    dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                
                
            else:
                dataHistostrue = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
                dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
               
        
        allHistos[ 'dataHisto' ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel ].Clone()
        allHistos[ 'dataHistoGenBin' ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()

        #fakes from nominal/signal MC
        fakeHistos = { signalLabel+'_fakereco'+k.split(('_fakereco'))[1]: v.Clone() for (k,v) in signalHistos.items()  if ('_fakereco' in k)}

        allHistos[ 'allBkgHisto' ] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
        allHistos[ 'allBkgHisto' ].Reset()
        allHistos[ 'allBkgHistoGenBin' ] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
        allHistos[ 'allBkgHistoGenBin' ].Reset()

        #include fakes in histo containing all bkg. (including other physics processes in the case of W/top) from MC
        for ih in fakeHistos:
            
            if ih.endswith(ivar+'_nom'+sel): 
                #print("FakeHisto Check:", ih)
                allHistos[ 'allBkgHisto' ].Add( fakeHistos[ih] )
            elif ih.endswith('_genBin'): 
                #print("FakeHisto Check:", ih)
                allHistos[ 'allBkgHistoGenBin' ].Add( fakeHistos[ih] )

        allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
        allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
        
        allHistos[ 'dataMinusBkgs' ] = allHistos[ 'dataHisto' ].Clone() 
        allHistos[ 'dataMinusBkgsGenBin' ] = allHistos[ 'dataHistoGenBin' ].Clone() 
        
        if process.startswith('data'):
            # if verbose: print(bkgHistos.keys())
            for ibkg in bkgHistos:
                # if verbose: print(f"Adding in {ibkg}")

                if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )
                if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone() )

            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )

        
        
        allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Clone() ) #add only true since fakes already contained in the allMC histo via the allBkgHisto
       
        
        allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Clone() ) 
        
        
        
        ################################################################################################################## 
        # rescale dijet MC to data for the case of individual years; reuse these rescaled histos in the all years case, thus the if

        if sel.startswith('_dijet') and not(year.startswith('all')): 

            scalingDict = {}    
            ### For dijet, scale QCD to data, as per AGE's past work
            scaleFactor=1.
            scaleFactorGenBin=1.
            altscaleFactor=1.
            altscaleFactorGenBin=1.
            nomIntegral = allHistos[ 'allMCHisto'].Integral()
            if process.startswith('MC'): 
                
                scaleFactor = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto'].Integral()
                scaleFactorGenBin = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' ].Integral()
                
                if not np.isclose(scaleFactor,scaleFactorGenBin): 
                    print(f"WARNING (something weird): gen-/reco-binning mc-to-data SFs are not the same, {scaleFactor,scaleFactorGenBin}")
                    print ("SF genBin and recobin, data int., nominal reco int., respectively:",scaleFactor,scaleFactorGenBin,dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral(),allHistos[ 'allMCHisto' ].Integral())
                
                for ihsig in dataHistos:
                    if ihsig.endswith(sel):
                        #print (ihsig)
                        dataHistos[ihsig].Scale( scaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactor

                    elif 'genBin' in ihsig: 
                        dataHistos[ihsig].Scale( scaleFactorGenBin )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin

                for ihsig in dataHistostrue:
                    if ihsig.endswith(sel):
                        dataHistostrue[ihsig].Scale( scaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactor
                    elif 'genBin' in ihsig: 
                        dataHistostrue[ihsig].Scale( scaleFactorGenBin )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin

                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        signalHistos[ihsig].Scale( scaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactor

                    elif 'genBin' in ihsig: 
                        signalHistos[ihsig].Scale( scaleFactorGenBin )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin
                        
                allHistos[ 'dataHisto' ].Scale( scaleFactor )
                allHistos[ 'dataHistoGenBin' ].Scale( scaleFactorGenBin )
                allHistos[ 'dataMinusBkgs' ].Scale( scaleFactor )
                allHistos[ 'dataMinusBkgsGenBin' ].Scale( scaleFactorGenBin )
                allHistos[ 'allMCHisto' ].Scale( scaleFactor )

                allHistos[ 'allMCHistoGenBin' ].Scale( scaleFactorGenBin )

                allHistos[ 'allBkgHistoGenBin' ].Scale( scaleFactorGenBin )

            if process.startswith('data'):
                                            
                scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto' ].Integral()
                scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' ].Integral()
                if not np.isclose(scaleFactor,scaleFactorGenBin): 
                    print(f"WARNING (something weird): gen-/reco-binning mc-to-data SFs are not the same, {scaleFactor,scaleFactorGenBin}")
                print ("SF genBin and recobin, data int., nominal reco int., respectively:",scaleFactor,scaleFactorGenBin,dataHistos['data_reco'+ivar+'_nom'+sel].Integral(),allHistos[ 'allMCHisto' ].Integral())
                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        #print (ihsig)
                        signalHistos[ihsig].Scale( scaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactor

                    elif ihsig.endswith('genBin'): 
                        signalHistos[ihsig].Scale( scaleFactorGenBin )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin
            
                allHistos[ 'allMCHisto' ].Scale( scaleFactor )
                scalingDict[f'scaling_allMCHisto'] = scaleFactor
                allHistos[ 'allMCHistoGenBin' ].Scale( scaleFactorGenBin )
                scalingDict[f'scaling_allMCHistoGenBin'] = scaleFactorGenBin

                #if process=='data':
                allHistos[ 'allBkgHisto' ].Scale( scaleFactor )
                scalingDict[f'scaling_allBkgHisto'] = scaleFactor
                allHistos[ 'allBkgHistoGenBin' ].Scale( scaleFactorGenBin )
                scalingDict[f'scaling_allBkgHistoGenBin'] = scaleFactorGenBin

                if len(sysUncert)!=0:
                    #scaleFactor_sys = 1.
                    #scaleFactorGenBin_sys = 1.
                    already_scaled = []
                    if verbose: 
                        #print(sysSignalLabels)
                        print(sysUncert)
                    for sys in sysUncert:
                        if sys in already_scaled: continue
                            
                        if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                        s = [i for i in sysSignalLabels if sys in i]
                        if verbose: 
                            print (sys,s)

                        for upDown in ["Up","Down"]:
                            if len(s)>1:
                                if ('2016' in sys or '2017' in sys or '2018' in sys) and 'jes' in sys:

                                    if '2016' in sys: s=[s[1]]
                                    elif '2017' in sys: s=[s[1]]
                                    elif '2018' in sys: s=[s[1]]
                                else:
                                    s=[s[0]]
                                #ihsig = s[0]+'_reco'+ivar+sys+upDown+sel

                            scaleFactor_sys = scaleFactor #if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[ihsig].Integral() #scaleFactor
                            scaleFactorGenBin_sys = scaleFactorGenBin #dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[ihsig].Integral() #scaleFactorGenBin

                            #elif len(s)==1 and ('jes' in sys or 'jer' in sys):
                            #    scaleFactor_sys = scaleFactor if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral()
                            #    scaleFactorGenBin_sys = scaleFactorGenBin  if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral()

                            #elif len(s)==1 and ('Weight' in sys):
                            #    scaleFactor_sys = scaleFactor if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral()
                            #    scaleFactorGenBin_sys = scaleFactorGenBin  if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral()
                            #alreadyScaled = 
                            for ihsig in sysSignalHistos:
                                if sys+upDown in ihsig:
                                    if '_recoJet' in ihsig and not('genBin' in ihsig):
                                        print ("hist_label,sys, SF genBin and recobin, data int., sys int., respectively:\n",
                                               ihsig,
                                               sys+upDown,
                                               scaleFactor_sys,
                                               scaleFactorGenBin_sys,
                                               dataHistos['data_reco'+ivar+'_nom'+sel].Integral(),
                                               sysSignalHistos[ihsig].Integral(), 
                                               sysSignalHistos[ihsig].Integral()/nomIntegral, '\n'
                                              )

                                    if ihsig.endswith(sel):
                                        #print (ihsig)
                                        sysSignalHistos[ihsig].Scale( scaleFactor_sys )
                                        scalingDict[f'scaling_{ihsig}'] = scaleFactor_sys

                                    elif ihsig.endswith('genBin'): 
                                        sysSignalHistos[ihsig].Scale( scaleFactorGenBin_sys )
                                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin_sys
                        if not(sys in already_scaled):already_scaled.append(sys)
                        print("Scaling done for foll. sysUncs:", already_scaled)
                                
            if not(process.startswith('MCSelfClosure')):
                print("Rescaling alternate signal MC")
                altscaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                altscaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )

                if 'cross' in process.lower(): 
                    altscaleFactor = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral() / (altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    altscaleFactorGenBin = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )

                print ("Alt MC SF genBin and recobin, respectively:",altSignalLabel,altscaleFactor,altscaleFactorGenBin)

                for ihsig in altSignalHistos:
                    if ihsig.endswith(sel):
                        altSignalHistos[ihsig].Scale( altscaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = altscaleFactor
                    elif 'genBin' in ihsig: 
                        altSignalHistos[ihsig].Scale( altscaleFactorGenBin )#.endswith('genBin')
                        scalingDict[f'scaling_{ihsig}'] = altscaleFactorGenBin
                        
                if extraMC:
                    
                    print("Rescaling alternate signal MC 1")
                    alt1scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    alt1scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( alt1SignalHistos[alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                    
                    print ("Alt MC 1 SF genBin and recobin, respectively:",alt1SignalLabel,alt1scaleFactor,alt1scaleFactorGenBin)

                    
                    for ihsig in alt1SignalHistos:
                        if ihsig.endswith(sel):
                            alt1SignalHistos[ihsig].Scale( alt1scaleFactor )
                            scalingDict[f'scaling_{ihsig}'] = alt1scaleFactor
                        elif 'genBin' in ihsig: 
                            alt1SignalHistos[ihsig].Scale( alt1scaleFactorGenBin )#.endswith('genBin')
                            scalingDict[f'scaling_{ihsig}'] = alt1scaleFactorGenBin
                            
                            
                    print("Rescaling alternate signal MC 2")
                    alt2scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    alt2scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( alt2SignalHistos[alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                    
                    print ("Alt MC 2 SF genBin and recobin, respectively:",alt2SignalLabel,alt2scaleFactor,alt2scaleFactorGenBin)

                    
                    for ihsig in alt2SignalHistos:
                        if ihsig.endswith(sel):
                            alt2SignalHistos[ihsig].Scale( alt2scaleFactor )
                            scalingDict[f'scaling_{ihsig}'] = alt2scaleFactor
                        elif 'genBin' in ihsig: 
                            alt2SignalHistos[ihsig].Scale( alt2scaleFactorGenBin )#.endswith('genBin')
                            scalingDict[f'scaling_{ihsig}'] = alt2scaleFactorGenBin
                    
                    

        ################################################################################################################## 

        ######## Cross check: plotting data vs all MC (scaled to data for QCD, other normalised as per usual by lumi and xs and genweights)
        print ('|------> Cross check: plotting data vs all MC')
        #print(dataHistostrue,'data_reco'+ivar+'_nom'+sel)
        plotSimpleComparison(dataHistos['data_reco'+ivar+'_nom'+sel].Clone(), 
                             'data', allHistos[ 'allMCHisto' ].Clone(), 'allMC', 
                             ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_nom", 
                             rebinX=1, version=sel+'_'+version, outputDir=outputDir )


        
        #################################### Make diagnostics/misc. plots prior to unfolding #################################

        print ('|------> Pre-unfolding cross-check plots '+ivar)


        ####### Cross check response matrix
        tmpGenHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionX()
        #if 'dijet' in sel: tmpGenHisto.Scale( scaleFactorGenBin )

        plotSimpleComparison( tmpGenHisto, 'projection', signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(), 'Regular Gen', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionGen", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        tmpRecoHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
        #if 'dijet' in sel: tmpRecoHisto.Scale( scaleFactor )

        plotSimpleComparison( tmpRecoHisto, 'projection', signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Clone(), 'Regular TrueReco', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionReco", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        ####### Plotting bkg subtracted reco vs. data and including fakes in bkgHistos

        
            

        tmpHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
        plotSimpleComparison( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone(), 'data', allHistos[ 'allBkgHisto' ].Clone(), 
                             'Bkg+fakes', ivar+'_from'+('Data' if process.startswith('data') else 'MC') + '_' + signalLabel + "_TestDataBkgFakes",
                             rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal true reco', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )
        
                                            
        
        getAndPlotPurity(signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone().RebinY(2),
                         reco=signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone().Rebin(2),
                         gen=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                         gen_bins=genBin,variables=variables,var=ivar, outputDir=outputDir,year=year,sel=sel)
        
        ######## Cross check: plotting response matrix
        print ('|------> Cross check: plotting response matrix for signal')
        ROOT.gStyle.SetPadRightMargin(0.15)
        #ROOT.gStyle.SetPalette(ROOT.kGistEarth)
        #ROOT.TColor.InvertPalette()
        can2D = ROOT.TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].GetXaxis().SetTitle('Gen '+variables[ivar]['label'])
        signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].GetYaxis().SetTitle('Reco '+variables[ivar]['label'])
        signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].GetYaxis().SetTitleOffset( 0.8 )
        signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
        
        CMS_lumi.relPosX = 0.12
        CMS_lumi.CMS_lumi(can2D, 4, 0)
        can2D.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+sel+'_responseMatrix'+version+'.'+ext)

        

        ######## TUnfold part
        print ('|------> TUnfolding starts:')

        ##### Defining options for TUnfold
        tunfolder = ROOT.TUnfoldDensity(
                                            signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel], ### response matrix. According to TUnfold, this distribution does NOT have to be normalized
                                            ROOT.TUnfold.kHistMapOutputHoriz,  #### kHistMapOutputVert if x->reco and y->gen, kHistMapOutputHoriz if x->gen and y->reco
                                            ROOT.TUnfold.kRegModeNone, #Curvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                            ROOT.TUnfold.kEConstraintNone if areaConstraint==None else ROOT.TUnfold.kEConstraintArea,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. 
                                            #ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                            )

        ##### Defining input (data recoJet )
        print ('|------> TUnfolding adding input:')

        tunfolder.SetInput( allHistos[ 'dataHisto' ])
        
        if process.startswith('MCCrossClosure'): 
            tunfolder_cross = ROOT.TUnfoldDensity(
                                            altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel], 
                                            ROOT.TUnfold.kHistMapOutputHoriz,  
                                            ROOT.TUnfold.kRegModeNone,#Curvature,   
                                            ROOT.TUnfold.kEConstraintNone if areaConstraint==None else ROOT.TUnfold.kEConstraintArea,    
                                            #ROOT.TUnfoldDensity.kDensityModeBinWidth  
                                            )

            ##### Defining input (data recoJet )
            print ('|------> TUnfolding adding input:')
            tunfolder_cross.SetInput( allHistos[ 'dataHisto' ])
        
        if process.startswith('data'):
            print ("Subtracting backgrounds")
            dummy=0
            bkgSources = []
            #if sel.startswith(('_W','_top')):
            dataMinusbkg_counter = allHistos[ 'dataHisto' ].Integral()
            for ibkg in bkgHistos:
                if ibkg.endswith('_reco'+ivar+'_nom'+sel):
                    #if verbose: 
                    print (ibkg,ibkg.split('_')[0]+ '%d'%dummy,
                           bkgHistos[ibkg].Integral(), dataMinusbkg_counter)

                    dataMinusbkg_counter-=bkgHistos[ibkg].Integral()
                    tunfolder.SubtractBackground( bkgHistos[ibkg].Clone(), ibkg.split('_')[0]+ '%d'%dummy )
                    dummy=dummy+1
                    bkgSources.append(ibkg.split('_')[0]+ '%d'%dummy)
                #if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone() )
            dataMinusbkg_counter-=signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral()
            print('Subtracted fakes', signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral(), dataMinusbkg_counter, signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),allHistos[ 'dataMinusBkgs' ].Integral())

            #if verbose: print('subtracting fakes')
            tunfolder.SubtractBackground( signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel], 'fakes')
            
            bkgSources.append('fakes')
        
        jes_uncorr_list = [
                            '_jesAbsolute_2016', '_jesBBEC1_2016', '_jesEC2_2016', '_jesHF_2016', '_jesRelativeSample_2016',
                            '_jesAbsolute_2017', '_jesBBEC1_2017', '_jesEC2_2017', '_jesHF_2017', '_jesRelativeSample_2017',
                            '_jesAbsolute_2018', '_jesBBEC1_2018', '_jesEC2_2018', '_jesHF_2018', '_jesRelativeSample_2018'
                          ]
        jes_corr_list = ['_jesAbsolute', '_jesBBEC1', '_jesEC2', '_jesFlavorQCD', '_jesHF', '_jesRelativeBal']
        jes_sources_names_uncorr = []
        jes_sources_names_corr = []
        ###### Adding SYS unc
        
        if len(sysUncert)>0 and process.startswith('data'):
            print ('|------> TUnfolding adding uncert:')
            dictUncHistos = {}
            #print (sysSignalHistos)
            for sys in sysUncert:
                #print (sys)
                
                if sys.startswith(('_jer', '_isrWeight', '_l1prefiringWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_jes', '_leptonWeight', '_btagWeight')):
                    
                    s = [i for i in sysSignalLabels if sys in i]
                    #if verbose: 
                    #print (s)
                    if ('2016' in sys or '2017' in sys or '2018' in sys) and len(s)>1:
                        if '2016' in s and 'VFP' in year: s=[s[1]]
                        elif '2016' in s and year.endswith('2016'): s=[s[1]]
                        elif '2017' in s: s=[s[2]]
                        elif '2018' in s: s=[s[3]]
                    else:
                        s=[s[0]]
                    print(f'|------> TUnfolding adding {sys} unc. with prefix-checked {s}')
                    #if verbose: 
                    #    print(sysSignalHistos)
                    dictUncHistos[sys+'Up'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Up'+sel].Clone()
                    dictUncHistos[sys+'Down'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Down'+sel].Clone()
                    for upDown in [ 'Up', 'Down' ]:
                        if verbose: 
                            print(f'|------> TUnfolding adding {sys+upDown}')
                        tunfolder.AddSysError(
                                            sysSignalHistos[s[0]+'_respWithMiss'+ivar+sys+upDown+sel],
                                            sys+upDown,
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                            #### 
                                            # kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. 
                                            # kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. 
                                            # kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+sys+upDown, ivar+'can2DNorm'+sys+upDown, 750, 500 )
                        sysSignalHistos[s[0]+'_respWithMiss'+ivar+sys+upDown+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+s[0]+sel+upDown+'Normalized_responseMatrix'+version+'.'+ext)

                #### adding model uncertainty
                elif sys.startswith(('_model')):
                    if verbose:
                        print('|------> TUnfolding adding modelUnc')
                    tunfolder.AddSysError(
                                        altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel],
                                        'modelUncTotal',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix,
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNormAltSignal', ivar+'can2DNormAltSignal', 750, 500 )
                    altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+altSignalLabel+sel+'Normalized_alt_responseMatrix'+version+'.'+ext)
                    dictUncHistos[sys] = altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel].Clone()
            
                #below ifs relevant only to W/top seln. 
                elif sys.startswith('_hdamp'): 
                    if verbose:
                        print('|------> TUnfolding adding hdampUnc')
                    
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_hdampUp_TuneCP5'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_hdampUp',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix,
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_hdampUp', ivar+'can2DNorm'+'_hdampUp', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_hdampUp_TuneCP5'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampUp'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_hdampDown_TuneCP5'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_hdampDown',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_hdampDown', ivar+'can2DNorm'+'_hdampDown', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_hdampDown_TuneCP5'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampDown'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_hdampUp'] = varSignalHistos['varTTToSemileptonic_hdampUp_TuneCP5'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_hdampDown'] = varSignalHistos['varTTToSemileptonic_hdampDown_TuneCP5'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_Tune'):
                    if verbose:
                        print('|------> TUnfolding adding TuneCP5Unc')
                    
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_TuneCP5Up',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_TuneCP5Up', ivar+'can2DNorm'+'_TuneCP5Up', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Up'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_TuneCP5Down',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_TuneCP5Down', ivar+'can2DNorm'+'_TuneCP5Down', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Down'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_TuneCP5Up'] = varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_TuneCP5Down'] = varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_erdON'): #elif sys.startswith('_CR'):#
                    if verbose:
                        print('|------> TUnfolding adding erdONUnc')
                    
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_erdON',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5_erdON', ivar+'can2DNorm'+'TuneCP5_erdON', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5_erdON'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_erdON'] = varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_reco'+ivar+'_nom'+sel].Clone()
                    
                    if verbose:
                        print('|------> TUnfolding adding Colour reconnection Unc')
                    
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_CR1',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR1', ivar+'can2DNorm'+'TuneCP5CR1', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+ivar+'_nom'+sel],
                                        '_CR2',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR2', ivar+'can2DNorm'+'TuneCP5CR2', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_CR1'] = varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_CR2'] = varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_reco'+ivar+'_nom'+sel].Clone()
                    
                elif sys.startswith('_mtop'): 
                    if verbose:
                        print('|------> TUnfolding adding mtopUnc')
                    mass_list = [ '171p5','173p5' ] #'166p5',
                    
                    for m in mass_list:
                        tunfolder.AddSysError(
                                             varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_respWithMiss'+ivar+'_nom'+sel],
                                             '_mtop%s'%m,
                                             ROOT.TUnfold.kHistMapOutputHoriz,
                                             ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_mtop%s_TuneCP5'%m, ivar+'can2DNorm'+'_mtop%s_TuneCP5'%m, 750, 500 )
                        varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_mtop%s_TuneCP5'%m +'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                        dictUncHistos['_mtop%s'%m] = varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_reco'+ivar+'_nom'+sel].Clone()

            #'''
            #!!!!!!!!!!!!!!FIXME!!!!!!!!!!!!!!!
            ### Making unc plot
            if not('all' in year):
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    dictUncHistos,
                                    ivar+'_'+signalLabel+'_JESAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES'
                                    )
            else:
                tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                for i in dictUncHistos.keys():
                    if i.startswith('_jes') and not('2016' in i or '2017' in i or '2018' in i):
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                    elif i.startswith('_jes') and ('2016' in i or '2017' in i or '2018' in i):
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()
                        
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,
                                    ivar+'_'+signalLabel+'_JESUncorrAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES'
                                    )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,
                                    ivar+'_'+signalLabel+'_JESCorrAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES'
                                    )
                
            if 'dijet' in sel:
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    dictUncHistos,
                                    ivar+'_'+signalLabel+'_NoJESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode=''
                                    )
            else:
                tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                for i in dictUncHistos.keys():
                    if i.startswith(('_model', '_jer', '_isrWeight', '_l1prefiringWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_leptonWeight', '_btagWeight')):
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                    elif i.startswith(('_mtop','_CR','_hdamp','_Tune','_erd')):
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()

                     
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,#dictUncHistos,
                                    ivar+'_'+signalLabel+'_NoJESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode=''
                                  )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,#dictUncHistos,
                                    ivar+'_'+signalLabel+'_TheorySys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode=''
                                  )
                del(tempDictUncHistos,tempDictUncHistos2)
                   
            #'''
        ###### Running the unfolding
        print ('|------> TUnfolding doUnfold:')
        tunfolder.DoUnfold(0)

        if process.startswith('MCCrossClosure'): tunfolder_cross.DoUnfold(0) 
            
        ##### Get output of unfolding 
        allHistos [ 'unfoldHisto'+ivar ] = tunfolder.GetOutput("unfoldHisto"+ivar).Clone()
        allHistos [ 'unfoldHisto'+ivar ].Sumw2()

                
        if process.startswith('MCCrossClosure'): 
            allHistos [ 'unfoldHistoCross'+ivar ] = tunfolder_cross.GetOutput("unfoldHistoCross"+ivar).Clone()
            
        unfoldingtot = allHistos [ 'unfoldHisto'+ivar ].Integral()
        print(f"For {ivar} in year (={year}), total unfolded event count = {unfoldingtot}")

        

        ############################ Get Probability matrix ########################################
        allHistos[ 'probaMatrix'+ivar ] = tunfolder.GetProbabilityMatrix('probaMatrix'+ivar).Clone()
        
        
                

        ##########################  Get various covariances  ###########################
        
        uncerUnfoldHisto = OrderedDict() 
        uncerUnfoldSystCov = OrderedDict()
        
        
        '''
        From TUnfold documentation:
        GetEmatrixSysUncorr(): uncorrelated errors on the input matrix histA, taken as the errors provided with the histogram. These are typically statistical errors from finite Monte Carlo samples.
        GetEmatrixSysSource()(GetDeltaSysSource()): correlated shifts of the input matrix histA. These shifts are taken as one-sigma effects when switchig on a given error soure. Several such error sources may be defined
        GetEmatrixSysBackgroundUncorr(): uncorrelated errors on background sources, originating from the errors provided with the background histograms
        GetEmatrixInput(): statistical uncertainty of the input (the measurement)
        
        '''
        print ('|------> TUnfolding: Obtaining various covariance matrices')
        
        allHistos[ 'cov'+ivar ] = tunfolder.GetEmatrixTotal("cov"+ivar, "Total Covariance Matrix")
        
        allHistos[ 'cov_uncorr_data_'+ivar ] = tunfolder.GetEmatrixInput("cov_uncorr_data"+ivar,
                                                                         "CM from Stat. Unc. of Input Distribution")
                
        allHistos[ 'cov_uncorr_'+ivar ] = tunfolder.GetEmatrixSysUncorr("cov_uncorr"+ivar, 
                                                                        "CM from uncorrelated uncertainties")
               
        allHistos[ 'cov_uncorr_bkg_'+ivar ] = tunfolder.GetEmatrixSysBackgroundUncorr('fakes', 
                                                                                      "CM from Uncorrelated Errors of Background Sources")
        
        if process.startswith('data'):
            if sel.startswith(('_W','_top')):
                for ibkg in bkgSources:
                    if 'fakes' not in ibkg: allHistos[ 'cov_uncorr_bkg_'+ivar ].Add(tunfolder.GetEmatrixSysBackgroundUncorr(ibkg, "CM from Uncorrelated Errors of Background Source "+ibkg))
        

        #### cov total = cov_uncorr + cov_uncorr_data + cov_uncorr_bkg
        #### cov stat = cov_uncorr + cov_uncorr_data
        
        

        
        ############### Build correlation matrix for unfolding#################################
        allHistos['correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone()
        allHistos['correlation_matrix_'+ivar].Reset()
        allHistos['correlation_matrix_'+ivar] = correlation_from_covariance(allHistos[ 'cov'+ivar ].Clone(),allHistos['correlation_matrix_'+ivar])
        
                
        ########################################################################################    
        
        # Create total uncertainty and sys uncertainty histos for plots 
        
        # first build up histos of systematic uncertainties from the background subtraction
        uncerUnfoldHisto[ivar+'_BkgTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_BkgTotal')
        uncerUnfoldHisto[ivar+'_BkgTotal'].Reset()

        for ibin in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):
            bkg_tot = np.sqrt(allHistos[ 'cov_uncorr_bkg_'+ivar ].GetBinContent(ibin,ibin))
            uncerUnfoldHisto[ivar+'_BkgTotal'].SetBinContent(ibin, bkg_tot)
            
        # second build up histos of total systematic uncertainties
        uncerUnfoldHisto[ivar+'_SystTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_SystTotal')
        uncerUnfoldHisto[ivar+'_SystTotal'].Reset()        
        
       
        ###### adding covariances from bkg subtraction and RM finite stats to the overall systematics covariance matrix
        allHistos['cov_systTotal'+ivar] = allHistos[ 'cov_uncorr_'+ivar ].Clone('cov_uncorr_bkg_'+ivar+'+cov_uncorr'+ivar)
        #allHistos['cov_systTotal'+ivar].Reset()
        allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_bkg_'+ivar ])      #array2hist(systcovTotal,allHistos['cov_systTotal'+ivar])
        
        tmp = OrderedDict()
        for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
            tmp[i] = 0
            for k in uncerUnfoldHisto:
                if k.endswith('Total') and not k.endswith(('SystTotal')):
                    tmp[i] = tmp[i] + ( uncerUnfoldHisto[k].GetBinContent( i )**2 )
                    #print(i, k, tmp[i], uncerUnfoldHisto[k].GetBinContent( i ), ( uncerUnfoldHisto[k].GetBinContent( i )**2 ))
        
        
        
                
        #allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_'+ivar ])
        #allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_bkg_'+ivar ])
        
        if len(sysUncert)>0: 
            for i,j in tmp.items():
                uncerUnfoldHisto[ivar+'_SystTotal'].SetBinContent( i, ROOT.TMath.Sqrt( j ) )    
                    
        #storing unnormalized unfoldings and uncs here, 
        #will normalize unfoldHisto and unc plot objects in the draw functions    
        
        
        uncerUnfoldHisto[ivar+'_StatTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_StatTotal')
        uncerUnfoldHisto[ivar+'_StatTotal'].Reset()
        uncerUnfoldHisto[ivar+'_TotalUnc'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_TotalUnc')
        uncerUnfoldHisto[ivar+'_TotalUnc'].Reset()
        
        uncerUnfoldHisto[ivar+'_CMErrTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_CMErrTotal')
        uncerUnfoldHisto[ivar+'_CMErrTotal'].Reset()
        uncerUnfoldHisto[ivar+'_CMMCStatErrTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_CMMCStatErrTotal')
        uncerUnfoldHisto[ivar+'_CMMCStatErrTotal'].Reset()
        uncerUnfoldHisto[ivar+'_CMDataStatErrTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_CMDataStatErrTotal')
        uncerUnfoldHisto[ivar+'_CMDataStatErrTotal'].Reset()
        
        
        allHistos[ 'unfoldHistowoUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone()        # No unc
        allHistos[ 'unfoldHistoStatUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+"unfoldHistoStatUnc")     # Unfolding and stat unc
        allHistos[ 'unfoldHistoBkgUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+"unfoldHistoBkgStatUnc")   # Bkg subtraction unc.
        allHistos[ 'unfoldHistoRMUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+"unfoldHistoRMUnc")         # RM stat sys unc.
        allHistos[ 'unfoldHistoSystUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+"unfoldHistoSystUnc")     # Unc. from systematics and variations
        #allHistos[ 'unfoldHistoTotUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone("unfoldHistoTotUnc")       # Total uncertainty

        
        ratioHistos = OrderedDict()
        ratioHistos[ 'StatUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone('StatUnc'+ivar)   
        ratioHistos[ 'StatUnc'+ivar ].Reset()
        ratioHistos[ 'TotalUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone('TotalUnc'+ivar)   
        ratioHistos[ 'TotalUnc'+ivar ].Reset()
        ratioHistos[ 'SystUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone('SystUnc'+ivar)
        ratioHistos[ 'SystUnc'+ivar ].Reset()
        
        #print ("BC of unfHist", "systot+cov in quadrature", "tot from cov", "syst tot", 'stat+unf unc from unfhisto')
        
        for ibin in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):
            
            unc_tot = np.sqrt( allHistos[ 'cov'+ivar ].GetBinContent(ibin,ibin) ) #total error extracted from from total error matrix from TUnfold
            bkg_tot = np.sqrt(allHistos[ 'cov_uncorr_bkg_'+ivar ].GetBinContent(ibin,ibin))
            datastat_tot = np.sqrt( allHistos[ 'cov_uncorr_data_'+ivar ].GetBinContent(ibin,ibin))
            rmstat_tot = np.sqrt( allHistos[ 'cov_uncorr_'+ivar ].GetBinContent(ibin,ibin))
            
            norm_unc = abs(allHistos[ 'unfoldHisto'+ivar ].GetBinContent(ibin))
            if unc_tot<=0.: 
                unc_tot=0.
            if datastat_tot<=0.: 
                datastat_tot=0.
            if rmstat_tot<=0.: 
                rmstat_tot=0.
            if bkg_tot<=0.: 
                bkg_tot=0.
               
            
            allHistos[ 'unfoldHisto'+ivar ].SetBinError(ibin, unc_tot)
            allHistos[ 'unfoldHistoBkgUnc'+ivar ].SetBinError(ibin, bkg_tot )
            allHistos[ 'unfoldHistoRMUnc'+ivar ].SetBinError(ibin, rmstat_tot)
            allHistos[ 'unfoldHistowoUnc'+ivar ].SetBinError(ibin, 0. )        # No unc
        
            uncerUnfoldHisto[ivar+'_TotalUnc'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMErrTotal'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMMCStatErrTotal'].SetBinContent(ibin, rmstat_tot)
            uncerUnfoldHisto[ivar+'_CMDataStatErrTotal'].SetBinContent(ibin, datastat_tot)
            
            ratioHistos[ 'TotalUnc'+ivar ].SetBinContent( ibin, 1. )
            ratioHistos[ 'StatUnc'+ivar ].SetBinContent( ibin, 1. )

            if norm_unc!=0:
                ratioHistos[ 'TotalUnc'+ivar ].SetBinError( ibin, unc_tot)#+syst_tot**2
                ratioHistos[ 'StatUnc'+ivar ].SetBinError( ibin, datastat_tot)
                
        ########################### Get systematic shifts of output ################################
        #storing individual systematics totals (ie, up/down or other variations), background subtraction systematics, 
        #and then the overall systematic unc (individual systs + bkgs)
        if len(sysUncert)>0 and process.startswith('data'):
            print ('|------> TUnfolding uncertainties:')

            
            
            for sys in sysUncert:
                
                sys_cov_up = tunfolder.GetEmatrixSysUncorr("cov_%s_Up"%sys+ivar)
                sys_cov_down = tunfolder.GetEmatrixSysUncorr("cov_%s_Down"%sys+ivar)
                
                if not sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_TuneCP5','_lepton', '_btag')):
                    #if 'jes' in sys and year.startswith('all'): continue
                    
                    tunfolder.GetEmatrixSysSource(sys_cov_up, sys+'Up')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, sys+'Down')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Up'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'] = sys_cov_down.Clone()
                    
                    for upDown in [ 'Up', 'Down' ]:
                        if verbose:
                            print (sys+upDown)
                        uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "+1#sigma" if 'up' in upDown.lower() else "-1#sigma")

                        if uncerUnfoldHisto[ivar+sys+upDown]: uncerUnfoldHisto[ivar+sys+upDown+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+sys+upDown].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                        try: uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                        except ReferenceError: uncerUnfoldHisto.pop( ivar+sys+upDown, None )
                        

                    # Create total uncertainty and sys uncertainty histos for plotting/further calculations.
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone() #total shifts not used but kept in case
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                    for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                        try: yup = abs(  uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs(  uncerUnfoldHisto[ivar+sys+'Down'].GetBinContent(i) )
                        except KeyError: ydn = 0
                                                        
                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                        #dy = 0.5*( (yup + ydn) ) #being conservative and doing an envelope based on max of up/down shift of a systematic
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )
                    
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'+"_shiftHist"] = get_syst_shifted_hist( uncerUnfoldHisto[ivar+sys.upper()+'Total'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                    if verbose:
                            print ("Done with adding %s to tot syst unc"%sys )    

                elif sys.startswith('_model'):
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_modelUncTotal')
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Reset()
                    
                    tmpModelHisto = tunfolder.GetDeltaSysSource('modelUncTotal', "unfoldHisto_"+ivar+"_modelUncTotalshift", "modelUncShift")
                    tunfolder.GetEmatrixSysSource(sys_cov_up, "modelUncTotal")
                    uncerUnfoldSystCov['systcov_'+ivar+'_Physics ModelTotal'] = sys_cov_up.Clone()
                    
                    #for i in range( 1, tmpModelHisto.GetNbinsX() + 1):
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = tmpModelHisto.Clone()#.SetBinContent( i, (tmpModelHisto.GetBinContent(i)))#abs(
                    
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                
                
                
                elif sys.startswith('_btag') and 'Corr' in sys:#corr?
                    btagUncs = [s for s in sysUncert if 'btag' in s]
                    for s in  btagUncs:
                        for upDown in [ 'Up', 'Down' ]:
                            s_upDown=f'{s}{upDown}'
                            sys_cov = tunfolder.GetEmatrixSysUncorr(f"cov_{sys}_{upDown}"+ivar)

                            tunfolder.GetEmatrixSysSource(sys_cov, f'{s_upDown}')

                            uncerUnfoldSystCov['systcov_'+ivar+f'{s_upDown}'] = sys_cov.Clone()
                            
                            uncerUnfoldHisto[ivar+f'{s_upDown}'] = tunfolder.GetDeltaSysSource(f'{s_upDown}', "unfoldHisto_" + ivar + f'{s_upDown}shift', "+1#sigma" if 'up' in (f'{s_upDown}').lower() else "-1#sigma")
                            
                            if uncerUnfoldHisto[ivar+f'{s_upDown}']: 
                                
                                uncerUnfoldHisto[ivar+f'{s_upDown}'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+f'{s_upDown}'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                    uncerUnfoldHisto[ivar+'_btagWeightTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_btagWeightCorrelatedUp')
                    uncerUnfoldHisto[ivar+'_btagWeightTotal'].Reset()

                    for i in range( 1, uncerUnfoldHisto[ivar+'_btagWeightTotal'].GetNbinsX() + 1):
                        shift_list = []
                        dy_temp=0
                        for s in [u for u in btagUncs if 'Corr' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_temp=0
                        for s in [u for u in btagUncs if 'Uncorr' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_temp=0
                        for s in [u for u in btagUncs if 'Eff' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                      
                        
                        dy_max = np.max(shift_list)
                        
                        uncerUnfoldHisto[ivar+'_btagWeightTotal'].SetBinContent( i, dy_max )
                                           
                    uncerUnfoldHisto[ivar+'_btagWeightTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_btagWeightTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                
                elif sys.startswith('_lepton') and 'ISOStat' in sys:
                    leptUncs = [s for s in sysUncert if 'lepton' in s]
                    for s in  leptUncs:
                        for upDown in [ 'Up', 'Down' ]:
                            s_upDown=f'{s}{upDown}'
                            sys_cov = tunfolder.GetEmatrixSysUncorr(f"cov_{s_upDown}"+ivar)

                            tunfolder.GetEmatrixSysSource(sys_cov, f'{s_upDown}')

                            uncerUnfoldSystCov['systcov_'+ivar+f'{s_upDown}'] = sys_cov.Clone()
                            
                            uncerUnfoldHisto[ivar+f'{s_upDown}'] = tunfolder.GetDeltaSysSource(f'{s_upDown}', "unfoldHisto_" + ivar + f'{s_upDown}shift', "+1#sigma" if 'up' in (f'{s_upDown}').lower() else "-1#sigma")
                            if uncerUnfoldHisto[ivar+f'{s_upDown}']: uncerUnfoldHisto[ivar+f'{s_upDown}'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+f'{s_upDown}'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                    uncerUnfoldHisto[ivar+'_leptonWeightTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_leptonWeightISOStatUp')
                    uncerUnfoldHisto[ivar+'_leptonWeightTotal'].Reset()

                    for i in range( 1, uncerUnfoldHisto[ivar+'_leptonWeightTotal'].GetNbinsX() + 1):
                           
                        shift_list = []
                        dy_temp=0
                        for s in [u for u in leptUncs if 'ISO' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_temp=0
                        for s in [u for u in leptUncs if 'RecoEff' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_temp=0
                        for s in [u for u in leptUncs if 'ID' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_temp=0
                        for s in [u for u in leptUncs if 'Trig' in u]:
                            try: yup = abs(uncerUnfoldHisto[ivar+f'{s}Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+f'{s}Down'].GetBinContent(i) )
                            except KeyError: ydn = 0

                            dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                            dy_temp = dy if dy>dy_temp else dy_temp
                        shift_list.append(dy_temp)
                        
                        dy_max = np.max(shift_list)
                        
                        uncerUnfoldHisto[ivar+'_leptonWeightTotal'].SetBinContent( i, dy_max )
                    
                    uncerUnfoldHisto[ivar+'_leptonWeightTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_leptonWeightTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                
                elif sys.startswith('_Tune'):
                    tunfolder.GetEmatrixSysSource(sys_cov_up, '_TuneCP5Up')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, '_TuneCP5Down')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+'_TuneCP5Up'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+'_TuneCP5Down'] = sys_cov_down.Clone()
                    for s in ['_TuneCP5Up', '_TuneCP5Down' ]:
                        uncerUnfoldHisto[ivar+s] = tunfolder.GetDeltaSysSource(s, "unfoldHisto_" + ivar + "%sshift"%s, "+1#sigma" if 'up' in s.lower() else "-1#sigma")
                        if uncerUnfoldHisto[ivar+s]: uncerUnfoldHisto[ivar+s+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+s].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                    uncerUnfoldHisto[ivar+'_UE tuneTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_TuneCP5Up')
                    uncerUnfoldHisto[ivar+'_UE tuneTotal'].Reset()

                    for i in range( 1, uncerUnfoldHisto[ivar+'_UE tuneTotal'].GetNbinsX() + 1):
                        try: yup = abs(uncerUnfoldHisto[ivar+'_TuneCP5Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs( uncerUnfoldHisto[ivar+'_TuneCP5Down'].GetBinContent(i) )
                        except KeyError: ydn = 0
                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)#dy = np.max([yup,ydn])
                        
                        uncerUnfoldHisto[ivar+'_UE tuneTotal'].SetBinContent( i, dy )
                    
                    uncerUnfoldHisto[ivar+'_UE tuneTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_UE tuneTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                elif sys.startswith('_hdamp'):
                    tunfolder.GetEmatrixSysSource(sys_cov_up, '_hdampUp')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, '_hdampDown')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+'_hdampUp'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+'_hdampDown'] = sys_cov_down.Clone()
                    for s in ['_hdampUp', '_hdampDown' ]:
                        uncerUnfoldHisto[ivar+s] = tunfolder.GetDeltaSysSource(s, "unfoldHisto_" + ivar + "%sshift"%s, "+1#sigma" if 'up' in s.lower() else "-1#sigma")
                        if uncerUnfoldHisto[ivar+s]: uncerUnfoldHisto[ivar+s+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+s].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                    uncerUnfoldHisto[ivar+'_hdampTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_hdampUP')
                    uncerUnfoldHisto[ivar+'_hdampTotal'].Reset()

                    for i in range( 1, uncerUnfoldHisto[ivar+'_hdampTotal'].GetNbinsX() + 1):
                        try: yup = abs(uncerUnfoldHisto[ivar+'_hdampUp'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs(uncerUnfoldHisto[ivar+'_hdampDown'].GetBinContent(i) )
                        except KeyError: ydn = 0
                        #dy = np.max([yup,ydn])# 0.5*(yup + ydn) 

                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                        uncerUnfoldHisto[ivar+'_hdampTotal'].SetBinContent( i, dy )
                    uncerUnfoldHisto[ivar+'_hdampTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_hdampTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())


                elif sys.startswith('_erdON'):#'_CR',
                    tunfolder.GetEmatrixSysSource(sys_cov_up, '_CR1')
                    tunfolder.GetEmatrixSysSource(sys_cov_down,'_CR2')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+'_CR1'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+'_CR2'] = sys_cov_down.Clone()
                    
                    
                    tunfolder.GetEmatrixSysSource(sys_cov_down, '_erdON')
                    uncerUnfoldSystCov['systcov_'+ivar+'_erdON'] = sys_cov_down.Clone()
                    
                    for CR in ['_CR1', '_CR2', '_erdON' ]:
                        uncerUnfoldHisto[ivar+CR] = tunfolder.GetDeltaSysSource(CR, "unfoldHisto_" + ivar + "%sshift"%CR, f"{CR}UncShift")
                        if uncerUnfoldHisto[ivar+CR]: uncerUnfoldHisto[ivar+CR+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+CR].Clone(),
                                                                                       allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                    uncerUnfoldHisto[ivar+'_CRTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_CR1')
                    uncerUnfoldHisto[ivar+'_CRTotal'].Reset()

                    for i in range( 1, uncerUnfoldHisto[ivar+'_CRTotal'].GetNbinsX() + 1):
                        #unfHistoBC = allHistos['unfoldHisto'+ivar].GetBinContent(i)
                        yup = (abs(uncerUnfoldHisto[ivar+'_CR1'].GetBinContent(i) )) #mnemonics for fun
                        ydn = (abs(uncerUnfoldHisto[ivar+'_CR2'].GetBinContent(i) ))
                        ymid = (abs(uncerUnfoldHisto[ivar+'_erdON'].GetBinContent(i) ))
                        dy = np.max([yup,ydn,ymid]) #considering CR unc as an envelope from max contribution per bin from CR variation sources
                        uncerUnfoldHisto[ivar+'_CRTotal'].SetBinContent( i, dy ) # considering uncertainty as an envelope
                    uncerUnfoldHisto[ivar+'_CRTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_CRTotal'].Clone(),
                                                                                                    allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                elif sys.startswith('_mtop'):
                    tunfolder.GetEmatrixSysSource(sys_cov_up, '_mtop173p5')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, '_mtop171p5')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+'_mtop171p5'] = sys_cov_down.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+'_mtop173p5'] = sys_cov_up.Clone()
                    
                    for m,l in zip(mass_list,['Down','Up']):
                        print ("mtop:", m,l)
                        uncerUnfoldHisto[ivar+'_mtop%s'%m] = tunfolder.GetDeltaSysSource('_mtop%s'%m, "unfoldHisto_"+ivar+"_mtop%sshift"%m, "+1#sigma" if '173' in m.lower() else "-1#sigma")
                        if uncerUnfoldHisto[ivar+'_mtop%s'%m]: uncerUnfoldHisto[ivar+'_mtop%s'%l+'_shiftHist'] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_mtop%s'%m].Clone(),
                                                                                                  allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                    uncerUnfoldHisto[ivar+'_Top MassTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_mtop171p5')
                    uncerUnfoldHisto[ivar+'_Top MassTotal'].Reset()

                    for i in range( 0, uncerUnfoldHisto[ivar+'_Top MassTotal'].GetNbinsX() + 1):
                        yup = (abs( uncerUnfoldHisto[ivar+'_mtop173p5'].GetBinContent(i) ))
                        ydn = (abs( uncerUnfoldHisto[ivar+'_mtop171p5'].GetBinContent(i) ))
                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                        uncerUnfoldHisto[ivar+'_Top MassTotal'].SetBinContent( i, dy ) # considering uncertainty as an envelope
                    uncerUnfoldHisto[ivar+'_Top MassTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_Top MassTotal'].Clone(),
                                                                                                    allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
        ################################################################################

        # Get folded distribution from unfolded distribution
        allHistos [ 'foldHisto'+ivar ] = tunfolder.GetFoldedOutput("folded"+ivar).Clone() 
        #allHistos [ 'foldHisto'+ivar ] = get_folded_unfolded(
        #                                                
        #                                                folded=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
        #                                                unfolded=allHistos['unfoldHisto'+ivar].Clone(), 
        #                                                cov_tot=allHistos['cov'+ivar].Clone(), 
        #                                                probaM=allHistos[ 'probaMatrix'+ivar] 
        #                                                           
        #                                                            ).Clone('foldHistoByHand'+ivar )
        
        #if verbose:
        print("native folded output integral",", true-reco Integral", ", custom folded output integral", ", gen-level (miss+accep) integral", ", unfolded (data-bkg) integral")

        print(tunfolder.GetFoldedOutput("folded"+ivar).Integral(), 
              signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),
              allHistos[ 'dataMinusBkgs' ].Integral(),                  
              allHistos [ 'foldHisto'+ivar ].Integral(),
              signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Integral(), 
              allHistos['unfoldHisto'+ivar].Integral())

        if process.startswith('data'): 
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs',  allHistos [ 'foldHisto'+ivar ].Clone(), 'folded', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            
            plotSimpleComparison( allHistos[ 'unfoldHisto'+ ivar ].Clone(), 'unfold',
                                  signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(), 'gen', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            
            
        
        print ('|------> Drawing unfold plot:')
        if not 'Closure' in process:
            drawUnfold(ivar=ivar, 
                       selection=sel, year=year,lumi=lumi, process=process,
                       dataJetHisto=allHistos[ 'dataMinusBkgs' ].Clone(),#allHistos[ 'dataHistoGenBin' ].Clone(),
                       genJetHisto=signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                       unfoldHistoStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(),
                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                       foldHisto=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                       recoJetHisto=signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY().Clone(),#signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Clone(),
                       cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),#ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                       cov_tot=allHistos['cov'+ivar].Clone(),#ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                       #ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                       altMCHisto =  altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       labelX=variables[ivar]['label'],
                       maxX=variables[ivar]['bins'][-1],
                       tlegendAlignment=variables[ivar]['alignLeg'],
                       
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_'+version+'.'+ext,
                       altMC1Histo = alt1SignalHistos[alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Clone(), #if 'dijet' in sel else None, 
                       altMC2Histo = alt2SignalHistos[alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Clone() if 'dijet' in sel else None, 
                       altMC1Histo_label = alt1SignalLabel, #if 'dijet' in sel else None, 
                       altMC2Histo_label = alt2SignalLabel if 'dijet' in sel else None,
                       extraMC=extraMC,# if 'dijet' in sel else False,
                       includeFSR = True if include_FSR_in_unfolded_result else False,
                       fsrUpHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Clone() if include_FSR_in_unfolded_result else None,#.Clone('_fsrWeightUp')
                       fsrDownHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Clone() if include_FSR_in_unfolded_result else None,#.Clone('_fsrWeightDown')
                       )
            
            drawUnfold_unNormalised(ivar=ivar, 
                       selection=sel, year=year,lumi=lumi, process=process,
                       dataJetHisto=allHistos[ 'dataMinusBkgs' ].Clone(),#allHistos[ 'dataHistoGenBin' ].Clone(),
                       genJetHisto=signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                       unfoldHistoStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(),
                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                       foldHisto=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                       recoJetHisto=signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY().Clone(),#signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Clone(),
                       cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),#ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                       cov_tot=allHistos['cov'+ivar].Clone(),#ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                       #ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                       altMCHisto =  altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       labelX=variables[ivar]['label'],
                       maxX=variables[ivar]['bins'][-1],
                       tlegendAlignment=variables[ivar]['alignLeg'],
                       
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_NO_NORM'+version+'.'+ext,
                       altMC1Histo = alt1SignalHistos[alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Clone(), #if 'dijet' in sel else None, 
                       altMC2Histo = alt2SignalHistos[alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Clone() if 'dijet' in sel else None, 
                       altMC1Histo_label = alt1SignalLabel, #if 'dijet' in sel else None, 
                       altMC2Histo_label = alt2SignalLabel if 'dijet' in sel else None,
                       extraMC=extraMC,# if 'dijet' in sel else False,
                       includeFSR = True if include_FSR_in_unfolded_result else False,
                       fsrUpHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Clone() if include_FSR_in_unfolded_result else None,#.Clone('_fsrWeightUp')
                       fsrDownHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Clone() if include_FSR_in_unfolded_result else None,#.Clone('_fsrWeightDown')
                       )
        else: 
            if 'Cross' in process:
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHistoCross=allHistos['unfoldHistoCross'+ivar ].Clone(),
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext
                             )
                drawClosures_unnormalised(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHistoCross=allHistos['unfoldHistoCross'+ivar ].Clone(),
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_NO_NORM_'+version+'.'+ext
                             )
                
            else:
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=[], 
                             unfoldHistoCross=[] ,
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext
                             )
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=[], 
                             unfoldHistoCross=[] ,
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_NO_NORM_'+version+'.'+ext
                             )
                
                
        ######### Plotting Uncertainties
        print ('|------> Drawing unfolding uncertainty plot:')
        
        if not 'Closure' in process: 
            
            drawUncertainties_normalizedshifts(ivar=ivar, 
                                               unfoldHistoTotUnc=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                                               unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                                               unfoldHistoStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(), 
                                               unfoldHistoRMStatUnc=allHistos[ 'unfoldHistoRMUnc'+ivar ].Clone(),
                                               unfoldHistoBkgSubUnc=allHistos[ 'unfoldHistoBkgUnc'+ivar ].Clone(), 
                                               uncerUnfoldHisto=uncerUnfoldHisto, 
                                               cov_tot=allHistos['cov'+ivar].Clone(), 
                                               cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                               cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                               cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                               labelX=variables[ivar]['label'], 
                                               tlegendAlignment=variables[ivar]['alignLeg'],
                                               outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_UNC_'+version+'.'+ext,
                                               year=year, 
                                               unftot=unfoldingtot,selection=sel
                                                )

        ######### Plotting 2D matrices of various kinds
        print ('|------> Drawing various 2D matrices:')
        if process.startswith('data'):
            draw2D( ivar,  tunfolder.GetRhoItotal("rhoI"+ivar, "Global correlations").Clone(), variables[ivar], outputLabel='data_rhoI', outputDir=outputDir,selection=sel,version=version,year=year)
            draw2D( ivar,  allHistos[ 'correlation_matrix_'+ivar ].Clone(), variables[ivar], outputLabel='data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov'+ivar].Clone(), variables[ivar], outputLabel='dataTotal_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_'+ivar].Clone(), variables[ivar], outputLabel='uncorrUncRM_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_data_'+ivar].Clone(), variables[ivar], outputLabel='dataInpStats_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_bkg_'+ivar].Clone(), variables[ivar], outputLabel='BkgSubtractionSyst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_systTotal'+ivar].Clone(), variables[ivar], outputLabel='Syst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        
        draw2D( ivar,  allHistos[ 'probaMatrix'+ivar ].Clone(), variables[ivar], outputLabel='data_probaMatrix', outputDir=outputDir, addCorrelation=True, addCondition=True ,selection=sel,version=version,year=year)
        draw2D( ivar, signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone(), variables[ivar], outputLabel='data_respMatrix', outputDir=outputDir, addCorrelation=False ,selection=sel,version=version,year=year)
            
        normed_covs = OrderedDict()
        tot = allHistos[ 'unfoldHisto'+ivar ].Integral()
        for i in allHistos:
            if 'cov' in i: 
                
                if 'cov'+ivar in i:
                    ah = allHistos[ i ].Clone()
                    ah.Reset()
                    ah = correlation_from_covariance( allHistos[i].Clone(),
                                                      allHistos['correlation_matrix_'+ivar])

                    draw2D( ivar,  ah.Clone(), variables[ivar], outputLabel='Un-normed_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)

                
                
                draw2D( ivar, allHistos[i].Clone(), variables[ivar], outputLabel='Un-normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                _, normed_covs['Normed'+i] = GetNormalizedTMatrixandTH2( allHistos[i].Clone(), allHistos[i].GetTitle()+'Normed',
                                                            allHistos['unfoldHisto'+ivar].Clone() )
                draw2D( ivar, normed_covs['Normed'+i].Clone(), variables[ivar], outputLabel='Normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                
                jacobian_covTot = compute_jacobian(allHistos[ 'unfoldHisto'+ivar ].Clone(),tot)
                cov_arr = th2_to_np_arr(allHistos[i].Clone()) 
                transformed_cov = np.dot(jacobian_covTot, np.dot(cov_arr, jacobian_covTot.T))
                transformed_cov_matrix = numpy_to_hist2D(transformed_cov, allHistos[i].Clone())                
                normed_covs['Normed'+i] = transformed_cov_matrix.Clone()
                
                #normed_covs.append(normed_cov)
                draw2D( ivar, normed_covs['Normed'+i].Clone(), variables[ivar], outputLabel='V2Normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                
        ############### Build correlation matrix for unfolding#################################
        
        allHistos['Normed_correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone()
        allHistos['Normed_correlation_matrix_'+ivar].Reset()
        allHistos['Normed_correlation_matrix_'+ivar] = correlation_from_covariance(normed_covs[ 'Normed'+'cov'+ivar ].Clone(),allHistos['correlation_matrix_'+ivar])
        draw2D( ivar,  allHistos['Normed_correlation_matrix_'+ivar].Clone(), variables[ivar], outputLabel='Normed_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)

    
        '''
        draw2D( ivar, normed_cov['Normed'+'cov_uncorr_'+ivar].Clone(), variables[ivar], outputLabel='Normed_uncorrUncRM_covMatrix', outputDir=outputDir, addCorrelation=True)
        draw2D( ivar, normed_cov['Normed'+'cov_uncorr_data_'+ivar].Clone(), variables[ivar], outputLabel='Normed_dataInpStats_covMatrix', outputDir=outputDir, addCorrelation=True)
        draw2D( ivar, normed_cov['Normed'+'cov_uncorr_bkg_'+ivar].Clone(), variables[ivar], outputLabel='Normed_BkgSubtractionSyst_covMatrix', outputDir=outputDir, addCorrelation=True)
        draw2D( ivar, normed_cov['Normed'+'cov_systTotal'+ivar].Clone(), variables[ivar], outputLabel='Normed_Syst_covMatrix', outputDir=outputDir, addCorrelation=True)
        '''
        
        #try:
        #    print (allHistos,signalHistos,dataHistos, altSignalHistos)
        #except UnboundLocalError as exc:
        #    print (allHistos,signalHistos,dataHistos, '____', exc)
        ######### Saving Histograms
        def renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                
                #if 'fold' in isam: print (isam)
                
                ihis = hist.Clone()
                ihis.Sumw2()
                #scale=-1.
                #if 'dijet' in sel and not('all' in year):
                #    for key in scalingDict.keys():
                #        if key.split('scaling_')[1]==isam:
                #            scale = scalingDict[key]
                #            print(isam,key)
                #            break
                #    if not(scale==-1.): ihis.Scale(1./scale)
                #    else: print(isam, ihis, key, "No scaling found in scaling dict, leaving histo with whatever scaling it already has")
                
                ihis.SetName(isam)
                ihis.SetTitle(isam)
                #if 'TT_' in isam: 
                #    print(isam, ihis.GetName(), ihis.Integral())
                ihis.Write()

        outputRootName = outputDir+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root'
        print ('|------> Saving histograms in rootfile: ', outputRootName)
        outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
        
        
        renamingHistos( signalHistos )
        
        if not process.startswith('MC'):
            renamingHistos( sysSignalHistos )
            
        if not('self' in process.lower()):
            renamingHistos( altSignalHistos )
            
        if extraMC and process.startswith('data'):
            #print(alt1SignalHistos.keys())
            renamingHistos( alt1SignalHistos )
            if 'dijet' in sel: renamingHistos( alt2SignalHistos )

        if process.startswith('data') and  sel.startswith(('_W','_top')): 
            renamingHistos( varSignalHistos )
            
        renamingHistos(uncerUnfoldSystCov)
        renamingHistos(normed_covs)
        renamingHistos( dataHistos )
        #print ("adding varsighistos", varSignalHistos.items())
        #renamingHistos( dataHistos )
        renamingHistos( bkgHistos )
        #print(allHistos.keys())
        renamingHistos( allHistos )
        renamingHistos( uncerUnfoldHisto )
        
        tunfolder.Write()
        outputRoot.Close()

        if return_tunfolder_object:
            return tunfolder, allHistos, signalHistos, uncerUnfoldHisto, ratioHistos, uncerUnfoldSystCov
        
        print ('|------> Saving histograms in yodafile: ', outputRootName.replace('.root', '.yoda'))
        histToYoda = [  yoda.root.to_yoda( allHistos [ 'unfoldHisto'+ivar ] ) ]
        yoda.writeYODA( histToYoda, outputRootName.replace('.root', '.yoda') )

##############################################################################################





##############################################################################################
############ Histogram loader and rebinner for inputs to unfolding script below ##############
##############################################################################################

def loadHistograms(samples, var, sel, sysUnc=[],
                   isMC=True, addGenInfo=True, respOnly=False, lumi=1., noResp=False,
                   variables={}, year='2017', process='data', noRebin=False, outputFolder=None ):
    """docstring for loadHistograms"""
    
    import sys

    if not('dijet' in sel):
        sys.path.insert(0,'../../')
        from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict

    else:
        sys.path.insert(0,'../../')
        from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict

    
    
    if sysUnc==[]: SYSUNC = [ '_nom' ] 
    else: SYSUNC = [ s+u for u in ['Up', 'Down'] for s in sysUnc if not s.startswith(('_model', '_hdamp', '_Tune', '_CR', '_erdON', '_mtop')) ]
    flip = False
    tmpSYSUNC={}
    #print("Loading histos", sysUnc,SYSUNC,flip)
    allHistos = {}
    for isam in samples:
        if sysUnc!=[]:
            for i in (sysUnc):
                if i in isam:
                    flip=True
                    tmpSYSUNC = [i+u for u in ['Up','Down']]
                    continue
        
        if not flip: tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC]
        else: tmpList = ['reco'+var+syst+sel for syst in tmpSYSUNC]
        
        
        if isMC and addGenInfo and not flip:
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in SYSUNC if 'nom' in syst or 'fsr' in syst or 'isr' in syst or 'pdf' in syst ] 
            tmpList = tmpList + [ 'reco'+var+syst+sel for syst in SYSUNC if not(('reco'+var+syst+sel) in tmpList)]
            if not(noResp):
                tmpList = tmpList + [ 'accepgen'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'truereco'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'fakereco'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'missgen'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'respWithMiss'+var+syst+sel for syst in SYSUNC]
        
        elif isMC and addGenInfo and flip: 
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in tmpSYSUNC if 'nom' in syst or 'fsr' in syst or 'isr' in syst or 'pdf' in syst ]
            tmpList = tmpList + [ 'reco'+var+syst+sel for syst in tmpSYSUNC if not(('reco'+var+syst+sel) in tmpList) ]
            if not(noResp): 
                tmpList = tmpList + ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC]
        
        if respOnly and not flip: 
            if not(noResp): 
                tmpList = [ 'respWithMiss'+var+syst+sel for syst in SYSUNC] 
        elif respOnly and flip: 
            if not(noResp): 
                tmpList = [ 'respWithMiss'+var+syst+sel for syst in tmpSYSUNC ]
        
        set_tmpList = list(set(tmpList))
        if not(len(tmpList)==len(set_tmpList)):
            print(len(tmpList),len(set_tmpList))
            tmpList = set_tmpList
            print(len(tmpList),len(set_tmpList))
            
        for ih in tmpList:
            
            if isMC:
                try:
                    iFile = ROOT.TFile.Open(samples[isam][0],'r')
                    #print(isam)
                    allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone() #'jetObservables/'+
                    allHistos[isam+'_'+ih].SetDirectory(0)
                    allHistos[isam+'_'+ih].Sumw2()
                    iFile.Close()
                    MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1]['nGenWeights']
                    #if np.round(MCScale,4)!=np.round(samples[isam][1]['MCScaling'],4): 
                    #    #print (year,'MCScale mismatch alert', MCScale, samples[isam][1]['MCScaling'])
                    #    allHistos[isam+'_'+ih].Sumw2()
                    #    allHistos[isam+'_'+ih].Scale( MCScale )
                    #else: 
                    allHistos[isam+'_'+ih].Sumw2()
                    allHistos[isam+'_'+ih].Scale( MCScale )
                    #if 'QCD' in isam and ih.startswith('reco'): print(isam, MCScale, allHistos[isam+'_'+ih].Integral())
                except ReferenceError:
                    print(isam, samples[isam][0], "FILE NOT FOUND/CORRUPTED; skipping for now")
                
            else:
                if 'dijet' in sel and process.startswith('data'):
                    #add pre/unprescaled trigger histograms to get full spectrum
                    tmpdataHistos = {}
                    itempFile=ROOT.TFile.Open(samples[isam][0],'r')
                    for it in checkDict( 'JetHT', dictSamples )[year]['triggerList']:
                        
                        
                        #print(it,ih,ih.replace( sel, '_'+it+sel ))#,itempFile,samples[isam][0])
                        tmpdataHistos[ it ] = itempFile.Get( ih.replace( sel, '_'+it+sel )).Clone()
                        tmpdataHistos[ it ].Sumw2()
                        tmpdataHistos[it].SetDirectory(0)
                        tmpdataHistos[ it ].Scale( checkDict( 'JetHT', dictSamples )[year]['triggerList'][it] )

                    itempFile.Close() 
                    allHistos[ isam+'_'+ih ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
                    allHistos[ isam+'_'+ih ].Reset()
                    for i in tmpdataHistos: 
                        tmpdataHistos[i].SetDirectory(0)
                        allHistos[isam+'_'+ih].Add( tmpdataHistos[i].Clone() )
                        allHistos[isam+'_'+ih].SetDirectory(0)
                else:
                    iFile = ROOT.TFile.Open(samples[isam][0],'r')
                    allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone()
                    allHistos[isam+'_'+ih].SetDirectory(0)
                    allHistos[isam+'_'+ih].Sumw2()
                    iFile.Close()
    def renamingHistos( dictHistos ):
        for isam, hist in dictHistos.items():
            ihis = hist.Clone()
            ihis.SetName(isam)
            ihis.SetTitle(isam)
            ihis.Write()

        
    #print(allHistos) 
    if isMC:
        
        if 'dijet' in sel.lower():

            c=0
            allHistos_upd={}
            for sys in SYSUNC:

                #print (f"Adding together histograms for {var} in {sel} for {sys}")
                tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k and sys in k}
                #print(tmpHistos.keys(),allHistos.keys())
                for ih in tmpHistos:
                    #if 'fsr' in sys: print(ih)
                    for jh in allHistos:
                        #if 'fsr' in sys: print(jh)
                        goodflag=False
                        if ('2016' in ih and '2016' in jh and '2016' in year) or ('2016' in ih and '2016' in jh and '2016_preVFP' in year) or ('2017' in ih and '2017' in jh and '2017' in year) or ('2018' in ih and '2018' in jh and '2018' in year) and (sys in ih) and (sys in jh):
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral())

                        elif not('201' in ih) and not ('201'in jh) and (sys in ih) and (sys in jh):
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral())

                        if goodflag and (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ) and sys in ih and sys in jh :
                            if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                    #if 'fsr' in sys: print(goodflag,ih,jh,tmpHistos[ih].Integral())

                                    tmpHistos[ih].Add( allHistos[jh].Clone() )
                                    tmpHistos[ih].SetDirectory(0)
                if len(tmpHistos)>0:
                    if c==0:
                        allHistos_upd = copy.deepcopy(tmpHistos)
                        c+=1
                    else:
                        allHistos_upd.update(copy.deepcopy(tmpHistos))
                        c+=1

            allHistos=copy.deepcopy(allHistos_upd)
        
        else:
            
            #c=0
            tmpHistos = { k:v for (k,v) in allHistos.items() if ('Pt-1000' in k) and ('QCD_Pt' in k) and ('MuEnriched' in k)}
            qcdFlag = False if (len(tmpHistos.keys())==0) else True

            #print(qcdFlag)
            if qcdFlag: 
                allHistos_upd=copy.deepcopy(allHistos)

                for k in allHistos.keys():
                    if 'QCD' in k:
                        #print(k)
                        del(allHistos_upd[k])
                        qcdFlag = True

            
                #if len(tmpHistos.keys())>0:
                #    print (f"Adding together QCD histograms for {var} in {sel}",tmpHistos.keys())


                for ih in tmpHistos:
                    for jh in allHistos:

                        if ('MuEnriched' in jh and 'QCD' in jh) and not ('Pt-1000' in jh ):
                            if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                    #if 'recoJet' in ih and not('genBin' in ih): print(ih, jh,tmpHistos[ih].Integral() )
                                    tmpHistos[ih].Add( allHistos[jh].Clone() )
                                    tmpHistos[ih].SetDirectory(0)
                    #if 'recoJet' in ih and not('genBin' in ih): print(ih, tmpHistos[ih].Integral() )
                #print(allHistos_upd.keys())
                
                allHistos_upd.update(copy.deepcopy(tmpHistos))
                allHistos=copy.deepcopy(allHistos_upd)
                
                #print(allHistos.keys())
                
                ##ttbar signal
                nomSysList = [k.split('TTToSemileptonic_')[1] for k in allHistos.keys() if ('SemiLeptonic' in k) and ('TTo' in k) and ('powheg' in k.lower()) and not('var' in k)]
                print("nom,wt, hist list", nomSysList)

                for sNomWt in nomSysList:
                    print("sNomWt",sNomWt)
                    tmpHistos = { k:v for (k,v) in allHistos.items() if ('SemiLeptonic' in k) and ('TTTo' in k) and ('powheg' in k.lower()) and (sNomWt in k) and not('var' in k)}
                    ttSigFlag = False if (len(tmpHistos.keys())==0) else True

                    #print(qcdFlag)
                    if ttSigFlag: 
                        allHistos_upd=copy.deepcopy(allHistos)

                        for k in allHistos.keys():
                            if 'TTTo' in k and not(k.startswith('var')) and (sNomWt in k):
                                #print(k)
                                del(allHistos_upd[k])
                                ttSigFlag = True

                    
                        #if len(tmpHistos.keys())>0:
                        #    print (f"Adding together QCD histograms for {var} in {sel}",tmpHistos.keys())


                        for ih in tmpHistos:
                            for jh in allHistos:

                                if ('TTTo' in jh and not(jh.startswith('var'))) and not ('Semi' in jh ) and (sNomWt in jh):
                                    if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                        if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                            #if 'recoJet' in ih and not('genBin' in ih): print(ih, jh,tmpHistos[ih].Integral() )
                                            tmpHistos[ih].Add( allHistos[jh].Clone() )
                                            tmpHistos[ih].SetDirectory(0)
                            #if 'recoJet' in ih and not('genBin' in ih): print(ih, tmpHistos[ih].Integral() )
                        #print(allHistos_upd.keys())
                        
                        allHistos_upd.update(copy.deepcopy(tmpHistos))
                        allHistos=copy.deepcopy(allHistos_upd)
                

                ##ttbar signal variations
                altSysList = [k.split('varTTToSemileptonic_')[1] for k in allHistos.keys() if ('Semileptonic' in k) and ('varTTTo' in k) and ('powheg' in k.lower())]
                print("altSys, hist list", altSysList)

                for altSys in altSysList:
                    print("altSys",altSys)
                    tmpHistos = { k:v for (k,v) in allHistos.items() if ('Semileptonic' in k) and ('varTTTo' in k) and ('powheg' in k.lower()) and (altSys in k)}
                    ttSigFlag = False if (len(tmpHistos.keys())==0) else True

                    #print(qcdFlag)
                    if ttSigFlag: 
                        allHistos_upd=copy.deepcopy(allHistos)

                        for k in allHistos.keys():
                            if 'varTTo' in k and (altSys in k):
                                #print(k)
                                del(allHistos_upd[k])
                                ttSigFlag = True

                    
                        #if len(tmpHistos.keys())>0:
                        #    print (f"Adding together QCD histograms for {var} in {sel}",tmpHistos.keys())


                        for ih in tmpHistos:
                            for jh in allHistos:

                                if ('varTTTo' in jh) and not('Semi' in jh ) and (altSys in jh):
                                    if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                        if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                            #if 'recoJet' in ih and not('genBin' in ih): print(ih, jh,tmpHistos[ih].Integral() )
                                            tmpHistos[ih].Add( allHistos[jh].Clone() )
                                            tmpHistos[ih].SetDirectory(0)
                            #if 'recoJet' in ih and not('genBin' in ih): print(ih, tmpHistos[ih].Integral() )
                        #print(allHistos_upd.keys())
                        
                        allHistos_upd.update(copy.deepcopy(tmpHistos))
                        allHistos=copy.deepcopy(allHistos_upd)
                    
            else: 
                del(tmpHistos)
                pass
            
        

    if not noRebin:
        #print("About to rebin histos:")#,allHistos.keys())#,tmpHistos.keys())            
        #print("Proceeding to rebin histograms")
        keyList=copy.deepcopy(list(allHistos.keys()))
        for ih in keyList:
            #if 'resp' in ih: print(ih)
            if len(variables[var]['bins'])==1:
                genBin = variables[var]['bins'][0]
                recoBin = variables[var]['bins'][0]/2
            else:
                genBin = variables[var]['bins']
                recoBin = variables[var]['bins_reco']
            if not('respWithMiss' in ih):
                #print(genBin,recoBin)

                if len(variables[var]['bins'])==1:
                    if 'recoJet' in ih:
                        allHistos[ih+'_genBin'] = allHistos[ih].Clone()
                        allHistos[ih+'_genBin'].Rebin( genBin )
                        allHistos[ih].Rebin( recoBin )
                    elif 'genJet' in ih: allHistos[ih].Rebin( genBin )

                else:
                    if 'recoJet' in ih:
                        allHistos[ih+'_genBin'] = allHistos[ih].Clone()
                        allHistos[ih+'_genBin'] = allHistos[ih+'_genBin'].Rebin( len(genBin)-1, allHistos[ih].GetName()+"_Rebin_genBin", array( 'd', genBin ) )
                        allHistos[ih] = allHistos[ih].Rebin( len(recoBin)-1, allHistos[ih].GetName()+"_Rebin", array( 'd', recoBin ) )
                    elif 'genJet' in ih:
                        allHistos[ih] = allHistos[ih].Rebin( len(genBin)-1, allHistos[ih].GetName()+"_Rebin", array( 'd', genBin ) )

            else:
                
                if len(variables[var]['bins'])==1: 
                    allHistos[ih].Rebin2D( genBin,recoBin )
                else:

                    #### fancy way to create variable binning TH2D
                    tmpHisto = ROOT.TH2F( allHistos[ih].GetName()+"_Rebin", allHistos[ih].GetName()+"_Rebin", len(genBin)-1, array( 'd', genBin), len(recoBin)-1, array( 'd', recoBin) )
                    #if tmpHisto.GetNbinsY()>500: print(ih, tmpHisto.GetNbinsY())
                    
                    tmpHisto = rebin_RM_withUF(allHistos[ih],genBin,recoBin).Clone(allHistos[ih].GetName()+'_Rebin')
                    #make_rebinned_2d_hist(allHistos[ih].Clone(), new_bin_edge_pairs,True)
                    tmpHisto.Sumw2()                        

                    tmpHisto.SetDirectory(0)
                    allHistos[ih] = copy.deepcopy(tmpHisto.Clone())
        
        if samples and outputFolder: 
            outputRootName = outputFolder+ f"/{sel.replace('_','')}/{year}" + '/loadedHistograms_main_'+process+isam+var+year+'.root'
            print ('|------> Saving histograms in rootfile: ', outputRootName)
            outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
            renamingHistos( copy.deepcopy(allHistos) )
            #print(allHistos)
            outputRoot.Close()    
    
    return allHistos
##############################################################################################



##############################################################################################
######################## some basic testers and helpers for unfolding ########################
##############################################################################################

def DoUnfolding(Response,Reco):
    tunfolder = ROOT.TUnfoldDensity(fonse,
                                    ROOT.TUnfold.kHistMapOutputHoriz,
                                    ROOT.TUnfold.kRegModeCurvature, 
                                    ROOT.TUnfold.kEConstraintNone, 
                                    ROOT.TUnfoldDensity.kDensityModeBinWidth)
    tunfolder.SetInput(Reco.Clone())
    tunfolder.DoUnfold(0.)
    return tunfolder.GetOutput("MC_unfolded").Clone()

def CrossClosure(response1,reco1,response2,reco2):
    unf11=DoUnfolding(response1.Clone(),reco1.Clone())
    unf12=DoUnfolding(response2.Clone(),reco1.Clone())
    unf21=DoUnfolding(response1.Clone(),reco2.Clone())
    unf22=DoUnfolding(response2.Clone(),reco2.Clone())
    return unf11.Clone(),unf12.Clone(),unf21.Clone(),unf22.Clone()

def SelfClosure(response1,reco2,response2,reco1):
    unf21=DoUnfolding(response1.Clone(),reco2.Clone())
    unf12=DoUnfolding(response2.Clone(),reco1.Clone())
    return unf21.Clone(),unf12.Clone()

def correctEfficiency_Addition(unfolded_hist,miss):
    aTH1=unfolded_hist.Clone()
    aTH1.Reset()
    #print(1./miss.Integral(),miss.Integral())
    #miss.Scale(1./(miss.Integral()))
    
    for i in range(1,unfolded_hist.GetNbinsX()+1):
        bc=unfolded_hist.GetBinContent(i)
        
        if bc<0:
            be = unfolded_hist.GetBinError(i)
            print(f'###### Warning: bin {i} has {bc} bin contents with err={be}')
        
        
        #print(i,bc,miss.GetBinContent(i))
        
        bc=bc+miss.GetBinContent(i)#=multiplicand
        print("#",i,bc,miss.GetBinContent(i))
        
        aTH1.SetBinContent(i, bc)

    aTH1.SetDirectory(0)   
    return aTH1
        
# a la https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L270        
def getMissRate(h_gen,h_missgen):
    h_missRate = h_missgen.Clone()
    genSubtract = np.zeros(h_missRate.GetNbinsX()+1)
    for i in range(0,h_missRate.GetNbinsX()+1):

        genSubtract[i]=(h_gen.GetBinContent(i)-h_missgen.GetBinContent(i))
        h_missRate.SetBinContent(i,h_missgen.GetBinContent(i)/genSubtract[i] if not genSubtract[i]==0 else 0.)

        #print(i,h_gen.GetBinContent(i),
        #      h_missgen.GetBinContent(i),genSubtract[i],
        #      (h_missgen.GetBinContent(i)/genSubtract[i] if not genSubtract[i]==0 else 0.))
    h_missRate.SetDirectory(0)
    
    return h_missRate

def correctEfficiency_Rate(unfolded_hist,miss_hist,gen_hist):
    missrate=getMissRate(gen_hist,miss_hist).Clone()
    aTH1=unfolded_hist.Clone()
    aTH1.Reset()
    #print(miss.Integral()/(miss.Integral()+unfoldhist.Integral())
    #miss.Scale(1./(miss.Integral()))
    
    for i in range(1,unfolded_hist.GetNbinsX()+1):
        bc=unfolded_hist.GetBinContent(i)
        
        if bc<0:
            be = unfolded_hist.GetBinError(i)
            print(f'###### Warning: bin {i} has {bc} bin contents with err={be}')
        
        ## missing gen (reco inefficiency) correction
        scaling = 1. #+ 
        scaling += missrate.GetBinContent(i)
        
        
        #print(i,bc,missrate.GetBinContent(i),scaling)
        
        bc*=scaling

        #print("#",i,bc)
        
        aTH1.SetBinContent(i, bc)
        
        
    aTH1.SetDirectory(0)   
    return aTH1

def correctByMiss_Rate(unfolded_hist,miss_hist,gen_hist):
    missrate=getMissRate(gen_hist,miss_hist).Clone()
    aTH1=unfolded_hist.Clone()
    aTH1.Reset()
    #print(miss.Integral()/(miss.Integral()+unfoldhist.Integral())
    #miss.Scale(1./(miss.Integral()))
    
    for i in range(1,unfolded_hist.GetNbinsX()+1):
        bc=unfolded_hist.GetBinContent(i)
        
        if bc<0:
            be = unfolded_hist.GetBinError(i)
            print(f'###### Warning: bin {i} has {bc} bin contents with err={be}')
        
        ## missing gen (reco inefficiency) correction
        scaling = 1. #+ 
        scaling -= missrate.GetBinContent(i)
        
        
        #print(i,bc,missrate.GetBinContent(i),scaling)
        
        bc*=scaling

        #print("#",i,bc)
        
        aTH1.SetBinContent(i, bc)
        
        
    aTH1.SetDirectory(0)   
    return aTH1

def getAndPlotPurity(h_resp_rebinned,reco,gen,gen_bins,variables,var,sel='_dijetSel',outputDir='../Results/',year='2017'):
    
    rebinned=h_resp_rebinned.Clone()#make_rebinned_2d_hist(h_resp.Clone(),new_bin_edge_pairs,)#rebinning to gen-level bins
    
    arr_rebinned,_ = th2_to_np_arr(rebinned.Clone())
    rebinned_array2d_normX = renorm(arr_rebinned, axis=0) # normalise axis to 1, renormed per x/gen bin
    rebinned_array2d_normY = renorm(arr_rebinned, axis=1) # normalise axis to 1, renormed per y/reco bin
    
    p_list=[]
    s_list=[]
    accepGen = rebinned.ProjectionX('accepGen',1,rebinned.GetNbinsX()+1)
    accepGen.Sumw2()
    accepGen.Divide(gen.Clone())
    trueReco = rebinned.ProjectionY('trueReco')#,1,rebinned.GetNbinsX()+1)
    trueReco.Sumw2()
    fakeReco = reco.Clone()
    fakeReco.Sumw2()
    fakeReco.Add(trueReco,-1.)
    fakeReco.Divide(reco.Clone())
    
    for ibin in range(len(gen_bins)-1):
        #print (f"Calculating p/s per bin in final new binning for bin: {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]}")
        purity = rebinned_array2d_normY[ibin][ibin] #contains fraction in a reco bin that are actually from the same gen bin
        stability = rebinned_array2d_normX[ibin][ibin] #contains fraction in a gen bin that are actually from the same reco bin
        p_list.append(purity)
        s_list.append(stability)
        #print (f"Purity, stability in bin {ibin}({gen_bins[i],gen_bins[i+1]}): {purity,stability}")
    

    if 'dijet' in sel:
        signalLabelBegin = 'QCD_HT_MG5-MLM+P8'
    else:
        signalLabelBegin='TTToSemiLeptonic'
    accepGen.SetDirectory(0)
    fakeReco.SetDirectory(0)
    makePSplot_simple(purity=array('d',p_list),stability=array('d',s_list),
                      accepGen=accepGen,fakeReco=fakeReco,
                      dictHistos=OrderedDict(),
                      variables=variables,
                      var=var,outputDir=outputDir,bins=gen_bins,year=year,
                      sel=sel,
                      signalLabelBegin=signalLabelBegin,
                      ext='pdf'
                     )
    return 1

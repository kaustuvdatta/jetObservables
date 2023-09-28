from collections import OrderedDict
import copy
import pprint 
import ROOT
import numpy as np
import array
from array import array
import bisect
#from legend import *
#ROOT.PyConfig.IgnoreCommandLineOptions = True
#ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
#ROOT.gErrorIgnoreLevel = ROOT.kWarning
####gReset()
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

ROOT.TH1.StatOverflows(ROOT.kTRUE)
ROOT.TH2.StatOverflows(ROOT.kTRUE)


from root_numpy import array2hist, hist2array
#import histoHelpers
from histoHelpers import *
#import unfoldingPlottersAndHelpers
from unfoldingPlottersAndHelpers import *
import os
import glob
import sys
import math
sys.path.insert(0,'../test/')
from DrawHistogram_dijetSel import *
sys.path.insert(0,'../../')
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle

#ROOT.gROOT.ForceStyle()
#tdrstyle.setTDRStyle()

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
                version='_Sept23',verbose=False, return_tunfolder_object = False
              ):
    
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    ROOT.TH1.StatOverflows(ROOT.kTRUE)
    ROOT.TH2.StatOverflows(ROOT.kTRUE)

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 41, 45, 48, 49, 37, 38, 33, 17]
    
    
    signalLabel = 'TTToSemiLeptonic'
    sigPlotLabel = 'Powheg+Pythia8'
    signalLabelBegin = 'TTToSemiLeptonic'
    varSignalLabelBegin = 'varTTToSemileptonic'
    sysSignalLabelBegin = 'sysTTToSemileptonic'

    altSignalLabelBegin = 'TTJets'
    altSigPlotLabel = 'aMC@NLO-FXFX+Pythia8'
    altSignalLabel = 'TTJets'
    
    varSignalLabels = ['varTTToSemileptonic_TuneCP5Up',
                       'varTTToSemileptonic_TuneCP5Down',
                       'varTTToSemileptonic_hdampUp_TuneCP5', 
                       'varTTToSemileptonic_hdampDown_TuneCP5', 
                       'varTTToSemileptonic_TuneCP5_erdON', 
                       'varTTToSemileptonic_TuneCP5CR1', 
                       'varTTToSemileptonic_TuneCP5CR2', 
                       #'varTTToSemileptonic_mtop166p5_TuneCP5', 
                       #'varTTToSemileptonic_mtop171p5_TuneCP5', 
                       #'varTTToSemileptonic_mtop173p5_TuneCP5', 

                      ]
    
    print ("Labels:", signalLabelBegin,altSignalLabelBegin)

    
    selection=sel
    
    
    for ivar in variables:
        if not('tau' in ivar) : continue
        #print (ivar)
        outputDir=outputFolder+sel.split('_')[1]+'/'+year+'/Unfolding/'+ivar+'/'+process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        genBin = variables[ivar]['bins']
        recoBin = variables[ivar]['bins_reco']

    ###########################################################################################
        if year.startswith('all'):
            signalHistos = {
                    signalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel),
                    signalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel),
                    signalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel),
                    signalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel),
                    signalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_gen'+ivar+'_nom'+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin': dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin')
                    }
            signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )
            
            signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )
            
            signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )            
            
            if not process.startswith('MC'):
                sysSignalHistos={}

                print ("Loading up all syst. variations from the following:", sysUncert)
                for sys in sysUncert:
                    print (sys)
                    if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                    for upDown in [ 'Up', 'Down' ]:
    
                        s = [i for i in sysSignalLabels if sys in i]
                        #if verbose: 
                        #print (s)
                        if len(s)>1:
                            if ('2016' in sys or '2017' in sys or '2018' in sys):

                                if '2016' in s: s=[s[1]]
                                elif '2017' in s: s=[s[2]]
                                elif '2018' in s: s=[s[3]]
                            else:
                                s=[s[0]]
                            if verbose: print (s[0]+'_reco'+ivar+sys+upDown+sel)
                        elif len(s)==1:
                            s=[s[0]]
                            if verbose: print (s[0]+'_reco'+ivar+sys+upDown+sel)
                            
                        else:                            
                            continue
                            
                            
                        if not process.startswith('MCCrossClosure') and (not ('2017' in sys) and not ('2018' in sys) and not('2016' in sys)): 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )

                        
                        # dealing with uncorrelated jes unc sources below
                        elif not process.startswith('MCCrossClosure') and '2016' in sys and not '2017' in sys and not '2018' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            
                            #ensuring pre and post VFP periods are treated differently

                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            

                        elif not process.startswith('MCCrossClosure') and '2017' in sys and not '2016' in sys and not '2018' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                         
                        elif not process.startswith('MCCrossClosure') and '2018' in sys and not '2017' in sys and not '2016' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))
                         
                        ################ added recohistos & resp matrices from nominal_2018(/2017) and jesUncorrUnc_2017(/2018) ########################
                        
            if not('self' in process.lower()):
                altSignalHistos = {
                    altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel),
                    altSignalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel),
                    altSignalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel),
                    altSignalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel),
                    altSignalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel),
                    altSignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel),
                    }

                altSignalHistos[ altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel) )

                altSignalHistos[ altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel) )

                altSignalHistos[ altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_accepgen'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_gen'+ivar+'_nom'+sel) )
            
            if process.startswith('data') and selection.startswith(('_W','_top')):
                varSignalHistos={}
                s=[]

                '''
                for sys in sysUncert:
                    if sys.startswith(('_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): 
                        f = [i for i in varSignalLabels if sys in i and i not in s]
                        s.append()
                print (s)
                s = list(set(s))
                print ("sys",sysUncert, s)
                '''
                print ("Processing signal variations from amongst the foll.:", sysUncert)

                for sys in sysUncert:
                    print (sys)
                    s = [i for i in varSignalLabels if (sys.split('_')[1] in i)]
                    for j in s:
                        print (j+'_reco'+ivar+'_nom'+sel)
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_reco'+ivar+'_nom'+sel)
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_reco'+ivar+'_nom'+sel) )
                        
                        #varSignalHistos[ j+'_gen'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_gen'+ivar+'_nom'+sel)
                        #varSignalHistos[ j+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_gen'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_gen'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_gen'+ivar+'_nom'+sel) )
                        
                        #varSignalHistos[ j+'_accepgen'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_accepgen'+ivar+'_nom'+sel)
                        #varSignalHistos[ j+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_accepgen'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_accepgen'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_accepgen'+ivar+'_nom'+sel) )
                        
                        #varSignalHistos[ j+'_fakereco'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_fakereco'+ivar+'_nom'+sel)
                        #varSignalHistos[ j+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_fakereco'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_fakereco'+ivar+'_nom'+sel) )
                        #varSignalHistos[ j+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_fakereco'+ivar+'_nom'+sel) )
                        
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ] = dataFile[ivar+'_2016_preVFP'].Get(j+'_respWithMiss'+ivar+'_nom'+sel)
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )
                        varSignalHistos[ j+'_respWithMiss'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_respWithMiss'+ivar+'_nom'+sel) )

            if process.startswith('MC'):
                allHistos = {
                        'dataHisto' : dataFile[ivar+'_2016_preVFP'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel ),
                        'dataHistoGenBin' : dataFile[ivar+'_2016_preVFP'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' )
                        }
                allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2016'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel ) )
                allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2016'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ))
                allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2017'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel ) )
                allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2017'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ) )
                allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2018'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel ) )
                allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2018'].Get( signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ) )

            else:# process.startswith('data'):
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

            dataHistos = { }
            bkgHistos = { }

        else:
            print('|-------> Running single year '+year)
            ### Getting input histos
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
            
            tmp2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(altSignalLabelBegin)  }
            
            if not process.startswith('MCSelfClosure'): 
                altSignalHistos = loadHistograms( tmp2SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                  lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
            
            if process.startswith('data'): 
                sysSignalHistos = loadHistograms( sysSigFiles, ivar, sel, sysUnc=sysUncert, respOnly=False, isMC=True, lumi=lumi,
                                                  year=year, process=process, variables=variables,outputFolder=outputFolder)

                if selection.startswith(('_W','_top')): 
                    varSignalHistos = loadHistograms( varSigFiles, ivar, sel, sysUnc=[], respOnly=False, isMC=True, lumi=lumi,
                                                  year=year,process=process, variables=variables, outputFolder=outputFolder )
                    #print(varSignalHistos.keys())

            if selection.startswith(('_W','_top')): 

                bkgHistos = loadHistograms( bkgFiles, ivar, sel, sysUnc=[], respOnly=False, isMC=True, lumi=lumi, 
                                            year=year,process=process, variables=variables, outputFolder=outputFolder )
            else: bkgHistos = {}

            if verbose: print ("All signal histos:", signalHistos)

            ####### Fix fake and true reco ###nothing to fix here anymore :) 
            
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone()
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Add( signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel], -1 )

            
            signalHistos[signalLabel+'_missgen'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone()
            signalHistos[signalLabel+'_missgen'+ivar+'_nom'+sel].Add( signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel], -1 )
            
            if not process.startswith('MCSelfClosure'):
                
                altSignalHistos[altSignalLabel+'_fakereco'+ivar+'_nom'+sel] = altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel].Clone()
                altSignalHistos[altSignalLabel+'_fakereco'+ivar+'_nom'+sel].Add( altSignalHistos[altSignalLabel+'_truereco'+ivar+'_nom'+sel], -1 )

                
                altSignalHistos[altSignalLabel+'_missgen'+ivar+'_nom'+sel] = altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone()
                altSignalHistos[altSignalLabel+'_missgen'+ivar+'_nom'+sel].Add( altSignalHistos[altSignalLabel+'_accepgen'+ivar+'_nom'+sel], -1 )    
                
                    
              
            if sel.startswith('_dijet'): 
                if process.startswith("MC"):
                    dataHistostrue = { 'data_reco'+k.split(('_truereco'))[1] : v for (k,v) in signalHistos.items() if ('_truereco' in k and not ('genBin' in k ))}
                    dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistostrue['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                    
                    dataHistos = { 'data_reco'+k.split(('_reco'))[1] : v for (k,v) in signalHistos.items() if ('_reco' in k and not ('genBin' in k ))}
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                    
                else:
                    dataHistostrue = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
                    dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
                    
                ######## Cross check: plotting data vs all MC Scaled
                print ('|------> Cross check: plotting data vs all MC')
                #print(dataHistostrue,'data_reco'+ivar+'_nom'+sel)
                allHistos = {}
                allHistos[ 'allBkgHisto' ] = dataHistostrue['data_reco'+ivar+'_nom'+sel].Clone()
                allHistos[ 'allBkgHisto' ].Reset()
                allHistos[ 'allBkgHistoGenBin' ] = dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
                allHistos[ 'allBkgHistoGenBin' ].Reset()
                if process.startswith('data'):
                    for ibkg in bkgHistos:
                        if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )
                        if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone() )
                allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
                allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Clone() )
                allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
                allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone() )

            
            
            else:
                if process.startswith("MC"):
                    dataHistostrue = { 'data_reco'+k.split(('_truereco'))[1] : v for (k,v) in signalHistos.items() if ('_truereco' in k and not ('genBin' in k ))}
                    dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistostrue['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                    
                    dataHistos = { 'data_reco'+k.split(('_reco'))[1] : v for (k,v) in signalHistos.items() if ('_reco' in k and not ('genBin' in k ))}
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                    
                else:
                    dataHistostrue = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
                    dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
            
                print ('|------> Cross check: plotting data vs all MC')
                allHistos = {}
                allHistos[ 'allBkgHisto' ] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
                allHistos[ 'allBkgHisto' ].Reset()
                allHistos[ 'allBkgHistoGenBin' ] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
                allHistos[ 'allBkgHistoGenBin' ].Reset()
                if process.startswith('data'):
                    for ibkg in bkgHistos:
                        if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )                            
                        if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone().Rebin( len(genBin)-1, 'ibkg_rebin', array( 'd', genBin ) ) )
                            
                allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
                allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Clone() )
                allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
                allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone() )

            ### For dijet, scale QCD to data
            scaleFactor=1.
            scaleFactorGenBin=1.
            altscaleFactor=1.
            altscaleFactorGenBin=1.
            if sel.startswith('_dijet'):
                scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto' ].Integral()
                scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' ].Integral()
                print ("SF genBin and recobin, respectively, nominal reco:",scaleFactor,scaleFactorGenBin)
                allHistos[ 'allMCHisto' ].Scale( scaleFactor )
                allHistos[ 'allMCHistoGenBin' ].Scale( scaleFactorGenBin )
                #if process=='data':
                allHistos[ 'allBkgHisto' ].Scale( scaleFactor )
                allHistos[ 'allBkgHistoGenBin' ].Scale( scaleFactorGenBin )
                
                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        #print (ihsig)
                        signalHistos[ihsig].Scale( scaleFactor )
                    elif ihsig.endswith('genBin'): signalHistos[ihsig].Scale( scaleFactorGenBin )
                    
                if not (process.startswith('MC')):
                    scaleFactor_sys = 1.
                    scaleFactorGenBin_sys = 1.
                    
                    if verbose: 
                        print(sysSignalLabels)
                        print(sysUncert)
                    for sys in sysUncert:
                        if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                        s = [i for i in sysSignalLabels if sys in i]
                        if verbose: 
                            print (sys,s)
                        
                        for upDown in ["Up","Down"]:
                            if len(s)>1:
                                if ('2016' in sys or '2017' in sys or '2018' in sys):

                                    if '2016' in sys: s=[s[1]]
                                    elif '2017' in sys: s=[s[1]]
                                    elif '2018' in sys: s=[s[1]]
                                else:
                                    s=[s[0]]
                                scaleFactor_sys = scaleFactor if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral()
                                scaleFactorGenBin_sys = scaleFactorGenBin  if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral()

                            elif len(s)==1 and ('jes' in sys or 'jer' in sys):
                                scaleFactor_sys = scaleFactor if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral()
                                scaleFactorGenBin_sys = scaleFactorGenBin  if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral()

                            elif len(s)==1 and ('Weight' in sys):
                                scaleFactor_sys = scaleFactor if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral()
                                scaleFactorGenBin_sys = scaleFactorGenBin  if not(sys.startswith(('_jes','_jer'))) else dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() /  sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral()

                            for ihsig in sysSignalHistos:
                                if sys+upDown in ihsig:
                                    #print ("SF genBin and recobin, respectively:",sys+upDown,scaleFactor_sys,scaleFactorGenBin_sys)

                                    if ihsig.endswith(sel):
                                        #print (ihsig)
                                        sysSignalHistos[ihsig].Scale( scaleFactor_sys )
                                    elif ihsig.endswith('genBin'): sysSignalHistos[ihsig].Scale( scaleFactorGenBin_sys )
                
                
                plotSimpleComparison( dataHistos['data_reco'+ivar+'_nom'+sel].Clone(), 
                                     'data', allHistos[ 'allMCHisto' ].Clone(), 'allBkgs', 
                                     ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_nom", 
                                     rebinX=1, version=sel+'_'+version, outputDir=outputDir )
                
                if not(process.startswith('MCSelfClosure')):
                    print("Rescaling alternate signal MC")
                    altscaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    altscaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                    print ("Alt MC SF genBin and recobin, respectively:",altscaleFactor,altscaleFactorGenBin)

                    for ihsig in altSignalHistos:
                        if ihsig.endswith(sel):
                            altSignalHistos[ihsig].Scale( altscaleFactor )
                        elif ihsig.endswith('genBin'): altSignalHistos[ihsig].Scale( altscaleFactorGenBin )
                
                

            print ('|------> Unfolding '+ivar)

            
            ####### Cross check response matrix
            tmpGenHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionX()
            plotSimpleComparison( tmpGenHisto, 'projection', signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel].Clone(), 'Regular AccepGen', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionGen", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )

            tmpRecoHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
            plotSimpleComparison( tmpRecoHisto, 'projection', signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Clone(), 'Regular TrueReco', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionReco", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            tmpRecoHisto.Scale( scaleFactor )
            
            ####### plotting Removing of bkgs from data
            fakeHistos = { k:v for (k,v) in signalHistos.items()  if ('_fakereco' in k and not ('genBin' in k ))}
            for ih in fakeHistos:
                if ih.endswith(ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( fakeHistos[ih] )
                if ih.endswith(ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( fakeHistos[ih] )
            plotSimpleComparison( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone(), 'data', allHistos[ 'allBkgHisto' ].Clone(), 'Bkg+fakes', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataBkgFakes", rebinX=1, 
                                  version=sel+'_'+version, outputDir=outputDir )
            
            if process.startswith('MC'):
                allHistos[ 'dataHisto' ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel ].Clone()
                allHistos[ 'dataHistoGenBin' ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            else:
                allHistos[ 'dataHisto' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone()
                allHistos[ 'dataHistoGenBin' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            
        if process.startswith('data') and not year.startswith('all'): #just for plotting
            allHistos[ 'dataMinusBkgs' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone()
            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            allHistos[ 'dataMinusBkgsGenBin' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )

            tmpHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
            
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal true reco', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
        
        elif process.startswith('data') and year.startswith('all'):
            
            print ('|------> Cross check: plotting data vs all MC')
            allHistos[ 'allBkgHisto' ] = allHistos['dataHisto'].Clone()
            allHistos[ 'allBkgHisto' ].Reset()
            allHistos[ 'allBkgHistoGenBin' ] = allHistos['dataHistoGenBin'].Clone()
            allHistos[ 'allBkgHistoGenBin' ].Reset()
            if process.startswith('data'):
                for ibkg in bkgHistos:
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel): 
                        allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): 
                        allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone().Rebin( len(genBin)-1, 'ibkg_rebin', array( 'd', genBin ) ) )

            allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
            allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Clone() )
            allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
            allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone() )

            allHistos[ 'dataMinusBkgs' ] = allHistos[ 'dataHisto' ].Clone()
            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            allHistos[ 'dataMinusBkgsGenBin' ] = allHistos[ 'dataHistoGenBin'].Clone()
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )

            
            tmpHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal true reco', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
    
        
        ######## plotting purity and stability for chosen binning scheme
        print ('|------> plotting purity and stability for chosen binning scheme')
        ###############Hack, to fix later, move func below to helper scripts
        def makePSplot_simple(purity,stability,var,outputDir,ext='pdf',
                              dictHistos=OrderedDict(),bins=[0.,1.],sel='_WSel',signalLabelBegin='TTToSemiLeptonic'):

            if not os.path.exists(outputDir): os.makedirs(outputDir)

            ROOT.gStyle.SetPadRightMargin(0.05)
            canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)
            p=ROOT.TH1D("Purity",";;",len(bins)-1,array('d',bins))
            s=ROOT.TH1D("Stability",";;",len(bins)-1,array('d',bins))

            for i in range(len(purity)):
                p.SetBinContent(i+1,purity[i])
                s.SetBinContent(i+1,stability[i])

            p.SetLineWidth(2)
            p.SetLineColor(colorPallete[3])
            #purity.SetLineStyle(1)
            s.SetLineWidth(2)
            s.SetLineColor(colorPallete[4])

            legend=ROOT.TLegend(0.15,0.15,0.90,0.35)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetLineColor(0)
            legend.SetTextSize(0.04)
            #legend.SetNColumns( 4 )


            legend.AddEntry( p, 'Purity', 'l' )
            legend.AddEntry( s,  'Stability' , 'l' )

            p.SetMaximum( 1. )
            p.SetMinimum( 0.2 )
            p.SetLineColor( ROOT.kBlack )
            p.SetMarkerColor( ROOT.kBlack )
            p.SetMarkerSize( 0.5 )
            p.SetLineWidth( 2 )

            s.SetLineColor( ROOT.kMagenta )
            s.SetMarkerColor( ROOT.kMagenta )
            s.SetMarkerSize( 0.5 )
            s.SetLineWidth( 2 )
            p.GetXaxis().SetTitle( nSubVariables[var]['label'] )
            p.Draw('hist')
            s.Draw('hist same')

            dictHistos[ 'purityGraph_'+var ] = p.Clone()
            dictHistos[ 'stabilityGraph_'+var ] = s.Clone()

            
            legend.Draw()
            CMS_lumi.extraText = "Simulation"
            CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
            CMS_lumi.relPosX = 0.11
            CMS_lumi.CMS_lumi(canvas, 4, 0)
            
            canvas.Update()
            canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+version+year+'.'+ext)
            if ext.startswith('pdf'):
                canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+version+'_'+year+'.png')
            return 1#numBins
        getAndPlotPurity(signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone().RebinY(2),genBin,variables,ivar, outputDir=outputDir,year=year,sel=sel)
        
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
                                            ROOT.TUnfold.kRegModeCurvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                            ROOT.TUnfold.kEConstraintNone,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. (Need to look into this!!)
                                            ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                            )

        ##### Defining input (data recoJet )
        print ('|------> TUnfolding adding input:')

        tunfolder.SetInput( allHistos[ 'dataHisto' ])
        
        if process.startswith('MCCrossClosure'): 
            tunfolder_cross = ROOT.TUnfoldDensity(
                                            altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel], 
                                            ROOT.TUnfold.kHistMapOutputHoriz,  
                                            ROOT.TUnfold.kRegModeCurvature,   
                                            ROOT.TUnfold.kEConstraintNone,    
                                            ROOT.TUnfoldDensity.kDensityModeBinWidth  
                                            )

            ##### Defining input (data recoJet )
            print ('|------> TUnfolding adding input:')
            tunfolder_cross.SetInput( allHistos[ 'dataHisto' ])
        
        if process.startswith('data'):
            print ("Subtracting backgrounds")
            dummy=0
            bkgSources = []
            if sel.startswith(('_W','_top')):
            
                for ibkg in bkgHistos:
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel):
                        #print (ibkg.split('_')[0]+ '%d'%dummy)
                        tunfolder.SubtractBackground( bkgHistos[ibkg].Clone(), ibkg.split('_')[0]+ '%d'%dummy )
                        dummy=dummy+1
                        bkgSources.append(ibkg.split('_')[0]+ '%d'%dummy)
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone() )

            
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

                elif sys.startswith('_CR'):
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

                elif sys.startswith('_erdON'): 
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

            '''
            !!!!!!!!!!!!!!FIXME!!!!!!!!!!!!!!!
            ### Making unc plot
            plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                dictUncHistos,
                                ivar+'_'+signalLabel+'_allSys',
                                labelX=variables[ivar]['label'],
                                version=sel+'_'+version,
                                year= ( '2016+2017+2018' if year.startswith('all') else year ),
                                outputDir=outputDir
                                )
            '''
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

        allHistos [ 'foldHisto'+ivar ] = tunfolder.GetFoldedOutput("folded"+ivar).Clone()
        if process.startswith('data'): 
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs',  allHistos [ 'foldHisto'+ivar ].Clone(), 'folded', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            plotSimpleComparison( allHistos[ 'unfoldHisto'+ ivar ].Clone(), 'unfold',  signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel].Clone(), 'accepgen', 
                                  ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", 
                                  rebinX=1, version=sel+'_'+version, outputDir=outputDir )




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
        
        #tmp = OrderedDict()
        #dummyM = hist2array(allHistos[ 'cov'+ivar ].Clone())
        #if verbose: print (dummyM.shape)
        #systcovTotal = np.zeros(dummyM.shape)
        #print (systcovTotal.shape)
        #for k in uncerUnfoldHisto:
        #    corrflag = False
        #    sysErr=np.zeros(uncerUnfoldHisto[k].GetNbinsX())
        #    if not k.endswith('Total') or k.endswith(('SystTotal')): continue

        #    for x in jes_corr_list:
        #        if corrflag==True: break
        #        if x.upper() in k:
        #            corrflag=True
        #            continue
                    
        #    #print ("CORRFLAG",corrflag)
        #    for i in range( 0, uncerUnfoldHisto[k].GetNbinsX()):


        #        sysErr[i] = uncerUnfoldHisto[k].GetBinContent(i+1)
        #        if corrflag: 
        #            systcovTotal+=np.outer(sysErr,sysErr)
        #            #print ("CORRFLAG",corrflag)
        #        else: systcovTotal+=np.diag(sysErr*sysErr)
        #print (systcovTotal.shape,corrflag)
        
        
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
                
                sys_cov_up = tunfolder.GetEmatrixSysUncorr("cov_%s_up"%sys+ivar)
                sys_cov_down = tunfolder.GetEmatrixSysUncorr("cov_%s_down"%sys+ivar)
                
                if not sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_TuneCP5')):
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
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone()      
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
                    
                    for i in range( 1, tmpModelHisto.GetNbinsX() + 1):
                        uncerUnfoldHisto[ivar+'_Physics ModelTotal'].SetBinContent( i, abs(tmpModelHisto.GetBinContent(i)))
                    
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                    
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


                elif sys.startswith(('_CR','_erdON')):
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

                
        allHistos [ 'foldHisto'+ivar ] = get_folded_unfolded(folded=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                                                     unfolded=allHistos['unfoldHisto'+ivar].Clone(), 
                                                     cov_tot=allHistos['cov'+ivar].Clone(), 
                                                     probaM=allHistos[ 'probaMatrix'+ivar ]
                                                     )
        
        
        #####################################
        #print (allHistos[ 'unfoldHisto'+ivar ].Clone())
        ###### Plot unfolding results
        #print ('|------> Drawing Data/MC comparison:')
        #if process.startswith('data'):

        
        #               )
        if return_tunfolder_object:
            return tunfolder, allHistos, uncerUnfoldHisto, ratioHistos, uncerUnfoldSystCov
        
        print ('|------> Drawing unfold plot:')
        if not 'Closure' in process:
            drawUnfold(ivar=ivar, 
                       selection=sel, year=year,lumi=lumi, process=process,
                       dataJetHisto=allHistos[ 'dataHistoGenBin' ].Clone(),
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
                       
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_'+version+'.'+ext
                       )
        else: 
            if 'Cross' in process:
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel],
                             unfoldHistoCross=allHistos['unfoldHistoCross'+ivar ].Clone(),
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel],
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext
                             )
            else:
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=[], 
                             unfoldHistoCross=[] ,
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel],
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext
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
            draw2D( ivar,  tunfolder.GetRhoItotal("rhoI"+ivar, "Global correlations"), variables[ivar], outputLabel='data_rhoI', outputDir=outputDir,selection=sel,version=version,year=year)
            draw2D( ivar,  allHistos[ 'correlation_matrix_'+ivar ].Clone(), variables[ivar], outputLabel='data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov'+ivar].Clone(), variables[ivar], outputLabel='dataTotal_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_'+ivar].Clone(), variables[ivar], outputLabel='uncorrUncRM_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_data_'+ivar].Clone(), variables[ivar], outputLabel='dataInpStats_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_bkg_'+ivar].Clone(), variables[ivar], outputLabel='BkgSubtractionSyst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_systTotal'+ivar].Clone(), variables[ivar], outputLabel='Syst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        
        draw2D( ivar,  allHistos[ 'probaMatrix'+ivar ].Clone(), variables[ivar], outputLabel='data_probaMatrix', outputDir=outputDir, addCorrelation=True, addCondition=True ,selection=sel,version=version,year=year)
        draw2D( ivar, signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone(), variables[ivar], outputLabel='data_respMatrix', outputDir=outputDir, addCorrelation=True ,selection=sel,version=version,year=year)
            
        normed_covs = OrderedDict()
        for i in allHistos:
            if 'cov' in i: 
                _, normed_covs['Normed'+i] = GetNormalizedTMatrixandTH2( allHistos[i].Clone(), allHistos[i].GetTitle()+'Normed',
                                                            allHistos['unfoldHisto'+ivar].Clone() )
                #normed_covs.append(normed_cov)
                draw2D( ivar, normed_covs['Normed'+i].Clone(), variables[ivar], outputLabel='Normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
        
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
        
        #print (allHistos,signalHistos,dataHistos,altSignalHistos)
        ######### Saving Histograms
        def renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                ihis = hist.Clone()
                ihis.SetName(isam)
                ihis.SetTitle(isam)
                ihis.Write()

        outputRootName = outputDir+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root'
        print ('|------> Saving histograms in rootfile: ', outputRootName)
        outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
        renamingHistos( signalHistos )
        if not process.startswith('MC'):
            renamingHistos( sysSignalHistos )
        if not('self' in process.lower()):
            renamingHistos( altSignalHistos )
        if process.startswith('data') and  selection.startswith(('_W','_top')): 
            renamingHistos( varSignalHistos )
        renamingHistos(uncerUnfoldSystCov)
        renamingHistos(normed_covs)
        renamingHistos( dataHistos )
        #print ("adding varsighistos", varSignalHistos.items())
        #renamingHistos( dataHistos )
        renamingHistos( bkgHistos )
        renamingHistos( allHistos )
        renamingHistos( uncerUnfoldHisto )
        
        tunfolder.Write()
        outputRoot.Close()

        #print ('|------> Saving histograms in yodafile: ', outputRootName.replace('.root', '.yoda'))
        #histToYoda = [  yoda.root.to_yoda( allHistos [ 'unfoldHisto'+ivar ] ) ]
        #yoda.writeYODA( histToYoda, outputRootName.replace('.root', '.yoda') )

##############################################################################################





##############################################################################################
############ Histogram loader and rebinner for inputs to unfolding script below ##############
##############################################################################################

def loadHistograms(samples, var, sel, sysUnc=[],
                   isMC=True, addGenInfo=True, respOnly=False, lumi=1., noResp=False,
                   variables={}, year='2017', process='data', noRebin=False, outputFolder=None ):
    """docstring for loadHistograms"""

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
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in SYSUNC if 'nom' in syst] 
            if not(noResp):tmpList = tmpList + [ 'accepgen'+var+syst+sel for syst in SYSUNC ]
            if not(noResp):tmpList = tmpList + [ 'truereco'+var+syst+sel for syst in SYSUNC ]
            tmpList = tmpList + [ 'reco'+var+syst+sel for syst in SYSUNC ]
            if not(noResp): tmpList = tmpList + [ 'respWithMiss'+var+syst+sel for syst in SYSUNC]
        
        elif isMC and addGenInfo and flip: 
            tmpList = tmpList + [ 'gen'+var+'_nom'+sel]
            tmpList = tmpList + ['reco'+var+syst+sel for syst in tmpSYSUNC]
            if not(noResp): tmpList = tmpList + ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC]
        
        if respOnly and not flip: 
            if not(noResp): tmpList = [ 'respWithMiss'+var+syst+sel for syst in SYSUNC] 
        elif respOnly and flip: 
            if not(noResp): tmpList = ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC ]
                
        #print (f'Processing {isam}')#,f'filename={samples[isam][0]}')
        for ih in tmpList:
            #print(ih)
            #print (f'From the above, processing histo: {ih},{isam+"_"+ih}')
            #print (samples[isam][0])
            
            if isMC:
                iFile = ROOT.TFile.Open(samples[isam][0],'r')
                allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone() #'jetObservables/'+
                allHistos[isam+'_'+ih].SetDirectory(0)
                allHistos[isam+'_'+ih].Sumw2()
                iFile.Close()
                MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1]['nGenWeights']
                if np.round(MCScale,5)!=np.round(samples[isam][1]['MCScaling'],5): 
                    print (year,'MCScale mismatch alert', MCScale, samples[isam][1]['MCScaling'])
                    allHistos[isam+'_'+ih].Sumw2()
                    allHistos[isam+'_'+ih].Scale( MCScale )
                else: 
                    allHistos[isam+'_'+ih].Sumw2()
                    allHistos[isam+'_'+ih].Scale( MCScale )
            else:
                if 'dijet' in sel and process.startswith('data'):
                    #add pre/unprescaled trigger histograms to get full spectrum
                    tmpdataHistos = {}
                    itempFile=ROOT.TFile.Open(samples[isam][0],'r')
                    for it in checkDict( 'JetHT', dictSamples )[year]['triggerList']:
                        
                        
                        #print(it,ih,ih.replace( sel, '_'+it+sel ))#,itempFile,samples[isam][0])
                        tmpdataHistos[ it ] = itempFile.Get( ih.replace( sel, '_'+it+sel )).Clone()
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
    if 'dijet' in sel.lower():
        if isMC:
            c=0
            allHistos_upd={}
            for sys in SYSUNC:

                print (f"Adding together histograms for {var} in {sel} for {sys}")
                tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k and sys in k}
                #print(tmpHistos.keys(),allHistos.keys())
                for ih in tmpHistos:
                    for jh in allHistos:
                        goodflag=False
                        if ('2016' in ih and '2016' in jh and '2016' in year) or ('2016' in ih and '2016' in jh and '2016_preVFP' in year) or ('2017' in ih and '2017' in jh and '2017' in year) or ('2018' in ih and '2018' in jh and '2018' in year) and sys in ih and sys in jh:
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral())

                        elif not('201' in ih) and not ('201'in jh) and sys in ih and sys in jh:
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral())

                        if goodflag and (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ) and sys in ih and sys in jh :
                            if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                    #print(goodflag,ih,jh,tmpHistos[ih].Integral())

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
        #else:
            
        

    if not noRebin:
        #print("About to rebin histos:")#,allHistos.keys())#,tmpHistos.keys())            
        print("Proceeding to rebin histograms")
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
            outputRootName = outputFolder.split('jetObservables_histograms')[0] + f"/{sel.replace('_','')}/{year}" + '/loadedHistograms_main_'+process+isam+var+year+'.root'
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
        
        ## missing gen (reco inefficiency) correction
        #multiplicand = 1 #+ 
        #multiplicand += miss.GetBinContent(i)
        
        #bc_new = unfolded_hist.GetBinContent(i)+miss.GetBinContent(i)
        
        print(i,bc,miss.GetBinContent(i))
        
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

        print(i,h_gen.GetBinContent(i),
              h_missgen.GetBinContent(i),genSubtract[i],
              (h_missgen.GetBinContent(i)/genSubtract[i] if not genSubtract[i]==0 else 0.))
    h_missRate.SetDirectory(0)
    
    return h_missRate

def correctEfficiency_Rate(unfolded_hist,miss_hist,gen_hist):
    missrate=getMissRate(h_gen,h_missgen).Clone()
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
        
        
        print(i,bc,missrate.GetBinContent(i),scaling)
        
        bc*=scaling

        print("#",i,bc)
        
        aTH1.SetBinContent(i, bc)
        
        
    aTH1.SetDirectory(0)   
    return aTH1

def getAndPlotPurity(h_resp_rebinned,gen_bins,variables,var,sel='_dijetSel',outputDir='../Results/',year='2017'):
    
    rebinned=h_resp_rebinned.Clone()#make_rebinned_2d_hist(h_resp.Clone(),new_bin_edge_pairs,)#rebinning to gen-level bins
    
    arr_rebinned,_ = th2_to_np_arr(rebinned.Clone())
    rebinned_array2d_normX = renorm(arr_rebinned, axis=0) # normalise axis to 1, renormed per x/gen bin
    rebinned_array2d_normY = renorm(arr_rebinned, axis=1) # normalise axis to 1, renormed per y/reco bin
    
    p_list=[]
    s_list=[]
    
    for ibin in range(len(gen_bins)-1):
        #print (f"Calculating p/s per bin in final new binning for bin: {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]}")
        purity = rebinned_array2d_normY[ibin][ibin] #contains fraction in a reco bin that are actually from the same gen bin
        stability = rebinned_array2d_normX[ibin][ibin] #contains fraction in a gen bin that are actually from the same reco bin
        p_list.append(purity)
        s_list.append(stability)
        #print (f"Purity, stability in bin {ibin}({gen_bins[i],gen_bins[i+1]}): {purity,stability}")
    


    
    makePSplot_simple(purity=array('d',p_list),stability=array('d',s_list),
                      dictHistos=OrderedDict(),variables=variables,var=var,outputDir=outputDir,bins=gen_bins,year=year)
    return 1

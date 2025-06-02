from collections import OrderedDict
import copy
import pprint 
import ROOT
import numpy as np
import array
from array import array
import bisect
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from root_numpy import array2hist, hist2array
from histoHelpers import *
from unfoldingPlottersAndHelpers import *
#from drawPaperPlots import *
import os
import glob
import sys
import math
import yoda

sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
sys.path.insert(0,'../../')
from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict
import json

import gc
#lumi=0.
#canvas = {}
#textBox=ROOT.TLatex()
#textBox.SetTextSize(0.10)
#textBox.SetTextAlign(12)

##############################################################################################

#################################### MAIN Unfolding script ###################################

##############################################################################################
""" 
    Original script for 1-D single obs. unfoldings. The script is still used to load up histograms from outputHistograms saved-per-year
    for each unfolded observable in each selection when doing the all-years and combined/simultaneous unfolding. 
"""


def doSimultaneousUnfolding():
    pass

def runTUnfold(
                dataFile,variables, sel, sysUncert, process, ext, sigFiles=None, bkgFiles=None,  lumi=1., sysSignalLabels=[], 
                year='2017',runMLU=False, sysSigFiles={}, varSigFiles={}, outputFolder='../Results',
                version='_Sept23',verbose=False, return_tunfolder_object = False, mainMC='MLM_HTbin', altMC='Ptbin',extraMC=False,
                noScale_dijets=False, include_FSR_in_unfolded_result=True, areaConstraint=False, dict_t3_samples=None
              ):
    
    
    """ 
    
    Response matrix for unfoldings have miss corrections applied already a la UF bin on y-axis;
    
    a) this means that the unfolded distribution's integral should match the overall gen-count in closure tests (not just accepgen counts), 
    b) and, the folded unfolded distribution from tunfold will match the true-reco in closure tests, but not exactly the data-bkg in data unfolding since the unfolded contains miss-gen counts, however unfolded minus mmissgen folded forward will close with data-bkg reco;
    """    
    
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    ROOT.TH1.StatOverflows(ROOT.kTRUE)
    ROOT.TH2.StatOverflows(ROOT.kTRUE)

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 41, 45, 48, 49, 37, 38, 33, 17]
    dict_condition_numbers = OrderedDict()
    if 'dijet' in sel:
        dict_MCScaling = OrderedDict()
        
        
    #print(variables.keys())
    import sys
    #unf_var_dict = OrderedDict()
    if not('dijet' in sel):
        sys.path.insert(0,'../../')
        from datasets_WtopSel_RunIISummer20UL_SampleDictPrep_newXS import dictSamples, checkDict
        
        dictSamples = dict_t3_samples
        
        signalLabel = 'TTToSemiLeptonic'
        sigPlotLabel = 'PWHG+P8'
        signalLabelBegin = 'TTTo'
        varSignalLabelBegin = 'varTTToSemileptonic'
        sysSignalLabelBegin = 'sysTTToSemiLeptonic'
        fsrLabel = 'sysTTToSemiLeptonic_fsrWeight'
        if altMC.startswith('TT_TuneCH3'):
            
            altSignalLabelBegin = 'TT_TuneCH3'
            altSigPlotLabel = 'PWHG+H7'
            altSignalLabel = 'TT_TuneCH3'
            
            alt1SignalLabelBegin = 'TTJets'
            alt1SigPlotLabel = 'aMC@NLO-FXFX+P8'
            alt1SignalLabel = 'TTJets'
            alt1 = 'TTJets'
            alt1MC = 'TTJets'
            
        elif altMC.startswith('TTJets'):
            
            altSignalLabelBegin = 'TTJets'
            altSigPlotLabel = 'aMC@NLO-FXFX+P8'
            altSignalLabel = 'TTJets'
            
            alt1SignalLabelBegin = 'TT_TuneCH3'
            alt1SigPlotLabel = 'PWHG+H7'
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
                     #'TTToHadronic', 'TTTo2L2Nu', 
                     'WJetsToLNu', 
                     'ST_s-channel_4f_leptonDecays', 
                     'ST_t-channel_top_5f_InclusiveDecays', 'ST_t-channel_antitop_5f_InclusiveDecays', 
                     'ST_tW_top_5f_NoFullyHadronicDecays', 'ST_tW_antitop_5f_NoFullyHadronicDecays', 
                     'WW', 'ZZ', 'WZ', 
                     'DYJetsToLL_M-50',
                     'QCD_Pt-1000_MuEnrichedPt5'
                    ]
    else:
        sys.path.insert(0,'../../')
        from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
        
        dictSamples = dict_t3_samples
        
        if extraMC:
            alt1MC='HTbin'
            alt2MC='Ptbin'
        else:
            alt1MC=None
            alt2MC=None
        
        if mainMC.startswith('MLM_HTbin'):
            signalLabelBegin = 'MLMQCD_HT'
            signalLabel = 'MLMQCD_HT2000toInf'
            sigPlotLabel = 'MG5-MLM+P8'


        elif mainMC.startswith('HTbin'):
            signalLabelBegin = 'QCD_HT'
            signalLabel = 'QCD_HT2000toInf'
            sigPlotLabel = 'MG5+P8'

        elif mainMC.startswith('H7MLM_HTbin'):
            signalLabelBegin = 'H7MLMQCD_HT'
            signalLabel = 'H7MLMQCD_HT2000toInf'
            sigPlotLabel = 'MG5-MLM+H7'

        if altMC.startswith('MLM_HTbin'):
            altSignalLabelBegin = 'MLMQCD_HT'
            altSignalLabel = 'MLMQCD_HT2000toInf'
            altSigPlotLabel = 'MG5-MLM+P8'
            
        elif altMC.startswith('HTbin'):
            altSignalLabelBegin = 'QCD_HT'
            altSignalLabel = 'QCD_HT2000toInf'
            altSigPlotLabel = 'MG5+P8'
            
        elif altMC.startswith('H7MLM_HTbin'):
            altSignalLabelBegin = 'H7MLMQCD_HT'
            altSignalLabel = 'H7MLMQCD_HT2000toInf'
            altSigPlotLabel = 'MG5-MLM+H7'
            
        if extraMC:
            if alt1MC.startswith('HTbin'):
                alt1SignalLabelBegin = 'QCD_HT'
                alt1SignalLabel = 'QCD_HT2000toInf'
                alt1SigPlotLabel = 'MG5+P8'
                
            elif alt1MC.startswith('Ptbin'):
                alt1SignalLabelBegin = 'QCD_Pt_'
                alt1SignalLabel = 'QCD_Pt_3200toInf'
                alt1SigPlotLabel = 'P8+P8'
                
            elif alt1MC.startswith('H7MLM_HTbin'):
                alt1SignalLabelBegin = 'H7MLMQCD_HT'
                alt1SignalLabel = 'H7MLMQCD_HT2000toInf'
                alt1SigPlotLabel = 'MG5-MLM+H7'

            if alt2MC.startswith('Ptbin'):
                alt2SignalLabelBegin = 'QCD_Pt_'
                alt2SignalLabel = 'QCD_Pt_3200toInf'
                alt2SigPlotLabel = 'P8+P8'
                
            elif alt2MC.startswith('HTbin'):
                alt2SignalLabelBegin = 'QCD_HT'
                alt2SignalLabel = 'QCD_HT2000toInf'
                alt2SigPlotLabel = 'MG5+P8'
                
            elif alt2MC.startswith('H7MLM_HTbin'):
                alt2SignalLabelBegin = 'H7MLMQCD_HT'
                alt2SignalLabel = 'H7MLMQCD_HT2000toInf'
                alt2SigPlotLabel = 'MG5-MLM+H7'
                
            
        dict_MCScaling[signalLabel] = OrderedDict()
        if not('self' in process.lower()):
            dict_MCScaling[altSignalLabel] = OrderedDict()
            if extraMC:
                dict_MCScaling[alt1SignalLabel] = OrderedDict()
                dict_MCScaling[alt2SignalLabel] = OrderedDict()
        
        sysSignalLabelBegin = 'sysMLMQCD'  if mainMC.startswith('MLM') else 'sysQCD' 
        fsrLabel = 'sysMLMQCD_fsrWeight_HT2000toInf' if mainMC.startswith('MLM') else 'sysQCD_fsrWeight_HT2000toInf' 
        
        
    
    dict_of_dicts = OrderedDict()
    for ivar in variables:
        gc.collect()
        
        outputDir=outputFolder+sel.split('_')[1]+'/'+year+'/Unfolding/'+ivar+'/'+process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        genBin = variables[ivar]['bins']
        recoBin = variables[ivar]['bins_reco']
        
        
        #################################################################################################
        ####### using the build_all_years_histograms to accumulate histos for each 1-D unfolding ########
        ###### from output/rebinned histograms of inidividual years' 1D unfolding in a dictionary #######
        ##### then these are used for 1D unfolding of the combination of years (single observables) #####
        #################################################################################################
        if year.startswith('all'):
            dataHistos = { }
            dataHistostrue = {}
            bkgHistos = { }
            sysUncs_added = []
            years_list = ['2016_preVFP','2016','2017','2018']

            dict_of_dicts[ivar] = build_all_years_histograms(
                                                            process=process,
                                                            year=year,
                                                            years_list=years_list,
                                                            dataFile=dataFile,
                                                            signalLabel=signalLabel,
                                                            altSignalLabel=altSignalLabel,
                                                            alt1SignalLabel=alt1SignalLabel if not ("mc" in process.lower()) and extraMC else False,
                                                            alt2SignalLabelBegin=alt2SignalLabelBegin if not ("mc" in process.lower()) and extraMC else False,
                                                            alt2SignalLabel=alt2SignalLabel if not ("mc" in process.lower()) and extraMC else False,
                                                            sysSignalLabels=sysSignalLabels,
                                                            sysUncert=sysUncert,
                                                            varSignalLabels=varSignalLabels if not 'dijet' in sel else None,
                                                            bkgLabels=bkgLabels if not 'dijet' in sel else None,
                                                            fsrLabel=fsrLabel,
                                                            ivar=ivar,
                                                            sel=sel,
                                                            genBin=genBin,
                                                            extraMC=extraMC,
                                                            verbose=True
                                                        )
            for dictKey in dict_of_dicts[ivar].keys():
                try:
                    for key in dict_of_dicts[ivar][dictKey].keys():
                        dict_of_dicts[ivar][dictKey][key].SetDirectory(0)
                except AttributeError:
                    print(f"WARNING: L250 in unfold script: error, {dictKey,ivar}, has NoneType obj in dict_of_dicts")
                    pass
            
            dataHistos       = dict_of_dicts[ivar]["dataHistos"]
            dataHistostrue   = dict_of_dicts[ivar]["dataHistostrue"]
            signalHistos     = dict_of_dicts[ivar]["signalHistos"]
            sysSignalHistos  = dict_of_dicts[ivar]["sysSignalHistos"]
            altSignalHistos  = dict_of_dicts[ivar]["altSignalHistos"]
            
            if not ("mc" in process.lower()) and extraMC:
                alt1SignalHistos = dict_of_dicts[ivar]["alt1SignalHistos"]
                
            allHistos        = dict_of_dicts[ivar]["allHistos"]
            
            if extraMC==True and 'dijet' in sel:
                alt2SignalHistos = dict_of_dicts[ivar]["alt2SignalHistos"]
                
            if ('W' in sel or 'top' in sel) and not("mc" in process.lower()):
                bkgHistos = dict_of_dicts[ivar]["bkgHistos"]
                varSignalHistos  = dict_of_dicts[ivar]["varSignalHistos"]


        else:
            #print('|-------> Running single year '+year)
            ### Getting input histos
            allHistos = {}
            
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            #print(mainSigFiles.keys())
            
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
            #print(signalHistos.keys())
            
            tmp2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(altSignalLabelBegin)  }
            
            if not('self' in process.lower()): 
                altSignalHistos = loadHistograms( tmp2SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                  lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
                #for ih in altSignalHistos.keys(): print(altSignalHistos[ih].Integral(),altSignalHistos[ih].GetName())
                #print(altSignalHistos.keys())
                    
            if 'data' in process and extraMC:
                tmp3SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt1SignalLabelBegin)  }
                alt1SignalHistos = loadHistograms( tmp3SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                  lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
                #print(alt1SignalHistos.keys())

                if alt2SignalLabelBegin!=None:
                    tmp4SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt2SignalLabelBegin)  }
                    alt2SignalHistos = loadHistograms( tmp4SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, 
                                                      lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder )
                    #print(alt2SignalHistos.keys())

            
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
                
                
            else: 
                
                bkgHistos = {}

            
                
            if process.startswith("MC"):
                dataHistostrue = { 'data_reco'+k.split(('_truereco'))[1] : v.Clone() for (k,v) in signalHistos.items() if ('_truereco' in k)}
                
               

                dataHistos = { 'data_reco'+k.split(('_reco'))[1] : v.Clone() for (k,v) in signalHistos.items() if ('_reco' in k)}
                
                if 'dijet' in sel:
                    dataHistosForRescaling = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process='data', variables=variables,outputFolder=outputFolder)
                    
                    dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Clone('data_reco'+ivar+'_nom'+sel+'_genBin')
                    
                    dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                
                
            else:
                dataHistostrue = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
                dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder)
               
        
        allHistos[ 'dataHisto'+ivar ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel ].Clone()
        allHistos[ 'dataHistoGenBin'+ivar ] = dataHistostrue[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
        
        if "data" in process:

            print("VARIOUS INTEGRALS: data, data genBin, recoMC, genMC, fsrUp/Down MC")
            print("All years:")
            print(
                allHistos[f"dataHisto{ivar}"].Integral(),
                allHistos[f"dataHistoGenBin{ivar}"].Integral(),
                signalHistos[f"{signalLabel}_reco{ivar}_nom{sel}"].Integral(),
                signalHistos[f"{signalLabel}_gen{ivar}_nom{sel}"].Integral(),
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}" in sysSignalHistos
                else "N/A",
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}" in sysSignalHistos
                else "N/A",
            )
        
        #fakes from nominal/signal MC
        fakeHistos = { signalLabel+'_fakereco'+k.split(('_fakereco'))[1]: v.Clone() for (k,v) in signalHistos.items()  if ('_fakereco' in k)}
        #print(year, sel, ivar, fakeHistos.keys())
        
        allHistos[ 'allBkgHisto'+ivar ] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
        allHistos[ 'allBkgHisto'+ivar ].Reset()
        allHistos[ 'allBkgHistoGenBin'+ivar ] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
        allHistos[ 'allBkgHistoGenBin'+ivar ].Reset()
        
        print("Bkg. Histo collection integral:", allHistos[ 'allBkgHisto'+ivar ].Integral())
        #include fakes in histo containing all bkg. (including other physics processes in the case of W/top) from MC
        for ih in fakeHistos:
            
            if ih.endswith(ivar+'_nom'+sel): 
                #print("FakeHisto Check:", ih)
                allHistos[ 'allBkgHisto'+ivar ].Add( fakeHistos[ih] )
            elif ih.endswith('_genBin'): 
                #print("FakeHisto Check:", ih)
                allHistos[ 'allBkgHistoGenBin'+ivar ].Add( fakeHistos[ih] )
        print("Bkg. Histo collection integral:", allHistos[ 'allBkgHisto'+ivar ].Integral())
        allHistos[ 'allMCHisto'+ivar ] = allHistos[ 'allBkgHisto'+ivar ].Clone()
        allHistos[ 'allMCHistoGenBin'+ivar ] = allHistos[ 'allBkgHistoGenBin'+ivar ].Clone()
        
        if process.startswith('data'):# and not('dijet' in sel):
            
            if verbose: print(bkgHistos.keys())
            for ibkg in bkgHistos:
                # if verbose: 
                #print(f"Adding in {ibkg}")

                if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto'+ivar ].Add( bkgHistos[ibkg].Clone() )
                if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin'+ivar ].Add( bkgHistos[ibkg].Clone() )

        
        
        allHistos[ 'dataMinusBkgs'+ivar ] = allHistos[ 'dataHisto'+ivar ].Clone() 
        allHistos[ 'dataMinusBkgsGenBin'+ivar ] = allHistos[ 'dataHistoGenBin'+ivar ].Clone() 
        
        
        allHistos[ 'allMCHisto'+ivar ].Add( signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Clone() ) #add only true since fakes already contained in the allMC histo via the allBkgHisto
       

        allHistos[ 'allMCHistoGenBin'+ivar ].Add( signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel+'_genBin' ].Clone() ) 
        
        
        
        ################################################################################################################## 
        # rescale dijet MC to data for the case of individual years; reuse these rescaled histos in the all years case, thus the if
        # logic's exclusion of combo of years below; 
        #if 'data' in process:
        #    print("VARIOUS INTEGRALS: Data, Data genbin, recoMC, genMC, fsrUp/Down MC" )
        #    print(year)
        #    print(allHistos[ 'dataHisto'+ivar ].Integral(), allHistos[ 'dataHistoGenBin' +ivar].Integral(), signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Integral(), sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Integral(),        sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Integral() )
            
        if sel.startswith('_dijet') and not('all' in year):#.startswith('all')): 

            scalingDict = {}    
            ### For dijet, scale QCD to data, as per AGE's past work
            scaleFactor=1.
            scaleFactorGenBin=1.
            altscaleFactor=1.
            altscaleFactorGenBin=1.
            nomIntegral = allHistos[ 'allMCHisto'+ivar].Integral()
            if process.startswith('MC'): 
                
                scaleFactor = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto'+ivar].Integral()
                scaleFactorGenBin = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' +ivar].Integral()
                
                #scaleFactor = round(scaleFactor,15)
                #scaleFactorGenBin = round(scaleFactorGenBin,15)
                
                if not np.isclose(scaleFactor,scaleFactorGenBin): 
                    print(f"WARNING (something weird): gen-/reco-binning mc-to-data SFs are not the same, {scaleFactor,scaleFactorGenBin}")
                    print ("SF genBin and recobin, data int., nominal reco int., respectively:",scaleFactor,scaleFactorGenBin,f"{dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral():.16f}",f"{allHistos[ 'allMCHisto' +ivar].Integral():.16f}")
                
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
                        
                allHistos[ 'dataHisto' +ivar].Scale( scaleFactor )
                allHistos[ 'dataHistoGenBin' +ivar].Scale( scaleFactorGenBin )
                allHistos[ 'dataMinusBkgs'+ivar ].Scale( scaleFactor )
                allHistos[ 'dataMinusBkgsGenBin'+ivar ].Scale( scaleFactorGenBin )
                allHistos[ 'allMCHisto' +ivar].Scale( scaleFactor )

                allHistos[ 'allMCHistoGenBin'+ivar ].Scale( scaleFactorGenBin )

                allHistos[ 'allBkgHistoGenBin'+ivar ].Scale( scaleFactorGenBin )
                dict_MCScaling[signalLabel][ivar] = scaleFactor

            if process.startswith('data'):
                                            
                scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto' +ivar].Integral()
                scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' +ivar].Integral()
                #scaleFactor = round(scaleFactor,15)
                #scaleFactorGenBin = round(scaleFactorGenBin,15)
                
                dict_MCScaling[signalLabel][ivar] = scaleFactor

                if not np.isclose(scaleFactor,scaleFactorGenBin): 
                    print(f"WARNING (something weird): gen-/reco-binning mc-to-data SFs are not the same, {scaleFactor,scaleFactorGenBin}")
                #print ("SF genBin and recobin, data int., nominal reco int., respectively:",scaleFactor,scaleFactorGenBin,dataHistos['data_reco'+ivar+'_nom'+sel].Integral(),allHistos[ 'allMCHisto'+ivar ].Integral())
                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        
                        signalHistos[ihsig].Scale( scaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactor
                        
                        
                    elif ihsig.endswith('genBin'): 
                        signalHistos[ihsig].Scale( scaleFactorGenBin )
                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin
            
                allHistos[ 'allMCHisto'+ivar ].Scale( scaleFactor )
                scalingDict[f'scaling_allMCHisto'] = scaleFactor
                allHistos[ 'allMCHistoGenBin'+ivar ].Scale( scaleFactorGenBin )
                scalingDict[f'scaling_allMCHistoGenBin'+ivar] = scaleFactorGenBin

                #if process=='data':
                allHistos[ 'allBkgHisto' +ivar].Scale( scaleFactor )
                scalingDict[f'scaling_allBkgHisto'] = scaleFactor
                allHistos[ 'allBkgHistoGenBin' +ivar].Scale( scaleFactorGenBin )
                scalingDict[f'scaling_allBkgHistoGenBin'] = scaleFactorGenBin

                if len(sysUncert)!=0:
                    scaleFactor_sys = 1.
                    scaleFactorGenBin_sys = 1.
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
                            
                            if sys.startswith(('_je','_const','_unclust', '_pu', '_l1', '_isr', '_fsr')):#experimental sources, sources where XS unchanged
                                scaleFactor_sys = scaleFactor  
                                scaleFactorGenBin_sys = scaleFactorGenBin 
                                
                            elif ('pdfweight' in sys.lower()):#theory sources that change XS
                                #https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#Modelling_uncertainties_in_gener 
                                #"...In general modeling, uncertainties that change the total cross section should be normalized back to the reference cross section before any acceptance requirements."
                                #----> scaling for varns.: \alpha_s and pdf weights, which change XS, accordingly
                                #i.e., so they match the modified reference cross section a la the nominals
                                
                                scaleFactor_sys = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel].Integral() #scaleFactor #
                                scaleFactorGenBin_sys = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / sysSignalHistos[s[0]+'_reco'+ivar+sys+upDown+sel+'_genBin'].Integral() #scaleFactorGenBin #
                                
                                #scaleFactor_sys = round(scaleFactor_sys,15)
                                #scaleFactorGenBin_sys = round(scaleFactorGenBin_sys,15)    
                            
                            nom_gen_int = signalHistos['MLMQCD_HT2000toInf_gen'+ivar+'_nom'+sel].Integral()
                            
                            for ihsig in sysSignalHistos:
                                if sys+upDown in ihsig:
                                    if ('_genJet' in ihsig) and not('genBin' in ihsig) and (('fsr' in ihsig.lower() or 'pdf' in ihsig.lower())):# or 'Flavor' in ihsig ):#'_recoJet' in ihsig or 
                                        
                                        print ("hist_label,sys, SF genBin and recobin, data int., sys int., respectively:\n",
                                               ihsig,
                                               sys+upDown,
                                               scaleFactor_sys,
                                               scaleFactorGenBin_sys,
                                               f"{dataHistos['data_reco'+ivar+'_nom'+sel].Integral():.16f}",
                                               f"{sysSignalHistos[ihsig].Integral():.16f}", 
                                               f"{nom_gen_int:.16f}",
                                               sysSignalHistos[ihsig].Integral()/(nomIntegral if 'reco' in ihsig or 'resp' in ihsig else nom_gen_int), '\n'
                                               
                                              )
                                        
                                    if ihsig.endswith(sel):
                                        #print (ihsig)
                                        sysSignalHistos[ihsig].Scale( scaleFactor_sys )
                                        scalingDict[f'scaling_{ihsig}'] = scaleFactor_sys

                                    elif ihsig.endswith('genBin'): 
                                        sysSignalHistos[ihsig].Scale( scaleFactorGenBin_sys )
                                        scalingDict[f'scaling_{ihsig}'] = scaleFactorGenBin_sys
                                    
                                    if ('_genJet' in ihsig) and not('genBin' in ihsig) and ('fsr' in ihsig or 'Flavor' in ihsig ):#'_recoJet' in ihsig or 
                                        """
                                        print ("hist_label,sys, SF genBin and recobin, data int., sys int., respectively:\n",
                                        
                                               ihsig,
                                               sys+upDown,
                                               scaleFactor_sys,
                                               scaleFactorGenBin_sys,
                                               dataHistos['data_reco'+ivar+'_nom'+sel].Integral(),
                                               sysSignalHistos[ihsig].Integral(), 
                                               nom_gen_int,
                                               sysSignalHistos[ihsig].Integral()/(nomIntegral if 'reco' in ihsig or 'resp' in ihsig else nom_gen_int), '\n'                                               
                                              )
                                        """
                        if not(sys in already_scaled):already_scaled.append(sys)
                print("Scaling done for foll. sysUncs:", already_scaled)
                                
            if not(process.startswith('MCSelfClosure')):
                if 'data' in process.lower():
                    #print("Rescaling alternate signal MC")
                    altscaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    altscaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )

                elif 'cross' in process.lower(): 
                    altscaleFactor = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel].Integral() / (altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    altscaleFactorGenBin = dataHistosForRescaling['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                
                #altscaleFactor = round(altscaleFactor,15)
                #altscaleFactorGenBin = round(altscaleFactorGenBin,15)
                
                print ("Alt MC SF genBin and recobin, respectively:",altSignalLabel,altscaleFactor,altscaleFactorGenBin)

                for ihsig in altSignalHistos:
                    if ihsig.endswith(sel):
                        altSignalHistos[ihsig].Scale( altscaleFactor )
                        scalingDict[f'scaling_{ihsig}'] = altscaleFactor
                    elif 'genBin' in ihsig: 
                        altSignalHistos[ihsig].Scale( altscaleFactorGenBin )#.endswith('genBin')
                        scalingDict[f'scaling_{ihsig}'] = altscaleFactorGenBin
                
                dict_MCScaling[altSignalLabel][ivar] = altscaleFactor


                
                if extraMC:
                    
                    #print("Rescaling alternate signal MC 1")
                    alt1scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (alt1SignalHistos[ alt1SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    alt1scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( alt1SignalHistos[alt1SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                    
                    #alt1scaleFactor = round(alt1scaleFactor,15)
                    #alt1scaleFactorGenBin = round(alt1scaleFactorGenBin,15)
                    print ("Alt MC 1 SF genBin and recobin, respectively:",alt1SignalLabel,alt1scaleFactor,alt1scaleFactorGenBin)

                    
                    for ihsig in alt1SignalHistos:
                        if ihsig.endswith(sel):
                            alt1SignalHistos[ihsig].Scale( alt1scaleFactor )
                            scalingDict[f'scaling_{ihsig}'] = alt1scaleFactor
                        elif 'genBin' in ihsig: 
                            alt1SignalHistos[ihsig].Scale( alt1scaleFactorGenBin )#.endswith('genBin')
                            scalingDict[f'scaling_{ihsig}'] = alt1scaleFactorGenBin
                            
                            
                    #print("Rescaling alternate signal MC 2")
                    alt2scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / (alt2SignalHistos[ alt2SignalLabel+'_reco'+ivar+'_nom'+sel ].Integral()  )
                    alt2scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / ( alt2SignalHistos[alt2SignalLabel+'_reco'+ivar+'_nom'+sel+'_genBin'].Integral() )
                    
                    #alt2scaleFactor = round(alt2scaleFactor,15)
                    #alt2scaleFactorGenBin = round(alt2scaleFactorGenBin,15)
                    print ("Alt MC 2 SF genBin and recobin, respectively:",alt2SignalLabel,alt2scaleFactor,alt2scaleFactorGenBin)

                    
                    for ihsig in alt2SignalHistos:
                        if ihsig.endswith(sel):
                            alt2SignalHistos[ihsig].Scale( alt2scaleFactor )
                            scalingDict[f'scaling_{ihsig}'] = alt2scaleFactor
                        elif 'genBin' in ihsig: 
                            alt2SignalHistos[ihsig].Scale( alt2scaleFactorGenBin )#.endswith('genBin')
                            scalingDict[f'scaling_{ihsig}'] = alt2scaleFactorGenBin
                    
                    dict_MCScaling[alt1SignalLabel][ivar] = alt1scaleFactor
                    dict_MCScaling[alt2SignalLabel][ivar] = alt2scaleFactor

        ################################################################################################################## 
        if process.startswith('data'):


            #print('data_2016_preVFP', (dataFile[ivar+'_2016_preVFP'].Get( f"dataHisto{ivar}" )).Integral(),  
            #      'data_2016', (dataFile[ivar+'_2016'].Get( f"dataHisto{ivar}" )).Integral(), 
            #      'data_2017', (dataFile[ivar+'_2017'].Get( f"dataHisto{ivar}" )).Integral(), 
            #      'data_2018', (dataFile[ivar+'_2018'].Get( f"dataHisto{ivar}" )).Integral(), )

            print("VARIOUS INTEGRALS: data, data genBin, recoMC, genMC, fsrUp/Down MC")
            print("All years:")
            print(
                allHistos[f"dataHisto{ivar}"].Integral(),
                allHistos[f"dataHistoGenBin{ivar}"].Integral(),
                signalHistos[f"{signalLabel}_reco{ivar}_nom{sel}"].Integral(),
                signalHistos[f"{signalLabel}_gen{ivar}_nom{sel}"].Integral(),
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}" in sysSignalHistos
                else "N/A",
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}" in sysSignalHistos
                else "N/A",
            )
            
            print("DOING dataminusBkgs")
            print(allHistos[ 'dataMinusBkgs' +ivar].Integral(),allHistos[ 'dataMinusBkgsGenBin'+ivar ].Integral())

            allHistos[ 'dataMinusBkgs'+ivar ].Add( allHistos[ 'allBkgHisto'+ivar ].Clone(), -1 )
            allHistos[ 'dataMinusBkgsGenBin'+ivar ].Add( allHistos[ 'allBkgHistoGenBin' +ivar].Clone(), -1 )
            
            print(allHistos[ 'dataMinusBkgs'+ivar ].Integral(),allHistos[ 'dataMinusBkgsGenBin'+ivar ].Integral())
            
        ######## Cross check: plotting data vs all MC (scaled to data for QCD, after normalising as per usual by lumi and xs and genweights; only the latter of course applies to W/top which aren't scaled to data)
        
        #if 'data' in process:
        #    print("VARIOUS INTEGRALS: Data, Data genbin, recoMC, genMC, fsrUp/Down MC" )
        #    print(year)
        #    print(allHistos[ 'dataHisto'+ivar ].Integral(), allHistos[ 'dataHistoGenBin' +ivar].Integral(), signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Integral(), signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Integral(), sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Integral(),sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Integral() )
            
        print ('|------> Cross check: plotting data vs all MC')
        #print(dataHistostrue,'data_reco'+ivar+'_nom'+sel)
        plotSimpleComparison(dataHistos['data_reco'+ivar+'_nom'+sel].Clone(), 
                             'data', allHistos[ 'allMCHisto'+ivar ].Clone(), 'allMC', 
                             ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_nom", 
                             rebinX=1, version=sel+'_'+version, outputDir=outputDir )


        
        #################################### Make diagnostics/misc. plots prior to unfolding #################################

        #print ('|------> Pre-unfolding cross-check plots '+ivar)


        ####### Cross check response matrix
        tmpGenHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionX()
        #if 'dijet' in sel: tmpGenHisto.Scale( scaleFactorGenBin )

        plotSimpleComparison( tmpGenHisto, 'projection', signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(), 'RegularGen', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionGen", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        tmpRecoHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
        #if 'dijet' in sel: tmpRecoHisto.Scale( scaleFactor )

        plotSimpleComparison( tmpRecoHisto, 'projection', signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Clone(), 'Regular_TrueReco', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionReco", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        ####### Plotting bkg subtracted reco vs. data and including fakes in bkgHistos

        
            

        tmpHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY()
        plotSimpleComparison( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone(), 'data', allHistos[ 'allBkgHisto'+ivar ].Clone(), 
                             'Bkg+fakes', ivar+'_from'+('Data' if process.startswith('data') else 'MC') + '_' + signalLabel + "_TestDataBkgFakes",
                             rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        plotSimpleComparison( allHistos[ 'dataMinusBkgs'+ivar ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal_truereco', 
                              ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", 
                              rebinX=1, version=sel+'_'+version, outputDir=outputDir )
        
                                            
        
        getAndPlotPurity(signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone().RebinY(2),
                         reco=signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone().Rebin(2),
                         accepgen=signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel].Clone(),
                         gen_bins=genBin,variables=variables,var=ivar, 
                         lumi=lumi,
                         outputDir=outputDir,year=year,sel=sel)
        
        ######## Cross check: plotting response matrix
        #print ('|------> Cross check: plotting response matrix for signal')
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
        print ('|------> (T)Unfolding starts:')

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
        if not ('data' in process):
            tunfolder.SetInput(signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel])
            print ("Subtracting fakes for nom. unf.")
            tunfolder.SubtractBackground(signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel], 'fakes')
        
        else:
            tunfolder.SetInput(allHistos[ 'dataHisto'+ivar ])
        
        if process.startswith('MCCrossClosure'): 
            tunfolder_cross = ROOT.TUnfoldDensity(
                                            altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel], 
                                            ROOT.TUnfold.kHistMapOutputHoriz,  
                                            ROOT.TUnfold.kRegModeNone,#Curvature,   
                                            ROOT.TUnfold.kEConstraintNone if areaConstraint==None else ROOT.TUnfold.kEConstraintArea,    
                                            #ROOT.TUnfoldDensity.kDensityModeBinWidth  
                                            )

            ##### Defining input (data recoJet )
            #print ('|------> TUnfolding adding input:')
            tunfolder_cross.SetInput(signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel])#.SetInput( allHistos[ 'dataHisto'+ivar ])
            print ("Subtracting fakes for alt. unf.")

            tunfolder_cross.SubtractBackground(signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel], 'fakes')
        
        if process.startswith('data'):
            print ("Subtracting backgrounds")
            if 'all' in year:
                print( 'integrals of fakes, true, dataMinusBkg, initial data', signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral(), signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),allHistos[ 'dataMinusBkgs'+ivar ].Integral(),allHistos[ 'dataHisto' +ivar].Integral())
            
            #dummy=0
            bkgSources = []
            #if sel.startswith(('_W','_top')):
            dataMinusbkg_counter = allHistos[ 'dataHisto'+ivar ].Integral()
            if not('dijet' in sel):
                suff = '_reco'+ivar+'_nom'+sel
                temp_bkgHistos = OrderedDict()
                scaleUnc_bkgHistos = OrderedDict()
                temp_bkgHistos['ST'+suff] = None
                temp_bkgHistos['WJets'+suff] = None
                temp_bkgHistos['QCD'+suff] = None
                temp_bkgHistos['DY'+suff] = None
                temp_bkgHistos['VV'+suff] = None
                
                scaleUnc_bkgHistos['ST'+suff] = 0.23
                scaleUnc_bkgHistos['WJets'+suff] = 0.19
                scaleUnc_bkgHistos['QCD'+suff] = 1.
                scaleUnc_bkgHistos['DY'+suff] = 1.#setting to 100% as per top mass measurement paper  #Z+jets and VV scale unc. set to 50% in ATLAS 2019 substructure paper
                scaleUnc_bkgHistos['VV'+suff] = 1.#setting to 100% as per top mass measurement paper  #Z+jets and VV scale unc. set to 50% in ATLAS 2019 substructure paper
                
                for ibkg in bkgHistos:
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel):
                        print(ibkg,bkgHistos[ibkg].Integral())
                        
                        if 'ST' in ibkg:
                            if not(temp_bkgHistos['ST'+suff] == None):
                                temp_bkgHistos['ST'+suff].Add(bkgHistos[ibkg])
                            else:
                                temp_bkgHistos['ST'+suff] = bkgHistos[ibkg].Clone('ST'+suff+'_bkgSubComb')
                                
                        elif 'WJ' in ibkg:
                            if not(temp_bkgHistos['WJets'+suff] == None):
                                temp_bkgHistos['WJets'+suff].Add(bkgHistos[ibkg])
                            else:
                                temp_bkgHistos['WJets'+suff] = bkgHistos[ibkg].Clone('WJets'+suff+'_bkgSubComb')
                                
                        elif 'WW' in ibkg or 'WZ' in ibkg or 'ZZ' in ibkg:
                            if not(temp_bkgHistos['VV'+suff] == None):
                                temp_bkgHistos['VV'+suff].Add(bkgHistos[ibkg])
                            else:
                                temp_bkgHistos['VV'+suff] = bkgHistos[ibkg].Clone('VV'+suff+'_bkgSubComb')
                                
                        elif 'DY' in ibkg:
                            if not(temp_bkgHistos['DY'+suff] == None):
                                temp_bkgHistos['DY'+suff].Add(bkgHistos[ibkg])
                            else:
                                temp_bkgHistos['DY'+suff] = bkgHistos[ibkg].Clone('DY'+suff+'_bkgSubComb')
                        
                        elif 'QCD' in ibkg:
                            if not(temp_bkgHistos['QCD'+suff] == None):
                                temp_bkgHistos['QCD'+suff].Add(bkgHistos[ibkg])
                            else:
                                temp_bkgHistos['QCD'+suff] = bkgHistos[ibkg].Clone('QCD'+suff+'_bkgSubComb')
                                
                        else:
                            print(f"WARNING: unknown background (label) included in bkgHistos dict, hist name: {ibkg}" )
                print( f"Number of background events from ST: {temp_bkgHistos['ST'+suff].Integral()}")
                print( f"Number of background events from WJets: {temp_bkgHistos['WJets'+suff].Integral()}")
                print( f"Number of background events from QCD: {temp_bkgHistos['QCD'+suff].Integral()}")
                print( f"Number of background events from DY: {temp_bkgHistos['DY'+suff].Integral()}")
                print( f"Number of background events from VV: {temp_bkgHistos['VV'+suff].Integral()}")
                print( f"Data minus backgrounds counter: {dataMinusbkg_counter}")
            
            
                for ibkg in temp_bkgHistos:
                    
                    
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel):
                        dataMinusbkg_counter-=temp_bkgHistos[ibkg].Integral()
                        tunfolder.SubtractBackground( temp_bkgHistos[ibkg].Clone(), ibkg.split('_')[0],#+ '%d'%dummy,
                                                     1.,scaleUnc_bkgHistos[ibkg]  )
                        #dummy=dummy+1
                        if 'all' in year:

                            print(ibkg,temp_bkgHistos[ibkg].Integral())
                            print('Subtracted fakes and/or bkgs; integrals of fakes, true, dataMinusBkg, initial data', signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral(), signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),allHistos[ 'dataMinusBkgs'+ivar ].Integral(),allHistos[ 'dataHisto' +ivar].Integral())
                            print( f"Data minus backgrounds counter: {dataMinusbkg_counter}")
                        
                        bkgSources.append(ibkg.split('_')[0])#+ '%d'%dummy)
                
            # Subtract signal/nominal MC reco fakes
            dataMinusbkg_counter-=signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral()
            
            print('Subtracted fakes and/or bkgs; integrals of fakes, true, dataMinusBkg, initial data', signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Integral(), signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),allHistos[ 'dataMinusBkgs'+ivar ].Integral(),allHistos[ 'dataHisto' +ivar].Integral())
            print( f"Data minus backgrounds counter: {dataMinusbkg_counter}")
            #if verbose: print('subtracting fakes')
            tunfolder.SubtractBackground( signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel], 'fakes')
            
            bkgSources.append('fakes')
        
        
        
        ###### Adding SYS unc
        
        if len(sysUncert)>0 and process.startswith('data'):
            print ('|------> TUnfolding adding uncert:')
            dictUncHistos = {}
            #print (sysSignalHistos)
            for sys in sysUncert:
                #print (sys)
                
                if sys.startswith(('_jer', '_isrWeight', '_l1prefiringWeight', '_fsrWeight', '_puWeight', '_jes', '_leptonWeight', '_btagWeight', '_const', '_unclust' )) or ('pdfweight' in sys.lower()):
                    
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
                    #print(f'|------> TUnfolding adding {sys} unc. with prefix-checked {s}')
                    #if verbose: 
                    #    print(sys,s,s[0])
                    dictUncHistos[sys+'Up'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Up'+sel].Clone()
                    dictUncHistos[sys+'Down'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Down'+sel].Clone()
                    for upDown in [ 'Up', 'Down' ]:
                        
                       
                            
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
            
                #below if blocks for modelling systematics relevant only to ttbar (W/top) selections 
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

                elif (sys.startswith('_erdON') and not('_CR1' in sysUncert or '_CR2' in sysUncert)) or (sys.startswith('_CR1') and not('_erdON' in sysUncert or '_CR2' in sysUncert)) or (sys.startswith('_CR2') and not('_erdON' in sysUncert or '_CR1' in sysUncert)) or (sys.startswith('_erdON') and ('_CR1' in sysUncert and '_CR2' in sysUncert)): #elif sys.startswith('_CR'):#
                    if '_erdON' in sysUncert or '_erdON' in sys:
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
                    
                    #if verbose:
                    #    print('|------> TUnfolding adding Colour reconnection Unc')
                    
                    if '_CR1' in sysUncert or '_CR1' in sys:
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+ivar+'_nom'+sel],
                                            '_CR1',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR1', ivar+'can2DNorm'+'TuneCP5CR1', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR1'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)
                        dictUncHistos['_CR1'] = varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_reco'+ivar+'_nom'+sel].Clone()

                    if '_CR2' in sysUncert or '_CR2' in sys:
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+ivar+'_nom'+sel],
                                            '_CR2',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, 
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR2', ivar+'can2DNorm'+'TuneCP5CR2', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR2'+'_respWithMiss'+sel+'Normalized_responseMatrix'+version+'.'+ext)

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
                tempDictUncHistos = OrderedDict()
                for i in dictUncHistos.keys():
                    if (i.startswith('_jes') and not('const' in i)):
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos[i].Sumw2()
                        
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,
                                    ivar+'_'+signalLabel+'_JESAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=  year ,
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES',sysList=sysUncert
                                    )
            else:
                tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                for i in dictUncHistos.keys():
                    if (i.startswith('_jes') and not('2016' in i or '2017' in i or '2018' in i  or 'const' in i)):
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos[i].Sumw2()
                    elif i.startswith('_jes') and ('2016' in i or '2017' in i or '2018' in i):
                        #print(i)
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos2[i].Sumw2()
                
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,
                                    ivar+'_'+signalLabel+'_JESUncorrAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# '2016+2017+2018',# if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES',sysList=sysUncert
                                    )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,
                                    ivar+'_'+signalLabel+'_JESCorrAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# '2016+2017+2018',# if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJES',sysList=sysUncert
                                    )
                
                #tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                for i in dictUncHistos.keys():
                    if i.startswith('_jer'):# and ('2016' in i or '2017' in i or '2018' in i):
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos2[i].Sumw2()

                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,
                                    ivar+'_'+signalLabel+'_JERUncorrAllSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='', 
                                    mode='onlyJER',sysList=sysUncert
                                    )
               
                
            if 'dijet' in sel:
                tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                for i in dictUncHistos.keys():
                    if (i.startswith('_jes')): continue
                    
                    elif year=='all' and 'jer' in i: continue
                        
                    elif i.startswith('_const'):
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos2[i].Sumw2()
                    else:
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos[i].Sumw2()
                        
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,
                                    ivar+'_'+signalLabel+'_NoJESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode='',sysList=sysUncert
                                    )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,
                                    ivar+'_'+signalLabel+'_constESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode='',sysList=sysUncert
                                    )
            else:
                tempDictUncHistos = OrderedDict()
                tempDictUncHistos2 = OrderedDict()
                tempDictUncHistos3 = OrderedDict()
                for i in dictUncHistos.keys():
                    if (i.startswith('_jes')): continue
                    
                    elif year=='all' and 'jer' in i: continue
                    
                    elif i.startswith('_const'):
                        tempDictUncHistos3[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos3[i].Sumw2()
                        
                    elif i.startswith(('_model', '_unclust', '_isrWeight', '_l1prefiringWeight', '_fsrWeight', '_puWeight', '_leptonWeight', '_btagWeight')) or ('pdfweight' in sys.lower()):
                        tempDictUncHistos[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos[i].Sumw2()
                    elif i.startswith(('_mtop','_CR','_hdamp','_Tune','_erd')):
                        tempDictUncHistos2[i] = dictUncHistos[i].Clone()
                        tempDictUncHistos2[i].Sumw2()

                     
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos,#dictUncHistos,
                                    ivar+'_'+signalLabel+'_NoJESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode='',sysList=sysUncert
                                  )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos2,#dictUncHistos,
                                    ivar+'_'+signalLabel+'_TheorySys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode='',sysList=sysUncert
                                  )
                plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    tempDictUncHistos3,#dictUncHistos,
                                    ivar+'_'+signalLabel+'_constESSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+version,
                                    year=year,# ( '2016+2017+2018' if year.startswith('all') else year ),
                                    outputDir=outputDir,
                                    #sys_pref='_je', 
                                    mode='',sysList=sysUncert
                                  )
                del(tempDictUncHistos,tempDictUncHistos2,tempDictUncHistos3)
                   
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
                                                                         "CM from stat. unc. of Input Distribution")
        
        if 'cross' in process.lower():
            allHistos[ 'cov_cross'+ivar ] = tunfolder_cross.GetEmatrixTotal("cov_cross"+ivar, "Total Covariance Matrix Alt RM")
        
            allHistos[ 'cov_cross_uncorr_data_'+ivar ] = tunfolder_cross.GetEmatrixInput("cov_cross_uncorr_data"+ivar,
                                                                         "CM from stat. unc. of Input Distribution Alt RM")
                
        allHistos[ 'cov_uncorr_'+ivar ] = tunfolder.GetEmatrixSysUncorr("cov_uncorr"+ivar, 
                                                                        "CM from stat. unc. from response matrix")
               
        allHistos[ 'cov_uncorr_bkg_'+ivar ] = tunfolder.GetEmatrixSysBackgroundUncorr('fakes', 
                                                                                      "CM from Uncorrelated Errors of Background Sources")
        
        
        allHistos[f'cov_dataAndBkgs{ivar}' ] = allHistos[f'cov_uncorr_data_{ivar}'].Clone("cov_dataAndBkgs"+ivar)
        allHistos[f'cov_dataAndBkgs{ivar}' ].Reset()
        #if sel.startswith(('_W','_top')):
        tunfolder.GetEmatrix(allHistos[f'cov_dataAndBkgs{ivar}' ])
        allHistos[f'cov_dataAndBkgs{ivar}' ].SetTitle("CM from stat uncertainties on data and background stat.+scale")
      
        
        if process.startswith('data'):
            if sel.startswith(('_W','_top')):
                print("Adding bkg. uncorr.  to cov for fakes")
                
                for ibkg in bkgSources:
                    print(ibkg)
                    if 'fakes' not in ibkg: allHistos[ 'cov_uncorr_bkg_'+ivar ].Add(tunfolder.GetEmatrixSysBackgroundUncorr(ibkg, "CM from Uncorrelated Errors of Background Source "+ibkg))
        

        #### cov total = cov_uncorr + cov_uncorr_data + cov_uncorr_bkg
        #### cov stat = cov_uncorr + cov_uncorr_data
        
        allHistos[ 'unfoldHistowoUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone()        # No unc

        for ibin in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):

            allHistos[ 'unfoldHistowoUnc'+ivar ].SetBinError(ibin, 0. )        # No unc
        

        
        
                
        ########################### Get systematic shifts of output ################################
        #storing individual systematics totals (ie, up/down or other variations), background subtraction systematics, 
        #and then the overall systematic unc (individual systs + bkgs)
        if len(sysUncert)>0 and process.startswith('data'):
            print ('|------> TUnfolding uncertainties breakdown:', sysUncert)

            def handle_systematics_multiSource(ivar, sys_prefix, allHistos, sysUncert, uncerUnfoldHisto, uncerUnfoldSystCov):
                """
                Handle systematic variations for a given prefix (e.g., _btag, _jer).
                
                Args:
                    ivar (str): Variable identifier for histograms.
                    sys_prefix (str): Prefix for the systematic category (e.g., '_btag', '_jer').
                    tunfolder: TUnfold object to retrieve systematic uncertainties.
                    allHistos (dict): Dictionary of all histograms.
                    sysUncert (list of str): List of systematic uncertainty names.
                    uncerUnfoldHisto (dict): Dictionary to store unfolded histograms for uncertainties.
                    uncerUnfoldSystCov (dict): Dictionary to store covariance matrices for uncertainties.
                """
                
                sysUncs = [s for s in sysUncert if s.startswith(sys_prefix)]
                print("Adding syst. with multiple sources:", sysUncs)
                sysUncsUp = [s + 'Up' for s in sysUncs]
                sysUncsDown = [s + 'Down' for s in sysUncs]
                for sys in sysUncs:#sysUncsUp + sysUncsDown:
                    sys_cov_up = tunfolder.GetEmatrixSysUncorr("cov_%s_Up"%sys+ivar).Clone()
                    sys_cov_down = tunfolder.GetEmatrixSysUncorr("cov_%s_Down"%sys+ivar).Clone()
                    tunfolder.GetEmatrixSysSource(sys_cov_up, sys+'Up')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, sys+'Down')
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Up'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'] = sys_cov_down.Clone()
                    
                    for upDown in ['Up','Down']:
                        
                        uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "+1#sigma" if 'up' in upDown.lower() else "-1#sigma")
                        tmpHisto = allHistos [ 'unfoldHistowoUnc'+ivar ].Clone("unfoldHisto_"+ivar+sys+upDown+"shift_tmpClone")
                        tmpHisto.Reset()
                        
                        try: 
                            uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                        except ReferenceError: 
                            print(f"{sys+upDown} has no effect", uncerUnfoldHisto[ivar+sys+upDown], type(uncerUnfoldHisto[ivar+sys+upDown]))
                            #uncerUnfoldHisto.pop( ivar+sys+upDown, None )
                            uncerUnfoldHisto[ivar+sys+upDown] = tmpHisto.Clone()
                            del tmpHisto
                            gc.collect()

                        #if uncerUnfoldHisto[ivar+sys+upDown]: 
                        uncerUnfoldHisto[ivar+sys+upDown+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+sys+upDown].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                        

                    del sys_cov_up
                    del sys_cov_down
                    gc.collect()

                total_hist = allHistos[f'unfoldHisto{ivar}'].Clone(f'{ivar}{sys_prefix}Total')
                total_hist.Reset()
                up_hist = allHistos[f'unfoldHisto{ivar}'].Clone(f'{ivar}{sys_prefix}Up')
                up_hist.Reset()
                down_hist = allHistos[f'unfoldHisto{ivar}'].Clone(f'{ivar}{sys_prefix}Down')
                down_hist.Reset()

                for i in range(1, total_hist.GetNbinsX() + 1):
                    dy_up_temp = 0
                    dy_down_temp = 0
                    #print(i)
                    for s in sysUncsUp:
                        #print(s, i, uncerUnfoldHisto[f'{ivar}{s}'].GetBinContent(i))#,allHistos[f'unfoldHisto{ivar}'].GetBinContent(i) )
                        try:
                            dy_up_temp += uncerUnfoldHisto[f'{ivar}{s}'].GetBinContent(i)**2
                        except KeyError:
                            print(f'{ivar}{s} not found; setting up shift to 0')
                            pass
                    
                    for s in sysUncsDown:
                        #print(s, i, uncerUnfoldHisto[f'{ivar}{s}'].GetBinContent(i))#,allHistos[f'unfoldHisto{ivar}'].GetBinContent(i) )
                        try:
                            dy_down_temp += uncerUnfoldHisto[f'{ivar}{s}'].GetBinContent(i)**2
                        except KeyError:
                            print(f'{ivar}{s} not found; setting down shift to 0')
                            pass
                    
                    dy_up_temp = ROOT.TMath.Sqrt(dy_up_temp)
                    dy_down_temp = ROOT.TMath.Sqrt(dy_down_temp)
                    dy_tot = ROOT.TMath.Sqrt(dy_up_temp**2+dy_down_temp**2)#max(dy_up_temp, dy_down_temp)
                    
                    #print(sys_prefix,i, dy_up_temp, dy_down_temp, dy_tot)
                    
                    total_hist.SetBinContent(i, dy_tot)
                    up_hist.SetBinContent(i, dy_up_temp)
                    down_hist.SetBinContent(i, dy_down_temp)
                    

                uncerUnfoldHisto[f'{ivar}{sys_prefix}Total'] = total_hist.Clone()
                uncerUnfoldHisto[f'{ivar}{sys_prefix}Up'] = up_hist.Clone()
                uncerUnfoldHisto[f'{ivar}{sys_prefix}Down'] = down_hist.Clone()
                
                uncerUnfoldHisto[f'{ivar}{sys_prefix}Total_shiftHist'] = get_syst_shifted_hist(
                    total_hist.Clone(),
                    allHistos[f'unfoldHistowoUnc{ivar}'].Clone()
                )
                
                uncerUnfoldHisto[f'{ivar}{sys_prefix}Up_shiftHist'] = get_syst_shifted_hist(
                    up_hist.Clone(),
                    allHistos[f'unfoldHistowoUnc{ivar}'].Clone()
                )
                
                uncerUnfoldHisto[f'{ivar}{sys_prefix}Down_shiftHist'] = get_syst_shifted_hist(
                    down_hist.Clone(),
                    allHistos[f'unfoldHistowoUnc{ivar}'].Clone()
                )
            
            for sys_prefix in ['_jes'] + ([ '_btag' ] if not('dijet' in sel) else [] ) + ([ '_jer' ] if 'all' in year else [] ):
                #print(sys_prefix,sysUncert)
                handle_systematics_multiSource(ivar, sys_prefix, allHistos, sysUncert, uncerUnfoldHisto, uncerUnfoldSystCov)
            
            #### handle background rate uncertainties
            bkgRateCov = tunfolder.GetEmatrixSysUncorr(f"cov_bkgRates_{ivar}_dummy").Clone(f"cov_bkgRates_{ivar}")
            bkgRateCov.Reset()
            allHistos[f'cov_bkgRate_{ivar}' ] = bkgRateCov.Clone()

            bkgRateCovs = OrderedDict()
            if process.startswith('data'):
                if sel.startswith(('_W','_top')):
                    for ibkg in bkgSources:
                        if 'fake' not in ibkg: 
                            print(ibkg)
                            thisBkgCov = bkgRateCov.Clone(bkgRateCov.GetName()+f'_{ibkg}')
                            thisBkgCov.Reset()
                            #thisBkgDelta = allHistos[f'unfoldHisto{ivar}'].Clone('deltaSysBkgRate_'+ibkg)
                            #thisBkgDelta.Reset()

                            tunfolder.GetEmatrixSysBackgroundScale(thisBkgCov, ibkg)
                            thisBkgCov.SetTitle("CM from Rate Errors of Background Source "+ibkg)
                            allHistos[f'cov_bkgRate_{ivar}' ].Add(thisBkgCov)

                            uncerUnfoldHisto[ivar+ibkg+'_BkgRateUp'] = tunfolder.GetDeltaSysBackgroundScale(ibkg, 
                                                                                                                    "unfoldHisto_"+ivar+ibkg+"BkgRateUp_shift", 
                                                                                                                    ibkg+'BkgRateDeltaUp') #tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "+1#sigma" if 'up' in upDown.lower() else "-1#sigma")
                            uncerUnfoldHisto[ivar+ibkg+'_BkgRateUp'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+ibkg+'_BkgRateUp'].Clone(), 
                                                                                                                  allHistos[f'unfoldHistowoUnc{ivar}' ].Clone())

                            uncerUnfoldHisto[ivar+ibkg+'_BkgRateDown'] = tunfolder.GetDeltaSysBackgroundScale(ibkg, 
                                                                                                                    "unfoldHisto_"+ivar+ibkg+"BkgRateDown_shift", 
                                                                                                                    ibkg+'BkgRateDeltaDown') #tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "+1#sigma" if 'up' in upDown.lower() else "-1#sigma")
                            uncerUnfoldHisto[ivar+ibkg+'_BkgRateDown'].Scale(-1)
                            uncerUnfoldHisto[ivar+ibkg+'_BkgRateDown'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+ibkg+'_BkgRateDown'].Clone(), 
                                                                                                                  allHistos[f'unfoldHistowoUnc{ivar}' ].Clone())



            
            
            for sys in sysUncert:
                #print(sys,"cov_%s_Up"%sys+ivar,"cov_%s_Down"%sys+ivar)
                if sys.startswith(('_jes', '_btag', '_lepton')): 
                    continue
                if 'all' in year and ('_jer' in sys): 
                    continue
                #print(sys)
                sys_cov_up = tunfolder.GetEmatrixSysUncorr(f"cov_{sys}_Up{ivar}").Clone(f"cov_{sys}_Up{ivar}_dummy")#dummies
                sys_cov_down = tunfolder.GetEmatrixSysUncorr(f"cov_{sys}_Down{ivar}").Clone(f"cov_{sys}_Down{ivar}_dummy")#dummies
                
                if not sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')):
                    
                    tunfolder.GetEmatrixSysSource(sys_cov_up, sys+'Up')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, sys+'Down')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Up'] = sys_cov_up.Clone('systcov_'+ivar+sys+'Up_clone')
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'] = sys_cov_down.Clone('systcov_'+ivar+sys+'Down_clone')
                    
                    for upDown in [ 'Up', 'Down' ]:
                        
                        uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "+1#sigma" if 'up' in upDown.lower() else "-1#sigma")

                        uncerUnfoldHisto[ivar+sys+upDown+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+sys+upDown].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                        

                    # Create total uncertainty and sys uncertainty histos for plotting/further calculations.
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone() #total shifts not used but kept in case
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                    #print(f"Doing max of up/down for 'total' estimate of {sys} unc source")
                    for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                        try: yup = abs(  uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs(  uncerUnfoldHisto[ivar+sys+'Down'].GetBinContent(i) )
                        except KeyError: ydn = 0
                                                        
                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)#np.max([yup,ydn])#
                        #dy = 0.5*( (yup + ydn) ) #being conservative and doing an envelope based on max of up/down shift, per bin, of a systematic
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )
                    
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'+"_shiftHist"] = get_syst_shifted_hist( uncerUnfoldHisto[ivar+sys.upper()+'Total'].Clone(), allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                    #if verbose:
                    print ("Done with adding %s to tot syst unc"%sys )    

                elif sys.startswith('_model'):
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_modelUncTotal')
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Reset()
                    
                    tmpModelHisto = tunfolder.GetDeltaSysSource('modelUncTotal', "unfoldHisto_"+ivar+"_modelUncTotalshift", "modelUncShift")
                    tunfolder.GetEmatrixSysSource(sys_cov_up, "modelUncTotal")
                    uncerUnfoldSystCov['systcov_'+ivar+'_Physics ModelTotal'] = sys_cov_up.Clone('systcov_'+ivar+'_Physics ModelTotal_clone')
                    
                    #for i in range( 1, tmpModelHisto.GetNbinsX() + 1):
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = tmpModelHisto.Clone()#.SetBinContent( i, (tmpModelHisto.GetBinContent(i)))#abs(
                    
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


                elif sys.startswith('_erdON') and ('_CR1' in sysUncert and '_CR2' in sysUncert):#'_CR',
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
                        #dy = np.max([yup,ydn,ymid]) #considering CR unc as an envelope from max contribution per bin from CR variation sources
                        dy = np.max([yup,ydn,ymid]) #considering CR unc as an envelope from max contribution per bin from CR variation sources
                        #ROOT.TMath.Sqrt(yup**2+ydn**2+ymid**2)
                        uncerUnfoldHisto[ivar+'_CRTotal'].SetBinContent( i, dy ) # considering uncertainty as an envelope
                    uncerUnfoldHisto[ivar+'_CRTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_CRTotal'].Clone(),
                                                                                                    allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                
                elif (sys.startswith('_erdON') and not('_CR1' in sysUncert or '_CR2' in sysUncert)) or (sys.startswith('_CR1') and not('_erdON' in sysUncert or '_CR2' in sysUncert)) or (sys.startswith('_CR2') and not('_erdON' in sysUncert or '_CR1' in sysUncert)):#'_CR',
                                   
                    CR = sys
                    tunfolder.GetEmatrixSysSource(sys_cov_up, CR)
                    
                    uncerUnfoldHisto[ivar+'_CRTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_CRTotal')
                    uncerUnfoldHisto[ivar+'_CRTotal'].Reset()
                    uncerUnfoldHisto[ivar+CR+"Total"] = tunfolder.GetDeltaSysSource(CR, "unfoldHisto_" + ivar + "%sshift"%CR, f"{CR}UncShift")
                    
                    if uncerUnfoldHisto[ivar+CR+"Total"]: 
                        uncerUnfoldHisto[ivar+CR+"Total_shiftHist"] = get_syst_shifted_hist( 
                                                                                        uncerUnfoldHisto[ivar+CR+"Total"].Clone(),
                                                                                        allHistos [ 'unfoldHistowoUnc'+ivar ].Clone()
                                                                                      )

                        uncerUnfoldSystCov['systcov_'+ivar+'_CRTotal'] = sys_cov_up.Clone('systcov_'+ivar+CR+'Total_clone')
                    
                                    
                elif sys.startswith('_mtop'):
                    tunfolder.GetEmatrixSysSource(sys_cov_up, '_mtop173p5')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, '_mtop171p5')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+'_mtop171p5'] = sys_cov_down.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+'_mtop173p5'] = sys_cov_up.Clone()
                    
                    for m,l in zip(mass_list,['Down','Up']):
                        #print ("mtop:", m,l)
                        uncerUnfoldHisto[ivar+'_mtop%s'%m] = tunfolder.GetDeltaSysSource('_mtop%s'%m, "unfoldHisto_"+ivar+"_mtop%sshift"%m, "+1#sigma" if '173' in m.lower() else "-1#sigma")
                        if uncerUnfoldHisto[ivar+'_mtop%s'%m]: uncerUnfoldHisto[ivar+'_mtop%s'%l+'_shiftHist'] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_mtop%s'%m].Clone(),
                                                                                                  allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())
                    uncerUnfoldHisto[ivar+'_Top MassTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_mtop171p5')
                    uncerUnfoldHisto[ivar+'_Top MassTotal'].Reset()

                    for i in range( 0, uncerUnfoldHisto[ivar+'_Top MassTotal'].GetNbinsX() + 1):
                        yup = (abs( uncerUnfoldHisto[ivar+'_mtop173p5'].GetBinContent(i) ))
                        ydn = (abs( uncerUnfoldHisto[ivar+'_mtop171p5'].GetBinContent(i) ))
                        #dy = np.max([yup,ydn])#
                        dy = ROOT.TMath.Sqrt(yup**2+ydn**2)
                        uncerUnfoldHisto[ivar+'_Top MassTotal'].SetBinContent( i, dy ) # considering uncertainty as an envelope
                    uncerUnfoldHisto[ivar+'_Top MassTotal'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_Top MassTotal'].Clone(),
                                                                                                    allHistos [ 'unfoldHistowoUnc'+ivar ].Clone())

                del(sys_cov_up,sys_cov_down)
                gc.collect()
        
        ############### Build correlation matrix for unfolding#################################
        allHistos['tunf_rhoIJ_correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone('tunf_rhoIJ_correlation_matrix_'+ivar)  
        allHistos['tunf_rhoIJ_correlation_matrix_'+ivar].Reset()
        tunfolder.GetRhoIJ(allHistos['tunf_rhoIJ_correlation_matrix_'+ivar])
        allHistos['correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone('correlation_matrix_'+ivar)
        allHistos['correlation_matrix_'+ivar].Reset()
        allHistos['correlation_matrix_'+ivar] = correlation_from_covariance(allHistos[ 'cov'+ivar ].Clone(),allHistos['correlation_matrix_'+ivar])
        
        draw2D( ivar,  allHistos['correlation_matrix_'+ivar].Clone(), variables[ivar], outputLabel='Un-normed_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)
        
                
        ########################################################################################    
        
        # Create total uncertainty and sys uncertainty histos for plots 
        
        # first build up histos of systematic uncertainties from the background subtraction
        uncerUnfoldHisto[ivar+'_BkgTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_BkgTotal')
        uncerUnfoldHisto[ivar+'_BkgTotal'].Reset()

        for ibin in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):
            bkg_tot = np.sqrt(allHistos[ 'cov_uncorr_bkg_'+ivar ].GetBinContent(ibin,ibin) + (allHistos[ 'cov_bkgRate_'+ivar ].GetBinContent(ibin,ibin) if 'data' in process and not('dijet' in sel) else 0.) )
            uncerUnfoldHisto[ivar+'_BkgTotal'].SetBinContent(ibin, bkg_tot)
            
        # second build up histos of total systematic uncertainties
        uncerUnfoldHisto[ivar+'_SystTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_SystTotal')
        uncerUnfoldHisto[ivar+'_SystTotal'].Reset()        
        
       
        ###### adding covariances from bkg subtraction and RM finite stats to the overall systematics covariance matrix
        allHistos['cov_systTotal'+ivar] = allHistos[ 'cov_uncorr_'+ivar ].Clone('cov_uncorr_bkg_'+ivar+'+cov_uncorr'+ivar)
        #allHistos['cov_systTotal'+ivar].Reset()
        allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_bkg_'+ivar ])    
        if 'data' in process and not('dijet' in sel): allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_bkgRate_'+ivar ])
        #array2hist(systcovTotal,allHistos['cov_systTotal'+ivar])
                              
        
        tmp = OrderedDict()
        for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
            tmp[i] = 0
            for k in uncerUnfoldHisto:
                if k.endswith('Total') and not k.endswith(('SystTotal')):
                    tmp[i] = tmp[i] + ( uncerUnfoldHisto[k].GetBinContent( i )**2 )
                    
        
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
        
        #unnormalised unfolded histo with unnormalised cov. unc from relevant covs.
        #allHistos[ 'unfoldHistowoUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone()        # No unc
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
            bkg_tot = np.sqrt(allHistos[ 'cov_uncorr_bkg_'+ivar ].GetBinContent(ibin,ibin) + (allHistos[ 'cov_bkgRate_'+ivar ].GetBinContent(ibin,ibin) if 'data' in process and not('dijet' in sel) else 0.))
            datastat_tot = np.sqrt( allHistos[ 'cov_dataAndBkgs'+ivar ].GetBinContent(ibin,ibin))#np.sqrt( allHistos[ 'cov_uncorr_data_'+ivar ].GetBinContent(ibin,ibin))
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
            #allHistos[ 'unfoldHistowoUnc'+ivar ].SetBinError(ibin, 0. )        # No unc
        
            uncerUnfoldHisto[ivar+'_TotalUnc'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMErrTotal'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMMCStatErrTotal'].SetBinContent(ibin, rmstat_tot)
            uncerUnfoldHisto[ivar+'_CMDataStatErrTotal'].SetBinContent(ibin, datastat_tot)
            
            ratioHistos[ 'TotalUnc'+ivar ].SetBinContent( ibin, 1. )
            ratioHistos[ 'StatUnc'+ivar ].SetBinContent( ibin, 1. )

            if norm_unc!=0:
                ratioHistos[ 'TotalUnc'+ivar ].SetBinError( ibin, unc_tot)#+syst_tot**2
                ratioHistos[ 'StatUnc'+ivar ].SetBinError( ibin, datastat_tot)
        
        ################################################################################

        # Get folded distribution from unfolded distribution
        allHistos [ 'foldHisto'+ivar ] = tunfolder.GetFoldedOutput("folded"+ivar).Clone('foldHisto_fromTUnfold'+ivar ) 
        allHistos [ 'foldHisto2'+ivar ] = get_folded_unfolded(
                                                        
                                                        folded=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                                                        unfolded=allHistos['unfoldHisto'+ivar].Clone(), 
                                                        cov_tot=allHistos['cov'+ivar].Clone(), 
                                                        probaM=allHistos[ 'probaMatrix'+ivar] 
                                                                   
                                                                    ).Clone('foldHisto_ByHand'+ivar )
        
        #if verbose:
        print(f"native folded output integral, true-reco Integral, data-bkg Integral, custom folded output integral, gen-level (miss+accep) integral, unfolded (data-bkg) integral")

        print(tunfolder.GetFoldedOutput("folded"+ivar).Integral(), 
              signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Integral(),
              allHistos[ 'dataMinusBkgs' +ivar].Integral(),                  
              allHistos [ 'foldHisto2'+ivar ].Integral(),
              signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Integral(), 
              allHistos['unfoldHisto'+ivar].Integral())

        if process.startswith('data'): 
            plotSimpleComparison( allHistos[ 'dataMinusBkgs'+ivar ].Clone(), 'data-Bkgs',  allHistos [ 'foldHisto'+ivar ].Clone(), 'folded', 
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
                       dataJetHisto=allHistos[ 'dataMinusBkgs'+ivar ].Clone(),
                       genJetHisto=signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),#=signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionX('genJetHisto_fromProjX', 0, signalHistos[ signalLabel+'_respWithMiss'+ivar+'_nom'+sel].GetNbinsY()+1).Clone(),
                       #unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                       #unfoldHistoStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(),
                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                       foldHisto=tunfolder.GetFoldedOutput("folded"+ivar).Clone(), recoJetHisto=signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY().Clone(),
                       cov_datastat_tot= allHistos['cov_dataAndBkgs'+ivar].Clone(),  #allHistos['cov_uncorr_data_'+ivar].Clone(),
                       cov_tot=allHistos['cov'+ivar].Clone(),
                       altMCHisto =  altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),#altMCHisto =  altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionX('altMCHisto_fromProjX', 0, altSignalHistos[altSignalLabel+'_respWithMiss'+ivar+'_nom'+sel].GetNbinsY()+1).Clone(),
                       labelX=variables[ivar]['label'],
                       maxX=variables[ivar]['bins'][-1],
                       tlegendAlignment=variables[ivar]['alignLeg'],
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_'+version+'.'+ext,
                       altMC1Histo = alt1SignalHistos[alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Clone(), 
                       altMC2Histo = alt2SignalHistos[alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Clone() if 'dijet' in sel else None, 
                       altMC1Histo_label = alt1SigPlotLabel,  
                       altMC2Histo_label = alt2SigPlotLabel if 'dijet' in sel else None,
                       nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel,
                       extraMC=extraMC,
                       includeFSR = True if include_FSR_in_unfolded_result else False,
                       fsrUpHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Clone() if include_FSR_in_unfolded_result else None,
                       fsrDownHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Clone() if include_FSR_in_unfolded_result else None,
                       )
            
            drawUnfold(ivar=ivar, 
                       selection=sel, year=year,lumi=lumi, process=process,
                       dataJetHisto=allHistos[ 'dataMinusBkgs'+ivar ].Clone(),
                       genJetHisto=signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       #unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                       #unfoldHistoStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(),
                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                       foldHisto=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                       recoJetHisto=signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].ProjectionY().Clone(),
                       cov_datastat_tot= allHistos['cov_dataAndBkgs'+ivar].Clone(),  #allHistos['cov_uncorr_data_'+ivar].Clone(),
                       cov_tot=allHistos['cov'+ivar].Clone(),
                       altMCHisto=  altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                       labelX=variables[ivar]['label'],
                       maxX=variables[ivar]['bins'][-1],
                       tlegendAlignment=variables[ivar]['alignLeg'],
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_NO_NORM'+version+'.'+ext,
                       altMC1Histo = alt1SignalHistos[alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Clone(), 
                       altMC2Histo = alt2SignalHistos[alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Clone() if 'dijet' in sel else None, 
                       altMC1Histo_label = alt1SigPlotLabel,  
                       altMC2Histo_label = alt2SigPlotLabel if 'dijet' in sel else None,
                       nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel,
                       extraMC=extraMC,
                       includeFSR = True if include_FSR_in_unfolded_result else False,
                       fsrUpHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Clone() if include_FSR_in_unfolded_result else None,
                       fsrDownHisto = sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Clone() if include_FSR_in_unfolded_result else None,
                       noNorm=True
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
                             cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),
                             cov_tot=allHistos['cov'+ivar].Clone(),
                             cov_datastat_tot_cross=allHistos['cov_cross_uncorr_data_'+ivar].Clone(),
                             cov_tot_cross=allHistos['cov_cross'+ivar].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext,
                             nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel, noNorm=False

                             )
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHistoCross=allHistos['unfoldHistoCross'+ivar ].Clone(),
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),
                             cov_tot=allHistos['cov'+ivar].Clone(),
                             cov_datastat_tot_cross=allHistos['cov_cross_uncorr_data_'+ivar].Clone(),
                             cov_tot_cross=allHistos['cov_cross'+ivar].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_NO_NORM_'+version+'.'+ext,
                             nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel, noNorm=True
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
                             cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),
                             cov_tot=allHistos['cov'+ivar].Clone(),
                             cov_datastat_tot_cross=None,#allHistos['cov_cross_uncorr_data_'+ivar].Clone(),
                             cov_tot_cross=None,#allHistos['cov_cross'+ivar].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_'+version+'.'+ext,
                             nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel, noNorm=False
                             )
                drawClosures(ivar=ivar, selection=sel, year=year, lumi=lumi, process=process,
                             genJetHistoCross=[], 
                             unfoldHistoCross=[] ,
                             genJetHisto=signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone(),
                             unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                             ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                             ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                             ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                             cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),
                             cov_tot=allHistos['cov'+ivar].Clone(),
                             cov_datastat_tot_cross=None,#allHistos['cov_cross_uncorr_data_'+ivar].Clone(),
                             cov_tot_cross=None,#allHistos['cov_cross'+ivar].Clone(),
                             labelX=variables[ivar]['label'],
                             maxX=variables[ivar]['bins'][-1],
                             tlegendAlignment=variables[ivar]['alignLeg'],
                             outputName=outputDir+ivar+sel+'_from'+process+signalLabel+'_TUnfold_NO_NORM_'+version+'.'+ext,
                             nomMCHisto_label = sigPlotLabel, altMCHisto_label = altSigPlotLabel, noNorm=True
                             )
                
                
        

                
        ######### Plotting Uncertainties
        print ('|------> Drawing unfolding uncertainty plot:')
        tempuncerUnfoldHisto = OrderedDict()
        if not('dijet' in sel): 
            modelVarnUncerUnfoldHisto = OrderedDict()
            
        for u in uncerUnfoldHisto.keys():
            if ('_hdamp' in u or 'tune' in u.lower() or '_erdON' in u or '_CR1' in u or '_CR2' in u or '_mtop' in u) and not('dijet' in sel): 
                print(u)
                modelVarnUncerUnfoldHisto[u] = uncerUnfoldHisto[u].Clone(uncerUnfoldHisto[u].GetName()+'_compare_Normed')
            else:
                tempuncerUnfoldHisto[u] = uncerUnfoldHisto[u].Clone(uncerUnfoldHisto[u].GetName()+'_compare_Normed')
        
            
        if not 'Closure' in process and 'tau' in ivar: 
            
            doRelUncPlot(ivar,
                         year,
                         lumi,
                         sel,
                         variables,
                         allHistos, uncerUnfoldSystCov,
                         outputDir,
                         version=version,
                         ext='pdf'
                        )
            doRelUncPlot(ivar,
                         year,
                         lumi,
                         sel,
                         variables,
                         allHistos, uncerUnfoldSystCov,
                         outputDir,
                         version=version,
                         ext='png'
                        )
            
            
            """
            drawUncertainties_from_err_shifts( ivar=ivar, 
                                               unfoldHistoTotUnc=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                                               unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                                               unfoldHistoDataStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(), 
                                               unfoldHistoRMStatUnc=allHistos[ 'unfoldHistoRMUnc'+ivar ].Clone(),
                                               unfoldHistoBkgSubUnc=allHistos[ 'unfoldHistoBkgUnc'+ivar ].Clone(), 
                                               uncerUnfoldHisto=tempuncerUnfoldHisto, 
                                               cov_tot=allHistos['cov'+ivar].Clone(), 
                                               cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                               cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                               cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                               labelX=variables[ivar]['label'], 
                                               tlegendAlignment=variables[ivar]['alignLeg'],
                                               outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_UNC_'+version+'.'+ext,
                                               year=year, 
                                               unftot=unfoldingtot,selection=sel,
                                               norming=True
                                                )
            """
            drawUncertainties_from_err_shifts_unitNorm(ivar=ivar, 
                                                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                                                       uncerUnfoldHisto=tempuncerUnfoldHisto, 
                                                       cov_tot=allHistos['cov'+ivar].Clone(), 
                                                       cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                                       cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                                       cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                                       labelX=variables[ivar]['label'], 
                                                       tlegendAlignment=variables[ivar]['alignLeg'],
                                                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_UNC_UnitNorm_'+version+'.'+ext,
                                                       year=year, 
                                                       selection=sel,
                                                       norming=True, lumi=lumi
                                )
            
            if not('dijet' in sel):
                """
                drawUncertainties_from_err_shifts_theoryVariations( ivar=ivar, 
                                                   unfoldHistoTotUnc=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                                                   unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                                                   unfoldHistoDataStatUnc=allHistos[ 'unfoldHistoStatUnc'+ivar ].Clone(), 
                                                   unfoldHistoRMStatUnc=allHistos[ 'unfoldHistoRMUnc'+ivar ].Clone(),
                                                   unfoldHistoBkgSubUnc=allHistos[ 'unfoldHistoBkgUnc'+ivar ].Clone(), 
                                                   uncerUnfoldHisto=modelVarnUncerUnfoldHisto, 
                                                   cov_tot=allHistos['cov'+ivar].Clone(), 
                                                   cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                                   cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                                   cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                                   labelX=variables[ivar]['label'], 
                                                   tlegendAlignment=variables[ivar]['alignLeg'],
                                                   outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_TheoryVariationUNC_'+version+'.'+ext,
                                                   year=year, 
                                                   unftot=unfoldingtot,selection=sel,
                                                   norming=True
                                                    )
                """
                drawUncertainties_from_err_shifts_theoryVariations_unitNorm(ivar=ivar, 
                                                                            unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                                                                            uncerUnfoldHisto=modelVarnUncerUnfoldHisto, 
                                                                            cov_tot=allHistos['cov'+ivar].Clone(), 
                                                                            cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                                                            cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                                                            cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                                                            labelX=variables[ivar]['label'], 
                                                                            tlegendAlignment=variables[ivar]['alignLeg'],
                                                                            outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_TheoryVariationUNC_UnitNorm_'+version+'.'+ext,
                                                                            year=year, 
                                                                            selection=sel,
                                                                            norming=True, lumi=lumi
                                                    )
            
            

        ######### Plotting 2D matrices of various kinds
        print ('|------> Drawing various 2D matrices:')
        if process.startswith('data'):
            #draw2D( ivar,  tunfolder.GetRhoItotal("rhoI"+ivar, "Global correlations").Clone(), variables[ivar], outputLabel='data_rhoI', outputDir=outputDir,selection=sel,version=version,year=year)
            
            
            
            draw2D( ivar,  allHistos[ 'correlation_matrix_'+ivar ].Clone(), variables[ivar], outputLabel='data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year,pngToo=True)
            draw2D( ivar, allHistos[ 'cov'+ivar].Clone(), variables[ivar], outputLabel='dataTotal_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year,pngToo=True)
            draw2D( ivar, allHistos[ 'cov_uncorr_'+ivar].Clone(), variables[ivar], outputLabel='uncorrUncRM_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_uncorr_data_'+ivar].Clone(), variables[ivar], outputLabel='dataInpStats_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year)
            
            draw2D( ivar, allHistos[ 'cov_dataAndBkgs'+ivar].Clone(), variables[ivar], outputLabel='dataAndBkgStats_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year)
            
            draw2D( ivar, allHistos[ 'cov_uncorr_bkg_'+ivar].Clone(), variables[ivar], outputLabel='BkgSubtractionSyst_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year)
            draw2D( ivar, allHistos[ 'cov_systTotal'+ivar].Clone(), variables[ivar], outputLabel='Syst_covMatrix', outputDir=outputDir, addCorrelation=False,selection=sel,version=version,year=year)
        
            draw2D( ivar,  allHistos[ 'probaMatrix'+ivar ].Clone(), variables[ivar], outputLabel='data_probaMatrix', outputDir=outputDir, addCorrelation=True, addCondition=True ,selection=sel,version=version,year=year,pngToo=True)
            dict_condition_numbers[ivar] = get_condition_number(allHistos[ 'probaMatrix'+ivar ].Clone())

            draw2D( ivar, signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone(), variables[ivar], outputLabel='data_respMatrix', outputDir=outputDir, addCorrelation=False ,selection=sel,version=version,year=year,pngToo=True)
            
        normed_covs = OrderedDict()
        tot = allHistos[ 'unfoldHisto'+ivar ].Integral()
        
        unfTemp = allHistos[ 'unfoldHisto'+ivar ].Clone(allHistos[ 'unfoldHisto'+ivar ].GetName()+'_covNormingTemp')
        
        for i in allHistos:
            if 'cov' in i: 
                if 'dijet' in sel:
                    if ('bkgrate' in i.lower()):
                        continue
                cov_normed_np, normed_cov = get_normalised_cov(unfTemp, 
                                                               allHistos[ i ].Clone())
                normed_covs['Normed'+i] = normed_cov.Clone(i+'_normed')

                normed_BW_cov = scale_th2_by_bin_width(normed_covs['Normed'+i].Clone())
                
                if 'cov'+ivar in i:
                    ah = allHistos[ i ].Clone(i+ 'ah')
                    ah.Reset()
                    
                    ah_divBy_BW = allHistos[ i ].Clone(i+ 'ah_divByBW')
                    ah_divBy_BW.Reset()
                    
                    ah = correlation_from_covariance( normed_cov.Clone(),
                                                      ah)
                    ah_divBy_BW = correlation_from_covariance( normed_BW_cov.Clone(),
                                                      ah_divBy_BW)
                    
                    if process.startswith('data'):

                        draw2D( ivar,  ah.Clone(), variables[ivar], outputLabel='Normed_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)
                        draw2D( ivar,  ah.Clone(), variables[ivar], outputLabel='Normed_divByBW_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)

                normed_covs['Normed_divByBW'+i] = normed_BW_cov.Clone()
                
                
                if process.startswith('data'):
                    
                    draw2D( ivar, allHistos[i].Clone(), variables[ivar], outputLabel='Un-normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                
                    draw2D( ivar, normed_covs['Normed'+i].Clone(), variables[ivar], outputLabel='Normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                    
                    draw2D( ivar, normed_covs['Normed_divByBW'+i].Clone(), variables[ivar], outputLabel='Normed_divByBW_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                
                jacobian_covTot = compute_jacobian(allHistos[ 'unfoldHisto'+ivar ].Clone())
                cov_arr = th2_to_np_arr(allHistos[i].Clone()) 
                transformed_cov = np.dot(jacobian_covTot, np.dot(cov_arr, jacobian_covTot.T))
                transformed_cov_matrix = numpy_to_hist2D(transformed_cov, allHistos[i].Clone())                
                normed_covs['Normed'+i] = transformed_cov_matrix.Clone()
                
                #normed_covs.append(normed_cov)
                if process.startswith('data'):
                    draw2D( ivar, normed_covs['Normed'+i].Clone(), variables[ivar], outputLabel='V2Normed_'+i, outputDir=outputDir,selection=sel,version=version,year=year)
                
        ############### Compute correlation matrix for unfolding from nomralised total covariance ######################
        
        allHistos['Normed_correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone()
        allHistos['Normed_correlation_matrix_'+ivar].Reset()
        allHistos['Normed_correlation_matrix_'+ivar] = correlation_from_covariance(normed_covs[ 'Normed'+'cov'+ivar ].Clone(),allHistos['correlation_matrix_'+ivar])
        if process.startswith('data'):

            draw2D( ivar,  allHistos['Normed_correlation_matrix_'+ivar].Clone(), variables[ivar], outputLabel='Normed_data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year,pngToo=True)

    
        def unf_renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                
                #if 'fold' in isam: print (isam)
                try:
                    histName = hist.GetName()
                    sampleName = isam
                except AttributeError:
                    print(f"WARNING: exception raised: histo {isam} does not exists in dictHistos; continuing nonetheless")
                    continue
                if not('all' in year):
                    if 'HEM' in histName:
                        histName=histName.replace('jesHEMIssue', 'jesHEMIssue_2018')
                        sampleName=sampleName.replace('jesHEMIssue', 'jesHEMIssue_2018')
                        #print(sampleName,histName)

                    elif 'jer' in histName:
                        y = year.replace('_preVFP','') if 'VFP' in year else year
                        histName=histName.replace('jer', f'jer_{y}')
                        sampleName=sampleName.replace('jer', f'jer_{y}')
                        #print(sampleName,histName)
                    elif 'btag' in histName.lower() and ('uncorr' in histName.lower() or 'eff' in histName.lower()):
                        y = year.replace('_preVFP','') if 'VFP' in year else year
                        histName=histName.replace('Uncorrelated', f'Uncorrelated_{y}')#.replace('Efficiency', f'Efficiency_{y}')
                        sampleName=sampleName.replace('Uncorrelated', f'Uncorrelated_{y}')#.replace('Efficiency', f'Efficiency_{y}')
                        #print(sampleName,histName)
                    
                
                ihis = hist.Clone(histName+'_renaming_clone')
                ihis.Sumw2()
                
                ihis.SetName(sampleName)#)
                ihis.SetTitle(sampleName)#)
                
                
                
                ihis.Write()
        print("#################################################################")
        print("Dictionary of condition numbers")
        pprint.pprint(dict_condition_numbers)
        if 'dijet' in sel :
            print(f"Dictionary of data-to-MC XS scaling factor for dijets in eras {year} ")
            pprint.pprint(dict_MCScaling)
            with open(f'dijets_dataToMCSFDict_April25_{year}.json','w') as dictToSave:
                json.dump(dict_MCScaling,dictToSave)
        print("#################################################################")
        outputRootName = outputDir+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root'
        if os.path.exists(outputRootName):
            os.remove(outputRootName)
        
        print ('|------> Saving histograms in rootfile: ', outputRootName)
        
        outputRoot = ROOT.TFile.Open( outputRootName, 'RECREATE' )
        
        
        unf_renamingHistos( signalHistos )
        
        if not process.startswith('MC'):
            unf_renamingHistos( sysSignalHistos )
        
            #print ("sys, data int., fsr gen ints., respectively:\n",
            #
            #       dataHistos['data_reco'+ivar+'_nom'+sel].Integral(),
            #      '_fsrWeightUp',
            #       sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightUp'+sel].Integral(), '\n'
            #      '_fsrWeightDown',
            #       sysSignalHistos[f'{fsrLabel}'+'_gen'+ivar+'_fsrWeightDown'+sel].Integral(), '\n'
            #
            #      )    
        if not('self' in process.lower()):
            unf_renamingHistos( altSignalHistos )
            
        if extraMC and process.startswith('data'):
            #print(alt1SignalHistos.keys())
            unf_renamingHistos( alt1SignalHistos )
            if 'dijet' in sel: unf_renamingHistos( alt2SignalHistos )

        if process.startswith('data') and  sel.startswith(('_W','_top')): 
            unf_renamingHistos( varSignalHistos )
            
        unf_renamingHistos(uncerUnfoldSystCov)
        unf_renamingHistos(normed_covs)
        unf_renamingHistos( dataHistos )
        #print ("adding varsighistos", varSignalHistos.items())
        #renamingHistos( dataHistos )
        unf_renamingHistos( bkgHistos )
        #print(allHistos.keys())
        unf_renamingHistos( allHistos )
        unf_renamingHistos( uncerUnfoldHisto )
        
        tunfolder.Write()
        outputRoot.Close()
        #if 'all' in year:
        #    dataFile[ivar+'_2016_preVFP'].Close()
        #    dataFile[ivar+'_2016'].Close()
        #    dataFile[ivar+'_2017'].Close()
        #    dataFile[ivar+'_2018'].Close()
        
        
        if return_tunfolder_object:
            if 'data' in process and 'dijet' in sel:        
                return tunfolder,allHistos, signalHistos, uncerUnfoldHisto, ratioHistos, uncerUnfoldSystCov, normed_covs
            elif 'data' in process and not('dijet' in sel):        
                return tunfolder,allHistos, dataHistos, signalHistos, uncerUnfoldHisto, ratioHistos, uncerUnfoldSystCov, normed_covs, altSignalHistos, alt1SignalHistos, varSignalHistos, bkgHistos, sysSignalHistos
            else:
                if not('self' in process.lower()):
                    return tunfolder,allHistos, signalHistos,altSignalHistos
                else:
                    return tunfolder,allHistos, signalHistos
        #else:
        #    print ('|------> Saving histograms in yodafile: ', outputRootName.replace('.root', '.yoda'))
        #    histToYoda = [  yoda.root.to_yoda( allHistos [ 'unfoldHisto'+ivar ] ) ]
        #    yoda.writeYODA( histToYoda, outputRootName.replace('.root', '.yoda') )
        else:
            if 'MC' in process: 
                del(tunfolder,allHistos, signalHistos, uncerUnfoldHisto)#, ratioHistos, uncerUnfoldSystCov)
            else:
                
                if 'dijet' in sel: 
                    del(tunfolder,allHistos, signalHistos, uncerUnfoldHisto, ratioHistos, sysSignalHistos, uncerUnfoldSystCov)
                else: 
                    del(tunfolder,allHistos, signalHistos, uncerUnfoldHisto, ratioHistos, sysSignalHistos, bkgHistos, varSignalHistos, uncerUnfoldSystCov)

        gc.collect()
        #sys.stdout.flush()
    return 1
    #if return_tunfolder_object: 
    #    return unf_var_dict

##############################################################################################


def build_all_years_histograms(
                                process,
                                year,
                                years_list,
                                dataFile,
                                signalLabel,
                                altSignalLabel,
                                alt1SignalLabel,
                                alt2SignalLabelBegin,
                                alt2SignalLabel,
                                sysSignalLabels,
                                sysUncert,
                                varSignalLabels,
                                bkgLabels,
                                fsrLabel,
                                ivar,
                                sel,
                                genBin,
                                extraMC=False,
                                verbose=False
                            ):
    
    dataHistos = {}
    dataHistostrue = {}
    
    signalHistos = {}
    sysSignalHistos = {} if not process.startswith('MC') else None
    altSignalHistos = {}
    alt1SignalHistos = {}
    alt2SignalHistos = {}
    allHistos = {}
    
    bkgHistos = {} if 'W' in sel or 'top' in sel else None
    varSignalHistos = {} if 'W' in sel or 'top' in sel else None
        
    sysUncs_added = []
    
    if not year.startswith("all"):
        if verbose:
            print(f"[build_all_years_histograms] year = {year} does not start with 'all', "
                  "this function is specialized for the 'all' case.")
        return {
            "dataHistos": dataHistos,
            "dataHistostrue": dataHistostrue,
            "bkgHistos": bkgHistos,
            "signalHistos": signalHistos,
            "sysSignalHistos": sysSignalHistos,
            "altSignalHistos": altSignalHistos,
            "alt1SignalHistos": alt1SignalHistos,
            "alt2SignalHistos": alt2SignalHistos,
            "varSignalHistos": varSignalHistos,
            "allHistos": allHistos,
        }
    
    else:
        ### load nominal MC(/signal MC for W/top)
        signal_keys = [
                        f"{signalLabel}_respWithMiss{ivar}_nom{sel}",
                        f"{signalLabel}_reco{ivar}_nom{sel}",
                        f"{signalLabel}_truereco{ivar}_nom{sel}",
                        f"{signalLabel}_fakereco{ivar}_nom{sel}",
                        f"{signalLabel}_accepgen{ivar}_nom{sel}",
                        f"{signalLabel}_gen{ivar}_nom{sel}",
                        f"{signalLabel}_reco{ivar}_nom{sel}_genBin",
                        f"{signalLabel}_truereco{ivar}_nom{sel}_genBin",
                        f"{signalLabel}_fakereco{ivar}_nom{sel}_genBin",
                      ]

        for y in years_list:
            var_year = f"{ivar}_{y}"
            for key in signal_keys:
                fill_or_add_hist(signalHistos, key, dataFile, var_year, key)
                
        if not process.startswith('MC'):

            #print ("Loading up all syst. variations from the following:", sysUncert)
            for isys,sys in enumerate(sysUncert):
                #print (isys,sys)
                if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                for upDown in [ 'Up', 'Down' ]:
                    #print (isys,sys,upDown)
                    if sys+upDown in sysUncs_added: 
                        print("unc already added in, moving to next unc.",sys+upDown)
                        continue
                    s = [i for i in sysSignalLabels if sys in i]
                    s=[s[0]]    

                    if (not ('2017' in sys) and not ( '2018' in sys ) and not ( '2016' in sys) ): 
                        
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone(s[0]+'_reco'+ivar+sys+upDown+sel+'_clone')
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel).Clone(s[0]+'_respWithMiss'+ivar+sys+upDown+sel+'_clone')
                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel) )
                        if 'fsr' in sys and not(s[0]+'_gen'+ivar+sys+upDown+sel in list(sysSignalHistos.keys())):
                            
                            
                            sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_gen'+ivar+sys+upDown+sel).Clone(s[0]+'_gen'+ivar+sys+upDown+sel+'_clone')
                            
                            #print(
                            #        sys + upDown, s[0],
                            #        f"Input int. 2016_preVFP: {dataFile[ivar + '_2016_preVFP'].Get(s[0] + '_gen' + ivar + sys + upDown + sel).Integral()}",
                            #        f"Nom gen 2016_preVFP: {dataFile[ivar + '_2016_preVFP'].Get('MLMQCD_HT2000toInf' + '_gen' + ivar + '_nom' + sel).Integral()}",
                            #        f"Just filled {s[0]+'_gen'+ivar+sys+upDown+sel}, int.: {sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Integral()}",
                            #        f'Cumul. int. fsrUp : {sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}", None).Integral() if f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}" in sysSignalHistos else "N/A"}',
                            #        f'Cumul. int. fsrDown : {sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}", None).Integral() if f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}" in sysSignalHistos else "N/A"}',
                            #    )
                            
                            
                            sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )
                            
                            
                            
                            
                            sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )
                            
                                                        

                            sysSignalHistos[ s[0]+'_gen'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_gen'+ivar+sys+upDown+sel) )
                            
                            
                        sysUncs_added.append(sys+upDown)
                        
                        
                    
                    # dealing with uncorrelated syst. unc sources below
                    elif '2016' in sys and not( '2017' in sys) and not('2018' in sys):# and s[0] in sys:#.endswith(sys): 


                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone(s[0]+'_reco'+ivar+sys+upDown+sel+'_clone') 

                        #add 2016 systematic to 2016_preVFP systematic
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel))  

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))


                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel).Clone(s[0]+'_respWithMiss'+ivar+sys+upDown+sel+'_clone')

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))


                    elif '2017' in sys and not( '2016' in sys) and not('2018' in sys):

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone(s[0]+'_reco'+ivar+sys+upDown+sel+'_clone') 

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel).Clone(s[0]+'_respWithMiss'+ivar+sys+upDown+sel+'_clone')

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                    elif '2018' in sys and not( '2016' in sys) and not('2017' in sys):# and s[0] in sys:#s[0].endswith(sys): 

                        #print("Making reco/resp histos for all years for:", sys,s[0],s[0]+'_reco'+ivar+sys+upDown+sel)
                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel).Clone(s[0]+'_reco'+ivar+sys+upDown+sel+'_clone') 

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_respWithMiss'+ivar+sys+upDown+sel).Clone(s[0]+'_respWithMiss'+ivar+sys+upDown+sel+'_clone')

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                        sysSignalHistos[ s[0]+'_respWithMiss'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_respWithMiss'+ivar+'_nom'+sel))

                    sysUncs_added.append(sys+upDown)
            
            #print("Sys. uncs. added",sysUncs_added)
        
        # load alt MC's
        if "self" not in process.lower():
            alt_signal_keys = [
                f"{altSignalLabel}_respWithMiss{ivar}_nom{sel}",
                f"{altSignalLabel}_reco{ivar}_nom{sel}",
                f"{altSignalLabel}_reco{ivar}_nom{sel}_genBin",
                f"{altSignalLabel}_truereco{ivar}_nom{sel}",
                f"{altSignalLabel}_truereco{ivar}_nom{sel}_genBin",
                f"{altSignalLabel}_fakereco{ivar}_nom{sel}",
                f"{altSignalLabel}_accepgen{ivar}_nom{sel}",
                f"{altSignalLabel}_missgen{ivar}_nom{sel}",
                f"{altSignalLabel}_gen{ivar}_nom{sel}",
            ]
            for y in years_list:
                var_year = f"{ivar}_{y}"
                for k in alt_signal_keys:
                    fill_or_add_hist(altSignalHistos, k, dataFile, var_year, k)

        
        if not ("mc" in process.lower()) and extraMC:
            alt1_signal_keys = [
                f"{alt1SignalLabel}_respWithMiss{ivar}_nom{sel}",
                f"{alt1SignalLabel}_reco{ivar}_nom{sel}",
                f"{alt1SignalLabel}_reco{ivar}_nom{sel}_genBin",
                f"{alt1SignalLabel}_truereco{ivar}_nom{sel}",
                f"{alt1SignalLabel}_truereco{ivar}_nom{sel}_genBin",
                f"{alt1SignalLabel}_fakereco{ivar}_nom{sel}",
                f"{alt1SignalLabel}_accepgen{ivar}_nom{sel}",
                f"{alt1SignalLabel}_missgen{ivar}_nom{sel}",
                f"{alt1SignalLabel}_gen{ivar}_nom{sel}",
            ]
            for y in years_list:
                var_year = f"{ivar}_{y}"
                for k in alt1_signal_keys:
                    fill_or_add_hist(alt1SignalHistos, k, dataFile, var_year, k)

            if alt2SignalLabelBegin is not None:
                alt2_signal_keys = [
                    f"{alt2SignalLabel}_respWithMiss{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_reco{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_reco{ivar}_nom{sel}_genBin",
                    f"{alt2SignalLabel}_truereco{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_truereco{ivar}_nom{sel}_genBin",
                    f"{alt2SignalLabel}_fakereco{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_accepgen{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_missgen{ivar}_nom{sel}",
                    f"{alt2SignalLabel}_gen{ivar}_nom{sel}",
                ]
                for y in years_list:
                    var_year = f"{ivar}_{y}"
                    for k in alt2_signal_keys:
                        fill_or_add_hist(alt2SignalHistos, k, dataFile, var_year, k)

        if process.startswith('data') and sel.startswith(('_W','_top')):
            
            if verbose: 
                print("Processing bkgs from amongst the following labels:", bkgLabels)
            
            if process.startswith('data') and sel.startswith(('_W','_top')):
                for ibkg in bkgLabels:
                    bkg_key = f"{ibkg}_reco{ivar}_nom{sel}"
                    bkg_key_genbin = f"{ibkg}_reco{ivar}_nom{sel}_genBin"
                    for y in years_list:
                        var_year = f"{ivar}_{y}"
                        fill_or_add_hist(bkgHistos, bkg_key, file_dict=dataFile, var_year=var_year, hist_name=bkg_key)
                        fill_or_add_hist(bkgHistos, bkg_key_genbin, dataFile, var_year, bkg_key_genbin)

            
                s=[]


                if verbose: print ("Processing signal variations from amongst the following uncertainty sources: ", sysUncert)

                for sys in sysUncert:
                    #print (sys)
                    s = [i for i in varSignalLabels if (sys.split('_')[1] in i)]
                    for j in s:
                        if 'Tune' in sys and not('TuneCP5Up' in j or 'TuneCP5Down' in j): continue

                        var_signal_keys = [
                                            j+'_reco'+ivar+'_nom'+sel,
                                            j+'_respWithMiss'+ivar+'_nom'+sel
                                          ]
                        for y in years_list:
                            var_year = f"{ivar}_{y}"
                            for key in var_signal_keys:
                                fill_or_add_hist(varSignalHistos, key, dataFile, var_year, key)
                    
        ### DATA ###
        data_key_recobin = f"dataHisto{ivar}"
        data_key_genbin  = f"dataHistoGenBin{ivar}"
        for y in years_list:
            var_year = f"{ivar}_{y}"
            fill_or_add_hist(allHistos, data_key_recobin, dataFile, var_year, data_key_recobin)
            fill_or_add_hist(allHistos, data_key_genbin, dataFile, var_year, data_key_genbin)

        dataHistos[f"data_reco{ivar}_nom{sel}"] = allHistos[data_key_recobin].Clone( f"data_reco{ivar}_nom{sel}_copy" )
        dataHistos[f"data_reco{ivar}_nom{sel}_genBin"] = allHistos[data_key_genbin].Clone( f"data_reco{ivar}_nom{sel}_genBin_copy" )
        

        if verbose and "data" in process:


            print('data_2016_preVFP', (dataFile[ivar+'_2016_preVFP'].Get( f"dataHisto{ivar}" )).Integral(),  
                  'data_2016', (dataFile[ivar+'_2016'].Get( f"dataHisto{ivar}" )).Integral(), 
                  'data_2017', (dataFile[ivar+'_2017'].Get( f"dataHisto{ivar}" )).Integral(), 
                  'data_2018', (dataFile[ivar+'_2018'].Get( f"dataHisto{ivar}" )).Integral(), )

            print("VARIOUS INTEGRALS: data, data genBin, recoMC, genMC, fsrUp/Down MC")
            print("All years:")
            print(
                allHistos[f"dataHisto{ivar}"].Integral(),
                allHistos[f"dataHistoGenBin{ivar}"].Integral(),
                signalHistos[f"{signalLabel}_reco{ivar}_nom{sel}"].Integral(),
                signalHistos[f"{signalLabel}_gen{ivar}_nom{sel}"].Integral(),
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightUp{sel}" in sysSignalHistos
                else "N/A",
                sysSignalHistos.get(f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}", None).Integral()
                if f"{fsrLabel}_gen{ivar}_fsrWeightDown{sel}" in sysSignalHistos
                else "N/A",
            )



        if "MC" in process:
            dataHistostrue[f"data_reco{ivar}_nom{sel}"] = signalHistos[
                f"{signalLabel}_truereco{ivar}_nom{sel}"
            ].Clone(f"{signalLabel}_data_truereco{ivar}_nom{sel}")

            dataHistostrue[f"data_reco{ivar}_nom{sel}_genBin"] = signalHistos[
                f"{signalLabel}_truereco{ivar}_nom{sel}"
            ].Clone(f"{signalLabel}_data_truereco{ivar}_nom{sel}_genBin")

            dataHistostrue[f"data_reco{ivar}_nom{sel}_genBin"].Rebin(
                len(genBin) - 1,
                f"{signalLabel}_data_truereco{ivar}_nom{sel}_genBin",
                array("d", genBin),
            )
        else:
            dataHistostrue['data_reco'+ivar+'_nom'+sel] = allHistos[ 'dataHisto'+ivar ].Clone('data_reco'+ivar+'_nom'+sel+'_copy')
            dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'] = allHistos[ 'dataHistoGenBin'+ivar].Clone('data_reco'+ivar+'_nom'+sel+'_genBin'+'_copy')
        
        if sel.startswith(('_W','_top')):

            return {
                        "dataHistos": dataHistos,
                        "dataHistostrue": dataHistostrue,
                        "bkgHistos": bkgHistos,
                        "signalHistos": signalHistos,
                        "sysSignalHistos": sysSignalHistos,
                        "altSignalHistos": altSignalHistos,
                        "alt1SignalHistos": alt1SignalHistos,
                        "varSignalHistos": varSignalHistos,
                        "allHistos": allHistos,
                    }
        
        else:
            return {
                        "dataHistos": dataHistos,
                        "dataHistostrue": dataHistostrue,
                        "signalHistos": signalHistos,
                        "sysSignalHistos": sysSignalHistos,
                        "altSignalHistos": altSignalHistos,
                        "alt1SignalHistos": alt1SignalHistos,
                        "alt2SignalHistos": alt2SignalHistos,
                        "allHistos": allHistos,
                    }
def fill_or_add_hist(hist_dict, dict_key, file_dict, var_year, hist_name, clone_title=None):
    """
    Helper: If 'dict_key' not in 'hist_dict', clone from file. 
    Else, add the file's histogram to the existing one.
    
    hist_dict    : dictionary in which we store histograms by key
    dict_key     : the key under which we store the histogram in hist_dict
    file_dict    : your dataFile dict (maps e.g. 'pt_2016' to a TFile object)
    var_year     : string like 'pt_2016' or 'eta_2017', i.e. ivar+'_'+year
    hist_name    : the name of the histogram inside the TFile
    clone_title  : optional new name or title for the cloned histogram 
                   (useful if you want them unique)
    """
    h_infile = file_dict[var_year].Get(hist_name)
    if h_infile is None:
        print(f"[fill_or_add_hist] WARNING: Histogram {hist_name} not found in {var_year}")
        return
    
    if dict_key not in hist_dict:
        if clone_title is None:
            clone_title = dict_key + "_clone"
        hist_dict[dict_key] = h_infile.Clone(clone_title)
    else:
        hist_dict[dict_key].Add(h_infile)


##############################################################################################
############ Histogram loader and rebinner for inputs to unfolding script below ##############
##############################################################################################

def loadHistograms(samples, var, sel, sysUnc=[],
                   isMC=True, addGenInfo=True, respOnly=False, lumi=1., noResp=False,
                   variables={}, year='2017', process='data', noRebin=False, outputFolder=None,
                   recoOnly=False, respName='respWithMiss'
                  ):
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
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in SYSUNC if 'nom' in syst or 'fsr' in syst or 'isr' in syst or 'pdf' in syst.lower() ] 
            tmpList = tmpList + [ 'reco'+var+syst+sel for syst in SYSUNC if not(('reco'+var+syst+sel) in tmpList)]
            if not(noResp):
                tmpList = tmpList + [ 'accepgen'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'truereco'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'fakereco'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ 'missgen'+var+syst+sel for syst in SYSUNC ]
                tmpList = tmpList + [ f'{respName}'+var+syst+sel for syst in SYSUNC]
        
        elif isMC and addGenInfo and flip: 
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in tmpSYSUNC if 'nom' in syst or 'fsr' in syst or 'isr' in syst or 'pdf' in syst.lower() ]
            tmpList = tmpList + [ 'reco'+var+syst+sel for syst in tmpSYSUNC if not(('reco'+var+syst+sel) in tmpList) ]
            if not(noResp): 
                tmpList = tmpList + [f'{respName}'+var+syst+sel for syst in tmpSYSUNC]
        
        if respOnly and not flip: 
            if not(noResp): 
                tmpList = [ f'{respName}'+var+syst+sel for syst in SYSUNC] 
        elif respOnly and flip: 
            if not(noResp): 
                tmpList = [ f'{respName}'+var+syst+sel for syst in tmpSYSUNC ]
       
        elif recoOnly:
            tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC] if not flip else [ 'reco'+var+syst+sel for syst in tmpSYSUNC ]
        
        
        set_tmpList = list(set(tmpList))
        if not(len(tmpList)==len(set_tmpList)):
            print(len(tmpList),len(set_tmpList))
            tmpList = set_tmpList
            print(len(tmpList),len(set_tmpList))
        #print(tmpList)
        for ih in tmpList:
            
            if isMC:
                #print(ih, samples[isam][0])
                try:
                    #if 'pdf' in ih.lower():
                    #    print(isam,samples[isam][0],isam+'_'+ih )
                    iFile = ROOT.TFile.Open(samples[isam][0],'r')
                    #print(isam)
                    allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone(isam+'_'+ih) #'jetObservables/'+
                    hist = allHistos[isam+'_'+ih]

                    if variables[var]['bins'][-1] < allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0)):
    
                        print("########################################")
                        print(f"WARNING: loadHistograms()" )
                        print(f"(allHistos[isam+'_'+ih]) for {isam,ih} has events beyond last gen bin ({variables[var]['bins'][-2],variables[var]['bins'][-1]}) in bin ({allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0))}, {allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0)+1)} " )
                        print("########################################")
                    if isinstance(hist, ROOT.TH2):
                        nbX, nbY = hist.GetNbinsX(), hist.GetNbinsY()
                        ofX = sum(hist.GetBinContent(nbX+1, y) for y in range(1, nbY+1))
                        ofY = sum(hist.GetBinContent(x, nbY+1) for x in range(1, nbX+1))
                        ofXY = hist.GetBinContent(nbX+1, nbY+1)
                        if ofX or ofY or ofXY:
                            print("########################################")
                            print(f"WARNING: loadHistograms() — TH2 {ih!r} has overflow:")
                            print(f"  x-overflow total = {ofX:.1f}, y-overflow total = {ofY:.1f}, corner = {ofXY:.1f}")
                            print("########################################")
                    elif isinstance(hist, ROOT.TH1):
                        of = hist.GetBinContent(hist.GetNbinsX()+1)
                        if of:
                            print("########################################")
                            print(f"WARNING: loadHistograms() — TH1 {ih!r} has {of:.1f} entries in the x overflow bin")
                            print("########################################")

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
                    
                    #if 'jer' in isam:
                    #    print("##################################################################")
                    #    print(isam, samples[isam][0], ih, isam+'_'+ih,  "FILE opened")
                    #    print("##################################################################")
                    #if 'QCD' in isam and ih.startswith('reco'): print(isam, MCScale, allHistos[isam+'_'+ih].Integral())
                except ReferenceError:
                    print("##################################################################")
                    print(isam, samples[isam][0], ih, isam+'_'+ih,  "FILE NOT FOUND/CORRUPTED; skipping for now")
                    print("##################################################################")
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
                        hist = tmpdataHistos[ it ]

                        if variables[var]['bins'][-1] < tmpdataHistos[ it ].GetBinLowEdge(tmpdataHistos[ it ].FindLastBinAbove(0)):
                            print("########################################")
                            print(f"WARNING: loadHistograms()" )
                            print(f"(tmpdataHistos[ it ]) for {it} has events beyond last reco bin ({variables[var]['bins_reco'][-1],variables[var]['bins_reco'][-1]}) in bin ({tmpdataHistos[ it ].GetBinLowEdge(tmpdataHistos[ it ].FindLastBinAbove(0))}, {tmpdataHistos[ it ].GetBinLowEdge(tmpdataHistos[ it ].FindLastBinAbove(0)+1)} " )
                            print("########################################")
                        if isinstance(hist, ROOT.TH2):
                            nbX, nbY = hist.GetNbinsX(), hist.GetNbinsY()
                            ofX = sum(hist.GetBinContent(nbX+1, y) for y in range(1, nbY+1))
                            ofY = sum(hist.GetBinContent(x, nbY+1) for x in range(1, nbX+1))
                            ofXY = hist.GetBinContent(nbX+1, nbY+1)
                            if ofX or ofY or ofXY:
                                print("########################################")
                                print(f"WARNING: loadHistograms() — TH2 {ih!r} has overflow:")
                                print(f"  x-overflow total = {ofX:.1f}, y-overflow total = {ofY:.1f}, corner = {ofXY:.1f}")
                                print("########################################")
                        elif isinstance(hist, ROOT.TH1):
                            of = hist.GetBinContent(hist.GetNbinsX()+1)
                            if of:
                                print("########################################")
                                print(f"WARNING: loadHistograms() — TH1 {ih!r} has {of:.1f} entries in the x overflow bin")
                                print("########################################")
                        
                    itempFile.Close() 
                    allHistos[ isam+'_'+ih ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
                    allHistos[ isam+'_'+ih ].Reset()
                    for i in tmpdataHistos: 
                        tmpdataHistos[i].SetDirectory(0)
                        allHistos[isam+'_'+ih].Add( tmpdataHistos[i].Clone() )
                        allHistos[isam+'_'+ih].SetDirectory(0)
                else:
                    iFile = ROOT.TFile.Open(samples[isam][0],'r')
                    allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone(isam+'_'+ih)
                    hist = allHistos[isam+'_'+ih]
                    allHistos[isam+'_'+ih].SetDirectory(0)
                    allHistos[isam+'_'+ih].Sumw2()
                    iFile.Close()
                    
                    if variables[var]['bins'][-1] < allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0)):
                        print("########################################")
                        print(f"WARNING: loadHistograms()" )
                        print(f"(allHistos[isam+'_'+ih]) for {isam,ih} has events beyond last reco bin ({variables[var]['bins_reco'][-1],variables[var]['bins_reco'][-1]}) in bin ({allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0))}, {allHistos[isam+'_'+ih].GetBinLowEdge(allHistos[isam+'_'+ih].FindLastBinAbove(0)+1)} " )
                        print("########################################")
                        
                    if isinstance(hist, ROOT.TH2):
                        nbX, nbY = hist.GetNbinsX(), hist.GetNbinsY()
                        ofX = sum(hist.GetBinContent(nbX+1, y) for y in range(1, nbY+1))
                        ofY = sum(hist.GetBinContent(x, nbY+1) for x in range(1, nbX+1))
                        ofXY = hist.GetBinContent(nbX+1, nbY+1)
                        if ofX or ofY or ofXY:
                            print("########################################")
                            print(f"WARNING: loadHistograms() — TH2 {ih!r} has overflow:")
                            print(f"  x-overflow total = {ofX:.1f}, y-overflow total = {ofY:.1f}, corner = {ofXY:.1f}")
                            print("########################################")
                    elif isinstance(hist, ROOT.TH1):
                        of = hist.GetBinContent(hist.GetNbinsX()+1)
                        if of:
                            print("########################################")
                            print(f"WARNING: loadHistograms() — TH1 {ih!r} has {of:.1f} entries in the x overflow bin")
                            print("########################################")
                    
                    
    
    def renamingHistos( dictHistos ):
        for isam, hist in dictHistos.items():
            ihis = hist.Clone(hist.GetName()+'_clone')
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
                    #print(ih)
                    #print(tmpHistos[ih].Integral())
                    
                    #if 'fsr' in sys: print(ih)
                    for jh in allHistos:
                        #if 'fsr' in sys: print(jh)
                        #print(jh)
                        goodflag=False
                        if ('2016' in ih and '2016' in jh and '2016' in year) or ('2016' in ih and '2016' in jh and '2016_preVFP' in year) or ('2017' in ih and '2017' in jh and '2017' in year) or ('2018' in ih and '2018' in jh and '2018' in year) and (sys in ih) and (sys in jh):
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral())

                        elif not('201' in ih) and not ('201'in jh) and (sys in ih) and (sys in jh):
                            goodflag=True
                            #print(ih,jh,tmpHistos[ih].Integral(),allHistos[jh].Integral())

                        if goodflag and (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ) and sys in ih and sys in jh :
                            if ('_recoJet' in ih and '_recoJet' in jh) or ('truerecoJet' in ih and 'truerecoJet' in jh) or ('fakerecoJet' in ih and 'fakerecoJet' in jh) or ('_genJet' in ih and '_genJet' in jh) or  ('accepgenJet' in ih and 'accepgenJet' in jh) or ('missgenJet' in ih and 'missgenJet' in jh) or ('respWithMissJet' in ih and 'respWithMissJet' in jh) or ('good' in ih and 'good' in jh):
                                if ('genBin' in ih and 'genBin' in  jh) or (not('genBin' in ih) and not('genBin' in  jh)):
                                    #if 'fsr' in sys: print(goodflag,ih,jh,tmpHistos[ih].Integral())
                                    #print(ih,jh,tmpHistos[ih].Integral(),allHistos[jh].Integral())
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
            nomSysList = [k.split('TTToSemiLeptonic')[1] for k in allHistos.keys() if ('TTToSemiLeptonic' in k) and not('var' in k or 'jes' in k or 'jer' in k)]
            #print("nom,wt, hist list", nomSysList,allHistos.keys())

            for sNomWt in nomSysList:
                #print("sNomWt",sNomWt)
                tmpHistos = { k:v for (k,v) in allHistos.items() if ('TTToSemiLeptonic' in k) and (sNomWt in k) and not('var' in k or 'jes' in k or 'jer' in k)}
                ttSigFlag = False if (len(tmpHistos.keys())==0) else True

                if ttSigFlag: 
                    allHistos_upd=copy.deepcopy(allHistos)

                    for k in allHistos.keys():
                        if 'TTTo' in k and not(k.startswith('var')) and not('var' in k or 'jes' in k or 'jer' in k) and (sNomWt in k):
                            #print(k)
                            del(allHistos_upd[k])
                            ttSigFlag = True


                    #if len(tmpHistos.keys())>0:
                    #    print (f"Adding together ttbar histograms for {var} in {sel}",tmpHistos.keys())


                    for ih in tmpHistos:
                        for jh in allHistos:

                            if ('TTTo' in jh and not(jh.startswith('var')) and not('jes' in jh or 'jer' in jh)) and not ('Semi' in jh ) and (sNomWt in jh):
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
            altSysList = [k.split('varTTToSemileptonic')[1] for k in allHistos.keys() if ('varTTToSemileptonic' in k)]
            #print("altSys, hist list", altSysList)

            for altSys in altSysList:
                #print("altSys",altSys)
                tmpHistos = { k:v for (k,v) in allHistos.items() if ('varTTToSemileptonic' in k) and (altSys in k)}
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
                    #    print (f"Adding together alt TTbar histograms for {var} in {sel}",tmpHistos.keys())


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
                    
            ##ttbar systematics
            altSysList = [k.split('sysTTToSemiLeptonic')[1] for k in allHistos.keys() if ('sysTTToSemiLeptonic' in k) and ('jes' in k or 'jer' in k)]
            #print("altSys, hist list", altSysList)

            for altSys in altSysList:
                #print("altSys",altSys)
                tmpHistos = { k:v for (k,v) in allHistos.items() if ('sysTTToSemiLeptonic' in k) and (altSys in k) and ('jes' in k or 'jer' in k)}
                ttSigFlag = False if (len(tmpHistos.keys())==0) else True

                #print(qcdFlag)
                if ttSigFlag: 
                    allHistos_upd=copy.deepcopy(allHistos)

                    for k in allHistos.keys():
                        if 'sysTTo' in k and (altSys in k) and ('jes' in k or 'jer' in k):
                            #print(k)
                            del(allHistos_upd[k])
                            ttSigFlag = True


                    #if len(tmpHistos.keys())>0:
                    #    print (f"Adding together sys TTbar histograms for {var} in {sel}",tmpHistos.keys())


                    for ih in tmpHistos:
                        for jh in allHistos:

                            if ('sysTTTo' in jh) and not('Semi' in jh ) and (altSys in jh) and ('jes' in jh or 'jer' in jh):
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
            if not(f'{respName}' in ih):
                #print(var,genBin,recoBin)

                #if len(variables[var]['bins'])==1:
                #    if 'recoJet' in ih:
                #        allHistos[ih+'_genBin'] = allHistos[ih].Clone(allHistos[ih].GetName()+'_genBin')
                #        allHistos[ih+'_genBin'].Rebin( genBin )
                #        allHistos[ih].Rebin( recoBin )
                #    elif 'genJet' in ih: allHistos[ih].Rebin( genBin )

                #else:
                if 'recoJet' in ih:
                    #if ('MLM' in ih and not 'H7' in ih):
                    #    print(ih, genBin,allHistos[ih].GetName(),allHistos[ih].Integral())
                    allHistos[ih+'_genBin'] = allHistos[ih].Clone(allHistos[ih].GetName()+'_genBin')
                    allHistos[ih+'_genBin'] = allHistos[ih+'_genBin'].Rebin( len(genBin)-1, allHistos[ih].GetName()+"_Rebin_genBin", array( 'd', genBin ) )
                    allHistos[ih] = allHistos[ih].Rebin( len(recoBin)-1, allHistos[ih].GetName()+"_Rebin", array( 'd', recoBin ) )
                elif 'genJet' in ih:
                    allHistos[ih] = allHistos[ih].Rebin( len(genBin)-1, allHistos[ih].GetName()+"_Rebin", array( 'd', genBin ) )

            else:
                
                #if len(variables[var]['bins'])==1: 
                #    allHistos[ih].Rebin2D( genBin,recoBin )
                #else:

                #### fancy way to create variable binning TH2D
                #tmpHisto = ROOT.TH2F( allHistos[ih].GetName()+"_Rebin", allHistos[ih].GetName()+"_Rebin", len(genBin)-1, array( 'd', genBin), len(recoBin)-1, array( 'd', recoBin) )
                #tmpHisto.Sumw2()
                #if tmpHisto.GetNbinsY()>500: print(ih, tmpHisto.GetNbinsY())

                #tmpHisto = rebin_RM_withUFandOF2(allHistos[ih],genBin,recoBin).Clone(allHistos[ih].GetName()+"_Rebin")
                tmpHisto = rebin_RM_withUF(allHistos[ih],genBin,recoBin).Clone(allHistos[ih].GetName()+"_Rebin")
                #make_rebinned_2d_hist(allHistos[ih].Clone(), new_bin_edge_pairs,True)
                tmpHisto.Sumw2()                        

                tmpHisto.SetDirectory(0)
                allHistos[ih] = copy.deepcopy(tmpHisto.Clone(tmpHisto.GetName()+'_clone'))
        
        if samples and outputFolder: 
            if not(os.path.exists(outputFolder+ f"/{sel.replace('_','')}/{year}")):
                os.makedirs(outputFolder+ f"/{sel.replace('_','')}/{year}" )
            outputRootName = outputFolder+ f"/{sel.replace('_','')}/{year}" + '/loadedHistograms_main_'+process+isam+var+year+'.root'
            print ('|------> Saving histograms in rootfile: ', outputRootName)
            outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
            renamingHistos( copy.deepcopy(allHistos) )
            #print(allHistos.keys())
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

def getAndPlotPurity(h_resp_rebinned,reco,accepgen,gen_bins,variables,var,lumi,sel='_dijetSel',outputDir='../Results/',year='2017'):
    
    rebinned=h_resp_rebinned.Clone()#make_rebinned_2d_hist(h_resp.Clone(),new_bin_edge_pairs,)#rebinning to gen-level bins
    
    arr_rebinned,_ = th2_to_np_arr(rebinned.Clone())
    rebinned_array2d_normX = renorm(arr_rebinned, axis=0) # normalise axis to 1, renormed per x/gen bin
    rebinned_array2d_normY = renorm(arr_rebinned, axis=1) # normalise axis to 1, renormed per y/reco bin
    
    p_list=[]
    s_list=[]
    #accep_list=[]
    #fake_list=[]
    gen = rebinned.ProjectionX('gen'+var+year,0,rebinned.GetNbinsX())
    accepGen = accepgen.Clone()#rebinned.ProjectionX('accepGen',0,rebinned.GetNbinsX())#+1)
    accepGen.Sumw2()
    accepRate = accepGen.Clone('acceptance'+var+year)
    accepRate.Reset()
    accepRate.Sumw2()
    accepRate.Divide(accepGen, gen, 1, 1,'B')
    
    trueReco = rebinned.ProjectionY('trueReco'+var+year)#,1,rebinned.GetNbinsX()+1)
    trueReco.Sumw2()
    fakeReco = reco.Clone('fakeReco'+var+year)
    fakeReco.Sumw2()
    fakeReco.Add(trueReco,-1.)
    fakeRate = fakeReco.Clone('fakeRate'+var+year)
    fakeRate.Reset()
    fakeRate.Sumw2()
    
    fakeRate.Divide(fakeReco, reco, 1, 1, 'B')
    
    
    for ibin in range(len(gen_bins)-1):
        #print (f"Calculating p/s per bin in final new binning for bin: {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]}")
        purity = rebinned_array2d_normY[ibin][ibin] #contains fraction in a reco bin that are actually from the same gen bin
        stability = rebinned_array2d_normX[ibin][ibin] #contains fraction in a gen bin that are actually from the same reco bin
        p_list.append(purity)
        s_list.append(stability)
        #accep_list.append(accepRate.GetBinContent(ibin+1))
        #fake_list.append(fakeRate.GetBinContent(ibin+1))
        #print (f"Purity, stability in bin {ibin}({gen_bins[i],gen_bins[i+1]}): {purity,stability}")
    

    if 'dijet' in sel:
        signalLabelBegin = 'QCD_HT_MG5-MLM+P8'
    else:
        signalLabelBegin='TTToSemiLeptonic'
    accepGen.SetDirectory(0)
    fakeReco.SetDirectory(0)
    makePSplot_simple(purity=array('d',p_list),stability=array('d',s_list),
                      accepGen=accepRate,fakeReco=fakeRate,
                      dictHistos=OrderedDict(),
                      variables=variables,
                      var=var,lumi=lumi,
                      outputDir=outputDir,bins=gen_bins,year=year,
                      sel=sel,
                      signalLabelBegin=signalLabelBegin,
                      ext='pdf'
                     )
    return 1


    
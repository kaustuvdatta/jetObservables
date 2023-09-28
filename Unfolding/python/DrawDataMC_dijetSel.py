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
sys.path.insert(0,'../../')
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
#import pdb
lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)



def runPlots_DijetSel(
                        dataFile, sigFiles, bkgFiles, variables, sel, sysUncert, process, ext, lumi=1., sysSignalLabels=[],
                        year='2017',runMLU=False, sysSigFiles=[], varSigFiles=[], jetType='Central', outputFolder='../Results',
                        version='_Feb23', main='MLM_HTbin', alt0='H7MLM_HTbin', alt1='HTbin', alt2='Ptbin',
                        verbose=False
                      ):
    
    
    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 41, 45, 48, 49, 37, 38, 33, 17]

    

    if main.startswith('MLM_HTbin'):
        signalLabelBegin = 'MLMQCD_HT'
        signalLabel = 'MLMQCD_HT2000toInf'

    elif main.startswith('HTbin'):
        signalLabelBegin = 'QCD_HT'
        signalLabel = 'QCD_HT2000toInf'

    elif main.startswith('Ptbin'):
        signalLabelBegin = 'QCD_Pt_'
        signalLabel = 'QCD_Pt_3200toInf'

    elif main.startswith('H7MLM_HTbin'):
        signalLabelBegin = 'H7MLMQCD_HT'
        signalLabel = 'H7MLMQCD_HT2000toInf'

    sysSignalLabelBegin = 'sysMLMQCD' 

    if alt0.startswith('MLM_HTbin'):
        alt0SignalLabelBegin = 'MLMQCD_HT'
        alt0SignalLabel = 'MLMQCD_HT2000toInf'

    elif alt0.startswith('HTbin'):
        alt0SignalLabelBegin = 'QCD_HT'
        alt0SignalLabel = 'QCD_HT2000toInf'

    elif alt0.startswith('Ptbin'):
        alt0SignalLabelBegin = 'QCD_Pt_'
        alt0SignalLabel = 'QCD_Pt_3200toInf'

    elif alt0.startswith('H7MLM_HTbin'):
        alt0SignalLabelBegin = 'H7MLMQCD_HT'
        alt0SignalLabel = 'H7MLMQCD_HT2000toInf'

    if alt1.startswith('MLM_HTbin'):
        alt1SignalLabelBegin = 'MLMQCD_HT'
        alt1SignalLabel = 'MLMQCD_HT2000toInf'

    elif alt1.startswith('HTbin'):
        alt1SignalLabelBegin = 'QCD_HT'
        alt1SignalLabel = 'QCD_HT2000toInf'

    elif alt1.startswith('Ptbin'):
        alt1SignalLabelBegin = 'QCD_Pt_'
        alt1SignalLabel = 'QCD_Pt_3200toInf'

    elif alt1.startswith('H7MLM_HTbin'):
        alt1SignalLabelBegin = 'H7MLMQCD_HT'
        alt1SignalLabel = 'H7MLMQCD_HT2000toInf'
        
    if alt2.startswith('MLM_HTbin'):
        alt2SignalLabelBegin = 'MLMQCD_HT'
        alt2SignalLabel = 'MLMQCD_HT2000toInf'

    elif alt2.startswith('HTbin'):
        alt2SignalLabelBegin = 'QCD_HT'
        alt2SignalLabel = 'QCD_HT2000toInf'

    elif alt2.startswith('Ptbin'):
        alt2SignalLabelBegin = 'QCD_Pt_'
        alt2SignalLabel = 'QCD_Pt_3200toInf'

    elif alt2.startswith('H7MLM_HTbin'):
        alt2SignalLabelBegin = 'H7MLMQCD_HT'
        alt2SignalLabel = 'H7MLMQCD_HT2000toInf'
    
    
    print ("Labels:", signalLabelBegin,alt0SignalLabelBegin,alt1SignalLabelBegin,alt2SignalLabelBegin)
    
    
    #sysSignalLabelBegin = 'sysMLMQCD_' 
    selection=sel
    
    
    for ivar in variables:
        print (ivar)
        dataLabel = 'data_reco'+ivar+'_nom'+sel if not('nPV' in ivar) else 'data_'+ivar+'_nom'+sel
        sigLabel = signalLabel+'_reco'+ivar+'_nom'+sel if not('nPV' in ivar) else signalLabel+'_'+ivar+'_nom'+sel
        alt0sigLabel = alt0SignalLabel+'_reco'+ivar+'_nom'+sel if not('nPV' in ivar) else alt0SignalLabel+'_'+ivar+'_nom'+sel 
        alt1sigLabel = alt1SignalLabel+'_reco'+ivar+'_nom'+sel if not('nPV' in ivar) else alt1SignalLabel+'_'+ivar+'_nom'+sel
        alt2sigLabel = alt2SignalLabel+'_reco'+ivar+'_nom'+sel if not('nPV' in ivar) else alt2SignalLabel+'_'+ivar+'_nom'+sel
        
        
        outputDir=outputFolder+sel.split('_')[1]+'/'+year+'/DataMC/'+ivar+'/'+process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        genBin = variables[ivar]['bins']
        recoBin = variables[ivar]['bins_reco']

    ###########################################################################################
        if year.startswith('all'):
            signalHistos = {
                    sigLabel : dataFile[ivar+'_2016_preVFP'].Get(sigLabel),
                    #signalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_gen'+ivar+'_nom'+sel),
                    sigLabel+'_genBin': dataFile[ivar+'_2016_preVFP'].Get(sigLabel+'_genBin') if ('tau' in ivar) else None
                    }
            signalHistos[ sigLabel ].Add( dataFile[ivar+'_2016'].Get(sigLabel) )
            #signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): signalHistos[ sigLabel+'_genBin' ].Add( dataFile[ivar+'_2016'].Get(sigLabel+'_genBin') ) 
            
            signalHistos[ sigLabel ].Add( dataFile[ivar+'_2017'].Get(sigLabel) )
            #signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): signalHistos[ sigLabel+'_genBin' ].Add( dataFile[ivar+'_2017'].Get(sigLabel+'_genBin') )
            
            signalHistos[ sigLabel ].Add( dataFile[ivar+'_2018'].Get(sigLabel) )
            #signalHistos[ signalLabel+'_gen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): signalHistos[ sigLabel+'_genBin' ].Add( dataFile[ivar+'_2018'].Get(sigLabel+'_genBin') )            
            

            alt0SignalHistos = {
                    alt0sigLabel : dataFile[ivar+'_2016_preVFP'].Get(alt0sigLabel),
                    #alt0SignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt0SignalLabel+'_gen'+ivar+'_nom'+sel),
                    alt0sigLabel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt0sigLabel+'_genBin') if ('tau' in ivar) else None,
                    }

            alt0SignalHistos[ alt0sigLabel ].Add( dataFile[ivar+'_2016'].Get(alt0sigLabel) )
            #alt0SignalHistos[ alt0SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2016'].Get(alt0SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt0SignalHistos[alt0sigLabel+'_genBin'].Add( dataFile[ivar+'_2016'].Get(alt0sigLabel+'_genBin') )

            alt0SignalHistos[ alt0sigLabel ].Add( dataFile[ivar+'_2017'].Get(alt0sigLabel) )
            #alt0SignalHistos[ alt0SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2017'].Get(alt0SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt0SignalHistos[alt0sigLabel+'_genBin'].Add( dataFile[ivar+'_2017'].Get(alt0sigLabel+'_genBin') )

            alt0SignalHistos[ alt0sigLabel ].Add( dataFile[ivar+'_2018'].Get(alt0sigLabel) )
            #alt0SignalHistos[ alt0SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2018'].Get(alt0SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt0SignalHistos[alt0sigLabel+'_genBin'].Add( dataFile[ivar+'_2018'].Get(alt0sigLabel+'_genBin') )
            
            alt1SignalHistos = {
                alt1sigLabel : dataFile[ivar+'_2016_preVFP'].Get(alt1sigLabel),
                #alt1SignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel),
                alt1sigLabel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt1sigLabel+'_genBin') if ('tau' in ivar) else None,
                }

            alt1SignalHistos[ alt1sigLabel ].Add( dataFile[ivar+'_2016'].Get(alt1sigLabel) )
            #alt1SignalHistos[ alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2016'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt1SignalHistos[alt1sigLabel+'_genBin'].Add( dataFile[ivar+'_2016'].Get(alt1sigLabel+'_genBin') )

            alt1SignalHistos[ alt1sigLabel ].Add( dataFile[ivar+'_2017'].Get(alt1sigLabel) )
            #alt1SignalHistos[ alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2017'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt1SignalHistos[alt1sigLabel+'_genBin'].Add( dataFile[ivar+'_2017'].Get(alt1sigLabel+'_genBin') )

            alt1SignalHistos[ alt1sigLabel ].Add( dataFile[ivar+'_2018'].Get(alt1sigLabel) )
            #alt1SignalHistos[ alt1SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2018'].Get(alt1SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt1SignalHistos[alt1sigLabel+'_genBin'].Add( dataFile[ivar+'_2018'].Get(alt1sigLabel+'_genBin') )
            
            
            alt2SignalHistos = {
                alt2sigLabel : dataFile[ivar+'_2016_preVFP'].Get(alt2sigLabel),
                #alt2SignalLabel+'_gen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel),
                alt2sigLabel+'_genBin' : dataFile[ivar+'_2016_preVFP'].Get(alt2sigLabel+'_genBin') if ('tau' in ivar) else None,
                }

            alt2SignalHistos[ alt2sigLabel ].Add( dataFile[ivar+'_2016'].Get(alt2sigLabel) )
            #alt2SignalHistos[ alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2016'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt2SignalHistos[alt2sigLabel+'_genBin'].Add( dataFile[ivar+'_2016'].Get(alt2sigLabel+'_genBin') )

            alt2SignalHistos[ alt2sigLabel ].Add( dataFile[ivar+'_2017'].Get(alt2sigLabel) )
            #alt2SignalHistos[ alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2017'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt2SignalHistos[alt2sigLabel+'_genBin'].Add( dataFile[ivar+'_2017'].Get(alt2sigLabel+'_genBin') )

            alt2SignalHistos[ alt2sigLabel ].Add( dataFile[ivar+'_2018'].Get(alt2sigLabel) )
            #alt2SignalHistos[ alt2SignalLabel+'_gen'+ivar+'_nom'+sel].Add( dataFile[ivar+'_2018'].Get(alt2SignalLabel+'_gen'+ivar+'_nom'+sel) )
            if ('tau' in ivar): alt2SignalHistos[alt2sigLabel+'_genBin'].Add( dataFile[ivar+'_2018'].Get(alt2sigLabel+'_genBin') )
            

            #if process.startswith('data'):
            allHistos = {
                        'dataHisto' : dataFile[ivar+'_2016_preVFP'].Get( 'dataHisto' ),
                        'dataHistoGenBin' : dataFile[ivar+'_2016_preVFP'].Get( 'dataHistoGenBin' ) if ('tau' in ivar) else None
                        }
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2016'].Get( 'dataHisto' ) )
            if ('tau' in ivar): allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2016'].Get( 'dataHistoGenBin' ) )
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2017'].Get( 'dataHisto' ) )
            if ('tau' in ivar): allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2017'].Get( 'dataHistoGenBin' ) )
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2018'].Get( 'dataHisto' ) )
            if ('tau' in ivar): allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2018'].Get( 'dataHistoGenBin' ) )


            dataHistos = { }
            bkgHistos = { }

        else:
            print('|-------> Running single year '+year,lumi)
            ### Getting input histos
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[], respOnly=False, lumi=lumi, noResp=True, year=year, process=process, variables=variables,outputFolder=outputFolder,noRebin=False, jetType=jetType)
            alt0SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt0SignalLabelBegin)  }
            alt1SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt1SignalLabelBegin)  }
            alt2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(alt2SignalLabelBegin)  }
            bkgHistos = {}
            #if not process.startswith('MC'): 
            alt0SignalHistos = loadHistograms( alt0SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, lumi=lumi, noResp=True, year=year, process=process, variables=variables,outputFolder=outputFolder,noRebin=False, jetType=jetType)#True if not('tau'in ivar ) else False)
            alt1SignalHistos = loadHistograms( alt1SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, lumi=lumi, noResp=True, year=year, process=process, variables=variables,outputFolder=outputFolder,noRebin=False, jetType=jetType)#True if not('tau'in ivar ) else False)
            alt2SignalHistos = loadHistograms( alt2SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, lumi=lumi, noResp=True, year=year, process=process, variables=variables,outputFolder=outputFolder,noRebin=False, jetType=jetType)#True if not('tau'in ivar ) else False)

        
            dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder=outputFolder,noResp=True,noRebin=False, jetType=jetType)#True if not('tau'in ivar ) else False)
            #print(dataHistos,dataHistos[dataLabel].Integral(),signalHistos, signalHistos[sigLabel].Integral())
            
            allHistos = {}
            
            print(dataHistos,signalHistos)
            
            ### For dijet, scale QCD to data
            scaleFactor=1.
            scaleFactorGenBin=1.
            altscaleFactor=1.
            altscaleFactorGenBin=1.
            
            #print('data_reco'+ivar+'_nom'+sel,sigLabel)
            #print(dataHistos['data_reco'+ivar+'_nom'+sel].Integral(), signalHistos[sigLabel].Integral())
            scaleFactor = dataHistos[dataLabel].Integral() / signalHistos[sigLabel].Integral()
            scaleFactorGenBin = dataHistos[dataLabel+'_genBin'].Integral() / signalHistos[sigLabel+'_genBin'].Integral() if ('tau' in ivar) else 1.
            print ("Nominal MC, MG5-MLM+P8, SF genBin and recobin, respectively:",scaleFactor,scaleFactorGenBin)
            
            
            for ihsig in signalHistos:
                if ihsig.endswith(sel):
                    #print (ihsig)
                    signalHistos[ihsig].Scale( scaleFactor )
                elif ihsig.endswith('genBin'): signalHistos[ihsig].Scale( scaleFactorGenBin )
                
            
            #if not(process.startswith('MC')):
            print("Rescaling alternate signal MC")
            
            ###### H7
            altscaleFactor = dataHistos[dataLabel].Integral() / (alt0SignalHistos[ alt0sigLabel ].Integral()  )
            altscaleFactorGenBin = dataHistos[dataLabel+'_genBin'].Integral() / ( alt0SignalHistos[alt0sigLabel+'_genBin'].Integral() ) if ('tau' in ivar) else 1.
            print ("Alt MC, MG5-MLM+H7, SF genBin and recobin, respectively:",altscaleFactor,altscaleFactorGenBin)

            for ihsig in alt0SignalHistos:
                if ihsig.endswith(sel):
                    #print (ihsig)
                    alt0SignalHistos[ihsig].Scale( altscaleFactor )
                elif ihsig.endswith('genBin'): alt0SignalHistos[ihsig].Scale( altscaleFactorGenBin )

            ###### MG5(no MLM)
            altscaleFactor = dataHistos[dataLabel].Integral() / (alt1SignalHistos[ alt1sigLabel ].Integral()  )
            altscaleFactorGenBin = dataHistos[dataLabel+'_genBin'].Integral() / ( alt1SignalHistos[alt1sigLabel+'_genBin'].Integral() ) if ('tau' in ivar) else 1.
            print ("Alt MC, MG5+P8, SF genBin and recobin, respectively:",altscaleFactor,altscaleFactorGenBin)

            for ihsig in alt1SignalHistos:
                if ihsig.endswith(sel):
                    #print (ihsig)
                    alt1SignalHistos[ihsig].Scale( altscaleFactor )
                elif ihsig.endswith('genBin'): alt1SignalHistos[ihsig].Scale( altscaleFactorGenBin )
                    
            ###### P8
            altscaleFactor = dataHistos[dataLabel].Integral() / (alt2SignalHistos[ alt2sigLabel ].Integral()  )
            altscaleFactorGenBin = dataHistos[dataLabel+'_genBin'].Integral() / ( alt2SignalHistos[alt2sigLabel+'_genBin'].Integral() ) if ('tau' in ivar) else 1.
            print ("Alt MC, P8+P8, SF genBin and recobin, respectively:",altscaleFactor,altscaleFactorGenBin)

            for ihsig in alt2SignalHistos:
                if ihsig.endswith(sel):
                    #print (ihsig)
                    alt2SignalHistos[ihsig].Scale( altscaleFactor )
                elif ihsig.endswith('genBin'): alt2SignalHistos[ihsig].Scale( altscaleFactorGenBin )
                
                                            
            
            allHistos[ 'dataHisto' ] = dataHistos[ dataLabel ].Clone()
            allHistos[ 'dataHistoGenBin' ] = dataHistos[ dataLabel+'_genBin'].Clone() if ('tau' in ivar) else None
            
        
        print ('|------> Drawing Data/MC comparison:')
        if process.startswith('data'):

            drawDataMCReco(ivar=ivar, 
                           selection=sel, year=year,lumi=lumi, process=process,
                           dataJetHisto=allHistos[ 'dataHisto' ].Clone(),# if not('tau'in ivar ) else allHistos[ 'dataHistoGenBin' ].Clone(),
                           nominal_recoJetHisto = signalHistos[ sigLabel ].Clone(),# if not('tau'in ivar ) else signalHistos[ sigLabel+'_genBin' ].Clone(),
                           alt0_recoJetHisto =  alt0SignalHistos[ alt0sigLabel ].Clone(),# if not('tau'in ivar ) else alt0SignalHistos[ alt0sigLabel+'_genBin' ].Clone(),
                           alt1_recoJetHisto =  alt1SignalHistos[ alt1sigLabel ].Clone(),# if not('tau'in ivar ) else alt1SignalHistos[ alt1sigLabel+'_genBin' ].Clone(),
                           alt2_recoJetHisto =  alt2SignalHistos[ alt2sigLabel ].Clone(),# if not('tau'in ivar ) else alt2SignalHistos[ alt2sigLabel+'_genBin' ].Clone(),
                           labelX=variables[ivar]['label'],
                           jetType=jetType,
                           maxX=variables[ivar]['bins'][-1],
                           tlegendAlignment=variables[ivar]['alignLeg'],
                           outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_Reco_DataMC_'+version+'.'+ext,
                           log=True if 'pt' in ivar else False
                       )
        
        ######### Saving Histograms
        def renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                if isam==None or hist==None: 
                    print(isam,hist)
                    continue
                
                ihis = hist.Clone()
                ihis.SetName(isam)
                ihis.SetTitle(isam)
                ihis.Write()

        outputRootName = outputDir+'/DataMC_main_'+signalLabel+'_alt0_'+alt0SignalLabel+f'_{jetType}.root'
        print ('|------> Saving histograms in rootfile: ', outputRootName)
        outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
        renamingHistos( signalHistos )
        if not('self' in process.lower()):
            renamingHistos( alt0SignalHistos )
            renamingHistos( alt1SignalHistos )
            renamingHistos( alt2SignalHistos )
        
        renamingHistos( dataHistos )
        renamingHistos( allHistos )
        
        outputRoot.Close()

        #print ('|------> Saving histograms in yodafile: ', outputRootName.replace('.root', '.yoda'))
        #histToYoda = [  yoda.root.to_yoda( allHistos [ 'unfoldHisto'+ivar ] ) ]
        #yoda.writeYODA( histToYoda, outputRootName.replace('.root', '.yoda') )

        
        
        
def loadHistograms(samples, var, sel, sysUnc=[],
                   isMC=True, addGenInfo=True, respOnly=False, lumi=1., noResp=False,
                   variables={}, year='2017', process='data', noRebin=False, outputFolder=None,jetType='Central' ):
    """docstring for loadHistograms"""

    if sysUnc==[]: SYSUNC = [ '_nom' ] 
    else: SYSUNC = [ s+u for u in ['Up', 'Down'] for s in sysUnc if not s.startswith(('_model', '_hdamp', '_Tune', '_CR', '_erdON', '_mtop')) ]
    flip = False
    tmpSYSUNC={}
    #print(sysUnc,SYSUNC,flip)
    allHistos = {}
    for isam in samples:
        if sysUnc!=[]:
            for i in (sysUnc):#[ '_jer', '_isrWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_jes' ]):
                if i in isam:
                    flip=True
                    tmpSYSUNC = [i+u for u in ['Up','Down']]
                    continue
        #print(SYSUNC,tmpSYSUNC)
        if not ('nPV'  in var):
            if not flip: tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC]
            else: tmpList = ['reco'+var+syst+sel for syst in tmpSYSUNC]
            #print (tmpList)
            #if isMC and addGenInfo: tmpList = tmpList + [ 'gen'+var+sel ] + [ 'respWithMiss'+var+syst+sel for syst in SYSUNC ]

            if isMC and addGenInfo and not flip:
                tmpList = tmpList + [ 'gen'+var+syst+sel for syst in SYSUNC if 'nom' in syst] 
                if not(noResp):tmpList = tmpList + [ 'accepgen'+var+syst+sel for syst in SYSUNC ]
                if not(noResp):tmpList = tmpList + [ 'truereco'+var+syst+sel for syst in SYSUNC ]
                #tmpList = tmpList + [ 'reco'+var+syst+sel for syst in SYSUNC ]
                if not(noResp): tmpList = tmpList + [ 'respWithMiss'+var+syst+sel for syst in SYSUNC]

            elif isMC and addGenInfo and flip: 
                tmpList = tmpList + [ 'gen'+var+'_nom'+sel]
                #tmpList = tmpList + ['reco'+var+syst+sel for syst in tmpSYSUNC]
                if not(noResp): tmpList = tmpList + ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC]
                #tmpList = tmpList + ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC]

            if respOnly and not flip: 
                if not(noResp): tmpList = [ 'respWithMiss'+var+syst+sel for syst in SYSUNC] 
                #print (SYSUNC,tmpList)
            elif respOnly and flip: 
                if not(noResp): tmpList = ['respWithMiss'+var+syst+sel for syst in tmpSYSUNC ] #+ [ 'missgen'+var+sel ]
        elif 'nPV' in var: 
            tmpList = [var+'_nom'+sel]
            print(var,[syst for syst in tmpSYSUNC],sel,isam,tmpList,samples[isam][0])
           
        print (f'Processing {var,isam+":",tmpList},filename={samples[isam][0]}')
        for ih in tmpList:
            #print(ih,samples[isam][0])
            #print (f'Processing {isam} {ih},{isam+"_"+ih}')
            print (samples[isam][0],ih)            
            if isMC:
                iFile = ROOT.TFile.Open(samples[isam][0],'r')
                allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone() #'jetObservables/'+
                allHistos[isam+'_'+ih].SetDirectory(0)
                allHistos[isam+'_'+ih].Sumw2()
                
                iFile.Close()
                #tmpIsam = 'TT' if isam.startswith('data') else isam
                MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1]['nGenWeights']
                #print (samples[isam])
                if MCScale!=samples[isam][1]['MCScaling']: 
                    print (year,'MCScale mismatch alert', MCScale, samples[isam][1]['MCScaling'])
                else: allHistos[isam+'_'+ih].Scale( MCScale )
            else:
                if sel.startswith('_dijet'):
                    tmpdataHistos = {}
                    #print("Working on trigger separated histos for data")
                    itempFile=ROOT.TFile.Open(samples[isam][0],'r')
                    for it in checkDict( 'JetHT', dictSamples )[year]['triggerList']:
                        
                        print(it,ih,ih.replace( sel, '_'+it+sel ))#,itempFile,samples[isam][0])
                        tmpdataHistos[ it ] = itempFile.Get( ih.replace( sel, '_'+it+sel )).Clone()
                        tmpdataHistos[it].SetDirectory(0)
                        tmpdataHistos[ it ].Scale( checkDict( 'JetHT', dictSamples )[year]['triggerList'][it] )
                        #print(it,ih,checkDict( 'JetHT', dictSamples )[year]['triggerList'][it] )
                        
                    itempFile.Close() 
                    #print(tmpdataHistos,isam+'_'+ih,next(iter(tmpdataHistos)))
                    allHistos[ isam+'_'+ih ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
                    allHistos[ isam+'_'+ih ].Reset()
                    for i in tmpdataHistos: 
                        tmpdataHistos[i].SetDirectory(0)
                        allHistos[isam+'_'+ih].Add( tmpdataHistos[i].Clone() )
                        allHistos[isam+'_'+ih].SetDirectory(0)
                        
                        
                else:
                    allHistos[isam+'_'+ih] = samples[isam][0].Get( ih )


        
    #print(allHistos) 
    if sel.startswith('_dijet'):
        if isMC:
            c=0
            allHistos_upd={}
            for sys in SYSUNC:

                print (f"Adding together histograms for {var} in {sel} for {sys}")
                tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k and sys in k}
                #print(tmpHistos.keys())#,allHistos.keys())
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
                #if 'resp' in ih: 
                #print(ih)
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
                        elif 'genJet' in ih: 
                            allHistos[ih].Rebin( genBin )
                        elif 'nPV' in ih:
                            allHistos[ih].Rebin(recoBin)
                            allHistos[ih+'_genBin']=allHistos[ih].Clone()
                            allHistos[ih+'_genBin'].Rebin(genBin)

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
                        tmpHisto = rebin_RM_withUF(allHistos[ih],genBin,recoBin).Clone(allHistos[ih].GetName()+'_Rebin')
                        #make_rebinned_2d_hist(allHistos[ih].Clone(), new_bin_edge_pairs,True)
                        tmpHisto.Sumw2()                        

                        tmpHisto.SetDirectory(0)
                        allHistos[ih] = copy.deepcopy(tmpHisto.Clone())
                        
        
        #if noRebin: 
        #    allHistosUpd = copy.deepcopy(allHistos)
        #    for ih in allHistos:
        #        if 'recoJet' in ih:
        #            allHistosUpd[ih+'_genBin'] = allHistos[ih].Clone()
        #            allHistosUpd[ih+'_genBin'].Rebin(2)
        #    allHistos.update(allHistosUpd)
        if samples and outputFolder: 
            outputRootName = outputFolder.split('jetObservables_histograms')[0] + f"/{sel.replace('_','')}/{year}/DataMC/" + '/loadedHistograms_'+isam+var+year+f'_{jetType}.root' 
            print ('|------> Saving loaded and rebinned histograms in rootfile: ', outputRootName)#, allHistos)
            def renamingHistos( dictHistos ):
                for isam, hist in dictHistos.items():
                    if isam==None or hist==None: 
                        print(isam,hist)
                        continue
                    ihis = hist.Clone()
                    ihis.SetName(isam)
                    ihis.SetTitle(isam)
                    ihis.Write()
                #print("finished renaming")
            outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
            renamingHistos( allHistos) 
        
            #tunfolder.Write()
            outputRoot.Close()    
    #for isam in samples: samples[isam][0].Close()
    #pprint.pprint(allHistos)
    
    return allHistos
##############################################################################################
        
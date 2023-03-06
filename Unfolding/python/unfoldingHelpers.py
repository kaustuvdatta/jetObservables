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
from root_numpy import array2hist, hist2array
import histoHelpers
from histoHelpers import *
import os
import glob
import sys
import math
#sys.path.insert(0,'../test/')
#from Drawhistogram_dijetSel import *
from datasets_dijetSel_RunIISummer20UL_nomWts import dictSamples, checkDict
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
#import pdb
lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)
from histoHelpers import *
##########################################################################

######################## MAIN Unfolding script ###########################
#1 load histos, #2 unfold 
##########################################################################

def loadHistograms(samples, var, sel, sysUnc=[], 
                   isMC=True, addGenInfo=True, respOnly=False, lumi=1., 
                       variables={}, year='2017', process='data', noRebin=False, outputFolder=None ):
    """docstring for loadHistograms"""

    if sysUnc==[]: SYSUNC = [ '_nom' ] 
    else: SYSUNC = [ s+u for u in ['Up', 'Down'] for s in sysUnc if not s.startswith(('_model', '_hdamp', '_Tune', '_CR', '_erdON', '_mtop')) ]
    flip = False
    
    allHistos = {}
    for isam in samples:
        if sysUnc!=[]:
            for i in ([ '_jer', '_isrWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_jes' ]):
                if i in isam:
                    flip=True
                    tmpSYSUNC = [i+u for u in ['Up','Down']]
                    continue
                    
        if not flip: tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC]
        else: tmpList = ['reco'+var+syst+sel for syst in tmpSYSUNC]
        #print (tmpList)
        #if isMC and addGenInfo: tmpList = tmpList + [ 'gen'+var+sel ] + [ 'resp'+var+syst+sel for syst in SYSUNC ]
        
        if isMC and addGenInfo and not flip:
            tmpList = tmpList + [ 'gen'+var+syst+sel for syst in SYSUNC if 'nom' in syst] 
            tmpList = tmpList + [ 'accepgen'+var+syst+sel for syst in SYSUNC ]
            tmpList = tmpList + [ 'truereco'+var+syst+sel for syst in SYSUNC ]
            tmpList = tmpList + [ 'resp'+var+syst+sel for syst in SYSUNC]
        
        elif isMC and addGenInfo and flip: 
            tmpList = tmpList + [ 'gen'+var+'_nom'+sel]
            for syst in tmpSYSUNC:
                tmpList = tmpList + ['resp'+var+syst+sel]
        
        if respOnly and not flip: 
            tmpList = [ 'resp'+var+syst+sel for syst in SYSUNC] 
            #print (SYSUNC,tmpList)
        elif respOnly and flip: 
            tmpList = ['resp'+var+syst+sel for syst in tmpSYSUNC ] #+ [ 'missgen'+var+sel ]
            #print (tmpSYSUNC,tmpList)
        #print (tmpList)
        #tmpList.sort(reverse=True)
        print (tmpList)
        print (f'Processing {isam}')
        for ih in tmpList:
            #print(ih)
            #print (f'Processing {isam} {ih},{isam+"_"+ih}')
            #print (samples[isam][0])
            if isMC:
                iFile = ROOT.TFile.Open(samples[isam][0],'r')
                allHistos[isam+'_'+ih] = iFile.Get( ih ).Clone() #'jetObservables/'+
                allHistos[isam+'_'+ih].SetDirectory(0)
                iFile.Close()
                #tmpIsam = 'TT' if isam.startswith('data') else isam
                MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1]['nGenWeights']
                #print (samples[isam])
                if MCScale!=samples[isam][1]['MCScaling']: 
                    print (year,'MCScale mismatch alert', MCScale, samples[isam][1]['MCScaling'])
                else: allHistos[isam+'_'+ih].Scale( MCScale )
            else:
                if sel.startswith('_dijet'):# and process.startswith('data'):
                    tmpdataHistos = {}
                    for it in checkDict( 'JetHT', dictSamples )[year]['triggerList']:
                        tmpdataHistos[ it ] = samples[isam][0].Get( ih.replace( sel, '_'+it+sel ) )
                        tmpdataHistos[ it ].Scale( checkDict( 'JetHT', dictSamples )[year]['triggerList'][it] )
                    allHistos[ isam+'_'+ih ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
                    allHistos[ isam+'_'+ih ].Reset()
                    for i in tmpdataHistos: allHistos[isam+'_'+ih].Add( tmpdataHistos[i] )
                else:
                    allHistos[isam+'_'+ih] = samples[isam][0].Get( ih )

            if not noRebin:
                if len(variables[var]['bins'])==1:
                    genBin = variables[var]['bins'][0]
                    recoBin = variables[var]['bins'][0]/2
                else:
                    genBin = variables[var]['bins']
                    recoBin = variables[var]['bins_reco']#np.sort(np.append( variables[var]['bins'], np.array([ (variables[var]['bins'][i]+variables[var]['bins'][i+1])/2 for i in range(len(variables[var]['bins'])-1) ]  ) ))
                    #print(genBin,recoBin)
                if not ih.startswith('resp'):
                    if len(variables[var]['bins'])==1:
                        if ih.startswith(('reco','fake', 'true')):
                            allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih].Clone()
                            allHistos[isam+'_'+ih+'_genBin'].Rebin( genBin )
                            allHistos[isam+'_'+ih].Rebin( recoBin )
                        else: allHistos[isam+'_'+ih].Rebin( genBin )
                    else:
                        if ih.startswith(('reco','fake','true')):
                            allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih].Clone()
                            allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih+'_genBin'].Rebin( len(genBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin_genBin", array( 'd', genBin ) )
                            allHistos[isam+'_'+ih] = allHistos[isam+'_'+ih].Rebin( len(recoBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin", array( 'd', recoBin ) )
                        elif ih.startswith(('accep','miss','gen')):
                            allHistos[isam+'_'+ih] = allHistos[isam+'_'+ih].Rebin( len(genBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin", array( 'd', genBin ) )
                else:
                    if len(variables[var]['bins'])==1: allHistos[isam+'_'+ih].Rebin2D( genBin,recoBin )
                    else:
                        
                        #print(genBin,recoBin)
                        genBin = [(i*1000)/1000 for i in genBin]#*1000array('d',genBin)
                        recoBin = [(i*1000)/1000 for i in recoBin]##array('d',recoBin)
                        
                        print(genBin,recoBin)
                        
                        axis = allHistos[isam+'_'+ih].GetYaxis()
                        old_reco_bin_edges = [axis.GetBinLowEdge(i) for i  in range(1,axis.GetNbins()+2)]
                        new_reco_bin_edges = recoBin
                        new_bin_edge_pairs = get_bin_edge_pairs(new_bins=new_reco_bin_edges,old_reco_bin_edges=old_reco_bin_edges)

                        print(new_bin_edge_pairs,recoBin[len(recoBin)-1],recoBin[len(recoBin)-2])
                        
                        tmpHisto = make_rebinned_2d_hist(allHistos[isam+'_'+ih].Clone(), new_bin_edge_pairs,True)
                        
                        for i in range(tmpHisto.GetNbinsX()+2):
                            for j in range(tmpHisto.GetNbinsY()+2):
                                if i == 0 or j == 0:
                                    tmpHisto.SetBinContent(i, j, allHistos[isam+'_'+ih].GetBinContent(i, j))
                                    tmpHisto.SetBinError(i, j, allHistos[isam+'_'+ih].GetBinError(i, j))
                                elif i == tmpHisto.GetNbinsX()+1 or j == tmpHisto.GetNbinsY()+1:
                                    tmpHisto.SetBinContent(i, j, allHistos[isam+'_'+ih].GetBinContent((allHistos[isam+'_'+ih].GetNbinsX()+1), (allHistos[isam+'_'+ih].GetNbinsY()+1)))
                                    tmpHisto.SetBinError(i, j, allHistos[isam+'_'+ih].GetBinError((allHistos[isam+'_'+ih].GetNbinsX()+1), (allHistos[isam+'_'+ih].GetNbinsY()+1)))
                                    
                        tmpHisto.Sumw2()
                        tmpHisto.SetDirectory(0)
                        allHistos[isam+'_'+ih] = tmpHisto
                        
                        #print(f"done rebinning {allHistos[isam+'_'+ih].GetName()}")
                        
                    #if isMC: allHistos[isam+'_'+ih].Scale( MCScale )
                    ##### For tests, projections directly from 2D
                    #if not respOnly:
                    #   
                    #    if not genBin: genBin = variables[var]['bins']

                    #    allHistos[isam+'_'+ih.replace('resp', 'accepgen')] = copy.deepcopy(allHistos[isam+'_'+ih].ProjectionX('px',firstybin=1))
                    #    allHistos[isam+'_'+ih.replace('resp', 'truereco')] = copy.deepcopy(allHistos[isam+'_'+ih].ProjectionY('py',firstxbin=1))
                    #    genBin=array('d',genBin)
                    #    recoBin=array('d',recoBin)
                        
                    #allHistos[isam+'_'+ih.replace('resp', 'truereco')+'_genBin']=allHistos[isam+'_'+ih.replace('resp', 'truereco')].Rebin(len(genBin)-1,'_rebinned', genBin)
                    #print (allHistos)
                    #allHistos[isam+'_'+ih.replace('resp', 'missgen')] = allHistos[isam+'_'+ih.replace('resp', 'gen')].Clone()
                    #allHistos[isam+'_'+ih.replace('resp', 'missgen')].Add(allHistos[isam+'_'+ih.replace('resp', 'accepgen')],-1)
                    #allHistos[isam+'_'+ih.replace('resp', 'fakereco')+'_genBin'] = allHistos[isam+'_'+ih.replace('resp', 'reco')].Clone()
                    #allHistos[isam+'_'+ih.replace('resp', 'fakereco')+'_genBin'].Add(allHistos[isam+'_'+ih.replace('resp', 'truereco')+'_genBin'],-1)
        
        def renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                ihis = hist.Clone()
                ihis.SetName(isam)
                ihis.SetTitle(isam)
                ihis.Write()

        if outputFolder: 
            outputRootName = outputFolder.split('jetObservables_histograms')[0]+'/loadedHistograms_main_'+isam+year+'.root'
            #print ('|------> Saving histograms in rootfile: ', outputRootName)
            outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
            renamingHistos( allHistos )
        
            #tunfolder.Write()
            outputRoot.Close()
        
    if sel.startswith('_dijet'):
        print (f"Adding together histograms for {var} in {sel}")
        tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k }
        for ih in tmpHistos:
            for jh in allHistos:
                if (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ):
                    tmpHistos[ih].Add( allHistos[jh].Clone() )
                    tmpHistos[ih].SetDirectory(0)
        if len(tmpHistos)>0: allHistos = copy.deepcopy(tmpHistos)
            
    #for isam in samples: samples[isam][0].Close()
    pprint.pprint(allHistos)
    
    return allHistos

############Unfolder below#################



def runTUnfold(
                dataFile, sigFiles, bkgFiles, variables, sel, sysUncert, process, ext, lumi=1.,
                year='2017',runMLU=False, sysSigFiles=[], varSigFiles=[], outputFolder='../Results',
                version='_DFeb23', mainMC='HTbin', altMC='Ptbin'
              ):
    
    
    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 41, 45, 48, 49, 37, 38, 33, 17]

    if mainMC.startswith('Ptbin'):
        signalLabelBegin = 'QCD_Pt_'
        signalLabel = 'QCD_Pt_3200toInf'
    elif mainMC.startswith('HTbin'):
        signalLabelBegin = 'QCD_HT'
        signalLabel = 'QCD_HT2000toInf'
    elif mainMC.startswith('herwig'):
        signalLabelBegin = 'QCD_Pt-'
        signalLabel = 'QCD_Pt-150to3000'
        
    if altMC.startswith('Ptbin'):
        altSignalLabelBegin = 'QCD_Pt_'
        altSignalLabel = 'QCD_Pt_3200toInf'
    elif altMC.startswith('HTbin'):
        altSignalLabelBegin = 'QCD_HT'
        altSignalLabel = 'QCD_HT2000toInf'
    elif altMC.startswith('herwig'):
        altSignalLabelBegin = 'QCD_Pt-'
        altSignalLabel = 'QCD_Pt150to3000'
        #altSignalLabel = 'QCD_Pt-15to7000'

    sysSignalLabelBegin = 'sysQCD_' 
    selection=sel
    for ivar in variables:
        #print (ivar)
        outputDir=outputFolder+sel.split('_')[1]+'/'+year+'/Unfolding/'+ivar+'/'+process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)
    ###########################################################################################
        if year.startswith('all'):
            signalHistos = {
                    signalLabel+'_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_resp'+ivar+'_nom'+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel),
                    signalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel),
                    signalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel),
                    signalLabel+'_accepgen'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin': dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin')
                    }
            signalHistos[ signalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_resp'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )
            
            signalHistos[ signalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )
            
            signalHistos[ signalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_truereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_truereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_fakereco'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_accepgen'+ivar+'_nom'+sel) )
            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )            
            if not process.startswith('MCSelfClosure'):
                sysSignalHistos={}

                #print ("Loading up all syst. variations from the following:", sysUncert)
                for sys in sysUncert:
                    for upDown in [ 'Up', 'Down' ]:
                        #print (sys)
    
                        if sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_Tune')): continue
                        s = [i for i in sysSignalLabels if sys in i]
                        #print (s)
                        if ('2016' in sys or '2017' in sys or '2018' in sys) and len(s)>1:
                            if '2016' in s: s=[s[1]]
                            elif '2017' in s: s=[s[2]]
                            elif '2018' in s: s=[s[3]]
                        else:
                            s=[s[0]]
                        #print (s[0]+'_reco'+ivar+sys+upDown+sel)
                        
                        if not process.startswith('MCCrossClosure') and not '2017' in sys and not '2018' in sys and not '2016' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_resp'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(s[0]+'_resp'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(s[0]+'_resp'+ivar+sys+upDown+sel) )
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(s[0]+'_resp'+ivar+sys+upDown+sel) )

                        
                        # dealing with uncorrelated jes unc sources below
                        elif not process.startswith('MCCrossClosure') and '2016' in sys and not '2017' in sys and not '2018' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016_preVFP'].Get(s[0]+'_resp'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            
                            #ensuring pre and post VFP periods are treated differently

                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2016'].Get(s[0]+'_resp'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            

                        elif not process.startswith('MCCrossClosure') and '2017' in sys and not '2016' in sys and not '2018' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(s[0]+'_resp'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                         
                         elif not process.startswith('MCCrossClosure') and '2018' in sys and not '2017' in sys and not '2016' in sys: 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_reco'+ivar+sys+upDown+sel) 
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2018'].Get(s[0]+'_resp'+ivar+sys+upDown+sel)
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016_preVFP'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2016'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                            sysSignalHistos[ s[0]+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+'_nom'+sel))
                         
                        ################ added recohistos & resp matrices from nominal_2018(/2017) and jesUncorrUnc_2017(/2018) ########################
                        

                altSignalHistos = {
                    altSignalLabel+'_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel),
                    altSignalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel),
                    altSignalLabel+'_truereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_truereco'+ivar+'_nom'+sel),
                    altSignalLabel+'_fakereco'+ivar+'_nom'+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel),
                    altSignalLabel+'_accepgen'+ivar+sel : dataFile[ivar+'_2016_preVFP'].Get(altSignalLabel+'_accepgen'+ivar+sel),
                    }

                altSignalHistos[ altSignalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+sel].Add( dataFile[ivar+'_2016'].Get(altSignalLabel+'_accepgen'+ivar+sel) )
                
                altSignalHistos[ altSignalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+sel].Add( dataFile[ivar+'_2017'].Get(altSignalLabel+'_accepgen'+ivar+sel) )
                
                altSignalHistos[ altSignalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_reco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_fakereco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_fakereco'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_accepgen'+ivar+sel].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_accepgen'+ivar+sel) )
                
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
                    #print ("Processing signal variations from amongst the foll.:", sysUncert)

                    for sys in sysUncert:
                        #print (sys)
                        s = [i for i in varSignalLabels if (sys.split('_')[1] in i)]
                        for j in s:
                            #print (j+'_reco'+ivar+'_nom'+sel)
                            varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ] = dataFile[ivar+'_2017'].Get(j+'_reco'+ivar+'_nom'+sel)
                            varSignalHistos[ j+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_reco'+ivar+'_nom'+sel) )
                            varSignalHistos[ j+'_resp'+ivar+'_nom'+sel ] = dataFile[ivar+'_2017'].Get(j+'_resp'+ivar+'_nom'+sel)
                            varSignalHistos[ j+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(j+'_resp'+ivar+'_nom'+sel) )

                
                
            allHistos = {
                    'dataHisto' : dataFile[ivar+'_2016_preVFP'].Get( 'dataHisto' ),
                    'dataHistoGenBin' : dataFile[ivar+'_2016'].Get( 'dataHistoGenBin' )
                    }
            allHistos[ 'dataHisto' ].Add( dataFile[ivar+'_2018'].Get( 'dataHisto' ) )
            allHistos[ 'dataHistoGenBin' ].Add( dataFile[ivar+'_2018'].Get( 'dataHistoGenBin' ) )

            dataHistos = { }
            bkgHistos = { }

        else:
            print('|-------> Running single year '+year)
            ### Getting input histos
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            print (ivar,sel)
            #pprint.pprint(mainSigFiles)
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[], respOnly=False, lumi=lumi, year=year, process=process, variables=variables,outputFolder='../Results' )
            tmp2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(altSignalLabelBegin)  }
            if not process.startswith('MCSelfClosure'): altSignalHistos = loadHistograms( tmp2SigFiles, ivar, sel, sysUnc=[], isMC=True, respOnly=False, lumi=lumi, year=year, process=process, variables=variables)
            if not process.startswith('MCSelfClosure'): sysSignalHistos = loadHistograms( sysSigFiles, ivar, sel, sysUnc=sysUncert, respOnly=False, isMC=True, lumi=lumi, year=year, process=process, variables=variables)
            if process.startswith('data') and selection.startswith(('_W','_top')): varSignalHistos = loadHistograms( varSigFiles, ivar, sel, sysUnc=[], respOnly=False, isMC=True, lumi=lumi, year=year, process=process, variables=variables)
            bkgHistos = loadHistograms( bkgFiles, ivar, sel, sysUnc=[], lumi=lumi, year=year, process=process, variables=variables) if process.startswith('data') else {}

            #print ("Alt histos:", bkgHistos)
            
            ####### Fix fake and true reco
            signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone()
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Add( signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel], -1 )

            signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionX()
            signalHistos[signalLabel+'_missgen'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_gen'+ivar+'_nom'+sel].Clone()
            signalHistos[signalLabel+'_missgen'+ivar+'_nom'+sel].Add( signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel], -1 )
            
            if not process.startswith('MCSelfClosure'):
                #print(altSignalHistos)
                altSignalHistos[altSignalLabel+'_truereco'+ivar+'_nom'+sel]= altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
                altSignalHistos[altSignalLabel+'_fakereco'+ivar+'_nom'+sel] = altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel].Clone()
                altSignalHistos[altSignalLabel+'_fakereco'+ivar+'_nom'+sel].Add( altSignalHistos[altSignalLabel+'_truereco'+ivar+'_nom'+sel], -1 )

                altSignalHistos[altSignalLabel+'_accepgen'+ivar+'_nom'+sel] = altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].ProjectionX()
                altSignalHistos[altSignalLabel+'_missgen'+ivar+'_nom'+sel] = altSignalHistos[altSignalLabel+'_gen'+ivar+'_nom'+sel].Clone()
                altSignalHistos[altSignalLabel+'_missgen'+ivar+'_nom'+sel].Add( altSignalHistos[altSignalLabel+'_accepgen'+ivar+'_nom'+sel], -1 )          
              
            if sel.startswith('_dijet'): 
                dataHistostrue = loadHistograms( dataFile, ivar, sel, isMC= False, lumi=lumi, year=year, process=process, variables=variables )
                ######## Cross check: plotting data vs all MC Scaled
                print ('|------> Cross check: plotting data vs all MC')
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
                    genBin = variables[ivar]['bins']
                    dataHistos = { 'data_reco'+k.split(('_truereco'))[1] : v for (k,v) in signalHistos.items() if ('_truereco' in k and not ('genBin' in k ))}
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                else:
                    genBin = variables[ivar]['bins']
                    dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, lumi=lumi, year=year, process=process, variables=variables )
            
                #print ('|------> Cross check: plotting data vs all MC')
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
            scaleFactor =1
            if sel.startswith('_dijet'):
                scaleFactor = dataHistostrue['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto' ].Integral()
                scaleFactorGenBin = dataHistostrue['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' ].Integral()
                allHistos[ 'allMCHisto' ].Scale( scaleFactor )
                if process=='data':
                    allHistos[ 'allBkgHisto' ].Scale( scaleFactor )
                    allHistos[ 'allBkgHistoGenBin' ].Scale( scaleFactorGenBin )
                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        #print (ihsig)
                        signalHistos[ihsig].Scale( scaleFactor )
                    if ihsig.endswith('genBin'): signalHistos[ihsig].Scale( scaleFactorGenBin )
                    #if 'resp' in ihsig: signalHistos[ihsig].Scale( scaleFactor )   ### dont notice the difference
                if not process.startswith('MC'):
                    for ihsig in sysSignalHistos:
                        if ihsig.endswith(sel):
                            #print (ihsig)
                            sysSignalHistos[ihsig].Scale( scaleFactor )
                        if ihsig.endswith('genBin'): sysSignalHistos[ihsig].Scale( scaleFactorGenBin )
            
                if process.startswith("MC"):
                    genBin = variables[ivar]['bins']
                    dataHistos = { 'data_reco'+k.split(('_truereco'))[1] : v for (k,v) in signalHistos.items()  if ('_truereco' in k and not ('genBin' in k ))}
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
                    dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Rebin( len(genBin)-1, 'data_reco'+ivar+'_nom'+sel+'_genBin', array( 'd', genBin ) )
                else: 
                    dataHistos = dataHistostrue
                plotSimpleComparison( dataHistos['data_reco'+ivar+'_nom'+sel].Clone(), 'data', allHistos[ 'allMCHisto' ].Clone(), 'allBkgs', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_nom", rebinX=1, version=sel+'_'+version, outputDir=outputDir )

            print ('|------> Unfolding '+ivar)

            
            ####### Cross check response matrix
            tmpGenHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionX()
            plotSimpleComparison( tmpGenHisto, 'projection', signalHistos[signalLabel+'_accepgen'+ivar+'_nom'+sel].Clone(), 'Regular AccepGen', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionGen", rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            tmpRecoHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            plotSimpleComparison( tmpRecoHisto, 'projection', signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Clone(), 'Regular TrueReco', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionReco", rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            tmpRecoHisto.Scale( scaleFactor )
            
            ####### plotting Removing of bkgs from data
            fakeHistos = { k:v for (k,v) in signalHistos.items()  if ('_fakereco' in k and not ('genBin' in k ))}
            for ih in fakeHistos:
                if ih.endswith(ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( fakeHistos[ih] )
                if ih.endswith(ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( fakeHistos[ih] )
            plotSimpleComparison( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone(), 'data', allHistos[ 'allBkgHisto' ].Clone(), 'Bkg+fakes', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataBkgFakes", rebinX=1, version=sel+'_'+version, outputDir=outputDir )

            allHistos[ 'dataHisto' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone()
            #allHistos[ 'dataHisto' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            #allHistos[ 'dataHisto' ].Scale( 1/allHistos[ 'dataHisto' ].Integral() )
            allHistos[ 'dataHistoGenBin' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            #allHistos[ 'dataHistoGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )
            #allHistos[ 'dataHistoGenBin' ].Scale( 1/allHistos[ 'dataHistoGenBin' ].Integral() )

        if process.startswith('data') and not year.startswith('all'):
            allHistos[ 'dataMinusBkgs' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone()
            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            #allHistos[ 'dataMinusBkgs' ].Scale( 1/allHistos[ 'dataMinusBkgs' ].Integral() )
            allHistos[ 'dataMinusBkgsGenBin' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )
            #allHistos[ 'dataMinusBkgsGenBin' ].Scale( 1/allHistos[ 'dataMinusBkgsGenBin' ].Integral() )

            tmpHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            #tmpHisto.Scale( 1/tmpHisto.Integral() )

            tmpHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            #tmpHisto.Scale( 1/tmpHisto.Integral() )
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal true reco', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", rebinX=1, version=sel+'_'+version, outputDir=outputDir )
        
        elif process.startswith('data'):# and year.startswith('all'):
            
            print ('|------> Cross check: plotting data vs all MC')
            allHistos[ 'allBkgHisto' ] = allHistos['dataHisto'].Clone()
            allHistos[ 'allBkgHisto' ].Reset()
            allHistos[ 'allBkgHistoGenBin' ] = allHistos['dataHistoGenBin'].Clone()
            allHistos[ 'allBkgHistoGenBin' ].Reset()
            if process.startswith('data'):
                for ibkg in bkgHistos:
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )
                    if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone().Rebin( len(genBin)-1, 'ibkg_rebin', array( 'd', genBin ) ) )
            allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
            allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Clone() )
            allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
            allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone() )

            allHistos[ 'dataMinusBkgs' ] = allHistos[ 'dataHisto' ].Clone()
            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            #allHistos[ 'dataMinusBkgs' ].Scale( 1/allHistos[ 'dataMinusBkgs' ].Integral() )
            allHistos[ 'dataMinusBkgsGenBin' ] = allHistos[ 'dataHistoGenBin'].Clone()
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )
            #allHistos[ 'dataMinusBkgsGenBin' ].Scale( 1/allHistos[ 'dataMinusBkgsGenBin' ].Integral() )

            
            tmpHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            #tmpHisto.Scale( 1/tmpHisto.Integral() )
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs', tmpHisto.Clone(), 'signal true reco', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        #print(sysSignalHistos)
            
        ######## Cross check: plotting response matrix
        print ('|------> Cross check: plotting response matrix for signal')
        ROOT.gStyle.SetPadRightMargin(0.15)
        #ROOT.gStyle.SetPalette(ROOT.kGistEarth)
        #ROOT.TColor.InvertPalette()
        can2D = ROOT.TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetXaxis().SetTitle('Accepted Gen '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitle('True Reco '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitleOffset( 0.8 )
        #signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetRange( variables[ivar]['bins'][0], variables[ivar]['bins'][-1] )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].Draw("colz")
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
        
        CMS_lumi.relPosX = 0.12
        CMS_lumi.CMS_lumi(can2D, 4, 0)
        can2D.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+sel+'_responseMatrix'+version+'.'+ext)

        

        ######## TUnfold part
        print ('|------> TUnfolding starts:')

        ##### Defining options for TUnfold
        tunfolder = ROOT.TUnfoldDensity(
                                            signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel], ### response matrix. According to TUnfold, this distribution does NOT have to be normalized
                                            ROOT.TUnfold.kHistMapOutputHoriz,  #### kHistMapOutputVert if x->reco and y->gen, kHistMapOutputHoriz if x->gen and y->reco
                                            ROOT.TUnfold.kRegModeCurvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                            ROOT.TUnfold.kEConstraintNone,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. (Need to look into this!!)
                                            ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                            )

        ##### Defining input (data recoJet )
        print ('|------> TUnfolding adding input:')
        #tunfolder.SetInput( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone() )
        #print (allHistos.keys())
        tunfolder.SetInput( allHistos[ 'dataHisto' ])
        
        if process.startswith('MCCrossClosure'): 
            tunfolder_cross = ROOT.TUnfoldDensity(
                                            altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel], ### response matrix. According to TUnfold, this distribution does NOT have to be normalized
                                            ROOT.TUnfold.kHistMapOutputHoriz,  #### kHistMapOutputVert if x->reco and y->gen, kHistMapOutputHoriz if x->gen and y->reco
                                            ROOT.TUnfold.kRegModeCurvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                            ROOT.TUnfold.kEConstraintNone,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. (Need to look into this!!)
                                            ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                            )

            ##### Defining input (data recoJet )
            print ('|------> TUnfolding adding input:')
            #tunfolder.SetInput( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone() )
            #print (allHistos.keys())
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
                
                if sys.startswith(('_jer', '_isrWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_jes')):
                    #print('|------> TUnfolding adding %sUnc'%sys)
                    
                    s = [i for i in sysSignalLabels if sys in i]
                    print (s)
                    if ('2016' in sys or '2017' in sys or '2018' in sys) and len(s)>1:
                        if '2016' in s and 'VFP' in year: s=[s[1]]
                        elif '2016' in s and year.endswith('2016'): s=[s[1]]
                        elif '2017' in s: s=[s[2]]
                        elif '2018' in s: s=[s[3]]
                    else:
                        s=[s[0]]
                    
                    dictUncHistos[sys+'Up'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Up'+sel].Clone()
                    dictUncHistos[sys+'Down'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Down'+sel].Clone()
                    for upDown in [ 'Up', 'Down' ]:
                        #print (sys+upDown)
                        tunfolder.AddSysError(
                                            sysSignalHistos[s[0]+'_resp'+ivar+sys+upDown+sel],
                                            sys+upDown,
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+sys+upDown, ivar+'can2DNorm'+sys+upDown, 750, 500 )
                        sysSignalHistos[s[0]+'_resp'+ivar+sys+upDown+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+s[0]+sel+upDown+'Normalized_responseMatrix'+version+'.'+ext)

                #### adding model uncertainty
                elif sys.startswith(('_model')):
                    #if not process.startswith('MCSelfClosure'):
                    #print('|------> TUnfolding adding modelUnc')
                    tunfolder.AddSysError(
                                        altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel],
                                        'modelUncTotal',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNormAltSignal', ivar+'can2DNormAltSignal', 750, 500 )
                    altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+altSignalLabel+sel+'Normalized_alt_responseMatrix'+version+'.'+ext)
                    dictUncHistos[sys] = altSignalHistos[altSignalLabel+'_reco'+ivar+'_nom'+sel].Clone()
            
            
                elif sys.startswith('_hdamp'): 
                    #if process.startswith('MC'): continue
                    #print('|------> TUnfolding adding hdampUnc')
                    #dictUncHistos[sys+'Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                    #dictUncHistos[sys+'Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_hdampUP_TuneCP5'+'_resp'+ivar+'_nom'+sel],
                                        '_hdampUP',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_hdampUP', ivar+'can2DNorm'+'_hdampUP', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_hdampUP_TuneCP5'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampUP'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_hdampDOWN_TuneCP5'+'_resp'+ivar+'_nom'+sel],
                                        '_hdampDOWN',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_hdampDOWN', ivar+'can2DNorm'+'_hdampDOWN', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_hdampDOWN_TuneCP5'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampDOWN'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_hdampUP'] = varSignalHistos['varTTToSemileptonic_hdampUP_TuneCP5'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_hdampDOWN'] = varSignalHistos['varTTToSemileptonic_hdampDOWN_TuneCP5'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_Tune'):
                    #if process.startswith('MC'): continue
                    #print('|------> TUnfolding adding TuneCP5Unc')
                    #dictUncHistos[sys+'Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                    #dictUncHistos[sys+'Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_resp'+ivar+'_nom'+sel],
                                        '_TuneCP5Up',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_TuneCP5Up', ivar+'can2DNorm'+'_TuneCP5Up', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Up'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_resp'+ivar+'_nom'+sel],
                                        '_TuneCP5Down',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_TuneCP5Down', ivar+'can2DNorm'+'_TuneCP5Down', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Down'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_TuneCP5Up'] = varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_TuneCP5Down'] = varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_CR'):
                    #if process.startswith('MC'): continue
                    #print('|------> TUnfolding adding Colour reconnection Unc')
                    #dictUncHistos[sys+'Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                    #dictUncHistos[sys+'Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_resp'+ivar+'_nom'+sel],
                                        '_CR1',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR1', ivar+'can2DNorm'+'TuneCP5CR1', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR1'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_resp'+ivar+'_nom'+sel],
                                        '_CR2',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5CR2', ivar+'can2DNorm'+'TuneCP5CR2', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5CR2'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_CR1'] = varSignalHistos['varTTToSemileptonic_TuneCP5CR1'+'_reco'+ivar+'_nom'+sel].Clone()
                    dictUncHistos['_CR2'] = varSignalHistos['varTTToSemileptonic_TuneCP5CR2'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_erdON'): 
                    #if process.startswith('MC'): continue
                    #print('|------> TUnfolding adding erdONUnc')
                    #dictUncHistos[sys+'Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                    #dictUncHistos[sys+'Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                    tunfolder.AddSysError(
                                        varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_resp'+ivar+'_nom'+sel],
                                        '_erdON',
                                        ROOT.TUnfold.kHistMapOutputHoriz,
                                        ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                    can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5_erdON', ivar+'can2DNorm'+'TuneCP5_erdON', 750, 500 )
                    varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                    can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5_erdON'+'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                    dictUncHistos['_erdON'] = varSignalHistos['varTTToSemileptonic_TuneCP5_erdON'+'_reco'+ivar+'_nom'+sel].Clone()

                elif sys.startswith('_mtop'): 
                    #if process.startswith('MC'): continue
                    #print('|------> TUnfolding adding mtopUnc')
                    mass_list = [ '171p5','173p5' ] #'166p5',
                    #dictUncHistos[sys+'Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                    #dictUncHistos[sys+'Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()

                    for m in mass_list:
                        tunfolder.AddSysError(
                                             varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_resp'+ivar+'_nom'+sel],
                                             '_mtop%s'%m,
                                         ROOT.TUnfold.kHistMapOutputHoriz,
                                         ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                        )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'_mtop%s_TuneCP5'%m, ivar+'can2DNorm'+'_mtop%s_TuneCP5'%m, 750, 500 )
                        varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_mtop%s_TuneCP5'%m +'_resp'+sel+'Normalized_responseMatrix'+version+'.'+ext)

                        dictUncHistos['_mtop%s'%m] = varSignalHistos['varTTToSemileptonic_mtop%s_TuneCP5'%m+'_reco'+ivar+'_nom'+sel].Clone()

            '''
            ### Making unc plot
            plotSysComparison2( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                dictUncHistos,
                                ivar+'_'+signalLabel+'_allSys',
                                labelX=variables[ivar]['label'],
                                version=sel+'_'+version,
                                year= ( '2017+2018' if year.startswith('all') else year ),
                                outputDir=outputDir
                                )
            '''
        ###### Running the unfolding
        print ('|------> TUnfolding doUnfold:')
        tunfolder.DoUnfold(0)
        if process.startswith('MCCrossClosure'): tunfolder_cross.DoUnfold(0) 
            
        ##### Get output of unfolding ###############################################################
        allHistos [ 'unfoldHisto'+ivar ] = tunfolder.GetOutput("unfoldHisto"+ivar).Clone()
        allHistos [ 'unfoldHisto'+ivar ].Sumw2()
        if process.startswith('MCCrossClosure'): allHistos [ 'unfoldHistoCross'+ivar ] = tunfolder_cross.GetOutput("unfoldHistoCross"+ivar).Clone()
        unfoldingtot = allHistos [ 'unfoldHisto'+ivar ].Integral()

        allHistos [ 'foldHisto'+ivar ] = tunfolder.GetFoldedOutput("folded"+ivar).Clone()
        if process.startswith('data'): 
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ].Clone(), 'data-Bkgs',  allHistos [ 'foldHisto'+ivar ].Clone(), 'folded', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", rebinX=1, version=sel+'_'+version, outputDir=outputDir )
            plotSimpleComparison( allHistos[ 'unfoldHisto'+ ivar ].Clone(), 'unfold',  signalHistos[signalLabel+'_accepgen'+ivar+sel].Clone(), 'accepgen', ivar+'_from'+('Data' if process.startswith('data') else 'MC')+'_'+signalLabel+"_Test", rebinX=1, version=sel+'_'+version, outputDir=outputDir )

        ############################ Get Probability matrix ########################################
        allHistos[ 'probaMatrix'+ivar ] = tunfolder.GetProbabilityMatrix('probaMatrix'+ivar).Clone()
        
        
        ########################### Get systematic shifts of output ################################
        uncerUnfoldHisto = OrderedDict() #storing individual systematics totals (ie, up/down or other variations), background subtraction systematics, and then the overall systematic unc (individual systs + bkgs)
        uncerUnfoldSystCov = OrderedDict()
        if len(sysUncert)>0 and process.startswith('data'):
            print ('|------> TUnfolding uncertainties:')

            
            
            for sys in sysUncert:
                
                sys_cov_up = tunfolder.GetEmatrixSysUncorr("cov_%s_up"%sys+ivar)
                sys_cov_down = tunfolder.GetEmatrixSysUncorr("cov_%s_down"%sys+ivar)
                '''
                if sys in jes_uncorr_list and year.startswith('all') and sys not in jes_corr_list: 
                    
                    print('|------> TUnfolding extracting %sUnc info to root file but from per year unfoldings'%sys)
                    s = [i for i in sysSignalLabels if sys in i]
                    #print (s)
                        
                    if ('2017' in sys or '2018' in sys) and len(s)>1:
                        print ("this list should have only one item, 2.0")
                    else:
                        s=[s[0]]
                    
                    dictUncHistos[sys+'Up'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Up'+sel].Clone()
                    dictUncHistos[sys+'Down'] = sysSignalHistos[s[0]+'_reco'+ivar+sys+'Down'+sel].Clone()
                        
                    for upDown in [ 'Up', 'Down' ]:
                        #print (sys+upDown)
                        
                        uncerUnfoldHisto[ivar+sys+upDown] = dictUncHistos[sys+upDown]#tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "-1#sigma")
                        try: uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                        except ReferenceError: uncerUnfoldHisto.pop( ivar+sys+upDown, None )
                        

                    # Create total uncertainty and sys uncertainty plots.
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone()       # Syst uncertainty
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                    for i in range( 0, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                        try: yup = abs(  uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs( uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'].GetBinContent(i,i) )
                        except KeyError: ydn = 0
                                                        
                        dy = 0.5*(yup + ydn) 
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )
                    print ("Done with adding %s to tot syst unc"%sys )    
                    
                
                if sys in jes_corr_list and year.startswith('all') and sys not in jes_uncorr_list:
                    
                    tunfolder.GetEmatrixSysSource(sys_cov_up, "systcov_"+ivar+sys+'Up'+"shift")
                    tunfolder.GetEmatrixSysSource(sys_cov_down, "systcov_"+ivar+sys+'Down'+"shift")
                    
                    #print (tunfolder.GetEmatrixSysSource(sys_cov_up, "systcov_"+ivar+sys+'Up'+"shift"))
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Up'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'] = sys_cov_down.Clone()

                    
                    for upDown in [ 'Up', 'Down' ]:
                        #print (sys+upDown)
                        uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "-1#sigma")
                        try: uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                        except ReferenceError: uncerUnfoldHisto.pop( ivar+sys+upDown, None )
                        

                    # Create total uncertainty and sys uncertainty plots.
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone()       # Syst uncertainty
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                    for i in range( 0, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                        try: yup = abs(  uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs( uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'].GetBinContent(i,i) )
                        except KeyError: ydn = 0
                                                        
                        dy = 0.5*(yup + ydn) 
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )
                    print ("Done with adding %s to tot syst unc"%sys )    
                '''
                if not sys.startswith(('_model', '_CR', '_erdON', '_mtop', '_hdamp', '_TuneCP5')):
                    #if 'jes' in sys and year.startswith('all'): continue
                    
                    tunfolder.GetEmatrixSysSource(sys_cov_up, sys+'Up')
                    tunfolder.GetEmatrixSysSource(sys_cov_down, sys+'Down')
                    
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Up'] = sys_cov_up.Clone()
                    uncerUnfoldSystCov['systcov_'+ivar+sys+'Down'] = sys_cov_down.Clone()
                    for upDown in [ 'Up', 'Down' ]:
                        #print (sys+upDown)
                        uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "-1#sigma")
                        #print (uncerUnfoldHisto[ivar+sys+upDown])
                        if uncerUnfoldHisto[ivar+sys+upDown]: uncerUnfoldHisto[ivar+sys+upDown+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+sys+upDown].Clone(),
                                                                                               allHistos [ 'unfoldHisto'+ivar ].Clone())
                        #uncerUnfoldSystCov['systcov_'+ivar+sys+upDown] = tunfolder.GetEmatrixSysSource(sys+upDown, "systcov_"+ivar+sys+upDown+"shift")
                        try: uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                        except ReferenceError: uncerUnfoldHisto.pop( ivar+sys+upDown, None )
                        

                    # Create total uncertainty and sys uncertainty plots.
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone()       # Syst uncertainty
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                    for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                        try: yup = abs(  uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                        except KeyError: yup = 0
                        try: ydn = abs(  uncerUnfoldHisto[ivar+sys+'Down'].GetBinContent(i) )
                        except KeyError: ydn = 0
                                                        
                        dy = 0.5*( (yup + ydn) )
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )
                    uncerUnfoldHisto[ivar+sys.upper()+'Total'+"_shiftHist"] = get_syst_shifted_hist( uncerUnfoldHisto[ivar+sys.upper()+'Total'].Clone(),
                                                                                                    allHistos [ 'unfoldHisto'+ivar ].Clone())
                    #print ("Done with adding %s to tot syst unc"%sys )    

                elif sys.startswith('_model'):
                    #if not process.startswith('MCSelfClosure'):
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_modelUncTotal')
                    uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Reset()
                    tmpModelHisto = tunfolder.GetDeltaSysSource('modelUncTotal', "unfoldHisto_"+ivar+"_modelUncTotalshift", "-1#sigma")
                    #uncerUnfoldHisto[ivar+'_Physics Model'+"_shiftHist"] = get_syst_shifted_hist(tmpModelHisto.Clone().Clone(),
                    #                                                                                allHistos [ 'unfoldHisto'+ivar ].Clone())
                    tunfolder.GetEmatrixSysSource(sys_cov_up, "modelUncTotal")
                    uncerUnfoldSystCov['systcov_'+ivar+'_Physics ModelTotal'] = sys_cov_up.Clone()
                    for i in range( 1, tmpModelHisto.GetNbinsX() + 1):
                        uncerUnfoldHisto[ivar+'_Physics ModelTotal'].SetBinContent( i, abs(tmpModelHisto.GetBinContent(i)))
                    uncerUnfoldHisto[ivar+'_Physics Model'+"_shiftHist"] = get_syst_shifted_hist(uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Clone(),allHistos [ 'unfoldHisto'+ivar ].Clone())
                    #print('_modelUncTotal')
                    
                
        ################################################################################
        
        #print (uncerUnfoldSystCov)
        
        ##########################  Get various covariances  ###########################
        '''
        GetEmatrixSysUncorr(): uncorrelated errors on the input matrix histA, taken as the errors provided with the histogram. These are typically statistical errors from finite Monte Carlo samples.
        GetEmatrixSysSource()(GetDeltaSysSource()): correlated shifts of the input matrix histA. These shifts are taken as one-sigma effects when switchig on a given error soure. Several such error sources may be defined
        GetEmatrixSysBackgroundUncorr(): uncorrelated errors on background sources, originating from the errors provided with the background histograms
        GetEmatrixInput(): statistical uncertainty of the input (the measurement)
        
        '''
        #print ("Here's a list of all available uncerUnfoldHistos:", uncerUnfoldHisto)
        print ('|------> TUnfolding covariances')
        allHistos[ 'cov'+ivar ] = tunfolder.GetEmatrixTotal("cov"+ivar, "Total Covariance Matrix")
        allHistos[ 'cov_uncorr_data_'+ivar ] = tunfolder.GetEmatrixInput("cov_uncorr_data"+ivar, "Covariance Matrix coming from input data distributions")
        allHistos[ 'cov_uncorr_'+ivar ] = tunfolder.GetEmatrixSysUncorr("cov_uncorr"+ivar, "Covariance Matrix at hadron level coming from RM")
        allHistos[ 'cov_uncorr_bkg_'+ivar ] = tunfolder.GetEmatrixSysBackgroundUncorr('fakes', "Covariance Matrix from uncorrelated errors of Background sources")
        
        if process.startswith('data'):
            if sel.startswith(('_W','_top')):
                for ibkg in bkgSources:
                    if 'fakes' not in ibkg: allHistos[ 'cov_uncorr_bkg_'+ivar ].Add(tunfolder.GetEmatrixSysBackgroundUncorr(ibkg, "Covariance Matrix from uncorrelated errors of Background sources"))
        
        #### cov = cov_uncorr + cov_uncorr_data + other uncertaitnies
        #### stat = cov_uncorr + cov_uncorr_data
        
        

        
        ############### Build correlation matrix for unfolding#################################
        allHistos['correlation_matrix_'+ivar] = allHistos[ 'cov'+ivar ].Clone()
        allHistos['correlation_matrix_'+ivar].Reset()
        allHistos['correlation_matrix_'+ivar] = correlation_from_covariance(allHistos[ 'cov'+ivar ].Clone(),allHistos['correlation_matrix_'+ivar])
        
        
        '''
        for ibinx in range( 1, allHistos[ 'cov'+ivar ].GetNbinsX()):
            variance = allHistos[ 'cov'+ivar ].GetBinContent(ibinx,ibinx)

            for ibiny in range( 1, allHistos[ 'cov'+ivar ].GetNbinsY()):
                cov_xy = allHistos[ 'cov'+ivar ].GetBinContent(ibinx,ibiny)
                if variance!=0: 
                    allHistos['correlation_matrix_'+ivar].SetBinContent(ibinx,ibiny,cov_xy/variance)
        '''
        ########################################################################################    
        
        # Create total uncertainty and sys uncertainty histos for plots 
        
        # first build up histos of systematic uncertainties from the BACKGROUND SUBTRACTION 
        uncerUnfoldHisto[ivar+'_BkgTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_BkgTotal')
        uncerUnfoldHisto[ivar+'_BkgTotal'].Reset()
        for ibin in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):
            bkg_tot = np.sqrt(allHistos[ 'cov_uncorr_bkg_'+ivar ].GetBinContent(ibin,ibin))
            uncerUnfoldHisto[ivar+'_BkgTotal'].SetBinContent(ibin, bkg_tot)
            
        # second build up histos of total systematic uncertainties
        uncerUnfoldHisto[ivar+'_SystTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_SystTotal')
        uncerUnfoldHisto[ivar+'_SystTotal'].Reset()        
        
        tmp = OrderedDict()
        dummyM = hist2array(allHistos[ 'cov'+ivar ].Clone())
        print (dummyM.shape)
        systcovTotal = np.zeros(dummyM.shape)
        #print (systcovTotal.shape)
        for k in uncerUnfoldHisto:
            corrflag = False
            sysErr=np.zeros(uncerUnfoldHisto[k].GetNbinsX())
            if not k.endswith('Total') or k.endswith(('SystTotal')): continue

            for x in jes_corr_list:
                if corrflag==True: break
                if x.upper() in k:
                    corrflag=True
                    continue
                    
            #print ("CORRFLAG",corrflag)
            for i in range( 0, uncerUnfoldHisto[k].GetNbinsX()):


                sysErr[i] = uncerUnfoldHisto[k].GetBinContent(i+1)
                if corrflag: 
                    systcovTotal+=np.outer(sysErr,sysErr)
                    #print ("CORRFLAG",corrflag)
                else: systcovTotal+=np.diag(sysErr*sysErr)

        #print (systcovTotal.shape,corrflag)

        allHistos['cov_systTotal'+ivar] = allHistos[ 'cov'+ivar ].Clone()
        allHistos['cov_systTotal'+ivar].Reset()
        allHistos['cov_systTotal'+ivar] = array2hist(systcovTotal,allHistos['cov_systTotal'+ivar])
        
        tmp = OrderedDict()
        for i in range( 1, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
            tmp[i] = 0
            for k in uncerUnfoldHisto:
                if k.endswith('Total') and not k.endswith(('SystTotal')):
                    tmp[i] = tmp[i] + ( uncerUnfoldHisto[k].GetBinContent( i )**2 )
                    #print(i, k, tmp[i], uncerUnfoldHisto[k].GetBinContent( i ), ( uncerUnfoldHisto[k].GetBinContent( i )**2 ))
        
        # adding covariances from bkg subtraction and RM finite stats to the overall systematics covariance matrix
        
                
        allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_'+ivar ])
        allHistos['cov_systTotal'+ivar].Add(allHistos[ 'cov_uncorr_bkg_'+ivar ])
        
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
            #if len(sysUncert)>0: syst_tot = uncerUnfoldHisto[ivar+'_SystTotal'].GetBinContent(ibin)
            #else: syst_tot=0
            orig_statunfunc = allHistos[ 'unfoldHisto'+ivar ].GetBinError(ibin) #obtain the stat + unf uncertatinties directly from output (which by default contains no systematics)
            
            norm_unc = abs(allHistos[ 'unfoldHisto'+ivar ].GetBinContent(ibin))
            if unc_tot<=0.: 
                unc_tot=0.
            if datastat_tot<=0.: 
                datastat_tot=0.
            if rmstat_tot<=0.: 
                rmstat_tot=0.
               
            
            allHistos[ 'unfoldHisto'+ivar ].SetBinError(ibin, unc_tot)
            allHistos[ 'unfoldHistoBkgUnc'+ivar ].SetBinError(ibin, bkg_tot )
            allHistos[ 'unfoldHistowoUnc'+ivar ].SetBinError(ibin, 0. )        # No unc
        
            uncerUnfoldHisto[ivar+'_TotalUnc'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMErrTotal'].SetBinContent(ibin, unc_tot )
            uncerUnfoldHisto[ivar+'_CMMCStatErrTotal'].SetBinContent(ibin, rmstat_tot)
            uncerUnfoldHisto[ivar+'_CMDataStatErrTotal'].SetBinContent(ibin, datastat_tot)
            #uncerUnfoldHisto[ivar+'_StatTotal'].SetBinContent(ibin, orig_statunfunc)
            
            ratioHistos[ 'TotalUnc'+ivar ].SetBinContent( ibin, 1. )
            ratioHistos[ 'StatUnc'+ivar ].SetBinContent( ibin, 1. )
            #ratioHistos[ 'SystUnc'+ivar ].SetBinContent( ibin, 1. )

            if norm_unc!=0:
                ratioHistos[ 'TotalUnc'+ivar ].SetBinError( ibin, unc_tot)#+syst_tot**2
                ratioHistos[ 'StatUnc'+ivar ].SetBinError( ibin, datastat_tot)
                #if len(sysUncert)>0: ratioHistos[ 'SystUnc'+ivar ].SetBinError( ibin, syst_tot/norm_unc)
            #print (norm_unc, np.sqrt(unc_tot**2 ), unc_tot, orig_statunfunc)

        allHistos [ 'foldHisto'+ivar ] = get_folded_unfolded(folded=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                                                            unfolded=allHistos['unfoldHisto'+ivar].Clone(), 
                                                            cov_tot=allHistos['cov'+ivar].Clone(), 
                                                             probaM=allHistos[ 'probaMatrix'+ivar ]
                                                            
                                                             )
        
        
        #####################################
        #print (allHistos[ 'unfoldHisto'+ivar ].Clone())
        ###### Plot unfolding results

        print ('|------> Drawing unfold plot:')
        if not 'Closure' in process:
            drawUnfold(ivar=ivar, 
                       selection=sel, year=year,
                       process=process,
                       dataJetHisto=allHistos[ 'dataHistoGenBin' ].Clone(),
                       genJetHisto=signalHistos[ signalLabel+'_accepgen'+ivar+sel ].Clone(),
                       unfoldHisto=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                       unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                       foldHisto=tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                       recoJetHisto=signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone(),
                       cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(),#ratiototUncHisto=ratioHistos[ 'TotalUnc'+ivar ].Clone(),
                       cov_tot=allHistos['cov'+ivar].Clone(),#ratioUncHisto=ratioHistos[ 'StatUnc'+ivar ].Clone(),
                       #ratiosystUncHisto =ratioHistos[ 'SystUnc'+ivar ].Clone(),
                       altMCHisto =  altSignalHistos[altSignalLabel+'_accepgen'+ivar+sel].Clone(),
                       labelX=variables[ivar]['label'],
                       maxX=variables[ivar]['bins'][-1],
                       tlegendAlignment=variables[ivar]['alignLeg'],
                       
                       outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+signalLabel+'_TUnfold_'+version+'.'+ext
                       )
        else: 
            if 'Cross' in process:
                drawClosures(ivar=ivar, selection=sel, year=year, process=process,
                             
                             genJetHistoCross=altSignalHistos[altSignalLabel+'_accepgen'+ivar+'_nom'+sel].Clone(), 
                             unfoldHistoCross=allHistos[ 'unfoldHistoCross'+ivar ].Clone() ,
                             
                             genJetHisto=signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Clone(),
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
                drawClosures(ivar=ivar, selection=sel, year=year, process=process,
                             genJetHistoCross=[], 
                             unfoldHistoCross=[] ,
                             genJetHisto=signalHistos[ signalLabel+'_accepgen'+ivar+'_nom'+sel ].Clone(),
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
        
        if not 'Closure' in process: drawUncertainties_normalizedshifts(ivar=ivar,unfoldHistoTotUnc=allHistos[ 'unfoldHisto'+ivar ].Clone(),
                                                                        unfoldHistowoUnc=allHistos[ 'unfoldHistowoUnc'+ivar ],
                                                                        uncerUnfoldHisto=uncerUnfoldHisto, 
                                                                        cov_tot=allHistos['cov'+ivar].Clone(), 
                                                                        cov_datastat_tot=allHistos['cov_uncorr_data_'+ivar].Clone(), 
                                                                        cov_rmstat_tot=allHistos['cov_uncorr_'+ivar].Clone(), 
                                                                        cov_bkg_tot=allHistos['cov_uncorr_bkg_'+ivar].Clone(),
                                                                        labelX=variables[ivar]['label'], 
                                                                        tlegendAlignment=variables[ivar]['alignLeg'],
                                                                        outputName=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_UNC_'+version+'.'+ext,
                                                                        year=year, unftot=unfoldingtot
                                                                        )
        
        ######### Plotting 2D matrices of various kinds
        print ('|------> Drawing various 2D matrices:')
        draw2D( ivar,  tunfolder.GetRhoItotal("rhoI"+ivar, "Global correlations"), variables[ivar], outputLabel='data_rhoI', outputDir=outputDir,selection=sel,version=version,year=year)
        draw2D( ivar,  allHistos[ 'correlation_matrix_'+ivar ].Clone(), variables[ivar], outputLabel='data_correlationMatrix', outputDir=outputDir,selection=sel,version=version,year=year)
        draw2D( ivar,  allHistos[ 'probaMatrix'+ivar ].Clone(), variables[ivar], outputLabel='data_probaMatrix', outputDir=outputDir, addCorrelation=True, addCondition=True ,selection=sel,version=version,year=year)
        draw2D( ivar, signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].Clone(), variables[ivar], outputLabel='data_respMatrix', outputDir=outputDir, addCorrelation=True ,selection=sel,version=version,year=year)
        draw2D( ivar, allHistos[ 'cov'+ivar].Clone(), variables[ivar], outputLabel='dataTotal_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        draw2D( ivar, allHistos[ 'cov_uncorr_'+ivar].Clone(), variables[ivar], outputLabel='uncorrUncRM_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        draw2D( ivar, allHistos[ 'cov_uncorr_data_'+ivar].Clone(), variables[ivar], outputLabel='dataInpStats_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        draw2D( ivar, allHistos[ 'cov_uncorr_bkg_'+ivar].Clone(), variables[ivar], outputLabel='BkgSubtractionSyst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)
        draw2D( ivar, allHistos[ 'cov_systTotal'+ivar].Clone(), variables[ivar], outputLabel='Syst_covMatrix', outputDir=outputDir, addCorrelation=True,selection=sel,version=version,year=year)

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
        
        print (allHistos,signalHistos,dataHistos)
        ######### Saving Histograms
        def renamingHistos( dictHistos ):
            for isam, hist in dictHistos.items():
                ihis = hist.Clone()
                ihis.SetName(isam)
                ihis.SetTitle(isam)
                ihis.Write()

        outputRootName = outputDir+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root'
        #print ('|------> Saving histograms in rootfile: ', outputRootName)
        outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
        renamingHistos( signalHistos )
        if not process.startswith('MCSelfClosure'):
            renamingHistos( sysSignalHistos )
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
#####################################################################################################################################






############# Plotting functions #################

def plotSysComparison( nomHisto, dictUncHistos, outputName, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ):

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 219, 92, 48, 49, 37, 38, 33, 17, 50, 205, 225, 94, 221, 16,  225, 128]
    
    
    outputFileName = outputName+'_'+version+'.'+ext
    print ('Processing.......', outputFileName)

    binWidth = nomHisto.GetBinWidth(1)

    legend=ROOT.TLegend(0.35,0.6,0.80,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.02)
    legend.SetBorderSize(0)

    multiGraph = ROOT.TMultiGraph()
    gnom = ROOT.TGraphAsymmErrors()
    gnom.Divide( nomHisto.Clone(), nomHisto.Clone(), 'pois' )
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetMarkerStyle(1)
    gnom.SetLineWidth(2)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    multiGraph.Add( gnom )

    dictgraph = {}
    dummy=0
    
    #print (dictUncHistos)
    for ih in dictUncHistos:
        #print(ih, dummy, len(colors))
        dictgraph[ih] = ROOT.TGraphAsymmErrors()
        dictgraph[ih].Divide( dictUncHistos[ih], nomHisto, 'pois' )
        if not dummy==len(colors)-1:
            dictgraph[ih].SetLineColor( colors[dummy] )
            dictgraph[ih].SetLineStyle( 2 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        else:
            dummy = dummy-len(colors)+2
            dictgraph[ih].SetLineColor( colors[dummy] )
            dictgraph[ih].SetLineStyle( 3 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        if 'jes' in ih and ('2017' in ih or '2018' in ih): legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        else: legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        multiGraph.Add( dictgraph[ih] )
        dummy=dummy+1

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 750, 500 )
    if log: canUnc.SetLogy()
    multiGraph.GetYaxis().SetTitle( 'Ratio Unc/Nominal' )
    multiGraph.GetXaxis().SetTitle( labelX )
    multiGraph.SetMaximum( 3. )
    multiGraph.SetMinimum( -1.)
    multiGraph.Draw('ALP')

    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+year
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    if ext.startswith('pdf'):
        canUnc.SaveAs( outputDir + outputFileName.replace('pdf', 'png') )
    del canUnc
    
def drawClosures(ivar, selection, process, year, genJetHisto, genJetHistoCross, unfoldHisto, unfoldHistoCross,
                 ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, labelX, maxX, tlegendAlignment, outputName ):
    
    """docstring for drawUnfold"""
    print ("Drawing unfolding closure")
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.65,0.90,0.88)

    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)

    dataIntegral = genJetHisto.Integral()
    unfoldIntegral = unfoldHisto.Integral()
    if process.startswith('MCCrossClosure'): unfoldcrossIntegral = unfoldHistoCross.Integral()
    #tmpMax = float(unfoldHisto.GetMaximum())
    
    #unfoldHisto.Scale(1, 'width')  ### divide by bin width
    unfoldHisto.Scale(1/unfoldHisto.Integral(), 'width')  ### divide by bin width
    unfoldHisto.SetMarkerStyle(4)
    #unfoldHisto.SetMarkerSize(2)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    unfoldHisto.SetLineWidth(2)
    legend.AddEntry( unfoldHisto, ('Data (self-closure)' if process.startswith('MCSelfClosure') else 'MG5+P8 unf. w/ MG5+P8'), 'pe' )
    
    #genJetHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    genJetHisto.Scale(1/genJetHisto.Integral(), 'width')  ### divide by bin width
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(ROOT.kBlue)
    genJetHisto.SetMarkerStyle(0)
    genJetHisto.SetLineStyle(2)
    legend.AddEntry( genJetHisto, 'MG5+Pythia8 (gen)', 'lp' )
    unfoldHisto.GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( 1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )

    unfoldHisto.Draw( "E")
    genJetHisto.Draw( "histe same")
    
    if not process.startswith('MCSelfClosure'):
        #unfoldHisto.Scale(1, 'width')  ### divide by bin width
        unfoldHistoCross.Scale(1/unfoldHistoCross.Integral(), 'width')  ### divide by bin width
        unfoldHistoCross.SetMarkerStyle(26)
        #unfoldHistoCross.SetMarkerSize(2)
        unfoldHistoCross.SetMarkerColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineWidth(2)
        legend.AddEntry( unfoldHistoCross, 'MG5+P8 unf. w/ Herwig7', 'pe' )

        #genJetHisto.Scale(1, 'width')
        #genJetHisto.Scale(scaleFactor)
        genJetHistoCross.Scale(1/genJetHistoCross.Integral(), 'width')  ### divide by bin width
        genJetHistoCross.SetLineWidth(2)
        genJetHistoCross.SetLineColor(ROOT.kMagenta)
        genJetHistoCross.SetMarkerStyle(0)
        genJetHistoCross.SetLineStyle(2)
        legend.AddEntry( genJetHistoCross, 'Herwig7 (gen)', 'lp' )
        
        unfoldHistoCross.Draw( "E same")
        genJetHistoCross.Draw( "histe same")
        

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.045)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = 'Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
    elif selection.startswith("_W"): seltext = ' Boosted W region'
    elif selection.startswith("_top"): seltext = ' Boosted top region'
    selText.DrawLatex( ( 0.2 if tlegendAlignment.startswith('right') else 0.68 ), 0.87, seltext )

    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if selection.startswith('dijet') else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)

    pad2.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    tmpPad2.GetXaxis().SetTitle( labelX )
    tmpPad2.GetYaxis().SetTitle( "Sim./Data" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(0.4,1.6 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    
    if 'Self' in process:

        hRatioUp = ROOT.TGraphAsymmErrors()
        hRatioUp.Divide( genJetHisto, unfoldHisto, 'pois' )
        hRatioUp.SetLineColor(ROOT.kBlack)
        hRatioUp.SetMarkerColor(ROOT.kBlack)
        hRatioUp.SetLineWidth(2)
        hRatioUp.SetMarkerStyle(25)
        hRatioUp.Draw('P0')

    else:
        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfoldHistoCross, unfoldHisto, 'pois' )
        hRatioUp2.SetLineColor(ROOT.kBlack)
        hRatioUp2.SetMarkerColor(ROOT.kBlack)
        hRatioUp2.SetLineWidth(2)
        hRatioUp2.SetMarkerStyle(25)
        hRatioUp2.Draw('P0')
    
    
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)
    
def draw2D( ivar, histo, varInfo, outputDir, outputLabel='data', addCorrelation=False, addCondition=False, addInvertedMatrix=False,ext='pdf',selection='_dijetSel',version='vNew',year='2017' ):

    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_'+version+'.'+ext


    ROOT.gStyle.SetPadRightMargin(0.15)
    can2D = ROOT.TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
    histo.GetXaxis().SetTitle('Accepted Gen '+varInfo['label'])
    histo.GetYaxis().SetTitle('True Reco '+varInfo['label'])
    histo.GetYaxis().SetTitleOffset( 0.8 )
    #signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel]covHisto.GetYaxis().SetRange( variables[ivar]['bins'][0], variables[ivar]['bins'][-1] )
    histo.Draw("colz")
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(can2D, 4, 0)

    textBox.SetTextSize(0.04)
    if addCorrelation:
        #print('|-----> Correlation: ', histo.GetCorrelationFactor())
        textBoxCorr = textBox.Clone()
        textBoxCorr.DrawLatex( 0.05, varInfo['bins'][-1]-( .05*(varInfo['bins'][-1]-varInfo['bins'][0]) ), '#color[8]{Corr. Factor = '+str(round(histo.GetCorrelationFactor(),2))+'}' )

    if addCondition:   
        ## based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L41
        Nx = histo.GetNbinsX()
        Ny = histo.GetNbinsY()
        RMx = histo.ProjectionX( 'RMx', 0, -1 )

        m = ROOT.TMatrixD( Ny, Nx )   ### need to swap the axes
        for ibin in range(1, Nx+1):
            normalization = RMx.GetBinContent(ibin)
            if (normalization>0):
                for jbin in range( 1, Ny+1 ):
                    m[jbin-1][ibin-1] = histo.GetBinContent(ibin,jbin) / normalization
        svd = ROOT.TDecompSVD(m)
        v = ROOT.TVectorD( svd.GetSig() )
        Min = v[0]
        Max = v[0]
        for ibin in range( 0, Nx ):
            if (abs(v[ibin]) < 1e-5 ): break
            Min = v[ibin]
        conditionNumber = str(round( Max/Min, 2 )) if Min > 0 else 'NaN'
        print('|-----> Condition Number: ', conditionNumber)
        textBoxCond = textBox.Clone()
        textBoxCond.DrawLatex( 0.05, varInfo['bins'][-1]-( .1*(varInfo['bins'][-1]-varInfo['bins'][0]) ), '#color[8]{Cond. Number = '+conditionNumber+'}' )

    can2D.SaveAs(outputName)
    if ext.startswith('pdf'):
        can2D.SaveAs( outputName.replace('pdf', 'png') )

    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)


def getFilesInDictSamples_fromROOT(labelBegin='', year_list=['2017','2018','all']):
    Files={}
    for iy in year_list:
        for isam in dictSamples:
            if not checkDict( isam, dictSamples )[iy]['skimmerHisto'].endswith('root'): continue
            if isam.startswith(labelBegin):
                Files[isam] = [
                                ROOT.TFile.Open(inputFolder+checkDict( isam, dictSamples )[iy]['skimmerHisto'] ),
                                checkDict( isam, dictSamples )
                            ]
    return Files

def combinePlots( name, dictHistos, numBins, mainHistoLabel, otherHisto, otherHistoLabel, outputLabel, variables, ext, process, log, year, runMLU, version,  axisX='', outputDir='Plots/', ratioOnly=False):

    """docstring for combinePlots"""

    outputFileName = name+'_'+outputLabel+'_combinePlots_'+version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.10,0.80,0.60,0.90)
    legend.SetNColumns(3 if runMLU else 2)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.06)

    dictHistos['combData'] = ROOT.TH1F('combData', 'combData', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combData' ], mainHistoLabel, 'lep' )
    dictHistos['combUnfold'] = ROOT.TH1F('combUnfold', 'combUnfold', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combUnfold' ], otherHistoLabel, 'lep' )
    if runMLU:
        dictHistos['combMLU'] = dictHistos['combUnfold'].Clone()
        legend.AddEntry( dictHistos[ 'combMLU' ], 'MLU', 'lep' )

    tmpNbin = 0
    Xlabels = []
    for ivar,ih in dictHistos.items():
        if ivar.startswith('Jet') and not ivar.endswith(('21', '32')):
            Xlabels.append( '#'+variables[ivar]['label'].split('#')[1] )
            for ibin in range(1, ih['data'].GetNbinsX()+1):
                tmpNbin = tmpNbin+1
                dictHistos['combData'].SetBinContent( tmpNbin, ih[otherHisto].GetBinContent(ibin) )
                dictHistos['combData'].SetBinError( tmpNbin, ih[otherHisto].GetBinError(ibin) )
                dictHistos['combUnfold'].SetBinContent( tmpNbin, ih['unfold'].GetBinContent(ibin) )
                dictHistos['combUnfold'].SetBinError( tmpNbin, ih['unfold'].GetBinError(ibin) )
                if runMLU:
                    dictHistos['combMLU'].SetBinContent( tmpNbin, ih['MLU'].GetBinContent(ibin) )
                    dictHistos['combMLU'].SetBinError( tmpNbin, ih['MLU'].GetBinError(ibin) )


    if not ratioOnly:
        canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 750 )
        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetPadLeftMargin(0.08)
        ROOT.gStyle.SetPadTickX(0)
        pad1 = ROOT.TPad(ivar+'1', "Fit",0.,0.330,1.00,1.00,-1)
        pad2 = ROOT.TPad(ivar+'2', "Pull",0,0.00,1.00,0.40,-1);
        pad1.Draw()
        pad2.Draw()

        pad1.cd()
        #pad1.SetLogy()
        dictHistos['combData'].GetXaxis().SetNdivisions(100)
        dictHistos['combData'].GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#tau_{X}}' )
        dictHistos['combData'].GetYaxis().SetTitleSize( 0.06 )
        dictHistos['combData'].GetYaxis().SetTitleOffset( 0.6 )
        dictHistos['combData'].SetMaximum( dictHistos['combUnfold'].GetMaximum()*1.2 )
        dictHistos['combData'].SetMinimum( 0.001 )
        dictHistos['combData'].SetLineColor( ROOT.kBlack )
        dictHistos['combData'].SetLineWidth( 2 )

        dictHistos['combUnfold'].SetLineColor( ROOT.kMagenta )
        dictHistos['combUnfold'].SetLineWidth( 2 )
        dictHistos['combData'].Draw('E')
        dictHistos[ 'combUnfold' ].Draw('E same')
        if runMLU:
            dictHistos['combMLU'].SetLineColor( 8 )
            dictHistos['combMLU'].SetLineWidth( 2 )
            dictHistos[ 'combMLU' ].Draw('E same')

        ### division lines
        lines = {}
        for i in numBins[1:-1]:
            lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
            lines[i].SetLineColor(ROOT.kGray)
            lines[i].Draw('same')

        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.relPosX = 0.07
        CMS_lumi.CMS_lumi( pad1, 4, 0)
        legend.Draw()

        pad2.cd()
        pad2.SetGridy()
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.2)

        tmppad= pad2.DrawFrame(0,0.2,numBins[-1],1.8)
        tmppad.GetYaxis().SetTitle( "Data/Unfold" )
        tmppad.GetXaxis().SetTitle(  axisX )
        tmppad.GetYaxis().SetTitleOffset( 0.4 )
        tmppad.GetYaxis().CenterTitle()
        tmppad.SetLabelSize(0., 'x')
        tmppad.SetTitleSize(0., 'x')
        tmppad.SetLabelSize(0.10, 'y')
        tmppad.SetTitleSize(0.10, 'y')
        tmppad.SetNdivisions(100, 'x')
        tmppad.SetNdivisions(505, 'y')
        pad2.Modified()
        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dictHistos['combUnfold'], dictHistos[ 'combData' ], 'pois' )
        hRatio.SetMarkerStyle(8)
        hRatio.Draw('P0')
        if runMLU:
            hRatioMLU = ROOT.TGraphAsymmErrors()
            hRatioMLU.Divide( dictHistos['combMLU'], dictHistos[ 'combData' ], 'pois' )
            hRatioMLU.SetMarkerStyle(8)
            hRatioMLU.SetMarkerColor( 8 )
            hRatioMLU.Draw('P0 same')
            hRatio.Draw('P0 same')

        for i in lines: lines[i].Draw('same')

        textBox.SetTextSize(0.10)
        textBox.SetTextAlign(12)
        textBoxList = {}
        
        for i in range(1, len(numBins)):
            textBoxList[i] = textBox.Clone()
            textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., 0., Xlabels[i-1] )
        textBox.SetTextSize(0.04)
    else:

        legend=ROOT.TLegend(0.60,0.80,0.90,0.90)
        legend.SetNColumns(3 if runMLU else 2)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.06)

        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetPadLeftMargin(0.08)
        ROOT.gStyle.SetPadTickX(0)
        canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 500 )
        canvas[outputFileName].SetGridy()

        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dictHistos[ 'combData' ], dictHistos['combUnfold'], 'pois' )
        if process.startswith('MCCrossClosure'): labelLegend = 'MCSelfClosure Ind.Sample'
        elif process.startswith('MCClosure'): labelLegend = 'MCCrossClosure'
        else: labelLegend = process
        legend.AddEntry( hRatio, labelLegend, 'lep' )
        hRatio.SetMarkerStyle(8)

        hRatio.GetXaxis().SetNdivisions(100)
        hRatio.GetYaxis().SetTitle( 'Sim. / Data'+(' (MC)' if process.startswith('MC') else '') )
        hRatio.GetYaxis().SetTitleSize( 0.06 )
        hRatio.GetYaxis().SetTitleOffset( 0.6 )
        hRatio.GetXaxis().SetLimits( 0., numBins[-1] )
        hRatio.SetMaximum( 2. )
        hRatio.SetMinimum( 0. )

        hRatio.Draw('AP')

        ### division lines
        lines = {}
        for i in numBins[1:-1]:
            lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
            lines[i].SetLineColor(ROOT.kGray)
            lines[i].Draw('same')

        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.relPosX = 0.07
        CMS_lumi.CMS_lumi( canvas[outputFileName], 4, 0)
        legend.Draw()

        aBox = textBox.Clone()
        aBox.SetTextSize(0.06)
        #aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Central Jet' if name.startswith('recoJet2') else 'Outer Jet')+"}{Dijet Selection}}"  )
        if 'WSel' in selection: aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Leading jet')+"}{Boosted W selection}}"  )
        elif 'topSel'in selection: aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Leading jet')+"}{Boosted top selection}}"  )
        textBox.SetTextAlign(12)
        textBoxList = {}
        for i in range(1, len(numBins)):
            textBoxList[i] = textBox.Clone()
            textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., -0.1, Xlabels[i-1] )

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]
    ROOT.gStyle.SetPadTickX(1)

def combineRatioPlots( name, ratioDicts, numBins, outputLabel,process, ext, log, version, selection, axisX='', outputDir='Plots/'):
    """docstring for combineRatioPlots"""

    outputFileName = name+'_'+outputLabel+'_combineBLTPlots_'+version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.5,0.8,0.8,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.055)
    legend.SetBorderSize(0)
    
    tmpNbin = 0
    Xlabels = []
    yvaluesGen = np.array( [0]*numBins[-1], 'd' )
    yErrLowValuesGen = np.array( [0]*numBins[-1], 'd' )
    yErrHighValuesGen = np.array( [0]*numBins[-1], 'd' )
    yvaluesReco = np.array( [0]*numBins[-1], 'd' )
    yErrLowValuesReco = np.array( [0]*numBins[-1], 'd' )
    yErrHighValuesReco = np.array( [0]*numBins[-1], 'd' )
    for ivar,ih in ratioDicts.items():
        if ivar.startswith('Jet') and not ivar.endswith(('21', '32')):
            Xlabels.append( '#'+variables[ivar]['label'].split('#')[1] )
            for ibin in range( ih[0].GetN() ):
                a = ctypes.c_double(0.)
                b = ctypes.c_double(0.)
                ih[0].GetPoint(ibin, a, b)
                yvaluesGen[tmpNbin] = b.value
                yErrLowValuesGen[tmpNbin] = ih[0].GetErrorYlow(ibin)
                yErrHighValuesGen[tmpNbin] = ih[0].GetErrorYhigh(ibin)
                c = ctypes.c_double(0.)
                d = ctypes.c_double(0.)
                ih[1].GetPoint(ibin, c, d )
                yvaluesReco[tmpNbin] = d.value
                yErrLowValuesReco[tmpNbin] = ih[1].GetErrorYlow(ibin)
                yErrHighValuesReco[tmpNbin] = ih[1].GetErrorYhigh(ibin)
                tmpNbin = tmpNbin+1

    xvalues = np.array( range(1,numBins[-1]+1), 'd' )
    xvalues = xvalues - 0.5

    dictHistos = {}
    dictHistos['combRecoRatio'] = ROOT.TGraphAsymmErrors(len(xvalues), np.array(xvalues, 'd'), yvaluesReco, np.array([0.5]*len(xvalues), 'd' ), np.array([0.5]*len(xvalues), 'd' ), yErrLowValuesReco, yErrHighValuesReco )
    legend.AddEntry( dictHistos[ 'combRecoRatio' ], 'detector level (stat. unc. only)', 'lep' )
    dictHistos['combGenRatio'] = ROOT.TGraphAsymmErrors(len(xvalues), np.array(xvalues, 'd'), yvaluesGen, np.array([0]*len(xvalues), 'd' ), np.array([0]*len(xvalues), 'd' ), yErrLowValuesGen, yErrHighValuesGen )
    legend.AddEntry( dictHistos[ 'combGenRatio' ], 'hadron level (all uncorr. unc.)', 'lep' )


    canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 500 )
    canvas[outputFileName].SetGridy()
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.08)
    ROOT.gStyle.SetPadTickX(0)
    dictHistos['combRecoRatio'].SetLineColor( ROOT.kBlack )
    dictHistos['combRecoRatio'].SetMarkerColor( ROOT.kBlack )
    dictHistos['combRecoRatio'].SetLineWidth( 1 )

    dictHistos['combGenRatio'].SetLineColor( ROOT.kMagenta )
    dictHistos['combGenRatio'].SetMarkerColor( ROOT.kMagenta )
    dictHistos['combGenRatio'].SetLineWidth( 1 )

    multiGraph = ROOT.TMultiGraph()
    multiGraph.Add( dictHistos['combRecoRatio'] )
    multiGraph.Add( dictHistos['combGenRatio'] )

    multiGraph.GetXaxis().SetNdivisions(100)
    multiGraph.GetYaxis().SetTitle( 'Sim. / Data'+(' (MC)' if process.startswith('MC') else '') )
    multiGraph.GetYaxis().SetTitleSize( 0.06 )
    multiGraph.GetYaxis().SetTitleOffset( 0.6 )
    multiGraph.GetXaxis().SetLimits( 0., xvalues[-1]+1 )
    multiGraph.SetMaximum( 4.5 )
    multiGraph.SetMinimum( 0.5)
    multiGraph.Draw('AP')

    aBox = textBox.Clone()
    aBox.SetTextSize(0.06)
    #aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Central Jet' if name.startswith('recoJet2') else 'Outer Jet')+"}{Dijet Selection}}"  )
    if 'WSel' in selection: aBox.DrawLatex( 5, 4.10, '#bf{Boosted W selection}'  )
    elif 'topSel'in selection: aBox.DrawLatex( 5, 4.10, '#bf{Boosted top selection}'  )

    ### division lines
    lines = {}
    for i in numBins[1:-1]:
        lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
        lines[i].SetLineColor(ROOT.kGray)
        lines[i].Draw('same')

    CMS_lumi.lumiTextSize = 0.6
    CMS_lumi.relPosX = 0.07
    CMS_lumi.CMS_lumi( canvas[outputFileName], 4, 0)
    legend.Draw()

    textBox.SetTextAlign(12)
    textBoxList = {}
    
    print (len(numBins),len(Xlabels), numBins, Xlabels)

    for i in range(1, len(numBins)):
        textBoxList[i] = textBox.Clone()
        textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., -0.1, Xlabels[i-1] )
    
    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]
    ROOT.gStyle.SetPadTickX(1)
    
def bottomLineTest( ivar, dataHisto, mainHistoLabel, datarecoHisto, recoHisto, genHisto, datacovMatrix, hadcovMatrix, varInfo, outputLabel, outputDir, runMLU=False, ext='pdf',selection='_dijetSel',version='vNew' ):
    """based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L74"""

    outputDir=outputDir+'/'+ivar+'/'+process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)

    #recoHisto.Rebin( 2 )  ### because data and covMatrix have less number of bins
    #datarecoHisto.Rebin(2)
    #recoHisto.Scale( dataHisto.Integral()/recoHisto.Integral() )
    
    #genHisto.Scale( dataHisto.Integral()/genHisto.Integral() )

    #dataHisto.Scale(datarecoHisto.Integral()/dataHisto.Integral())

    #datarecoHisto.Scale( 1 )
    
    
    ##### computing chi2 and inverted matrix
    vector = []
    ndf = 0
    for ibin in range(1, datacovMatrix.GetNbinsX()+1):
        if (datacovMatrix.GetBinContent(ibin, ibin) > 0.):
            ndf = ndf + 1
            vector.append( datarecoHisto.GetBinContent( ibin ) - recoHisto.GetBinContent( ibin ) )
            
    print (ndf,datacovMatrix.GetNbinsX(),datacovMatrix.GetNbinsY())
    #assert ndf==covMatrix.GetNbinsX()
    matrix = np.eye( ndf, ndf )
    
    for ibin in range(1, datacovMatrix.GetNbinsX()+1):
        if (datacovMatrix.GetBinContent(ibin, ibin) > 0.):
            for jbin in range(1, datacovMatrix.GetNbinsY()+1):
                matrix[ibin-1][jbin-1] = datacovMatrix.GetBinContent( ibin, jbin )
    
    vector = np.array( vector )
    invMatrix = np.linalg.inv( matrix )
    chi2 = np.dot( vector, np.dot( invMatrix, vector ) )
    print('Detector level: chi2, ndf, chi2/ndf = ', chi2, ndf, chi2/ndf)
    
    
    vector = []
    ndf = 0
    for ibin in range(1, hadcovMatrix.GetNbinsX()+1):
        if (hadcovMatrix.GetBinContent(ibin, ibin) > 0.):
            ndf = ndf + 1
            vector.append( dataHisto.GetBinContent( ibin ) - genHisto.GetBinContent( ibin ) )
            
    print (ndf,hadcovMatrix.GetNbinsX(),hadcovMatrix.GetNbinsY())
    
    matrix = np.eye( ndf, ndf)
    
    for ibin in range(1, hadcovMatrix.GetNbinsX()+1):
        if (hadcovMatrix.GetBinContent(ibin, ibin) > 0.):
            for jbin in range(1, hadcovMatrix.GetNbinsY()+1):
                matrix[ibin-1][jbin-1] = hadcovMatrix.GetBinContent( ibin, jbin )

    #assert ndf==covMatrix.GetNbinsX()
    vector = np.array( vector )
    invMatrix = np.linalg.inv( matrix )
    chi2 = np.dot( vector, np.dot( invMatrix, vector ) )
    print('Hadron level: chi2, ndf, chi2/ndf = ', chi2, ndf, chi2/ndf)
    
    


    #### plotting inverted matrix
    invertedMatrix = hadcovMatrix.Clone()
    invertedMatrix.Reset()
    for ibin in range(1, hadcovMatrix.GetNbinsX()):
        for jbin in range(1, hadcovMatrix.GetNbinsY()):
            invertedMatrix.SetBinContent( ibin, jbin, invMatrix[ibin-1][jbin-1]  )
    draw2D( ivar, invertedMatrix, varInfo, outputLabel=outputLabel+'_invertedMatrix', outputDir=outputDir.split(ivar)[0],selection=selection,version=version )

    ##### plotting ratios together
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_bottomLineTest_'+version+'.'+ext

    canRatio = ROOT.TCanvas('canRatio'+ivar, 'canRatio'+ivar,  10, 10, 750, 500 )

    genRatio = ROOT.TGraphAsymmErrors()
    genRatio.Divide( dataHisto, genHisto, 'pois' )
    recoRatio = ROOT.TGraphAsymmErrors()
    recoRatio.Divide( datarecoHisto, recoHisto, 'pois' )

    legend=ROOT.TLegend(0.15,0.70,0.40,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.AddEntry( genRatio, 'Hadron level (w/ stat + unf. unc.)', 'pl' )
    legend.AddEntry( recoRatio, 'Detector level (stat unc.)', 'pl' )
    legend.SetBorderSize(0)
    
    genRatio.SetLineWidth(2)
    genRatio.GetXaxis().SetTitle(varInfo['label'])
    genRatio.GetXaxis().SetLimits( varInfo['bins'][0], varInfo['bins'][-1] )
    genRatio.GetYaxis().SetTitle('Data / Simulation')
    genRatio.GetYaxis().SetTitleOffset( 0.8 )
    genRatio.SetMarkerStyle(8)
    genRatio.GetYaxis().SetRangeUser(-5.,10.)
    genRatio.Draw('AP0')
    recoRatio.SetLineColor(ROOT.kRed)
    recoRatio.SetLineWidth(2)
    recoRatio.SetMarkerStyle(4)
    recoRatio.Draw('P0 same')

    lineOne = ROOT.TGraph(2, array('d', [0, 1]), array('d', [1, 1]))
    lineOne.Draw('same')

    legend.Draw()
    CMS_lumi.extraText = " Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canRatio, 4, 0)
    canRatio.SaveAs(outputName)

    return [genRatio, recoRatio]

def makeJacobian(aTH1, aJac):
    N = aTH1.Integral(0,aTH1.GetNbinsX())
    for i in range(0,aTH1.GetNbinsX()):
        for j in range(0,aTH1.GetNbinsX()):
            if i==j: aJac[i][j]=1.*(N-aTH1.GetBinContent(i))/N**2 
            else: aJac[i][j]=(-1.*aTH1.GetBinContent(i))/N**2 
                
                

def drawUnfold( ivar, selection, process, year, dataJetHisto, genJetHisto, unfoldHisto, unfoldHistowoUnc, altMCHisto,
               foldHisto, recoJetHisto, cov_tot, cov_datastat_tot, labelX, maxX, tlegendAlignment, outputName ):
    """docstring for drawUnfold"""
    print ("Drawing unfolding")
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad1.Draw()

    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.67,0.90,0.9)

    else: legend=ROOT.TLegend(0.20,0.67,0.40,0.9)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    
    bins = variables[ivar]['bins']

    unfoldHistoTot = unfoldHisto.Clone()
    #print (unfoldHisto)
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the covariance matrix, from the normalise -> the unnormalised space
    normed_cov_tot_matrix, normed_cov_tot = GetNormalizedTMatrixandTH2(cov_tot.Clone(),"normed_cov_tot", unfoldHisto.Clone())
    
    normed_cov_datastat_tot_matrix, normed_cov_datastat_tot = GetNormalizedTMatrixandTH2(cov_datastat_tot.Clone(),"normed_cov_dastat_tot", unfoldHisto.Clone())
    '''
    normed_cov_bkgsub_tot_matrix, normed_cov_bkgsub_tot = GetNormalizedTMatrixandTH2(cov_bkg_tot.Clone(),
                                                                                     "normed_cov_tot",
                                                                                     unfoldHistoTotUnc)
    
    normed_cov_rmstat_tot_matrix, normed_cov_rmstat_tot = GetNormalizedTMatrixandTH2(cov_rmstat_tot.Clone(),
                                                                                     "normed_cov_rmstat_tot", 
                                                                                     unfoldHistoTotUnc)
    '''
    
    #print (unfoldHisto)

    
    unfoldHistoDataStatErr = normalise_hist(unfoldHisto.Clone())    
    unfoldHisto = normalise_hist(unfoldHisto.Clone())    
    dataJetHisto = normalise_hist(dataJetHisto.Clone())
    genJetHisto = normalise_hist(genJetHisto.Clone())
    unfoldHistowoUnc = normalise_hist_divide_bin_width(unfoldHistowoUnc.Clone())
    altMCHisto = normalise_hist(altMCHisto.Clone())
    foldHisto = normalise_hist(foldHisto.Clone())
    recoJetHisto = normalise_hist(recoJetHisto.Clone())
    #print (unfoldHisto)

    #norm = unfoldHistoTot.Integral("width")/unfoldHisto.Integral("width")
    
    #normed_cov_tot = cov_tot.Clone()
    #normed_cov_tot.Scale(1./(norm*norm))
    #scale_th2_bin_widths(normed_cov_tot, variables[ivar]['bins'])
    #normed_cov_datastat_tot = cov_datastat_tot.Clone()
    #normed_cov_datastat_tot.Scale(1./(norm*norm))
    #scale_th2_bin_widths(normed_cov_datastat_tot, variables[ivar]['bins'])   
                                                         
    
    
    for ibin in range(1, unfoldHisto.GetNbinsX()+1):
        tot_err = normed_cov_tot.GetBinContent(ibin,ibin)
        tot_datastaterr = normed_cov_datastat_tot.GetBinContent(ibin,ibin)
        
        if tot_err<=0.: tot_err=0.
        else: tot_err = np.sqrt(tot_err)
        
        if tot_datastaterr<=0.: tot_datastaterr=0.
        else: tot_datastaterr = np.sqrt(tot_datastaterr)    
        
        unfoldHisto.SetBinError(ibin,tot_err)
        unfoldHistoDataStatErr.SetBinError(ibin,tot_datastaterr)
            
        
    unfoldHistoDataStatErr.Scale(1,'width')
    #dataIntegral = dataJetHisto.Integral()
    #unfoldIntegral = unfoldHisto.Integral()
    #tmpMax = float(unfoldHisto.GetMaximum())
    
    unfoldHisto.Scale(1, 'width')  ### divide by bin width
    #unfoldHisto.Scale(1/unfoldHisto.Integral(), 'width')  ### divide by bin width
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(2)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    genJetHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    #genJetHisto.Scale(1/genJetHisto.Integral(), 'width')  ### divide by bin width
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(ROOT.kRed)
    genJetHisto.SetMarkerColor(ROOT.kRed)
    genJetHisto.SetMarkerStyle(25)
    legend.AddEntry( genJetHisto, 'MG5+Pythia8', 'lp' )

   

    unfoldHisto.GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( 1.6*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )

    unfoldHisto.Draw( "E")
    genJetHisto.Draw( "histe same")

    altMCHisto.Scale(1, 'width')  ### divide by bin width
    altMCHisto.SetLineWidth(2)
    altMCHisto.SetLineColor(ROOT.kBlue)
    altMCHisto.SetMarkerColor(ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    legend.AddEntry( altMCHisto, 'Herwig7', 'lp' )
    altMCHisto.Draw("histe same")
    #unfoldHistowoUnc.Draw( "e1 same")
    #foldHisto.Draw( "histe same")
    #if selfClosure: foldHisto.Draw( "histe same")
    #dataJetHisto.Draw( "histe same")
    #recoJetHisto.Draw( "histe same")
    #if process.startswith('data'): recoJetHisto.Draw( "hist same")
    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.045)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = 'Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
    elif selection.startswith("_W"): seltext = 'Boosted W region'
    elif selection.startswith("_top"): seltext = 'Boosted top region'
    selText.DrawLatex( ( 0.21 if tlegendAlignment.startswith('right') else 0.55 ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65#leqm_{SD}<115 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>350 GeV, 140#leqm_{SD}<220 GeV'
    selText.DrawLatex( ( 0.21 if tlegendAlignment.startswith('right') else 0.55 ), 0.78, seltext )

    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if selection.startswith('dijet') else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistowoUnc.Clone()
    ratio_datastatUnc.Divide(unfoldHistoDataStatErr)
    #ratio_datastatUnc.Reset()
    ratio_totalUnc = unfoldHistowoUnc.Clone()
    ratio_totalUnc.Divide(unfoldHisto)
    #ratio_totalUnc.Reset()
    '''
    for ibin in range(0, unfoldHisto.GetNbinsX()):
        
        tot_err = unfoldHisto.GetBinError(ibin)
        tot_datastaterr = unfoldHistoDataStatErr.GetBinError(ibin)
        bc = unfoldHisto.GetBinContent(ibin)
        
        ratio_datastatUnc.SetBinContent(ibin,1.)
        ratio_totalUnc.SetBinContent(ibin,1.)
                                         
        if bc>0.:
            ratio_datastatUnc.SetBinError(ibin,tot_datastaterr/bc)
            ratio_totalUnc.SetBinError(ibin,tot_err/bc)
   
    '''
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    print (labelX)
    #tmpPad2.GetXaxis().SetTitle( labelX )
    tmpPad2.GetYaxis().SetTitle( "Sim./Data" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(-1., 3. )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    #hRatioDown.Draw('P same')
    ratio_datastatUnc.SetFillColor(ROOT.kAzure+7)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    ratio_datastatUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    ratio_datastatUnc.GetYaxis().SetTitle( "Sim./Data" )
    ratio_datastatUnc.GetYaxis().SetTitleOffset( 0.5 )
    ratio_datastatUnc.GetYaxis().SetRangeUser(-1., 3. )
    ratio_datastatUnc.GetYaxis().CenterTitle()
    ratio_datastatUnc.GetXaxis().SetLabelSize(0.12)
    ratio_datastatUnc.GetXaxis().SetTitleSize(0.12)
    ratio_datastatUnc.GetYaxis().SetLabelSize(0.12)
    ratio_datastatUnc.GetYaxis().SetTitleSize(0.12)
    ratio_datastatUnc.GetXaxis().SetNdivisions(505)
    ratio_datastatUnc.GetYaxis().SetNdivisions(505)
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)
    #ratio_datastatUnc.Scale(1,'width')
    ratio_datastatUnc.Draw('E5')

    ratio_totalUnc.SetFillColor(ROOT.kGray+3)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3254)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    #ratio_totalUnc.Scale(1,'width')    
    ratio_totalUnc.Draw('E5 SAME')
    
    #ratiosystUncHisto.SetFillColor(9)
    #ratiosystUncHisto.SetLineColor(0)
    #ratiosystUncHisto.SetFillStyle(3006)
    #ratiosystUncHisto.SetMarkerStyle(0)
    #ratiosystUncHisto.Scale(1,'width')    
    #ratiosystUncHisto.Draw('E3 same')

    hRatioUp = ROOT.TGraphAsymmErrors()
    hRatioUp.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatioUp.SetLineColor(ROOT.kRed)
    hRatioUp.SetMarkerColor(ROOT.kRed)
    #hRatioUp.SetLineWidth(2)
    hRatioUp.SetMarkerStyle(25)
    hRatioUp.Draw('P0 same')
    
    hRatioUp2 = ROOT.TGraphAsymmErrors()
    hRatioUp2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatioUp2.SetLineColor(ROOT.kBlue)
    hRatioUp2.SetMarkerColor(ROOT.kBlue)
    #hRatioUp.SetLineWidth(2)
    hRatioUp2.SetMarkerStyle(25)
    hRatioUp2.Draw('P0 same')
    
    ratioLegend=ROOT.TLegend(0.20,0.85,0.80,0.95)
    ratioLegend.SetTextSize(0.095)
    ratioLegend.SetNColumns(3)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_totalUnc, 'Data total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)
    
def drawUncertainties(ivar, unfoldHistoTotUnc,
                      uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot,
                      labelX, tlegendAlignment, outputName,year,unftot ):
    """docstring for drawUncUncertainties"""
    
    #Draw fractional uncertainties for all sources of systematics and statical errors
    #Systematics are obtained from the total of up/down shifts normalized per bin by the nominal unfolded 
    #histos respective bins'contents, this is with the exception of the bkg_tot_err and RM_staterr which are obtained 
    #from the variances (diagonal of the respective cov. matrix) and then normalized
    #Stat errors are obtained from the diagonals of the respective covariance matrices provided by TUnfold and normalized
    
    normed_cov_tot_matrix, normed_cov_tot = GetNormalizedTMatrixandTH2(cov_tot.Clone(),"normed_cov_tot", 
                                                                       unfoldHistoTotUnc.Clone())
    
    normed_cov_datastat_tot_matrix, normed_cov_datastat_tot = GetNormalizedTMatrixandTH2(cov_datastat_tot.Clone(),
                                                                                         "normed_cov_dastat_tot", 
                                                                                          unfoldHistoTotUnc.Clone())
   
    normed_cov_bkgsub_tot_matrix, normed_cov_bkgsub_tot = GetNormalizedTMatrixandTH2(cov_bkg_tot.Clone(),
                                                                                     "normed_cov_tot",
                                                                                     unfoldHistoTotUnc.Clone())
    
    normed_cov_rmstat_tot_matrix, normed_cov_rmstat_tot = GetNormalizedTMatrixandTH2(cov_rmstat_tot.Clone(),
                                                                                     "normed_cov_rmstat_tot", 
                                                                                     unfoldHistoTotUnc.Clone())
   
    
    unfoldHistoNoNorm=unfoldHistoTotUnc.Clone()
    
    unfoldHistoDataStatUnc=unfoldHistoTotUnc.Clone()
    unfoldHistoRMStatUnc=unfoldHistoTotUnc.Clone()
    unfoldHistoBkgSubUnc=unfoldHistoTotUnc.Clone()
    
    unfoldHistoTotUnc = normalise_hist_divide_bin_width(unfoldHistoTotUnc.Clone())
    unfoldHistoDataStatUnc = normalise_hist_divide_bin_width(unfoldHistoDataStatUnc.Clone())
    unfoldHistoRMStatUnc = normalise_hist_divide_bin_width(unfoldHistoRMStatUnc.Clone())
    unfoldHistoBkgSubUnc = normalise_hist_divide_bin_width(unfoldHistoBkgSubUnc.Clone())
    
    scale_th2_bin_widths(normed_cov_tot_matrix)
    scale_th2_bin_widths(normed_cov_datastat_tot_matrix)
    scale_th2_bin_widths(normed_cov_bkgsub_tot_matrix)
    scale_th2_bin_widths(normed_cov_rmstat_tot_matrix)
    
    for ibin in range(1, unfoldHistoNoNorm.GetNbinsX()+1):
        tot_err = normed_cov_tot.GetBinContent(ibin,ibin)
        tot_datastaterr = normed_cov_datastat_tot.GetBinContent(ibin,ibin)
        tot_rmstaterr = normed_cov_rmstat_tot.GetBinContent(ibin,ibin)
        tot_bkgsuberr = normed_cov_bkgsub_tot.GetBinContent(ibin,ibin)
        
        if tot_err<=0.: tot_err=0.
        else: tot_err = np.sqrt(tot_err)
        
        if tot_datastaterr<=0.: tot_datastaterr=0.
        else: tot_datastaterr = np.sqrt(tot_datastaterr)    
     
        if tot_rmstaterr<=0.: tot_rmstaterr=0.
        else: tot_rmstaterr = np.sqrt(tot_rmstaterr)
        
        if tot_bkgsuberr<=0.: tot_bkgsuberr=0.
        else: tot_bkgsuberr = np.sqrt(tot_bkgsuberr)
        
        unfoldHistoTotUnc.SetBinError(ibin,tot_err)
        unfoldHistoDataStatUnc.SetBinError(ibin,tot_datastaterr)
        unfoldHistoRMStatUnc.SetBinError(ibin,tot_datastaterr)
        unfoldHistoBkgSubUnc.SetBinError(ibin,tot_bkgsuberr)
        
        
    #unfoldHistoNoNorm.Scale(1,'width')
    #unfoldHistoTotUnc.Scale(1,'width')
    #unfoldHistoDataStatUnc.Scale(1,'width')
    #unfoldHistoRMStatUnc.Scale(1,'width')
    #unfoldHistoBkgSubUnc.Scale(1,'width')
    
    colors = [ 2, 3, 4, 6, 7, 8, 9, 50, 205, 225, 94, 221, 92, 16, 28, 219, 225, 128,]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 750, 500 )
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.7,0.90,0.92)
    else: legend=ROOT.TLegend(0.20,0.7,0.40,0.92)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.SetTextSize(0.02)
    legend.SetBorderSize(0)
    
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Reset()
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Reset()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Reset()
    bkgSubErrHist = unfoldHistoTotUnc.Clone()
    bkgSubErrHist.Reset()
    
    
    
    #print (uncerUnfoldHisto)
    
    #print (uncerUnfoldHisto.keys())
    normeduncerUnfoldHisto = OrderedDict()
    for k in uncerUnfoldHisto:
        if ('Total'in k ) and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            normeduncerUnfoldHisto[k]=uncerUnfoldHisto[k].Clone()
            #normeduncerUnfoldHisto[k].Scale(1,'width')
            #normeduncerUnfoldHisto[k].Divide(unfoldHistoTotUnc)
            for i in range(1,normeduncerUnfoldHisto[k].GetNbinsX()+1):
                if unfoldHistoNoNorm.GetBinContent(i)>0.:
                    if unfoldHistoTotUnc.GetBinContent(i)>0: 
                        perbin_scale = 1.0/(unfoldHistoNoNorm.GetBinContent(i)/unfoldHistoTotUnc.GetBinContent(i))
                    else: perbin_scale = 0
                    normeduncerUnfoldHisto[k].SetBinContent(i, 
                                                            normeduncerUnfoldHisto[k].GetBinContent(i)*perbin_scale)

    
            
            text = k.split('_')[-1].split('Total')[0].replace('TOTAL', '').replace('WEIGHT', '')
            if text.startswith(('ISR', 'FSR', 'JER', 'PU', 'PDF', 'Bkg')):
                uncerUnfoldHisto[k].SetLineStyle(2)
            uncerUnfoldHisto[k].SetLineColor(colors[dummy])
            uncerUnfoldHisto[k].SetLineWidth(2)
            #uncerUnfoldHisto[k].Scale( uncScaleFactor )
            uncerUnfoldHisto[k].Draw("hist same")
            if not 'damp' in k: legend.AddEntry( uncerUnfoldHisto[k], k.split('_')[-1].split('Total')[0].replace('TOTAL', '').replace('WEIGHT', ''), 'l' )
            else: legend.AddEntry( uncerUnfoldHisto[k], 'h_{damp}', 'l' )
            dummy=dummy+1
            
        
    #get fractionals for stat and total errors
    for i in range(1,dataStatErrHist.GetNbinsX()+1):
    
        if unfoldHistoTotUnc.GetBinContent(i)>0.:
            dataStatErrHist.SetBinContent(i,unfoldHistoDataStatUnc.GetBinError(i)/unfoldHistoTotUnc.GetBinContent(i))
            rmStatErrHist.SetBinContent(i,unfoldHistoRMStatUnc.GetBinError(i)/unfoldHistoTotUnc.GetBinContent(i))
            bkgSubErrHist.SetBinContent(i,unfoldHistoBkgSubUnc.GetBinError(i)/unfoldHistoTotUnc.GetBinContent(i))
            totalErrHist.SetBinContent(i,unfoldHistoTotUnc.GetBinError(i)/unfoldHistoTotUnc.GetBinContent(i))
        else:
            dataStatErrHist.SetBinContent(i,0)
            rmStatErrHist.SetBinContent(i,0)
            bkgSubErrHist.SetBinContent(i,0)
            totalErrHist.SetBinContent(i,0)
            
        dataStatErrHist.SetBinError(i,0)
        rmStatErrHist.SetBinError(i,0)
        bkgSubErrHist.SetBinError(i,0)
        totalErrHist.SetBinError(i,0)
            
    
    
    dataStatErrHist.SetLineWidth(2)
    dataStatErrHist.SetLineStyle(3)
    #dataStatErrHist.Scale( uncScaleFactor )
    dataStatErrHist.GetXaxis().SetTitle(labelX)
    dataStatErrHist.GetYaxis().SetTitle('Relative Uncertainty')
    dataStatErrHist.SetMaximum(5)
    dataStatErrHist.GetYaxis().SetRangeUser(5*(10**(-5)),15)
    dataStatErrHist.SetLineColor(88)
    dataStatErrHist.Draw('hist ')
    
    rmStatErrHist.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(3)
    #dataStatErrHist.Scale( uncScaleFactor )
    #rmStatErrHist.GetXaxis().SetTitle(labelX)
    #rmStatErrHist.GetYaxis().SetTitle('Relative Uncertainty')
    #rmStatErrHist.SetMaximum(5)
    #rmStatErrHist.GetYaxis().SetRangeUser(5*(10**(-4)),8)
    rmStatErrHist.SetLineColor(209)
    rmStatErrHist.Draw('hist same ')
    
    bkgSubErrHist.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(3)
    #dataStatErrHist.Scale( uncScaleFactor )
    #rmStatErrHist.GetXaxis().SetTitle(labelX)
    #rmStatErrHist.GetYaxis().SetTitle('Relative Uncertainty')
    #rmStatErrHist.SetMaximum(5)
    #rmStatErrHist.GetYaxis().SetRangeUser(5*(10**(-4)),8)
    bkgSubErrHist.SetLineColor(ROOT.kMagenta-4)
    bkgSubErrHist.Draw('hist same ')
    
    
    dummy=0
    #print (uncerUnfoldHisto.keys())
    
    for k in normeduncerUnfoldHisto:
        if ('Total'in k ) and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'Bkg' in k and not 'CM' in k and not 'jes' in k and not 'JES' in k:
            text = k.split('_')[-1].split('Total')[0].replace('TOTAL', '').replace('WEIGHT', '')
            if text.startswith(('ISR', 'FSR', 'JER', 'PU', 'PDF', 'Bkg')):
                normeduncerUnfoldHisto[k].SetLineStyle(2)
            normeduncerUnfoldHisto[k].SetLineColor(colors[dummy])
            normeduncerUnfoldHisto[k].SetLineWidth(2)
            #uncerUnfoldHisto[k].Scale( uncScaleFactor )
            normeduncerUnfoldHisto[k].Draw("hist same")
            if not 'damp' in k: legend.AddEntry( normeduncerUnfoldHisto[k], k.split('_')[-1].split('Total')[0].replace('TOTAL', '').replace('WEIGHT', ''), 'l' )
            else: legend.AddEntry( normeduncerUnfoldHisto[k], 'h_{damp}', 'l' )
            dummy=dummy+1
    #print ("adding jeshisto values")
    
   
    totalErrHist.SetLineWidth(3)
    totalErrHist.SetLineStyle(2)
    #totalErrHist.Scale( uncScaleFactor )
    #totalErrHist.GetXaxis().SetTitle(labelX)
    #totalErrHist.GetYaxis().SetTitle('Fractional Uncertainty')
    totalErrHist.SetMarkerStyle(4)
    totalErrHist.SetLineColor(1)
    #totalErrHist.GetYaxis().SetRangeUser(5*(10**(-4)),5)
    legend.AddEntry( rmStatErrHist, 'Bkg. Sub. Unc.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data Stat. Unc.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'RM Stat. Unc.', 'l' )    
    
    legend.AddEntry( totalErrHist, 'Total Unc.', 'l' )   
    
    totalErrHist.Draw('hist same')
    
    
    
    #uncerUnfoldHisto[ivar+'_SystTotal'].SetLineWidth(2)
    #uncerUnfoldHisto[ivar+'_SystTotal'].SetLineStyle(2)
    #uncerUnfoldHisto[ivar+'_SystTotal'].Scale( uncScaleFactor )
    #uncerUnfoldHisto[ivar+'_SystTotal'].SetLineColor(209)
    #legend.AddEntry( uncerUnfoldHisto[ivar+'_SystTotal'], 'Syst. Tot.', 'l' )    
    #uncerUnfoldHisto[ivar+'_SystTotal'].Draw('hist same')
    
    
    
    legend.Draw()
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    canUnc.SetLogy()
    canUnc.Update()

    canUnc.SaveAs(outputName)
    
    
def drawUncertainties_normalizedshifts(ivar, unfoldHistoTotUnc, unfoldHistowoUnc,
                      uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot,
                      labelX, tlegendAlignment, outputName,year,unftot ):
    """docstring for drawUncUncertainties"""
    #outputRoot = ROOT.TFile.Open( 'draw_unc'+ivar+'.root', 'recreate' )

    #Absolute shifts on the sbolute unfolded distribution are first created
    #this distribtuion is normalized taking into account the bin widths as for the nominal 
   
    colors = [ 95, 38, 6, 7, 8, 42, 50, 218, 225, 30, 16, 198, 190, 83, 167, 207, 209, 212, 216, 51, 61, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)
    #canUnc.SetBottomMargin(0.02)
    
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.20,0.65,0.55,0.9)
    else: legend=ROOT.TLegend(0.55,0.65,0.9,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm=unfoldHistoTotUnc.Clone()
    unfoldHistoDataStatUnc=unfoldHistoTotUnc.Clone()
    unfoldHistoRMStatUnc=unfoldHistoTotUnc.Clone()
    unfoldHistoBkgSubUnc=unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc = normalise_hist(unfoldHistoNoNorm.Clone())
    unfoldHistoTotUnc = normalise_hist(unfoldHistoNoNorm.Clone())
    unfoldHistoDataStatUnc = normalise_hist(unfoldHistoDataStatUnc.Clone())
    unfoldHistoRMStatUnc = normalise_hist(unfoldHistoRMStatUnc.Clone())
    unfoldHistoBkgSubUnc = normalise_hist(unfoldHistoBkgSubUnc.Clone())

    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    
    cov_tot.Sumw2()
    cov_datastat_tot.Sumw2()
    cov_bkg_tot.Sumw2()
    cov_rmstat_tot.Sumw2()
    
    normed_cov_tot_matrix, normed_cov_tot = GetNormalizedTMatrixandTH2(cov_tot.Clone(),"normed_cov_tot", 
                                                                       unfoldHistoNoNorm.Clone())
    
    normed_cov_datastat_tot_matrix, normed_cov_datastat_tot = GetNormalizedTMatrixandTH2( cov_datastat_tot.Clone(), "normed_cov_dastat_tot", unfoldHistoNoNorm.Clone() )
   
    normed_cov_bkgsub_tot_matrix, normed_cov_bkgsub_tot = GetNormalizedTMatrixandTH2(cov_bkg_tot.Clone(),
                                                                                     "normed_cov_tot",
                                                                                     unfoldHistoNoNorm.Clone())
    
    normed_cov_rmstat_tot_matrix, normed_cov_rmstat_tot = GetNormalizedTMatrixandTH2(cov_rmstat_tot.Clone(),
                                                                                     "normed_cov_rmstat_tot", 
                                                                                     unfoldHistoNoNorm.Clone())
    #add errors now to unf. histo without dealing with bins divided by width
    #then ask root to handle that part of it for you
    
     
    for ibin in range(1, unfoldHistoNoNorm.GetNbinsX()+1):
        tot_err = normed_cov_tot.GetBinContent(ibin,ibin)
        tot_datastaterr = normed_cov_datastat_tot.GetBinContent(ibin,ibin)
        tot_rmstaterr = normed_cov_rmstat_tot.GetBinContent(ibin,ibin)
        tot_bkgsuberr = normed_cov_bkgsub_tot.GetBinContent(ibin,ibin)
        
        if tot_err<=0.: tot_err=0.
        else: tot_err = np.sqrt(tot_err)
        
        if tot_datastaterr<=0.: tot_datastaterr=0.
        else: tot_datastaterr = np.sqrt(tot_datastaterr)    
     
        if tot_rmstaterr<=0.: tot_rmstaterr=0.
        else: tot_rmstaterr = np.sqrt(tot_rmstaterr)
        
        if tot_bkgsuberr<=0.: tot_bkgsuberr=0.
        else: tot_bkgsuberr = np.sqrt(tot_bkgsuberr)
        
        unfoldHistoTotUnc.SetBinError(ibin,tot_err)
        unfoldHistoDataStatUnc.SetBinError(ibin,tot_datastaterr)
        unfoldHistoRMStatUnc.SetBinError(ibin,tot_rmstaterr)
        unfoldHistoBkgSubUnc.SetBinError(ibin,tot_bkgsuberr)
    
    unfoldHistowoUnc.Scale(1, 'width')
    unfoldHistoNoNorm.Scale(1,'width')
    unfoldHistoTotUnc.Scale(1,'width')
    unfoldHistoDataStatUnc.Scale(1,'width')
    unfoldHistoRMStatUnc.Scale(1,'width')
    unfoldHistoBkgSubUnc.Scale(1,'width')
    
    dummy=0
    dummy_jes=0
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    up_counter=0
    down_counter=0
    
    upstyles=[20,22,21,23,29,34,47,33,43,39,41,]
    downstyles=[24,26,25,32,30,28,46,27,42,37,40]
    
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone()
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone()
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    for k in uncerUnfoldHisto:
        text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        text=text.upper()
        if 'crtotal' in text.lower() or 'model' in text.lower(): continue
        if ('_shifthist'in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1,'width')
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'FSR' in text or 'JER' in text or ('PU' in text and 'DAMP' not in text) or 'PDF' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                #print (k,text, upstyles[up_counter], up_counter, dummy)
                dummy=dummy+1    
                up_counter=up_counter+1
            elif 'DAMP' in text or 'MTOP' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(3)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                #print (k,text, upstyles[up_counter], up_counter, dummy)
                dummy=dummy+1    
                up_counter=up_counter+1
    dummy=0
    
    for k in uncerUnfoldHisto:
        text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        text=text.upper()
        if 'crtotal' in text.lower() or 'model' in text.lower(): continue
        if ('_shifthist' in k.lower() and 'down' in k.lower()):
            
            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1,'width')
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'FSR' in text or 'JER' in text or ('PU' in text and 'DAMP' not in text) or 'PDF' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                #print (k,text, downstyles[down_counter], down_counter, dummy)
                down_counter=down_counter+1
                dummy=dummy+1 
            elif 'DAMP' in text or 'MTOP' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(3)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                #print (k,text, downstyles[down_counter], down_counter, dummy)
                down_counter=down_counter+1
                dummy=dummy+1 

    dummy = dummy+2
    modelkey = 0
    '''
    cr_maxup = unfoldHistoTotUnc.Clone()
    cr_maxup.Reset()
    cr_maxdown = unfoldHistoTotUnc.Clone()
    cr_maxdown.Reset()
    
    cr_histos = OrderedDict()
    
    for k in uncerUnfoldHisto:
        text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        #print (text, dummy)
        if ('cr1' in text.lower() or 'cr2' in text.lower() or 'erd' in text.lower()) and '_shifthist' in k.lower():
            #print (text, dummy, k)
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            cr_histos[k] = normalise_hist(cr_histos[k].Clone())
            cr_histos[k].Scale(1,'width')
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
    '''
    
    for ibin in range(1,unfoldHistoTotUnc.GetNbinsX()+1):
        #cr_maxup.SetBinContent(ibin,1.)
        #cr_maxdown.SetBinContent(ibin,1.)
        #upmax_ibin = 1
        #downmax_ibin = 1
        
        for k in uncerUnfoldHisto:
            text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            #print (text, dummy)
            if 'model'in k.lower(): modelkey=k
            '''
            if ('cr1' in text.lower() or 'cr2' in text.lower() or 'erd' in text.lower()) and '_shifthist' in k.lower():
                #print (text, dummy, k)
                
                if upmax_ibin<cr_histos[k].GetBinContent(ibin):
                    upmax_ibin=cr_histos[k].GetBinContent(ibin)
                if downmax_ibin>cr_histos[k].GetBinContent(ibin):
                    downmax_ibin=cr_histos[k].GetBinContent(ibin)
               
                
                otherUncs[k] = uncerUnfoldHisto[k].Clone()
                otherUncs[k].Sumw2()
                otherUncs[k] = normalise_hist(otherUncs[k].Clone())
                otherUncs[k].Scale(1,'width')
                otherUncs[k].Divide(unfoldHistowoUnc.Clone())# = convert_syst_shift_to_error_ratio_hist(otherUncs[k].Clone(), unfoldHistoTotUnc.Clone())
                
                #if 'model'in text.lower(): legend.AddEntry( otherUncs[k], 'Physics Model', 'l' )
                #elif 'mass'in text.lower(): legend.AddEntry( otherUncs[k], 'Choice of m_{top}', 'l' )
                if 'cr1' in text.lower(): legend.AddEntry( otherUncs[k], 'CR1', 'l' )
                elif 'cr2' in text.lower(): legend.AddEntry( otherUncs[k], 'CR2', 'l' )
                elif 'erdon' in text.lower(): legend.AddEntry( otherUncs[k], 'erdOn', 'l' )
                dummy = dummy+1
                
            #print (upmax_ibin,downmax_ibin)
            cr_maxup.SetBinContent(ibin,upmax_ibin)
            cr_maxdown.SetBinContent(ibin,downmax_ibin)
            '''        
    #print ("OtherUncs", otherUncs)
    modelUnc = uncerUnfoldHisto[modelkey].Clone()
    modelUnc.Sumw2()
    modelUnc = normalise_hist(modelUnc.Clone())
    modelUnc.Scale(1,'width')
    modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
    modelUnc.SetLineStyle(1)
    modelUnc.SetLineWidth(1)
    modelUnc.SetMarkerSize(0)
    modelUnc.SetLineColor(colors[dummy])
    modelUnc.SetFillColor(0)

               
                
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    dataStatErrHist.Divide(unfoldHistowoUnc)
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    totalErrHist.GetYaxis().SetRangeUser(-0.5,2.5)
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColor(ROOT.kGray+3)
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3154)
    totalErrHist.SetLineColor(1)
    totalErrHist.Draw(' E5')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(3)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColor(ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(88)
    dataStatErrHist.Draw('E5 same')
    
    '''
    ###CR unc
    
    cr_maxup.SetLineColor(208)
    cr_maxdown.SetLineColor(208)
    cr_maxup.SetMarkerColor(208)
    cr_maxdown.SetMarkerColor(208)
    cr_maxup.SetMarkerStyle(upstyles[up_counter])
    cr_maxdown.SetMarkerStyle(downstyles[down_counter])
    cr_maxup.SetMarkerSize(2)
    cr_maxdown.SetMarkerSize(2)
    print (up_counter,down_counter)
    
    cr_maxup.SetMarkerSize(2)
    cr_maxdown.SetMarkerSize(2)
    cr_maxup.Draw('P same')
    cr_maxdown.Draw('P same')
    up_counter=up_counter+1
    down_counter=down_counter+1
    print (up_counter,down_counter)
    '''
    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    rmStatErrHist.Draw('L same ')
    h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineWidth(2)
    h2.SetLineStyle(7)
    bkgSubErrHist.SetLineColor(50)
    h2.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    bkgSubErrHist.Draw('L same ')
    h2.Draw("L same")
    
    #h3 = convert_error_bars_to_error_ratio_hist(modelUnc.Clone(),-1)
    #modelUnc = convert_error_bars_to_error_ratio_hist(modelUnc.Clone(),1)
    modelUnc.Draw('L same')
    #h3.Draw('L same')
    
    
    for k in otherUncs:
        text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print (text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        jesHistoUpMax.SetBinContent(ibin,1)
        jesHistoDownMax.SetBinContent(ibin,1)
        upmax_ibin = 1
        downmax_ibin = 1
        for k in normeduncerUnfoldHistoshiftsUp:
            text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()
            if 'JES' in text:
                if upmax_ibin<normeduncerUnfoldHistoshiftsUp[k].GetBinContent(ibin):
                    upmax_ibin=normeduncerUnfoldHistoshiftsUp[k].GetBinContent(ibin)
        
                
        for k in normeduncerUnfoldHistoshiftsDown:
            text = (k.split('_shifthist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()
            if 'JES' in text:
                if downmax_ibin>normeduncerUnfoldHistoshiftsDown[k].GetBinContent(ibin):
                    downmax_ibin=normeduncerUnfoldHistoshiftsDown[k].GetBinContent(ibin)
        #print ("JES",upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,downmax_ibin)
        
    #legend.AddEntry(cr_maxup,'CR model', 'p')
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(51)
    jesHistoDownMax.SetLineColor(51)
    jesHistoUpMax.SetMarkerColor(51)
    jesHistoDownMax.SetMarkerColor(51)
    jesHistoUpMax.SetMarkerStyle(upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(2)
    jesHistoDownMax.SetMarkerSize(2)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if 'jes' in k.lower() or 'model'in k.lower() or 'bkg' in k.lower() or 'crtotal' in k.lower(): 
            continue
        normeduncerUnfoldHistoshiftsUp[k].Draw("P same")
        
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        text=text.upper()
        
        if 'damp' in k: 
            normeduncerUnfoldHistoshiftsDown[k.replace('UP', 'DOWN')].Draw("P same")
            legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
        elif 'mtop' in k:
            normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
            legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'Choice of m_{top}', 'p' )
        else: 
            normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
            legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)

    
    
    legend.AddEntry( modelUnc, 'Physics Model', 'l' )    
    legend.AddEntry( bkgSubErrHist, 'Bkg. stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data Stat.', 'f' )    
    legend.AddEntry( rmStatErrHist, 'RM Stat.', 'l' )    
    legend.AddEntry( totalErrHist, 'Total Unc.', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
    
def plotSysComparison2( nomHisto, dictUncHistos, outputName, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ):
    """docstring for plot"""
    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 219, 92, 48, 49, 37, 38, 33, 17, 50, 205, 225, 94, 221, 16,  225, 128]
    
    
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 750, 500 )
    #if log: canUnc.SetLogy()
    

    outputFileName = outputName+'_'+version+'.'+ext
    #print ('Processing.......', outputFileName)

    binWidth = nomHisto.GetBinWidth(1)

    legend=ROOT.TLegend(0.35,0.6,0.80,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.02)
    legend.SetBorderSize(0)

    #multiGraph = ROOT.TMultiGraph()
    gnom = nomHisto.Clone()
    gnom.Divide( nomHisto.Clone() )
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetMarkerStyle(1)
    gnom.SetLineWidth(2)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    #multiGraph.Add( gnom )
    dictShifts = {}
    dummy=0
    
    gnom.GetYaxis().SetTitle( 'Ratio Unc/Nominal' )
    gnom.GetXaxis().SetTitle( labelX )
    gnom.SetMaximum( 3. )
    gnom.SetMinimum( -1.)
    gnom.Draw('L')

    #print (dictUncHistos)
    for ih in dictUncHistos:
        #print(ih, dummy, len(colors))
        dictShifts[ih] = dictUncHistos[ih].Clone()
        dictShifts[ih].Divide( nomHisto)
        
        if not dummy==len(colors)-1:
            dictShifts[ih].SetLineColor( colors[dummy] )
            dictShifts[ih].SetLineStyle( 2 )
        else:
            dummy = dummy-len(colors)+2
            dictShifts[ih].SetLineColor( colors[dummy] )
            dictShifts[ih].SetLineStyle( 3 )
        
        if 'jes' in ih: dictShifts[ih].SetLineStyle( 1 )
        dictShifts[ih].SetMarkerStyle(0)
        dictShifts[ih].SetLineWidth( 1 )
        
        if 'jes' in ih and ('2017' in ih or '2018' in ih): 
            print (ih)
            legend.AddEntry( dictShifts[ih], ih.split('_')[1] , 'l' )
        else: legend.AddEntry( dictShifts[ih], ih.split('_')[1] , 'l' )
        dummy=dummy+1
        dictShifts[ih].Draw("L SAME")
        
    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+year
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    if ext.startswith('pdf'):
        canUnc.SaveAs( outputDir + outputFileName.replace('pdf', 'png') )
    del canUnc

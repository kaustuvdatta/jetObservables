#!/usr/bin/env python

import ROOT
import time, os, math, sys, copy
from array import array
import argparse
from collections import OrderedDict
import numpy as np
from variables import nSubVariables
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrStyle
sys.path.insert(0,'../../Skimmer/test/')
from datasets import checkDict, dictSamples

####gReset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()
tdrStyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)

canvas = {}

colors = [ 2, 4, 9, 8, 28, 30, 42, 13 ]
##################################################
def makePlots( name, ext='png', outputDir='Plots/' ):
    """function to plot s and b histos"""

    numBins = 0
    numBinsList = [0]
    dictHistos = OrderedDict()
    uncertDictHistos = OrderedDict()
    for ivar in nSubVariables.keys():
        if args.only and not ivar.startswith( args.only ): continue
        if name.endswith( ivar.split('_')[0] ):
            ###### Loading histograms
            dictHistos[ivar] = {
                    'data' : dictFiles[ivar].Get( 'dataMinusBkgsGenBin' ),
                    'unfold' : dictFiles[ivar].Get( 'unfoldHisto'+ivar ),
                    'cov' : dictFiles[ivar].Get( 'cov'+ivar ),
                    'cov_uncorr_data' : dictFiles[ivar].Get( 'cov_uncorr_data_'+ivar ),
                    'fold' : dictFiles[ivar].Get( 'foldHisto'+ivar ),
                    'gen' : dictFiles[ivar].Get( signalLabel+'_accepgen'+ivar+'_'+args.selection ),
                    'reco' : dictFiles[ivar].Get( signalLabel+'_reco'+ivar+'_nom_'+args.selection ),
                    'gen'+altSignalLabel : dictFiles[ivar].Get( altSignalLabel+'_genJetfrom_resp'+ivar+'_nom_'+args.selection ),
                    'reco'+altSignalLabel : dictFiles[ivar].Get( altSignalLabel+'_recoJetfrom_resp'+ivar+'_nom_'+args.selection )
                    }

            ##### making ratio uncert histos
            uncertDictHistos[ivar] = { 'totalUnc' : dictHistos[ivar]['unfold'].Clone() }
            uncertDictHistos[ivar][ 'totalUnc' ].Reset()
            uncertDictHistos[ivar][ 'statUnc' ] =  uncertDictHistos[ivar][ 'totalUnc' ].Clone()
            for ibin in range( 0, dictHistos[ivar][ 'cov' ].GetNbinsX()+1 ):
                uncertDictHistos[ivar][ 'totalUnc' ].SetBinContent( ibin, 1 )
                uncertDictHistos[ivar][ 'totalUnc' ].SetBinError( ibin, ROOT.TMath.Sqrt( dictHistos[ivar][ 'cov' ].GetBinContent(ibin,ibin) ) )
                uncertDictHistos[ivar][ 'statUnc' ].SetBinContent( ibin, 1 )
                uncertDictHistos[ivar][ 'statUnc' ].SetBinError( ibin, ROOT.TMath.Sqrt( dictHistos[ivar][ 'cov_uncorr_data' ].GetBinContent(ibin,ibin) ) )

            ##### removing tau21 and tau32 from total plots
            if not ivar.endswith(('21', '32')):
                numBins = numBins + dictHistos[ivar]['data'].GetNbinsX()
                numBinsList.append(numBins)

            ##### defining what to include in each plot
            otherHistos = OrderedDict()
            mainHisto = dictHistos[ivar]['unfold']
            if args.process.startswith('MCClosure'):
                outputLabel = args.process+'_'+sigPlotLabel+'_'+altSigPlotLabel
                mainHistoLabel = 'Unfolded (Closure)'
                otherHisto = 'gen'
                otherHistoLabel=sigPlotLabel
                otherHistos[sigPlotLabel] = dictHistos[ivar]['gen']
                dictHistos[ivar]['gen'+altSignalLabel].Rebin(2)
                otherHistos[altSigPlotLabel] = dictHistos[ivar]['gen'+altSignalLabel]
            elif args.process.startswith('MCSelfClosure'):
                outputLabel = args.process+'_'+sigPlotLabel
                mainHistoLabel = 'Unfolded (Self-Closure)'
                otherHisto = 'gen'
                otherHistoLabel='Accepted Gen'
                otherHistos['Accepted Gen'] = dictHistos[ivar]['gen']
                dictHistos[ivar]['fold'].Rebin(2)
                otherHistos['Folded'] = dictHistos[ivar]['fold']
                dictHistos[ivar]['reco'].Rebin(2)
                otherHistos['True Reco '+sigPlotLabel] = dictHistos[ivar]['reco']
            elif args.process.startswith('data'):
                outputLabel = args.process+'_'+sigPlotLabel+'_'+altSigPlotLabel
                mainHistoLabel = 'Data'
                otherHisto = 'gen'
                otherHistoLabel=sigPlotLabel
                otherHistos[sigPlotLabel] = dictHistos[ivar]['gen']
                dictHistos[ivar]['gen'+altSignalLabel].Rebin(2)
                otherHistos[altSigPlotLabel] = dictHistos[ivar]['gen'+altSignalLabel]

            if args.runMLU:
                #Xvalues = []    ### might need it later
                #Yvalues = []
                #YvaluesCombine = []
                #YMaxError = []
                #YMinError = []
                combineTree = dictFiles[ivar+'Combine'].Get('limit')
                otherHistos['MLU'] = dictHistos[ivar]['gen'].Clone()
                otherHistos['MLU'].Reset()
                for ibin in range(1, dictHistos[ivar]['gen'].GetNbinsX()+1):
                    histoCont = dictHistos[ivar]['gen'].GetBinContent(ibin)
                    histoErr = dictHistos[ivar]['gen'].GetBinError(ibin)
                    #Xvalues.append( dictHistos[ivar]['gen'].GetBinCenter(ibin)  )

                    tmpHisto = ROOT.TH1F('tmp'+str(ibin), 'tmp'+str(ibin), 100, -5, 5)
                    combineTree.Draw("r_bin"+str(ibin)+">>tmp"+str(ibin))
                    combineCont = tmpHisto.GetMean()
                    combineErrUp = combineTree.GetMaximum("r_bin"+str(ibin)) - combineCont
                    combineErrDown = combineCont - combineTree.GetMinimum("r_bin"+str(ibin))
                    #print 'r_bin'+str(ibin), combineCont
                    #Yvalues.append( histoCont*combineCont )
                    otherHistos['MLU'].SetBinContent( ibin, histoCont*combineCont )

                    #YvaluesCombine.append( combineCont )
                    histError = 0 if histoCont==0 else ROOT.TMath.Power( histoErr/histoCont, 2 )
                    combineErrorUp = 0 if combineCont==0 else ROOT.TMath.Power( combineErrUp/combineCont, 2 )
                    combineErrorDown = 0 if combineCont==0 else ROOT.TMath.Power( combineErrDown/combineCont, 2 )
                    otherHistos['MLU'].SetBinError( ibin, (abs( histoCont*combineCont ) * ROOT.TMath.Sqrt( histError + combineErrorUp  ) ) )
                    #YMaxError.append( abs( histoCont*combineCont ) * ROOT.TMath.Sqrt( histError + combineErrorUp  ) )
                    #YMinError.append( abs( histoCont*combineCont ) * ROOT.TMath.Sqrt( histError + combineErrorDown  ) )

                #UnfoldGraph = ROOT.TGraphAsymmErrors( len(Xvalues), array( 'd', Xvalues), array( 'd', Yvalues), array( 'd', [0]*len(Xvalues)), array( 'd', [0]*len(Xvalues)), array( 'd', YMinError), array( 'd', YMaxError) )
                dictHistos[ivar]['MLU'] = otherHistos['MLU']

            #### plotting for each variable
            drawUnfold( ivar, mainHisto, mainHistoLabel, otherHistos, nSubVariables[ivar], uncertDictHistos[ivar], outputLabel=outputLabel, outputDir=outputDir )

    ##### making combine plots
    if not args.only: combinePlots( name+'_'+args.selection, dictHistos, numBinsList, mainHistoLabel, otherHisto, otherHistoLabel, outputLabel, outputDir=outputDir, axisX='' )

##########################################################################
def drawUnfold( ivar, mainHisto, mainHistoLabel, otherHistos, varInfo, uncertDictHistos, outputLabel, outputDir ):
    """docstring for drawUnfold"""

    outputDir=outputDir+'/'+ivar+'/'+args.process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+'_'+args.selection+'_'+outputLabel+'_TUnfold_'+args.version+'.'+args.ext


    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 750, 750 )
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    if varInfo['alignLeg'].startswith('right'): legend=ROOT.TLegend(0.60,0.65,0.90,0.88)
    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

    dataIntegral = mainHisto.Integral()
    mainHisto.Scale(1, 'width')  ### divide by bin width
    mainHisto.SetLineWidth(2)
    #mainHisto.SetLineStyle(2)
    mainHisto.SetLineColor(ROOT.kMagenta)
    legend.AddEntry( mainHisto, mainHistoLabel, 'l' )

    tmpMax = [ mainHisto.GetMaximum() ]
    dummy = 0
    for ilabel, ih in otherHistos.items():
        #ih.Scale(1, 'width')
        #ih.Scale(scaleFactor)
        ih.Scale(dataIntegral/ih.Integral(), 'width')  ### divide by bin width
        ih.SetLineWidth(2)
        ih.SetLineColor(colors[dummy])
        tmpMax.append( ih.GetMaximum() )
        legend.AddEntry( ih, ilabel, 'l' )
        dummy= dummy+1


    mainHisto.GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+varInfo['label'].split('#')[1]+'}' )
    #mainHisto.GetYaxis().SetTitleOffset(0.95)
    mainHisto.GetYaxis().SetTitleSize(0.05)
    mainHisto.SetMaximum( 1.2*max( tmpMax )  )

    mainHisto.Draw( "histe")
    for ih in otherHistos: otherHistos[ih].Draw("histe same")
    mainHisto.Draw( "histe same")

    legend.Draw()
    CMS_lumi.relPosX = 0.12
    CMS_lumi.lumiTextSize = 0.5
    CMS_lumi.CMS_lumi(pad1, 4, 0)

    pad2.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame( 0, 0., varInfo['bins'][-1], 1.9 )
    tmpPad2.GetXaxis().SetTitle( varInfo['label'] )
    tmpPad2.GetYaxis().SetTitle( "Gen/Unfolded" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    hRatioUp = ROOT.TGraphAsymmErrors()
    hRatioUp.Divide( mainHisto, otherHistos[next(iter(otherHistos))], 'pois' )
    #hRatioUp.SetLineColor(kRed-4)
    #hRatioUp.SetLineWidth(2)
    hRatioUp.SetMarkerStyle(8)
    hRatioUp.Draw('P0')
    #hRatioDown.Draw('P same')
    if args.runMLU:
        hRatioMLU = ROOT.TGraphAsymmErrors()
        hRatioMLU.Divide( otherHistos['MLU'], otherHistos[next(iter(otherHistos))], 'pois' )
        hRatioMLU.SetLineColor(8)
        #hRatioMLU.SetLineWidth(2)
        hRatioMLU.SetMarkerStyle(8)
        hRatioMLU.Draw('P same')

    uncertDictHistos['statUnc'].SetFillColor(ROOT.kBlue)
    uncertDictHistos['statUnc'].SetFillStyle(3004)
    uncertDictHistos['statUnc'].Draw('E2 same')
    uncertDictHistos['totalUnc'].SetFillColor(ROOT.kBlack)
    uncertDictHistos['totalUnc'].SetFillStyle(3004)
    uncertDictHistos['totalUnc'].Draw('E2 same')

    can.SaveAs(outputName)
    if args.ext.startswith('pdf'):
        can.SaveAs( outputName.replace('pdf', 'png') )
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)


##################################################
def combinePlots( name, dictHistos, numBins, mainHistoLabel, otherHisto, otherHistoLabel, outputLabel, axisX='', outputDir='Plots/'):
    """docstring for combinePlots"""

    outputFileName = name+'_'+outputLabel+'_combinePlots_'+args.version+'.'+args.ext
    if args.log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.10,0.80,0.60,0.90)
    legend.SetNColumns(3 if args.runMLU else 2)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.06)

    dictHistos['combData'] = ROOT.TH1F('combData', 'combData', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combData' ], mainHistoLabel, 'lep' )
    dictHistos['combUnfold'] = ROOT.TH1F('combUnfold', 'combUnfold', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combUnfold' ], otherHistoLabel, 'lep' )
    if args.runMLU:
        dictHistos['combMLU'] = dictHistos['combUnfold'].Clone()
        legend.AddEntry( dictHistos[ 'combMLU' ], 'MLU', 'lep' )

    tmpNbin = 0
    Xlabels = []
    for ivar,ih in dictHistos.items():
        if ivar.startswith('Jet') and not ivar.endswith(('21', '32')):
            Xlabels.append( '#'+nSubVariables[ivar]['label'].split('#')[1] )
            for ibin in range(1, ih['data'].GetNbinsX()+1):
                tmpNbin = tmpNbin+1
                dictHistos['combData'].SetBinContent( tmpNbin, ih[otherHisto].GetBinContent(ibin) )
                dictHistos['combData'].SetBinError( tmpNbin, ih[otherHisto].GetBinError(ibin) )
                dictHistos['combUnfold'].SetBinContent( tmpNbin, ih['unfold'].GetBinContent(ibin) )
                dictHistos['combUnfold'].SetBinError( tmpNbin, ih['unfold'].GetBinError(ibin) )
                if args.runMLU:
                    dictHistos['combMLU'].SetBinContent( tmpNbin, ih['MLU'].GetBinContent(ibin) )
                    dictHistos['combMLU'].SetBinError( tmpNbin, ih['MLU'].GetBinError(ibin) )


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
    if args.runMLU:
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
    if args.runMLU:
        hRatioMLU = ROOT.TGraphAsymmErrors()
        hRatioMLU.Divide( dictHistos['combMLU'], dictHistos[ 'combData' ], 'pois' )
        hRatioMLU.SetMarkerStyle(8)
        hRatioMLU.SetMarkerColor( 8 )
        hRatioMLU.Draw('P0 same')
        hRatio.Draw('P0 same')

    for i in lines: lines[i].Draw('same')

    textBox=ROOT.TLatex()
    #textBox.SetNDC()
    textBox.SetTextSize(0.10)
    textBox.SetTextAlign(12)
    textBoxList = {}
    for i in range(1, len(numBins)):
        textBoxList[i] = textBox.Clone()
        textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., 0., Xlabels[i-1] )

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if args.ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]
    ROOT.gStyle.SetPadTickX(1)

##########################################################
def plotUnfoldCombined( ivar, combineFile, dataHisto, mainHisto, mainHistoLabel, otherHistos, varInfo, uncertDictHistos, outputLabel, outputDir ):
    """"Very specific, takes rootfile from combine and compute unfolding"""

    outputDir=outputDir+'/'+ivar+'/'+args.process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+'_'+args.selection+'_'+outputLabel+'_MLU_'+args.version+'.'+args.ext

    Xvalues = []
    Yvalues = []
    YvaluesCombine = []
    YMaxError = []
    YMinError = []
    combineTree = combineFile.Get('limit')
    for ibin in range(1, dataHisto.GetNbinsX()+1):
        histoCont = dataHisto.GetBinContent(ibin)
        histoErr = dataHisto.GetBinError(ibin)
        Xvalues.append( dataHisto.GetBinCenter(ibin)  )

        tmpHisto = ROOT.TH1F('tmp'+str(ibin), 'tmp'+str(ibin), 100, -5, 5)
        combineTree.Draw("r_bin"+str(ibin)+">>tmp"+str(ibin))
        combineCont = tmpHisto.GetMean()
        combineErrUp = combineTree.GetMaximum("r_bin"+str(ibin)) - combineCont
        combineErrDown = combineCont - combineTree.GetMinimum("r_bin"+str(ibin))

        #print 'r_bin'+str(ibin), combineCont
        Yvalues.append( histoCont*combineCont )
        YvaluesCombine.append( combineCont )
        histError = 0 if histoCont==0 else ROOT.TMath.Power( histoErr/histoCont, 2 )
        combineErrorUp = 0 if combineCont==0 else ROOT.TMath.Power( combineErrUp/combineCont, 2 )
        combineErrorDown = 0 if combineCont==0 else ROOT.TMath.Power( combineErrDown/combineCont, 2 )
        YMaxError.append( abs( histoCont*combineCont ) * ROOT.TMath.Sqrt( histError + combineErrorUp  ) )
        YMinError.append( abs( histoCont*combineCont ) * ROOT.TMath.Sqrt( histError + combineErrorDown  ) )

    UnfoldGraph = ROOT.TGraphAsymmErrors( len(Xvalues), array( 'd', Xvalues), array( 'd', Yvalues), array( 'd', [0]*len(Xvalues)), array( 'd', [0]*len(Xvalues)), array( 'd', YMinError), array( 'd', YMaxError) )
    UnfoldGraph.SetMarkerStyle(8)
    UnfoldGraph.GetXaxis().SetTitle( mainHisto.GetXaxis().GetTitle() )
    UnfoldGraph.GetYaxis().SetTitle( 'Events / '+str(mainHisto.GetBinWidth(1)) )
    UnfoldGraph.GetYaxis().SetTitleOffset( 0.8 )

    if varInfo['alignLeg'].startswith('right'): legend=ROOT.TLegend(0.60,0.65,0.90,0.88)
    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.AddEntry( UnfoldGraph, 'maxLikelihood Unfolding', 'pe' )
    #legend.AddEntry( genMCHisto, 'Gen-level MC', 'l' )
    #legend.AddEntry( recoMCHisto, ('Data' if args.process.startswith('data') else 'Reco-level MC (Closure)'), 'l' )


    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can', 'can',  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1", "Main",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2", "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    tmpPad1= pad1.DrawFrame( 0, 0, 1, 1.5*dataHisto.GetMaximum() )
    pad1.Modified()
    UnfoldGraph.Draw("P0")
    #recoMCHisto.Draw("histe same")
    dataHisto.Draw("histe same")
    legend.Draw()
#    if args.process.startswith('data'):
#        CMS_lumi.extraText = "Preliminary"
#        CMS_lumi.lumi_13TeV = str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, 2016"
#    else:
#        CMS_lumi.extraText = "Simulation Preliminary"
#        CMS_lumi.lumi_13TeV = "13 TeV, 2016"
    CMS_lumi.relPosX = 0.13
    CMS_lumi.CMS_lumi(pad1, 4, 0)

    pad2.cd()
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame( 0, 0.5, 1, 1.5 )
    tmpPad2.GetXaxis().SetTitle( mainHisto.GetXaxis().GetTitle() )
    tmpPad2.GetYaxis().SetTitle( "True/Unfolded" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    hRatioUp = ROOT.TGraphAsymmErrors( len(Xvalues), array( 'd', Xvalues), array( 'd', YvaluesCombine), array( 'd', [0]*len(Xvalues)), array( 'd', [0]*len(Xvalues)), array( 'd', YMinError), array( 'd', YMaxError) )
    #hRatioUp.SetLineColor(kRed-4)
    #hRatioUp.SetLineWidth(2)
    hRatioUp.SetMarkerStyle(8)
    hRatioUp.Draw('P0')
    #hRatioDown.Draw('P same')

    can.SaveAs(outputName)
    if args.ext.startswith('pdf'):
        can.SaveAs( outputName.replace('pdf', 'png') )
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='MCClosure', dest='process', help='Process to draw, example: data, MCClosure, MCSelfClosure.' )
    parser.add_argument('-s', '--selection', action='store', default='dijetSel', help='Selection: dijetSel, WSel, topSel.' )
    parser.add_argument('-v', '--version', action='store', default='v02', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=0., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )
    parser.add_argument("--main", action='store', dest="main", choices=['Ptbin', 'HTbin', 'herwig'], default='Ptbin', help="For dijet sel, main signal QCD" )
    parser.add_argument("--alt", action='store', dest="alt", choices=['Ptbin', 'HTbin', 'herwig'], default='herwig', help="For dijet sel, alternative signal QCD" )
    parser.add_argument("--inputFolder", action='store', dest="inputFolder", default="", help="input folder" )
    parser.add_argument("--outputFolder", action='store', dest="outputFolder", default="", help="Output folder" )
    parser.add_argument("-r", "--runMLU", action='store_true', default=False, help="Run max likelihood unfold (combine)" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    for iy in ( ['2017', '2018'] if args.year.startswith('all') else [ args.year ] ):
        args.lumi = args.lumi + checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[iy]['lumi']

    #### Define main and alt signals
    if args.selection.startswith('dijet'):
        if args.main.startswith('Ptbin'):
            sigPlotLabel = 'Pythia8'
            signalLabelBegin = 'QCD_Pt_'
            signalLabel = 'QCD_Pt_3200toInf'
        elif args.main.startswith('HTbin'):
            sigPlotLabel = 'MG5+Pythia8'
            signalLabelBegin = 'QCD_HT'
            signalLabel = 'QCD_HT2000toInf'
        elif args.main.startswith('herwig'):
            sigPlotLabel = 'Herwig7'
            signalLabelBegin = 'QCD_Pt-'
            signalLabel = 'QCD_Pt-150to3000'
        if args.alt.startswith('Ptbin'):
            altSigPlotLabel = 'Pythia8'
            altSignalLabelBegin = 'QCD_Pt_'
            altSignalLabel = 'QCD_Pt_3200toInf'
        elif args.alt.startswith('HTbin'):
            altSigPlotLabel = 'MG5+Pythia8'
            altSignalLabelBegin = 'QCD_HT'
            altSignalLabel = 'QCD_HT2000toInf'
        elif args.alt.startswith('herwig'):
            altSigPlotLabel = 'Herwig7'
            altSignalLabelBegin = 'QCD_Pt-'
            altSignalLabel = 'QCD_Pt-150to3000'

    dictFiles = OrderedDict()
    for ivar, varInfo in nSubVariables.items():
        if args.selection.startswith('dijet') and ivar.startswith(('Jet1', 'Jet2')):
            if args.only and not ivar.startswith( args.only ): continue
            dictFiles[ivar] = ROOT.TFile.Open( args.inputFolder+'Plots/'+args.selection+'/Unfold/'+args.year+'/'+ivar+'/'+args.process+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root' )
            if args.runMLU:
                dictFiles[ivar+'Combine'] = ROOT.TFile.Open( 'Plots/'+args.selection+'/Unfold/'+args.year+'/'+ivar+'/'+args.process+'/higgsCombine'+ivar+'_from'+signalLabel+'.MultiDimFit.mH120.root' )

    outputDir = args.outputFolder+'Plots/'+args.selection+'/Unfold/'+args.year
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    CMS_lumi.extraText = 'Preliminary' if args.process.startswith('data') else "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+( '2017+2018' if args.year.startswith('all') else args.year )

    plotList = [ 'recoJet1', 'recoJet2' ] if args.selection.startswith('dijet') else [ 'recoJet' ]

    for i in plotList:
        makePlots( i, outputDir=outputDir)

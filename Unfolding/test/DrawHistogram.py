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

colors = [ 2, 4, 6, 8, 9, 28, 30, 42, 13, 46 ]
canvas = {}

##########################################################
def setSelection( listSel, xMin=0.65, yMax=0.65, align='right' ):

    for i in range( len( listSel ) ):
        textBox=TLatex()
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        if 'right' in align: textBox.SetTextAlign(31)
        textBox.SetTextFont(62) ### 62 is bold, 42 is normal
        textBox.DrawLatex(xMin, yMax, listSel[i])
        yMax = yMax -0.05


#####################################################
#####################################################
def plotSimpleComparison( inFile1, sample, inFile2, sample2, name, rebinX=1, xmin='', xmax='', labX=0.92, labY=0.50, axisX='', axisY='', log=False, ext='png', Norm=False, version='', outputDir='Plots/' ):
    """"Take two root files, make simple comparison plot"""

    outputFileName = name+'_'+sample+sample2+'_simpleComparisonPlot'+version+'.'+ext
    print('Processing.......', outputFileName)

    if isinstance( inFile1, ROOT.TTree ):
        histo = inFile1.Get( 'jetObservables/'+name )
        if rebinX!=1: histo.Rebin( rebinX )
        histo2 = inFile2.Get( 'jetObservables/'+name )
        if rebinX!=1: histo2.Rebin( rebinX )
    else:  ##inFile1 is a histogram
        histo = inFile1
        histo2 = inFile2

    binWidth = histo.GetBinWidth(1)

    legend=ROOT.TLegend(0.60,0.75,0.90,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    #histo.SetFillColor(48)
    histo.SetFillStyle(1001)

    #tdrStyle.SetPadRightMargin(0.05)
    canvas[name] = ROOT.TCanvas('c1'+name, 'c1'+name,  10, 10, 750, 500 )
    if log:
        canvas[name].SetLogy()
        outName = outputFileName.replace('_simplePlot','_Log_simplePlot')
    else: outName = outputFileName

    legend.AddEntry( histo, sample, 'f' )
    legend.AddEntry( histo2, sample2, 'f' )
    if xmax and xmin: histo.GetXaxis().SetRangeUser( xmin, xmax )
    histo.GetYaxis().SetTitleOffset(0.90)
    histo.SetLineColor(ROOT.kRed)
    histo2.SetLineColor(ROOT.kBlue)
    if Norm:
        histo.DrawNormalized('hist')
        histo2.DrawNormalized('hist same')
    else:
        histo.Draw('histe')
        histo2.Draw('hist same')
    if not axisY: histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: histo.GetXaxis().SetTitle( axisX )

    legend.Draw()

    canvas[name].SaveAs( outputDir+outName )
    if ext.startswith('pdf'):
        canvas[name].SaveAs( outputDir+outName.replace('pdf', 'png') )
    #del can

#####################################################
#####################################################
def plotData( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False,
                      addRatioFit=False, Norm=False, ext='png', outputDir='Plots/', legendAlignment='right' ):
    """"For dijets, make plots only for data"""

    outputFileName = name+'_dataPlots_'+args.version+'.'+args.ext
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.60,0.65,0.90,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    tmpdataHistos = OrderedDict()
    dummy = 1
    triggerList = [ 'AK8PFJet80', 'AK8PFJet140', 'AK8PFJet200', 'AK8PFJet260', 'AK8PFJet320', 'AK8PFJet400', 'AK8PFJet450', 'AK8PFJet500' ]
    for dataSamples in dataFile.keys():
        yearLabel = dataSamples.split('data')[1]
        if args.selection.startswith('dijet'):
            for it in list(reversed(triggerList)):
            #for it in checkDict( 'JetHT', dictSamples )[yearLabel]['triggerList']:
                tmpdataHistos[ it+dataSamples ] = dataFile[dataSamples].Get( 'jetObservables/'+name.replace( args.selection, it+'_'+args.selection ) )
                tmpdataHistos[ it+dataSamples ].Scale( checkDict( 'JetHT', dictSamples )[yearLabel]['triggerList'][it] )
                if isinstance(rebinX, int): tmpdataHistos[ it+dataSamples ] = tmpdataHistos[ it+dataSamples ].Rebin( rebinX )
                if name.startswith('nPV'): tmpdataHistos[ it+dataSamples ].Scale( 1./tmpdataHistos[ it+dataSamples ].Integral() )
                tmpdataHistos[ it+dataSamples ].SetLineColor( colors[dummy] )
                tmpdataHistos[ it+dataSamples ].SetLineWidth( 2 )
                legend.AddEntry( tmpdataHistos[ it+dataSamples ], it, 'lep' )
                dummy = dummy+1
        else:
            tmpdataHistos[ dataSamples ] = dataFile[yearLabel].Get( 'jetObservables/'+name )

    stackHisto = ROOT.THStack('stackHisto'+name, 'stack'+name)
    for i in tmpdataHistos: stackHisto.Add( tmpdataHistos[i] )
    binWidth = tmpdataHistos[next(iter(tmpdataHistos))].GetBinWidth(1)


    #tdrStyle.SetPadRightMargin(0.05)
    canvas[name] = ROOT.TCanvas('c1'+name, 'c1'+name,  10, 10, 750, 500 )
    if log:
        canvas[name].SetLogy()
        outName = outputFileName.replace('Plots','Plots_Log')
    else: outName = outputFileName
    if name.startswith('nPV'): stackHisto.Draw('histe nostack')
    else: stackHisto.Draw('histe')
    stackHisto.GetYaxis().SetTitleOffset(0.90)
    if log: stackHisto.SetMinimum(1)
    if xmax and xmin: stackHisto.GetXaxis().SetRangeUser( xmin, xmax )
    if not axisY: stackHisto.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: stackHisto.GetXaxis().SetTitle( axisX )
    legend.Draw()

    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.CMS_lumi(canvas[name], 4, 0)

    canvas[name].SaveAs( outputDir+'/'+outName )
    if ext.startswith('pdf'):
        canvas[name].SaveAs( outputDir+'/'+outName.replace('pdf', 'png') )
    #del can


######################################################
def plotSysComparison( nomHisto, dictUncHistos, outputName, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ):
    """docstring for plot"""

    outputFileName = outputName+'_'+version+'.'+ext
    print 'Processing.......', outputFileName

    binWidth = nomHisto.GetBinWidth(1)

    legend=ROOT.TLegend(0.70,0.65,0.90,0.87)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    multiGraph = ROOT.TMultiGraph()
    gnom = ROOT.TGraphAsymmErrors()
    gnom.Divide( nomHisto, nomHisto, 'pois' )
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetLineWidth(2)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    multiGraph.Add( gnom )

    dictgraph = {}
    dummy=0
    for ih in dictUncHistos:
        print(ih)
        dictgraph[ih] = ROOT.TGraphAsymmErrors()
        dictgraph[ih].Divide( dictUncHistos[ih], nomHisto, 'pois' )
        dictgraph[ih].SetLineColor( colors[dummy] )
        dictgraph[ih].SetLineStyle( 2 )
        dictgraph[ih].SetLineWidth( 2 )
        legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        multiGraph.Add( dictgraph[ih] )
        dummy=dummy+1

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 750, 500 )
    multiGraph.GetYaxis().SetTitle( 'Ratio Unc/Nominal' )
    multiGraph.GetXaxis().SetTitle( labelX )
    multiGraph.SetMaximum( 1.5 )
    multiGraph.SetMinimum( 0.5 )
    multiGraph.Draw('ALP')

    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+year
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    del canUnc

##################################################

##################################################
def plotSignalBkg( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False, Norm=False,
                      addRatioFit=False, ext='png', outputDir='Plots/', legendAlignment='right' ):
    """function to plot s and b histos"""

    outputFileName = name+'_AnalysisPlots_'+args.version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    if legendAlignment.startswith('right'): legend=ROOT.TLegend(0.60,(0.68 if args.selection.startswith('dijet') else 0.65),0.90,0.88)
    else: legend=ROOT.TLegend(0.17,(0.68 if args.selection.startswith('dijet') else 0.65),0.50,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    dataHistos = {}
    tmpdataHistos = {}
    for dataSamples in dataFile.keys():
        yearLabel = dataSamples.split('data')[1]
        if args.selection.startswith('dijet'):
            for it in checkDict( 'JetHT', dictSamples )[yearLabel]['triggerList']:
                tmpdataHistos[ it+dataSamples ] = dataFile[dataSamples].Get( 'jetObservables/'+name.replace( args.selection, it+'_'+args.selection ) )
                tmpdataHistos[ it+dataSamples ].Scale( checkDict( 'JetHT', dictSamples )[yearLabel]['triggerList'][it] )
        else:
            tmpdataHistos[ dataSamples ] = dataFile[yearLabel].Get( 'jetObservables/'+name )
    dataHistos[ 'DATA' ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
    dataHistos[ 'DATA' ].Reset()
    for i in tmpdataHistos: dataHistos['DATA'].Add( tmpdataHistos[i] )
    if isinstance(rebinX, int): dataHistos[ 'DATA' ] = dataHistos[ 'DATA' ].Rebin( rebinX )
    else:
        tmpDataHist = dataHistos[ 'DATA' ].Clone()
        dataHistos['DATA'] = tmpDataHist.Rebin( len(rebinX)-1, tmpDataHist.GetName()+"_rebinX", array( 'd', rebinX ) )

    legend.AddEntry( dataHistos[ 'DATA' ], 'Data', 'lep' )

    bkgHistos = OrderedDict()
    tmpbkgHistos = OrderedDict()
    if len(bkgFiles) > 0:
        for bkgSamples in bkgFiles:
            yearLabel = ''.join(bkgSamples[-4:])
            print(name, bkgSamples)
            bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( 'jetObservables/'+name )
            bkgHistos[ bkgSamples ].SetTitle(bkgSamples)
            bkgHistos[ bkgSamples ].Scale( args.lumi*bkgFiles[ bkgSamples ][1]['XS'] / bkgFiles[ bkgSamples ][1][yearLabel]['nGenWeights'] )
            if isinstance(rebinX, int): bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( rebinX )
            else:
                tmpBkgHist = bkgHistos[ bkgSamples ].Clone()
                bkgHistos[bkgSamples] = tmpBkgHist.Rebin( len(rebinX)-1, tmpBkgHist.GetName()+"_rebinX", rebinX )
            if bkgSamples.startswith(('ST', 'W', 'Diboson')):
                bkgHistos[ bkgSamples ].SetFillStyle( 1001 )
                bkgHistos[ bkgSamples ].SetFillColor( bkgFiles[ bkgSamples ][1]['color'] )
            else:
                bkgHistos[ bkgSamples ].SetFillColor( 0 )
            bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][1]['color'] )
            bkgHistos[ bkgSamples ].SetLineWidth( 2 )

    #### Merging samples
    yearLabel = '2018' if args.year.startswith('all') else args.year
    for bkg in bkgHistos:
        if args.selection.startswith('dijet'):
            if bkg.startswith('QCD_Pt_') and not bkg.endswith('Inf'+yearLabel):
                bkgHistos['QCD_Pt_3200toInf'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            elif bkg.startswith('QCD_HT') and not bkg.endswith('Inf'+yearLabel):
                bkgHistos['QCD_HT2000toInf'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            elif args.year.startswith('all') and bkg.startswith('QCD_Pt-150to30002017'):
                bkgHistos['QCD_Pt-150to30002018'].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            else:
                legend.AddEntry( bkgHistos[ bkg ], bkgFiles[ bkg ][1]['label'], 'le' )
        else:
            if bkg.startswith(('WZ','ZZ')):
                bkgHistos['WW'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            elif bkg.startswith('ST_t'):
                bkgHistos['ST_s-channel_4f_leptonDecays'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            elif bkg.startswith('TTTo2L2Nu'):
                bkgHistos['TTToSemiLeptonic'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            else:
                legend.AddEntry( bkgHistos[ bkg ], bkgFiles[ bkg ][1]['label'], ('l' if bkg.startswith('TT') else 'f' ) )

    hBkg = bkgHistos[next(iter(bkgHistos))].Clone()
    hBkg.Reset()
    ratioDict = {}
    binWidth = dataHistos['DATA'].GetBinWidth(1)

    if args.selection.startswith('dijet'): dataHistos['DATA'].Scale( 1/dataHistos['DATA'].Integral(), 'width' )
    stackHisto = ROOT.THStack('stackHisto'+name, 'stack'+name)
    bkgHistos = OrderedDict(reversed(bkgHistos.items()))
    for samples in bkgHistos:
        if not args.selection.startswith('dijetSel'):
            if samples.startswith('TTJets'): stackHisto.Add( bkgHistos[ samples ].Clone() )
            elif samples.startswith('TTTo'): hBkg.Add( bkgHistos[ samples ].Clone() )
            else:
                stackHisto.Add( bkgHistos[ samples ].Clone() )
                hBkg.Add( bkgHistos[ samples ].Clone() )
        else:
            bkgHistos[samples].Scale( dataHistos['DATA'].Integral()/bkgHistos[samples].Integral() )
            bkgHistos[samples].Scale( 1/bkgHistos[samples].Integral(), 'width' )
            stackHisto.Add( bkgHistos[samples].Clone() )
            ratioDict[samples] = ROOT.TGraphAsymmErrors()
            ratioDict[samples].Divide( dataHistos[ 'DATA' ].Clone(), bkgHistos[samples].Clone(), 'pois' )
            ratioDict[samples].SetLineColor( bkgFiles[samples][1]['color']  )
            ratioDict[samples].SetLineWidth(2)

    canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name,  10, 10, 750, 750 )
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    pad1 = ROOT.TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    #if log and not args.final: pad1.SetLogy()
    if log: pad1.SetLogy()
    if 'tau' in axisX:
        dataHistos['DATA'].GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+axisX.split('#')[1]+'}' )
    else:
        dataHistos['DATA'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    dataHistos[ 'DATA' ].SetMarkerStyle(8)
    dataHistos['DATA'].GetYaxis().SetTitleSize( 0.05 )
    dataHistos['DATA'].GetYaxis().SetTitleOffset( 1.2 )
    dataHistos['DATA'].GetXaxis().SetTitle( axisX )
    if xmax: dataHistos['DATA'].GetXaxis().SetRangeUser( xmin, xmax )
    dataHistos['DATA'].SetMaximum( dataHistos['DATA'].GetMaximum()*1.2 )
    #dataHistos['DATA'].SetMinimum( 2. )

    #stackHisto.SetMinimum( 0.1 )
    if not args.selection.startswith('dijet'):
        hBkg.SetLineColor( bkgFiles['TTToSemiLeptonic'][1]['color'] )
        hBkg.SetFillColor( 0 )
        hBkg.SetLineWidth(2)

    dataHistos[ 'DATA' ].Draw('E')
    if args.selection.startswith('dijet'): stackHisto.Draw('histe same nostack')
    else:
        hBkg.Draw("histe same")
        stackHisto.Draw('histe same')
    dataHistos[ 'DATA' ].Draw('E same')

    CMS_lumi.CMS_lumi( pad1, 4, 0)
    legend.Draw()

    pad2.cd()
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)

    tmpPad2= pad2.DrawFrame(xmin,0.5,xmax,1.5)
    tmpPad2.GetYaxis().SetTitle( "Data/Sim." )
    tmpPad2.GetXaxis().SetTitle(  axisX )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    if args.selection.startswith('dijet'):
        dummy2=0
        for ih in ratioDict:
            if dummy2==0: ratioDict[ih].Draw('P')
            else: ratioDict[ih].Draw('P same')
    else:
        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dataHistos[ 'DATA' ], hBkg, 'pois' )
        hRatio.SetLineColor( bkgFiles['TTToSemiLeptonic'][1]['color']  )
        hRatio.SetLineWidth(2)
        hBkg2 = stackHisto.GetHistogram()
        hRatio2 = ROOT.TGraphAsymmErrors()
        hRatio2.Divide( dataHistos[ 'DATA' ], hBkg2, 'pois' )
        hRatio2.SetLineColor( bkgFiles['TTJets'][1]['color']  )
        hRatio2.SetLineWidth(2)
        hRatio.Draw('P')
        hRatio2.Draw('P same')

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]

##################################################
def plotResolution( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False,
                      addRatioFit=False, Norm=False, ext='png', outputDir='Plots/', legendAlignment='right' ):
    """function to resolution"""

    outputFileName = name+'_Resolution_'+args.version+'.'+ext
    print('Processing.......', outputFileName)

    if legendAlignment.startswith('right'): legend=ROOT.TLegend(0.60,0.75,0.90,0.90)
    else: legend=ROOT.TLegend(0.20,0.75,0.50,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    bkgHistos = OrderedDict()
    maxList = []
    if len(bkgFiles) > 0:
        for bkgSamples in bkgFiles:
            yearLabel = ''.join(bkgSamples[-4:])
            bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( 'jetObservables/'+name )
            bkgHistos[ bkgSamples ].SetTitle(bkgSamples)
            bkgHistos[ bkgSamples ].Scale( args.lumi*bkgFiles[ bkgSamples ][1]['XS'] / bkgFiles[ bkgSamples ][1][yearLabel]['nGenWeights'] )
            print(bkgSamples, round(bkgHistos[ bkgSamples ].Integral(), 2) )
            if isinstance(rebinX, int): bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( rebinX )
            else:
                tmpBkgHist = bkgHistos[ bkgSamples ].Clone()
                bkgHistos[bkgSamples] = tmpBkgHist.Rebin( len(rebinX)-1, tmpBkgHist.GetName()+"_rebinX", rebinX )
            bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][1]['color'] )
            bkgHistos[ bkgSamples ].SetLineWidth( 2 )

    yearLabel = '2018' if args.year.startswith('all') else args.year
    for bkg in bkgFiles:
        if args.selection.startswith('dijet'):
            if bkg.startswith('QCD_Pt_') and not bkg.endswith('Inf'+yearLabel):
                bkgHistos['QCD_Pt_3200toInf'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            elif bkg.startswith('QCD_HT') and not bkg.endswith('Inf'+yearLabel):
                bkgHistos['QCD_HT2000toInf'+yearLabel].Add( bkgHistos[bkg] )
                bkgHistos.pop(bkg, None)
            else:
                legend.AddEntry( bkgHistos[ bkg ], bkgFiles[ bkg ][1]['label'], 'le' ) # if Norm else 'f' )

    stackHisto = ROOT.THStack('stackHisto'+name, 'stack'+name)
    for samples in bkgHistos:
        print(samples)
        bkgHistos[ samples ].Scale( 1/bkgHistos[ samples ].Integral() )
        stackHisto.Add( bkgHistos[ samples ].Clone() )
        binWidth = bkgHistos[ samples ].GetBinWidth(1)

    canvas[name] = ROOT.TCanvas('c1'+name, 'c1'+name,  10, 10, 750, 500 )
    if log:
        canvas[name].SetLogy()
        outName = outputFileName.replace('_simplePlot','_Log_simplePlot')
    else: outName = outputFileName

    if xmax and xmin: stackHisto.GetXaxis().SetRangeUser( xmin, xmax )
    #stackHisto.SetLineColor(ROOT.kRed)
    #stackHisto2.SetLineColor(ROOT.kBlue)
    stackHisto.Draw('nostack')
    if not axisY: stackHisto.GetYaxis().SetTitle( 'Normalized / '+str(binWidth) )
    stackHisto.GetYaxis().SetTitleOffset(0.90)
    if axisX: stackHisto.GetXaxis().SetTitle( axisX )

    CMS_lumi.lumi_13TeV = "13 TeV, "+args.year
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canvas[name], 4, 0)
    legend.Draw()

    canvas[name].SaveAs( outputDir+'/'+outName )
    if ext.startswith('pdf'):
        canvas[name].SaveAs( outputDir+'/'+outName.replace('pdf', 'png') )
    #del can

####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='bkgData', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
    parser.add_argument('-s', '--selection', action='store', default='dijetSel', help='Selection: dijetSel, WSel, topSel.' )
    parser.add_argument('-v', '--version', action='store', default='v0', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=0., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )
    parser.add_argument("--outputFolder", action='store', dest="outputFolder", default="", help="Output folder" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    args.inputFolder = os.environ['CMSSW_BASE']+'/src/jetObservables/Unfolding/test/Rootfiles/'
    dictSamples = OrderedDict(dictSamples)
    VER = args.version.split('_')[1] if '_' in args.version else args.version
    dataFile = {}
    bkgFiles = OrderedDict()
    for iy in ( ['2017', '2018'] if args.year.startswith('all') else [ args.year ] ):
        dataFile['data'+iy] = ROOT.TFile.Open(args.inputFolder+checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[iy]['skimmerHisto'])
        args.lumi = args.lumi + checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[iy]['lumi']

        for isam in dictSamples:
            if not checkDict( isam, dictSamples )[iy]['skimmerHisto'].endswith('root'): continue
            if isam.startswith(('JetHT', 'SingleMuon')): continue
            if args.selection.startswith('dijet') and not checkDict( isam, dictSamples )['selection'].startswith('dijet'): continue
            #if isam.startswith('QCD_HT'): continue ###### temporary
            bkgFiles[isam.split('_Tune')[0]+iy] = [
                                ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[iy]['skimmerHisto'] ),
                                checkDict( isam, dictSamples )
                            ]


    outputDir = args.outputFolder+'Plots/'+args.selection+'/'+('Resolution' if args.process.startswith('reso') else 'Basic')+'/'+args.year
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    CMS_lumi.extraText = "Preliminary" if args.process.startswith(('data', 'bkgData')) else "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if args.year.startswith('all') else ", "+args.year )

    plotList = [
            [ 'data', 'nPVs_only', 'Number of Primary Vertex', 0., 100., 2, 'right' ],
            [ 'data', 'LeadJetAK8_pt_only', 'Leading AK8 jet pt [GeV]', 100, 1500, 2, 'right' ] ,
            [ 'data', 'LeadJetAK8_eta_only', 'Leading AK8 jet eta [GeV]', 100, 1500, 2, 'right' ] ,
            [ 'data', 'recoJet1_sortedPt_pt_', 'Leading AK8 jet pt [GeV]', 100, 1000, 1, 'right' ] ,
            [ 'bkgData', 'nPVs', 'Number of Primary Vertex', 0., 100., 2, 'right' ],
            [ 'bkgData', 'recoPtAsym', 'pT Asymmetry', 0., .4, 2, 'right' ],
            [ 'bkgData', 'recoDeltaPhi', '#Delta #Phi(  j1, j2 )', 1., 3.5, 2, 'left' ],
            [ 'bkgData', 'recoJet1_sortedPt_pt', 'Leading AK8 jet pt [GeV]', 100, 1500, 5, 'right' ] ,
            [ 'bkgData', 'recoJet1_sortedPt_eta', 'Leading AK8 jet #eta', -3., 3., 2, 'right' ] ,
            [ 'bkgData', 'recoJet1_sortedPt_phi', 'Leading AK8 jet #Phi', -3.14, 3.14, 5, 'right' ] ,
            [ 'bkgData', 'recoJet1_sortedPt_mass', 'Leading AK8 jet mass [GeV]', 0, 200, 1, 'right' ] ,
            [ 'bkgData', 'recoJet2_sortedPt_pt', '2nd Leading AK8 jet pt [GeV]', 100, 1500, 5, 'right' ] ,
            [ 'bkgData', 'recoJet2_sortedPt_eta', '2nd Leading AK8 jet #eta', -3., 3., 2, 'right' ] ,
            [ 'bkgData', 'recoJet2_sortedPt_phi', '2nd Leading AK8 jet #Phi', -3.14, 3.14, 5, 'right' ] ,
            [ 'bkgData', 'recoJet2_sortedPt_mass', '2nd Leading AK8 jet mass [GeV]', 0, 200, 1, 'right' ] ,
            ]
    if not args.selection.startswith('dijet'):
            [ 'bkgData', 'recoJet_mass_nom', 'Leading AK8 jet softdrop mass [GeV]', ( 50. if args.selection.startswith('WSel') else 130. ), ( 140. if args.selection.startswith('WSel') else 300. ), 1, 'right' ] ,

    jetlabels = [ ('Jet1', 'Outer'), ('Jet2', 'Central') ] if args.selection.startswith('dijet') else [ ('Jet', 'Leading') ]
    for ijet in jetlabels:
        plotList.append( [ 'bkgData', 'reco'+ijet[0]+'_pt_nom', ijet[1]+' AK8 jet pt [GeV]', 100, 1500, 5, 'right' ] )
        plotList.append( [ 'bkgData', 'reco'+ijet[0]+'_eta_nom', ijet[1]+' AK8 jet eta', -3, 3, 2, 'left' ] )
        plotList.append( [ 'resol', 'resol'+ijet[0]+'_pt', ijet[1]+'AK8 jet reco/gen pt', 0., 2., 2, 'right' ] )
        plotList.append( [ 'resol', 'resol'+ijet[0]+'_sdmass', ijet[1]+'AK8 jet reco/gen sd mass', 0., 2., 2, 'right' ] )
    for ivar, varInfo in nSubVariables.items():
        if args.selection.startswith('dijet') and ivar.startswith('Jet_'): continue
        if not args.selection.startswith('dijet') and ivar.startswith(('Jet1_','Jet2_')): continue
        plotList.append( [ 'resol', 'resol'+ivar, varInfo['label']+' reco/gen', 0., 2., 2, 'right' ] )
        plotList.append( [ 'bkgData', 'reco'+ivar+'_nom', varInfo['label'], varInfo['bins'][0], varInfo['bins'][-1], varInfo['bins'], varInfo['alignLeg'] ] )

    if args.only: Plots = [ y[1:] for y in plotList if ( ( args.process in y[0] ) and ( args.only in y[1] ) )  ]
    else: Plots = [ x[1:] for x in plotList if ( ( args.process in x[0] ) )  ]
    if len(Plots)==0 :
        print('Variable not found. Have a nice day')
        sys.exit(0)

    for i in Plots:
        if ( 'bkgData' in args.process ):
            plotSignalBkg( i[0]+'_'+args.selection, i[2], i[3], i[4],
                            log=args.log, axisX=i[1], legendAlignment=i[5], outputDir=outputDir, ext=args.ext)
        elif ( 'resol' in args.process ):
            plotResolution( i[0]+'_'+args.selection, i[2], i[3], i[4],
                            log=args.log, axisX=i[1], Norm=args.norm, legendAlignment=i[5], outputDir=outputDir, ext=args.ext)
        elif ( 'data' in args.process ):
            plotData( i[0]+args.selection, i[2], i[3], i[4],
                            log=args.log, axisX=i[1], Norm=args.norm, legendAlignment=i[5], outputDir=outputDir, ext=args.ext)
        elif ( 'simple' in args.process ):
            plotSimpleComparison(
                    ###bkgFiles["TTToSemiLeptonic"][0], "TTToSemiLeptonic", signalFiles["ttHTobb"][0], "ttHTobb",
                    #ROOT.TFile('Rootfiles/'+VER+'/histograms_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_NOPUPPI_boosted.root'), "ttH_NOPUPPI",
                    ROOT.TFile('Rootfiles/'+VER+'/histograms_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_NOBTAG_boosted.root'), "ttH_NOBTAG",
                    ROOT.TFile('Rootfiles/'+VER+'/histograms_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_boosted.root'), "Nominal",
                    #ROOT.TFile('Rootfiles/'+VER+'/histograms_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_NOPUPPI_boosted.root'), "TTSemi_NOPUPPI",
                    ##ROOT.TFile('Rootfiles/'+VER+'/histograms_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_NOBTAG_boosted.root'), "TTSemi_NOBTAG",
                    #ROOT.TFile('Rootfiles/'+VER+'/histograms_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_boosted.root'), "Nominal",
                    i[0], xmin=i[2], xmax=i[3], rebinX=i[4], log=i[5], axisX=i[1] )

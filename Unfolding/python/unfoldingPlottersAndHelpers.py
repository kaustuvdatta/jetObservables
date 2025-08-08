import copy, pprint, array, bisect, scipy,os, sys, glob, math, gc
from array import array
from scipy import stats
from PIL import Image
from PyPDF2 import PdfMerger, PdfReader, PdfWriter

from collections import OrderedDict

import ROOT
import numpy as np


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()
ROOT.TH1.StatOverflows(ROOT.kTRUE)
ROOT.TH2.StatOverflows(ROOT.kTRUE)

from histoHelpers import *

sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
ROOT.gROOT.ForceStyle()

lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)

#############################################################################
###################### Helpers for plotting unfoldings ######################
#############################################################################

def plotSimpleComparison( inFile1, sample, inFile2, sample2, name, rebinX=1, xmin='', xmax='', labX=0.92, labY=0.50, axisX='', axisY='', log=False, ext='png', Norm=False, version='', outputDir='Plots/' ):
    """"Take two root files, make simple comparison plot"""

    outputFileName = name+'_'+sample+sample2+'_simpleComparisonPlot'+version+'.'+ext
    #print('Processing.......', outputFileName)

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
    
    if Norm:
        if 'data' in histo.GetName().lower():
            norm = histo.Integral()#Clone()
        elif 'data' in histo2.GetName().lower():
            norm = histo2.Integral()#Clone()
        else:
            norm=1.
    else:
        norm=1.
        

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
    histo.SetLineColorAlpha(ROOT.kRed,0.7)
    histo.SetMarkerColor(ROOT.kRed)
    histo.SetLineWidth(2)
    histo2.SetLineColor(ROOT.kBlue)
    #if Norm:
    histo.Scale(1./norm,'width')
    histo.Draw('histE')
    histo2.Scale(1./norm,'width')
    histo2.Draw('histE same')
    #else:
    ##    histo.Scale(1./norm,'width')
    #    histo2.Scale(1./norm,'width')

    #    histo.Draw('histe')
    #    histo2.Draw('histe same')
    if not axisY: histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: histo.GetXaxis().SetTitle( axisX )

    legend.Draw()

    canvas[name].SaveAs( outputDir+outName )
    if ext.startswith('pdf'):
        canvas[name].SaveAs( outputDir+outName.replace('pdf', 'png') )

def plotSysComparison( nomHisto, dictUncHistos, outputName, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ): #from Alejandro

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 219, 92, 48, 49, 37, 38, 33, 17, 50, 205, 225, 94, 221, 16,  225, 128]
    
    
    outputFileName = outputName+'_'+version+'.'+ext
    #print ('Processing.......', outputFileName)

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
    col_counter=0
    
    #print (dictUncHistos)
    for ih in dictUncHistos:
        #print(ih, col_counter, len(colors))
        dictgraph[ih] = ROOT.TGraphAsymmErrors()
        dictgraph[ih].Divide( dictUncHistos[ih], nomHisto, 'pois' )
        if not col_counter==len(colors)-1:
            dictgraph[ih].SetLineColor( colors[col_counter] )
            dictgraph[ih].SetLineStyle( 2 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        else:
            col_counter = col_counter-len(colors)+2
            dictgraph[ih].SetLineColor( colors[col_counter] )
            dictgraph[ih].SetLineStyle( 3 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        if 'jes' in ih and ('2017' in ih or '2018' in ih): legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        else: legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        multiGraph.Add( dictgraph[ih] )
        col_counter=col_counter+1

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

    
def drawDataMCReco( ivar, selection, year, lumi, process,
                    dataJetHisto, nominal_recoJetHisto, 
                    alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto,
                    labelX, jetType, maxX, tlegendAlignment, outputName,log=False ):
    """docstring for drawDataMCReco (dijets)"""
    print ("Drawing Data/MC")
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    #[ROOT.TColor.GetColor("#bd1f01"),ROOT.TColor.GetColor("#3f90da"),ROOT.TColor.GetColor("#ffa90e")]
    
    #ROOT.gStyle.SetPadRightMargin(0.05)
    #ROOT.gStyle.SetPadLeftMargin(0.15)
    #can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    #can.cd()
    
    if log==True:
        
        outputName=outputName.replace("DataMC_","DataMC_logPlot_")

    extraSpace = 0.02
    # Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    # Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    # Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    # Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    # references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw

    can = ROOT.TCanvas('canUnfolding'+ivar, 'canUnfolding'+ivar,  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
        
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0, Hdw / H, 1, 1, -1)
    
    
    #pad1.SetPad(0, Hdw / H, 1, 1)
    pad1.SetLeftMargin(L)
    pad1.SetRightMargin(R)
    pad1.SetTopMargin(Tup)
    pad1.SetBottomMargin(Bup)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1) 
    
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.63,0.55,0.90,0.89)
    else: legend=ROOT.TLegend(0.175,0.55,0.545,0.89) #legend=ROOT.TLegend(0.19,0.60,0.45,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.042)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    print(dataJetHisto,nominal_recoJetHisto,alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto)
    
    dataHisto = dataJetHisto.Clone()
    recoHisto = nominal_recoJetHisto.Clone()
    alt0recoHisto = alt0_recoJetHisto.Clone()
    alt1recoHisto = alt1_recoJetHisto.Clone()
    alt2recoHisto = alt2_recoJetHisto.Clone()
    
    if 'tau' in ivar: dataHisto.Scale(1, 'width')  ### divide by bin width
    #dataHisto.Scale(1/dataHisto.Integral(), 'width')  ### divide by bin width
    dataHisto.SetMarkerStyle(8)
    dataHisto.SetMarkerSize(1)
    dataHisto.SetMarkerColor(ROOT.kBlack)
    dataHisto.SetLineColor(ROOT.kBlack)
    #dataHisto.SetFillColorAlpha(16,0.7)
    #dataHisto.SetLineWidth(1)
    
    if 'tau' in ivar: recoHisto.Scale(1, 'width')
    
    recoHisto.SetLineWidth(2)
    recoHisto.SetLineColor(colors[0])
    recoHisto.SetMarkerColor(colors[0])
    recoHisto.SetMarkerStyle(25)
    recoHisto.SetMarkerSize(1)

    print(labelX)

    if 'tau' in ivar or '#' in labelX:
        dataHisto.GetYaxis().SetTitle( 'Events')#'#frac{dN}{d#'+labelX.split('#')[1]+'}' )##frac{1}{dN} , +'    [A.U.]'
        #print( '#frac{dN}{d#'+labelX.split('#')[1]+'}')
    else:
        bw = np.round(recoHisto.GetXaxis().GetBinLowEdge(3)-recoHisto.GetXaxis().GetBinLowEdge(2),3)
        label=None
        if 'pt' in ivar:
            label = 'p_{T}'
        elif 'mass'in ivar:
            label = 'm'
        elif 'softdrop' in ivar:
            label = 'm_{SD}'
        elif 'y' in ivar:
            label = 'y'
        elif 'npv' in ivar.lower():
            label= 'N_{PV}'
        else:
            pass
        #if label:
        #    dataHisto.GetYaxis().SetTitle( '#frac{dN}{d'+label+f'}} [Events/{bw} GeV]')# if 'pt' in ivar or 'mass' in ivar or 'softdrop' in ivar else '') ) 
        #else:
        dataHisto.GetYaxis().SetTitle( f' Events/[{bw} '+ ( 'GeV]' if 'pt' in ivar or 'mass' in ivar or 'softdrop' in ivar else 'A.U.]') ) 
    
    #dataHisto.GetYaxis().SetTitleOffset(1.15)
    dataHisto.GetYaxis().CenterTitle()
    #dataHisto.GetYaxis().SetTitleSize(0.05)
    #dataHisto.GetYaxis().SetLabelSize(0.05)
    #dataHisto.GetXaxis().SetTitleSize(0.0)
    #dataHisto.GetXaxis().SetLabelSize(0.0)
    #dataHisto.GetXaxis().SetTickLength(0.)
    
    dataHisto.GetXaxis().SetTitleOffset(999)    
    dataHisto.GetXaxis().SetLabelOffset(999)    
    dataHisto.GetYaxis().SetLabelOffset(0.011* H_ref / Hup)    
    dataHisto.GetYaxis().SetTitleOffset(extraSpace+1.1*Hup/H_ref)    
    dataHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    dataHisto.GetYaxis().SetLabelSize(0.045* H_ref / Hup)
    dataHisto.GetYaxis().SetTitleFont(42)
    dataHisto.SetMaximum( 1.7*max([ recoHisto.GetMaximum(), dataHisto.GetMaximum()] ) if not('pt') in ivar else 40.*max([ recoHisto.GetMaximum(), dataHisto.GetMaximum()] )  )
    dataHisto.SetMinimum(0. if not log else 0.01)
    
    dataHisto.SetTickLength(0.03, "XY")  # ?? ok if 1/3
    dataHisto.GetYaxis().SetMaxDigits(3)#,'y')
    ROOT.TGaxis.SetExponentOffset(-0.08, 0.014, "Y")        


    dataHisto.Draw( "AXIS")
    can.Update()

    legend.AddEntry( dataHisto, 'Data', 'pe' )

    
    dataHisto.Draw( "PE1 same")

    dataHisto.SetTitle('')
    can.SetTitle('')

    
    #alt1recoHisto.Scale(1, 'width')  ### divide by bin width
    #alt1recoHisto.SetLineWidth(3)
    #alt1recoHisto.SetLineColor(ROOT.kCyan+3)
    #alt1recoHisto.SetMarkerColor(ROOT.kCyan+3)
    #alt1recoHisto.SetMarkerStyle(25)
    #alt1recoHisto.SetMarkerSize(2)
    #legend.AddEntry( alt1recoHisto, 'MadGraph5+Pythia8', 'lp' )
    #alt1recoHisto.Draw("histe1 same")
    legend.AddEntry( recoHisto, 'MadGraph5+P8', 'lp' )

    recoHisto.Draw( "hist][ same")

    if 'tau' in ivar: alt0recoHisto.Scale(1, 'width')  ### divide by bin width
    alt0recoHisto.SetLineWidth(3)
    alt0recoHisto.SetLineStyle(2)
    alt0recoHisto.SetLineColor(colors[1])#ROOT.kBlue)
    alt0recoHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    alt0recoHisto.SetMarkerStyle(25)
    alt0recoHisto.SetMarkerSize(1)
    legend.AddEntry( alt0recoHisto, 'MadGraph5+H7', 'lp' )
    alt0recoHisto.Draw("hist][ same")
    
    
    if 'tau' in ivar: alt2recoHisto.Scale(1, 'width')  ### divide by bin width
    alt2recoHisto.SetLineWidth(3)
    alt2recoHisto.SetLineStyle(7)
    alt2recoHisto.SetLineColor(colors[2])#ROOT.kGray+4)
    alt2recoHisto.SetMarkerColor(colors[2])#ROOT.kGray+4)
    alt2recoHisto.SetMarkerStyle(25)
    alt2recoHisto.SetMarkerSize(1)
    legend.AddEntry( alt2recoHisto, 'P8+P8', 'lp' )
    alt2recoHisto.Draw("hist][ same")
    
    if log: ROOT.gPad.SetLogy()
    else: ROOT.gPad.SetLogy(0)
    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = ( 'Central Dijet' if 'Central' in jetType  else 'Forward Dijet' )#+' dijet region'
        dijetOffset = 0.0
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.63+dijetOffset ), 0.87, seltext )
    #selText.DrawLatex( 0.61 if 'dijet' in selection else 0.51, 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV, |y|<1.7' 
    elif selection.startswith("_W"): seltext = '#splitline{p_{T}>200 GeV, |y|<1.7}{65<m_{jet}<125 GeV}' 
    elif selection.startswith("_top"): seltext = '#splitline{p_{T}>400 GeV, |y|<1.7}{140<m_{jet}<300 GeV}'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.63+dijetOffset ), 0.77 if not('dijet' in selection) else 0.80, seltext )
    
    legend.Draw()
    
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        if year=='all': 
            #if 'dijet' in selection:
            CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
        else:
            CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0, 0, 1, Hdw / H,-1);
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
    
    pad2.SetGrid()
    
    pad2.SetTopMargin(Tdw)
    pad2.SetBottomMargin(Bdw)
    pad2.SetLeftMargin(L)
    pad2.SetRightMargin(R)
    
    
    pad2.Draw()
    pad2.Update()
    can.Update()
    pad2.cd()
    
    
    tmpPad2= pad2.DrawFrame( recoHisto.GetXaxis().GetBinLowEdge(1), 0., maxX, 1.9 )
    #print (labelX)
    #tmpPad2.GetYaxis().SetTitle( "Sim./Data." )
    #tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(0.25, 1.85 )
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()   
       
    if 'tau' in ivar: tmpPad2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: tmpPad2.GetXaxis().SetTitle(labelX)
        
    tmpPad2.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    tmpPad2.GetXaxis().SetTitleOffset(0.96)
    tmpPad2.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
    tmpPad2.GetYaxis().SetLabelSize(0.045 * H_ref / Hdw)#, "Y")
    tmpPad2.GetXaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "X")
    tmpPad2.GetXaxis().SetLabelSize(0.045 * H_ref / Hdw)#, "X")
    tmpPad2.GetXaxis().SetLabelOffset(0.011 * H_ref / Hdw)#, "X")
    #ratio_nominal.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    
    #if 'tau' in ivar: ratio_nominal.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #else: ratio_nominal.GetXaxis().SetTitle(labelX)
        
    tmpPad2.GetYaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetNdivisions(505)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    # Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  # ?? ok if 1/3
    tmpPad2.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

    # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    tmpPad2.GetYaxis().SetNdivisions(505)
    tmpPad2.GetYaxis().CenterTitle()
    
    ratio_nominal = ROOT.TGraphAsymmErrors()#len(l_bins)-1,x_bins,y_vals)
    ratio_nominal.SetName(ratio_nominal.GetName()+'_ratio')#+ivar)
    ratio_nominal.SetLineColor(colors[0])
    ratio_nominal.SetLineWidth(2)
    ratio_nominal.SetMarkerColor(colors[0])
    ratio_nominal.SetMarkerSize(1)
        
    ratio_nominal.Divide(recoHisto, dataHisto, 'pois' )
    #ratio_nominal.GetYaxis().SetRangeUser(0.1,2.1 )
    
    #ratio_nominal.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_nominal.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    ratio_nominal.GetXaxis().SetTitleOffset(0.96)
    ratio_nominal.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
    ratio_nominal.GetYaxis().SetLabelSize(0.045 * H_ref / Hdw)#, "Y")
    ratio_nominal.GetXaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "X")
    ratio_nominal.GetXaxis().SetLabelSize(0.045 * H_ref / Hdw)#, "X")
    ratio_nominal.GetXaxis().SetLabelOffset(0.011 * H_ref / Hdw)#, "X")
    #ratio_nominal.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    ratio_nominal.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    
    if 'tau' in ivar: ratio_nominal.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: ratio_nominal.GetXaxis().SetTitle(labelX)
        
    ratio_nominal.GetYaxis().SetTitleFont(42)
    ratio_nominal.GetXaxis().SetTitleFont(42)

    # Set tick lengths to match original (these are fractions of axis length)
    ratio_nominal.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  # ?? ok if 1/3
    ratio_nominal.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

    # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    ratio_nominal.GetYaxis().SetNdivisions(505)
    ratio_nominal.GetYaxis().CenterTitle()
    
    
    
    #ratio_nominal.SetLineColor(colors[0])
    #ratio_nominal.SetMarkerColor(colors[0])
    #ratio_nominal.SetMarkerSize(1)
    ratio_nominal.GetXaxis().SetNdivisions(505)
    ratio_nominal.GetYaxis().SetNdivisions(505)
    ratio_nominal.SetMarkerStyle(25)
    ratio_nominal.Draw('PE1 ')
    
    ratio_alt0MC = ROOT.TGraphAsymmErrors()
    ratio_alt0MC.SetName(ratio_alt0MC.GetName()+'_ratio')#+ivar)    
    ratio_alt0MC.Divide(  alt0recoHisto, dataHisto, 'pois' )
    ratio_alt0MC.SetLineStyle(2)
    ratio_alt0MC.SetLineWidth(3)
    ratio_alt0MC.SetLineColor(colors[1])
    ratio_alt0MC.SetMarkerColor(colors[1])
    ratio_alt0MC.SetMarkerStyle(25)
    ratio_alt0MC.SetMarkerSize(1)
    ratio_alt0MC.Draw('PE1 same')
    
    #ratio_alt1MC = ROOT.TGraphAsymmErrors()
    #ratio_alt1MC.Divide(  dataHisto,alt1recoHisto, 'pois' )
    #ratio_alt1MC.SetLineColor(ROOT.kCyan+3)
    #ratio_alt1MC.SetMarkerColor(ROOT.kCyan+3)
    #ratio_alt1MC.SetMarkerStyle(25)
    #ratio_alt1MC.SetMarkerSize(2)
    #ratio_alt1MC.Draw('PE1 same')
    
    ratio_alt2MC = ROOT.TGraphAsymmErrors()
    ratio_alt2MC.SetName(ratio_alt2MC.GetName()+'_ratio')#+ivar)    
    ratio_alt2MC.Divide(  alt2recoHisto, dataHisto,  'pois' )
    ratio_alt2MC.SetLineColor(colors[2])#ROOT.kGray+4)
    ratio_alt2MC.SetLineWidth(3)#ROOT.kGray+4)
    ratio_alt2MC.SetLineStyle(7)#ROOT.kGray+4)
    ratio_alt2MC.SetMarkerColor(colors[2])#ROOT.kGray+4)
    ratio_alt2MC.SetMarkerStyle(25)
    ratio_alt2MC.SetMarkerSize(1)
    ratio_alt2MC.Draw('PE1 same')
    
    
        
    pad2.Update()
    #ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    can.SaveAs(root_macro)
    gc.collect()
    
    plotHistos = [ 
                    dataHisto.Clone(dataHisto.GetName()+'_forHEPData'),
                    recoHisto.Clone(recoHisto.GetName()+'_forHEPData'),
                    alt0recoHisto.Clone(alt0recoHisto.GetName()+'_forHEPData'),
                    alt2recoHisto.Clone(alt2recoHisto.GetName()+'_forHEPData'),
                    ratio_nominal.Clone(ratio_nominal.GetName()+'_forHEPData'),
                    ratio_alt0MC.Clone(ratio_alt0MC.GetName()+'_forHEPData'),
                    ratio_alt2MC.Clone(ratio_alt2MC.GetName()+'_forHEPData'),
                  ]
    #hist1D_to_yoda = [yoda.root.to_yoda(  
    #                    unfoldHisto.Clone( 'normed_unfoldHisto'+ ivar + selection ) ) 
    #                 ]
    #yoda.writeYODA(hist1D_to_yoda, f"{outputDir}normed_1D_hists{ivar}{selection}.yoda" )
    print(outputName,png)

    outputFile = ROOT.TFile.Open(f"{outputName.split('.pdf')[0]}.root","RECREATE")
    for h in plotHistos:
        h.Write()
    can.Write()
    outputFile.Close()
    
    #ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)
                    


    
def drawUnfold(ivar, selection, process, year, lumi,
               dataJetHisto, genJetHisto, 
               unfoldHistowoUnc, altMCHisto, foldHisto, recoJetHisto,
               cov_tot, cov_datastat_tot, labelX, maxX, tlegendAlignment, outputName,
               altMC1Histo = None, altMC2Histo = None, altMC1Histo_label = None, altMC2Histo_label = None, 
               nomMCHisto_label = None, altMCHisto_label = None,
               extraMC=False, includeFSR = False, fsrUpHisto = None, fsrDownHisto=None, noNorm=False
              ):
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    # Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    # Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    # Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    # Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    # references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw

    can = ROOT.TCanvas('canUnfolding'+ivar, 'canUnfolding'+ivar,  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
    

    #can.SetFrameLineColor(1)
    #can.SetFrameLineStyle(1)
    #can.SetFrameLineWidth(1)
    
    
    #can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0, Hdw / H, 1, 1, -1)
    
    
    #pad1.SetPad(0, Hdw / H, 1, 1)
    pad1.SetLeftMargin(L)
    pad1.SetRightMargin(R)
    pad1.SetTopMargin(Tup)
    pad1.SetBottomMargin(Bup)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1) 
    
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.58,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.58,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.040)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    unfoldHisto = unfoldHistowoUnc.Clone('unfoldHisto'+ivar)
    unfoldHisto.Sumw2()
 
    unfoldHistoDataStatErr=unfoldHistowoUnc.Clone('unfoldHistoStatUnc'+ivar)
    unfoldHistoDataStatErr.Sumw2()
    
    dataJetHisto.SetTitle("")
    print("data(minus bkgs).Integral()",dataJetHisto.Integral())
    
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    unfoldHistoDataStatErr.SetTitle("")
    
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    
    recoJetHisto.SetTitle("")
    print("(RM proj.Y )recoJetHisto.Integral()",recoJetHisto.Integral())
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())
    
    
    #dataScaling = unfoldHisto.Integral()
    #print (dataScaling)
    
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the transformed covariance matrix, from the unnormalised -> the normalised space
    
    cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHisto, 
                                                       cov_tot.Clone())
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatErr, 
                                                                   cov_datastat_tot.Clone())
    #print(cov_normTot_np)
    #print(cov_norm_dataStat_np)
    

    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()

    unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')

    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.))
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_norm_dataStat_np)
    else:
        cov_abs_dataStat_np, _ = th2_to_ndarray(cov_datastat_tot.Clone())
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_abs_dataStat_np)
        
    unfoldHistoDataStatErr.Scale(1.,'width')

    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.))
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    else:
        cov_absTot_np, _ = th2_to_ndarray(cov_tot.Clone())
        get_th1_normedCovErrors(unfoldHisto, cov_absTot_np)
        
    unfoldHisto.Scale(1.,'width')

    
    dataJetHisto.Scale(1./(dataJetHisto.Integral() if not(noNorm) else 1.),'width')
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.),'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.),'width')
    foldHisto.Scale(1./(foldHisto.Integral() if not(noNorm) else 1.),'width')
    recoJetHisto.Scale(1./(recoJetHisto.Integral() if not(noNorm) else 1.),'width')
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        fsrUpHisto.Scale(1./(fsrUpHisto.Integral() if not(noNorm) else 1.),'width')
        fsrDownHisto.Sumw2()
        fsrDownHisto.Scale(1./(fsrDownHisto.Integral() if not(noNorm) else 1.),'width')
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        altMC1Histo.Scale(1./(altMC1Histo.Integral() if not(noNorm) else 1.),'width')
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo.Scale(1./(altMC2Histo.Integral() if not(noNorm) else 1.),'width')
            
            altMC2Histo.SetTitle("")

    
    
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(1)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(1)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(1)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(23)


        fsrDownHisto.SetMarkerSize(1)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(22)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )
    
    unfoldHisto.GetXaxis().SetTitleOffset(999)    
    unfoldHisto.GetXaxis().SetLabelOffset(999)    
    #unfoldHisto.GetYaxis().SetTitleOffset(0.012)    
    unfoldHisto.GetYaxis().SetLabelOffset(0.011* H_ref / Hup)    
    #unfoldHisto.GetYaxis().SetNdivisions(505)
    unfoldHisto.GetXaxis().SetNdivisions(505)
    unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+1.1*Hup/H_ref)    
    unfoldHisto.GetYaxis().SetTitleSize(0.055* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelSize(0.046* H_ref / Hup)
   
    if 'tau' in labelX: 
        
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    
    
    unfoldHisto.GetYaxis().SetTitleFont(42)
    unfoldHisto.SetMaximum( 1.7*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )
    unfoldHisto.SetTickLength(0.03, "XY")  # ?? ok if 1/3

    unfoldHisto.Draw( "AXIS")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    can.Update()
    
    altMCHisto.SetLineWidth(2)
    altMCHisto.SetMarkerSize(1)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 
        legend.AddEntry(fsrDownHisto, #('MadGraph5+P8, ' if 'dijet' in selection else 'POWHEG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MadGraph5+P8, ' if 'dijet' in selection else 'POWHEG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')

        
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'POWHEG+H7','lpe')#
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            altMC2Histo.SetLineWidth(2)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(1)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            altMC1Histo.SetLineWidth(2)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(1)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            
            legend.AddEntry( altMC1Histo, altMC1Histo_label.replace('-FXFX',''),'lpe')#'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1 same")

        
    genJetHisto.Draw( "histE1 same")
    altMCHisto.Draw("histE1 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1 same")
        fsrDownHisto.Draw("PE1 same")

    unfoldHisto.Draw( "E1 same")

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.2
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.80, seltext )
    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        if year=='all': 
            #if 'dijet' in selection:
            CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
        else:
            CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0, 0, 1, Hdw / H,-1);
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
    
    pad2.SetGrid()
    
    pad2.SetTopMargin(Tdw)
    pad2.SetBottomMargin(Bdw)
    pad2.SetLeftMargin(L)
    pad2.SetRightMargin(R)
    
    
    pad2.Draw()
    pad2.Update()
    can.Update()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistoDataStatErr.Clone('ratio_datastatUnc')
    ratio_datastatUnc.Divide(unfoldHistowoUnc)
    ratio_totalUnc = unfoldHisto.Clone('ratio_totalUnc')
    ratio_totalUnc.Divide(unfoldHistowoUnc)
    
    tmpPad2= pad2.DrawFrame( 0, 0.3, maxX, 1.9 )
    #print (labelX, label)
    
    #tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #tmpPad2.GetYaxis().CenterTitle()
    #tmpPad2.SetLabelSize(0.13, 'x')
    #tmpPad2.SetTitleSize(0.12, 'x')
    #tmpPad2.SetLabelSize(0.12, 'y')
    #tmpPad2.SetTitleSize(0.12, 'y')
    #tmpPad2.SetNdivisions(505, 'x')
    #tmpPad2.SetNdivisions(505, 'y')
    
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3345 if not('dijet' in selection) else 3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.94)
    ratio_totalUnc.SetTitleSize(0.055 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.055 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    if 'tau' in ivar: ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: ratio_totalUnc.GetXaxis().SetTitle(labelX)
        
    
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    # Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  # ?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    ratio_totalUnc.GetYaxis().SetNdivisions(505)
    ratio_totalUnc.GetYaxis().CenterTitle()
    #ratio_totalUnc.GetXaxis().SetLabelSize(0.12)
    #ratio_totalUnc.GetXaxis().SetTitleSize(0.13)

    #ratio_totalUnc.GetYaxis().SetLabelSize(0.12)
    #ratio_totalUnc.GetYaxis().SetTitleSize(0.12)
    ratio_totalUnc.GetXaxis().SetNdivisions(505)
    #ratio_totalUnc.GetYaxis().SetNdivisions(505)
    
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)

    ratio_totalUnc.SetFillColorAlpha(14,0.9)
    ratio_totalUnc.SetLineColor(14)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3354)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.5,0.5)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerStyle(23)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerStyle(22)
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
    
    hRatio.SetMarkerSize(1)
    hRatio.Draw('P0 same')
    
    hRatio2.SetMarkerSize(1)
    hRatio2.Draw('P0 same')
    
    hRatio5.SetMarkerSize(1)
    hRatio5.Draw('P0 same')
    
    if includeFSR:
        hRatio3.SetMarkerSize(1)
        hRatio3.Draw('P0 same')

        hRatio4.SetMarkerSize(1)
        hRatio4.Draw('P0 same')
    
    
    ratioLegend=ROOT.TLegend(0.19,0.78,0.69,0.88)
    ratioLegend.SetTextSize(0.09)
    ratioLegend.SetTextFont(42)
    ratioLegend.SetNColumns(2)
    ratioLegend.SetFillStyle(0)#ColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    ratioLegend.AddEntry( ratio_totalUnc, 'Total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    can.SaveAs(root_macro)

    #ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)    

    

def drawClosures(ivar, selection, process, year, lumi, genJetHisto, genJetHistoCross, unfoldHisto, unfoldHistoCross,
                 ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, 
                 cov_tot, cov_datastat_tot,
                 cov_tot_cross, cov_datastat_tot_cross,
                 labelX, maxX, tlegendAlignment, 
                 outputName, nomMCHisto_label = None, altMCHisto_label = None, noNorm = False ):
    
    if process.startswith('MCCrossClosure'):
    
        genJetHistoCross.SetTitle("") 
        unfoldHistoCross.SetTitle("")
        
    else:

        genJetHisto.SetTitle("") 
        unfoldHisto.SetTitle("")

    
    """docstring for drawClosures"""
    print ("Drawing unfolding closure")
    extraSpace = 0.01
    # Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    # Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    # Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    # Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    # references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw

    can = ROOT.TCanvas('canUnfolding'+ivar, 'canUnfolding'+ivar,  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
    

    #can.SetFrameLineColor(1)
    #can.SetFrameLineStyle(1)
    #can.SetFrameLineWidth(1)
    
    
    #can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0, Hdw / H, 1, 1, -1)
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    #pad1.SetPad(0, Hdw / H, 1, 1)
    pad1.SetLeftMargin(L)
    pad1.SetRightMargin(R)
    pad1.SetTopMargin(Tup)
    pad1.SetBottomMargin(Bup)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    
    #if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.56,0.61,0.88,0.89)
    #    else: legend=ROOT.TLegend(0.20,0.61,0.52,0.89)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.60,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.60,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    
        
    cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHisto.Clone(), 
                                                       cov_tot.Clone())
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHisto.Clone(), 
                                                                   cov_datastat_tot.Clone())
    #print(cov_normTot_np)
    #print(cov_norm_dataStat_np)
    
    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.))
    if not(noNorm): get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    unfoldHisto.Scale(1.,'width')
    
    
    if process.startswith('MCCrossClosure'): 
        print(genJetHisto.Integral(), unfoldHisto.Integral(),unfoldHistoCross.Integral(),genJetHistoCross.Integral())
        
        cov_normTot_np_cross, normed_covTot_cross = get_normalised_cov(unfoldHistoCross.Clone(), 
                                                                       cov_tot_cross.Clone())
        cov_norm_dataStat_np_cross, normed_cov_dataStat_cross = get_normalised_cov(unfoldHistoCross.Clone(), 
                                                                                   cov_datastat_tot_cross.Clone())
        
        unfoldHistoCross.Scale(1./(unfoldHistoCross.Integral() if not(noNorm) else 1.))
        if not(noNorm): get_th1_normedCovErrors(unfoldHistoCross, cov_normTot_np_cross)
        unfoldHistoCross.Scale(1.,'width')

    
    
    
    #unfoldHisto.Scale(1./(unfoldHisto.Integral() if not noNorm else 1.),'width')
    unfoldHisto.SetMarkerStyle(4)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    unfoldHisto.SetLineWidth(2)
    
    legend.AddEntry( unfoldHisto, (f'{nomMCHisto_label} (closure)' if process.startswith('MCSelfClosure') else f'#splitline{{{nomMCHisto_label} unf. with }}{{{nomMCHisto_label.replace("-FXFX","")} }}'), 'pe' )
    
    
    genJetHisto.Scale(1./(genJetHisto.Integral() if not noNorm else 1.),'width')
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(ROOT.kBlue)
    genJetHisto.SetMarkerStyle(0)
    genJetHisto.SetLineStyle(2)
    legend.AddEntry( genJetHisto, f'{nomMCHisto_label} (gen)', 'lp' )
    
    unfoldHisto.GetXaxis().SetTitleOffset(999)    
    unfoldHisto.GetXaxis().SetLabelOffset(999)    

    
    if 'tau' in labelX:
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
            
    unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+1.1*Hup/H_ref)    
    unfoldHisto.GetYaxis().SetTitleSize(0.055* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelSize(0.046* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelOffset(0.011* H_ref / Hup)
    #unfoldHisto.GetYaxis().SetNdivisions(505)
    unfoldHisto.GetXaxis().SetNdivisions(505)
    unfoldHisto.GetYaxis().SetTitleFont(42)

    unfoldHisto.SetTickLength(0.03, "XY")  # ?? ok if 1/3
    unfoldHisto.SetMaximum( 1.6*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum(), unfoldHistoCross.GetMaximum() if 'Cross' in process else unfoldHisto.GetMaximum() ] )  )
    
    unfoldHisto.Draw( "AXIS")
    can.Update()
    can.Modified()
    unfoldHisto.Draw( "E")    
    genJetHisto.Draw( "histe same")
    
    if not process.startswith('MCSelfClosure'):
        
        #unfoldHistoCross.Scale(1./(unfoldHistoCross.Integral() if not noNorm else 1.),'width')
        
        unfoldHistoCross.SetMarkerStyle(26)
        #unfoldHistoCross.SetMarkerSize(2)
        unfoldHistoCross.SetMarkerColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineWidth(2)
        legend.AddEntry( unfoldHistoCross, f'#splitline{{{nomMCHisto_label} unf. with }}{{{altMCHisto_label.replace("-FXFX","")} }}', 'pe')

        genJetHistoCross.Scale(1./(genJetHistoCross.Integral() if not noNorm else 1.),'width')
        
        genJetHistoCross.SetLineWidth(2)
        genJetHistoCross.SetLineColor(ROOT.kMagenta)
        genJetHistoCross.SetMarkerStyle(0)
        genJetHistoCross.SetLineStyle(2)
        legend.AddEntry( genJetHistoCross, f'{altMCHisto_label.replace("-FXFX","")} (gen)', 'lp')
        
        unfoldHistoCross.Draw( "E same")
        genJetHistoCross.Draw( "histe same")

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.20
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.80, seltext )

    
    legend.Draw()
    #if process.startswith('data'):
    CMS_lumi.extraText = "Simulation Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #    CMS_lumi.extraText = "Simulation Preliminary"
    #    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0, 0, 1, Hdw / H,-1);
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
        
    
    pad2.SetGrid()
    pad2.SetTopMargin(Tdw)
    pad2.SetBottomMargin(Bdw)
    pad2.SetLeftMargin(L)
    pad2.SetRightMargin(R)

    
    pad2.Draw()
    pad2.cd()
    
    tmpPad2= pad2.DrawFrame( 0, 0.3, maxX, 1.9 )
    print (labelX)
    
    tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    tmpPad2.GetXaxis().SetRangeUser( unfoldHisto.GetBinLowEdge(1),  unfoldHisto.GetBinLowEdge( unfoldHisto.GetNbinsX() + 1 ) )

   
    
    tmpPad2.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    tmpPad2.GetXaxis().SetTitleOffset(0.94)
    tmpPad2.SetTitleSize(0.055 * H_ref / Hdw, "Y")
    tmpPad2.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    tmpPad2.SetTitleSize(0.055 * H_ref / Hdw, "X")
    tmpPad2.SetLabelSize(0.046 * H_ref / Hdw, "X")
    tmpPad2.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    if 'tau' in labelX: 
        tmpPad2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else:
        tmpPad2.GetXaxis().SetTitle( labelX)
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Unf.}" if 'Self' in process else "#frac{Unf. alt.}{Unf. nom.}" )
    tmpPad2.GetYaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    # Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  # ?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    tmpPad2.GetYaxis().SetNdivisions(505)
    tmpPad2.GetXaxis().SetNdivisions(505)
    tmpPad2.GetYaxis().CenterTitle()
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    pad2.cd()
    
    if 'Self' in process:

        hRatioUp = ROOT.TGraphAsymmErrors()
        hRatioUp.Divide( genJetHisto, unfoldHisto, 'pois' )
        hRatioUp.SetLineColor(ROOT.kBlack)
        hRatioUp.SetMarkerColor(ROOT.kBlack)
        hRatioUp.SetLineWidth(2)
        hRatioUp.SetMarkerStyle(25)
        
        #hRatioUp.GetXaxis().SetLimits(0.,maxX)#unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        
        #hRatioUp.GetYaxis().SetRangeUser(0.3,1.9 )
    
        #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
        hRatioUp.GetYaxis().SetTitleOffset(extraSpace + (1.15 ) * Hdw / H_ref)
        hRatioUp.GetXaxis().SetTitleOffset(0.92)
        hRatioUp.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
        hRatioUp.GetYaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "Y")
        hRatioUp.GetXaxis().SetTitleSize(0.056 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelOffset(0.012 * H_ref / Hdw)#, "X")
        #hRatioUp.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        if 'tau' in labelX: 
            hRatioUp.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        else:
            hRatioUp.GetXaxis().SetTitle( labelX)
        
        
        hRatioUp.GetYaxis().SetTitle( "#frac{Sim.}{Unf.}" )
        hRatioUp.GetYaxis().SetTitleFont(42)
        hRatioUp.GetXaxis().SetTitleFont(42)

        #tmpPad2.GetXaxis().SetTitle(nameXaxis)
        #tmpPad2.GetYaxis().SetTitle(nameRatio)

        # Set tick lengths to match original (these are fractions of axis length)
        hRatioUp.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  # ?? ok if 1/3
        hRatioUp.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp.GetYaxis().SetNdivisions(505)
        hRatioUp.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp,1.5,0.5)
        
        
        hRatioUp.Draw('P0')

    else:
        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfoldHistoCross, unfoldHisto, 'pois' )
        hRatioUp2.SetLineColor(ROOT.kBlack)
        hRatioUp2.SetMarkerColor(ROOT.kBlack)
        hRatioUp2.SetLineWidth(2)
        hRatioUp2.SetMarkerStyle(25)
        
        #hRatioUp2.GetXaxis().SetLimits(0.,maxX)#unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        
        #hRatioUp2.GetYaxis().SetRangeUser(0.3,1.9 )
    
        #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
        hRatioUp2.GetYaxis().SetTitleOffset(extraSpace + (1.15 ) * Hdw / H_ref)
        hRatioUp2.GetXaxis().SetTitleOffset(0.92)
        hRatioUp2.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
        hRatioUp2.GetYaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "Y")
        hRatioUp2.GetXaxis().SetTitleSize(0.056 * H_ref / Hdw)#, "X")
        hRatioUp2.GetXaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "X")
        hRatioUp2.GetXaxis().SetLabelOffset(0.012 * H_ref / Hdw)#, "X")
        #hRatioUp2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        if 'tau' in labelX: 
            hRatioUp2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        else:
            hRatioUp2.GetXaxis().SetTitle( labelX)
        hRatioUp2.GetYaxis().SetTitle( "#frac{Unf. alt.}{Unf. nom.}" )
        hRatioUp2.GetYaxis().SetTitleFont(42)
        hRatioUp2.GetXaxis().SetTitleFont(42)

        #tmpPad2.GetXaxis().SetTitle(nameXaxis)
        #tmpPad2.GetYaxis().SetTitle(nameRatio)

        # Set tick lengths to match original (these are fractions of axis length)
        hRatioUp2.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  # ?? ok if 1/3
        hRatioUp2.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp2.GetYaxis().SetNdivisions(505)
        hRatioUp2.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp2,1.5,0.5)
        
        
        hRatioUp2.Draw('P0')
    
    
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    #ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)
    
def drawUnfoldFromNdim(ivar, selection, process, year, lumi,
               genJetHisto, 
               unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoTotUnc,
               altMCHisto, 
               labelX, maxX, tlegendAlignment, outputName,
               altMC1Histo = None, altMC2Histo = None, altMC1Histo_label = None, altMC2Histo_label = None, 
               nomMCHisto_label = None, altMCHisto_label = None,
               extraMC=False, includeFSR = False, fsrUpHisto = None, fsrDownHisto=None, noNorm=False
              ):
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    # Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    # Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    # Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    # Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    # references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw

    can = ROOT.TCanvas('canUnfolding'+ivar, 'canUnfolding'+ivar,  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
    

    #can.SetFrameLineColor(1)
    #can.SetFrameLineStyle(1)
    #can.SetFrameLineWidth(1)
    
    
    #can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0, Hdw / H, 1, 1, -1)
    
    
    #pad1.SetPad(0, Hdw / H, 1, 1)
    pad1.SetLeftMargin(L)
    pad1.SetRightMargin(R)
    pad1.SetTopMargin(Tup)
    pad1.SetBottomMargin(Bup)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1) 
    
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.60,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.60,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.040)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    unfoldHisto = unfoldHistoTotUnc.Clone('unfoldHisto'+ivar)
    unfoldHisto.Sumw2()
 
    unfoldHistoDataStatErr=unfoldHistoDataStatUnc.Clone('unfoldHistoStatUnc'+ivar)
    unfoldHistoDataStatErr.Sumw2()
    
    
    
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    unfoldHistoDataStatErr.SetTitle("")
    
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    
    recoJetHisto.SetTitle("")
    print("(RM proj.Y )recoJetHisto.Integral()",recoJetHisto.Integral())
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())
    
    
    #dataScaling = unfoldHisto.Integral()
    #print (dataScaling)
    
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the transformed covariance matrix, from the unnormalised -> the normalised space
    
    #cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHisto, 
    #                                                   cov_tot.Clone())
    #cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatErr, 
    #                                                               cov_datastat_tot.Clone())
    #print(cov_normTot_np)
    #print(cov_norm_dataStat_np)
    

    genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()

    #unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.))#,'width')

    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.))
    #if not(noNorm): 
    #    #get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_norm_dataStat_np)
    #else:
    #    cov_abs_dataStat_np, _ = th2_to_ndarray(cov_datastat_tot.Clone())
    #    #get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_abs_dataStat_np)
        
    unfoldHistoDataStatErr.Scale(1.)#,'width')

    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.))
    #if not(noNorm): 
    #    #get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    #else:
    #    cov_absTot_np, _ = th2_to_ndarray(cov_tot.Clone())
    #    #get_th1_normedCovErrors(unfoldHisto, cov_absTot_np)
        
    unfoldHisto.Scale(1.)#,'width')
    
    
    dataJetHisto.Scale(1./(dataJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.))#,'width')
    foldHisto.Scale(1./(foldHisto.Integral() if not(noNorm) else 1.))#,'width')
    recoJetHisto.Scale(1./(recoJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        fsrUpHisto.Scale(1./(fsrUpHisto.Integral() if not(noNorm) else 1.))#,'width')
        fsrDownHisto.Sumw2()
        fsrDownHisto.Scale(1./(fsrDownHisto.Integral() if not(noNorm) else 1.))#,'width')
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        altMC1Histo.Scale(1./(altMC1Histo.Integral() if not(noNorm) else 1.),'width')
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo.Scale(1./(altMC2Histo.Integral() if not(noNorm) else 1.),'width')
            
            altMC2Histo.SetTitle("")

    
    
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(1)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(1)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(1)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(23)


        fsrDownHisto.SetMarkerSize(1)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(22)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )
    
    unfoldHisto.GetXaxis().SetTitleOffset(999)    
    unfoldHisto.GetXaxis().SetLabelOffset(999)    

    unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+1.15*Hup/H_ref)    
    unfoldHisto.GetYaxis().SetTitleSize(0.056* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelSize(0.047* H_ref / Hup)
   
    if 'tau' in labelX: 
        
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    
    
    unfoldHisto.GetYaxis().SetTitleFont(42)
    unfoldHisto.SetMaximum( (1.7 if '21' in ivar or '32' in ivar else 1.65)*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )
    unfoldHisto.SetTickLength(0.03, "XY")  # ?? ok if 1/3

    unfoldHisto.Draw( "AXIS")
    can.Update()
    
    altMCHisto.SetLineWidth(2)
    altMCHisto.SetMarkerSize(1)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 

        legend.AddEntry(fsrDownHisto, #('MadGraph5+P8, ' if 'dijet' in selection else 'POWHEG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MadGraph5+P8, ' if 'dijet' in selection else 'POWHEG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'POWHEG+H7','lpe')#
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            altMC2Histo.SetLineWidth(2)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(1)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            altMC1Histo.SetLineWidth(2)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(1)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            
            legend.AddEntry( altMC1Histo, altMC1Histo_label.replace('-FXFX',''),'lpe')#'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1 same")

        
    genJetHisto.Draw( "histE1 same")
    altMCHisto.Draw("histE1 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1 same")
        fsrDownHisto.Draw("PE1 same")

    unfoldHisto.Draw( "E1 same")

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.80, seltext )
    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        if year=='all': 
            #if 'dijet' in selection:
            CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
        else:
            CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0, 0, 1, Hdw / H,-1);
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)
    
    pad2.SetGrid()
    
    pad2.SetTopMargin(Tdw)
    pad2.SetBottomMargin(Bdw)
    pad2.SetLeftMargin(L)
    pad2.SetRightMargin(R)
    
    
    pad2.Draw()
    pad2.Update()
    can.Update()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistoDataStatErr.Clone('ratio_datastatUnc')
    ratio_datastatUnc.Divide(unfoldHistowoUnc)
    ratio_totalUnc = unfoldHisto.Clone('ratio_totalUnc')
    ratio_totalUnc.Divide(unfoldHistowoUnc)
    
    tmpPad2= pad2.DrawFrame( 0, 0.3, maxX, 1.9 )
    print (labelX)
    
    #tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #tmpPad2.GetYaxis().CenterTitle()
    #tmpPad2.SetLabelSize(0.13, 'x')
    #tmpPad2.SetTitleSize(0.12, 'x')
    #tmpPad2.SetLabelSize(0.12, 'y')
    #tmpPad2.SetTitleSize(0.12, 'y')
    #tmpPad2.SetNdivisions(505, 'x')
    #tmpPad2.SetNdivisions(505, 'y')
    
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (1.155 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.96)
    ratio_totalUnc.SetTitleSize(0.056 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.047 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.056 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.047 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    if 'tau' in labelX: 
        ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else:
        ratio_totalUnc.GetXaxis().SetTitle( labelX )
    
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    # Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  # ?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    # Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    ratio_totalUnc.GetYaxis().SetNdivisions(505)
    ratio_totalUnc.GetYaxis().CenterTitle()
    #ratio_totalUnc.GetXaxis().SetLabelSize(0.12)
    #ratio_totalUnc.GetXaxis().SetTitleSize(0.13)

    #ratio_totalUnc.GetYaxis().SetLabelSize(0.12)
    #ratio_totalUnc.GetYaxis().SetTitleSize(0.12)
    #ratio_totalUnc.GetXaxis().SetNdivisions(505)
    #ratio_totalUnc.GetYaxis().SetNdivisions(505)
    
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)

    ratio_totalUnc.SetFillColorAlpha(14,0.8)
    ratio_totalUnc.SetLineColor(14)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3354)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.5,0.5)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerStyle(22)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerStyle(23)
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
    
    hRatio.SetMarkerSize(1)
    hRatio.Draw('P0 same')
    
    hRatio2.SetMarkerSize(1)
    hRatio2.Draw('P0 same')
    
    hRatio5.SetMarkerSize(1)
    hRatio5.Draw('P0 same')
    
    if includeFSR:
        hRatio3.SetMarkerSize(1)
        hRatio3.Draw('P0 same')

        hRatio4.SetMarkerSize(1)
        hRatio4.Draw('P0 same')
    
    
    ratioLegend=ROOT.TLegend(0.19,0.78,0.69,0.88)
    ratioLegend.SetTextSize(0.09)
    ratioLegend.SetTextFont(42)
    ratioLegend.SetNColumns(2)
    ratioLegend.SetFillStyle(0)#ColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    ratioLegend.AddEntry( ratio_totalUnc, 'Total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    #ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)    
       
    
def get_condition_number(histo):
    
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
    conditionNumber = round( Max/Min, 2 ) if Min > 0 else 1000000
    
    return conditionNumber
    
def draw2D( ivar, histo, varInfo, outputDir, outputLabel='data', addCorrelation=False, addCondition=False, addInvertedMatrix=False,ext='pdf',selection='_dijetSel',version='vNew',year='2017',pngToo=False, outputName=None,  ):

    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_'+version+'.'+ext

    
    #ROOT.gStyle.SetPadRightMargin(0.15)
    #ROOT.gStyle.SetPadTopMargin(0.08)
    
    W_ref = 800
    H_ref = 600
    y_offset = 1.02
    x_offset = 1.
    extraSpace = 0.02
    W = W_ref
    H = H_ref
    T = 0.08 * H_ref
    B = 0.10 * H_ref
    L = 0.13 * H_ref
    R = 0.03 * H_ref
    
    can2D = ROOT.TCanvas(ivar+'can2D'+histo.GetName(), ivar+'can2D'+histo.GetName(), 50, 50, W, H )# if not('body' in ivar) else ROOT.TCanvas(ivar+'can2D'+histo.GetName(), ivar+'can2D'+histo.GetName(), 50, 50, W, H )
    can2D.SetRightMargin(0.15)
    can2D.SetFillColor(0)
    can2D.SetBorderMode(0)
    can2D.SetFrameFillStyle(0)
    can2D.SetFrameBorderMode(0)
    can2D.SetLeftMargin(L / W + extraSpace)
    #can2D.SetRightMargin(R / W)
    #if with_z_axis:
    #if 'cov' in outputName.lower() or 'resp' in outputName.lower():
    #    can2D.SetRightMargin(B / W + 0.06)
    #else:
    can2D.SetRightMargin(B / W +  0.07) 
        
    can2D.SetTopMargin(T / H)
    can2D.SetBottomMargin(B / H + 0.02)
    x_proj = histo.ProjectionX('x_proj'+histo.GetName())
    x_proj_ax = histo.GetXaxis()#ProjectionX('x_proj'+histo.GetName())

    x_min = -0.5 if 'body' in ivar else x_proj.GetBinLowEdge(1)#0.#histo.GetBinLowEdge(1)
    x_max = x_proj_ax.GetBinLowEdge(histo.GetNbinsX()+1)
    y_min = -0.5 if 'body' in ivar else x_proj.GetBinLowEdge(1)#0. 
    y_max = x_max#max([y_max]+[10+i*100. for i in total_unc])
    #if 'tau' in ivar: 
    #    x_axis_title = '#'+x_axis_title.split('#')[1] 
        
        
    histo.GetYaxis().SetTitleFont(42)
    histo.GetYaxis().SetTitleOffset(y_offset)
    histo.GetXaxis().SetTitleFont(42)
    histo.GetXaxis().SetTitleOffset(x_offset)
    
    
    histo.GetZaxis().SetTitleOffset(1.)
    
    histo.GetYaxis().SetTitleSize(0.06-0.006)
    histo.GetXaxis().SetTitleSize(0.06-0.006)
    histo.GetZaxis().SetTitleSize(0.045)
    
    histo.GetYaxis().SetLabelSize(0.05-0.004)
    histo.GetXaxis().SetLabelSize(0.05-0.004)
    histo.GetZaxis().SetLabelSize(0.04)
    histo.GetYaxis().SetLabelFont(42)
    histo.GetXaxis().SetLabelFont(42)
    histo.GetZaxis().SetLabelFont(42)
    
    
    if 'resp' in outputName.lower() or 'prob' in outputName.lower(): 
        histo.GetXaxis().SetTitle('AK8 gen jet '+varInfo['label'].replace('AK8 jet ','') +f" {'basis' if 'body' in outputName.lower() else ''}")   
    elif 'cov' in outputName.lower() or 'corr' in outputName.lower() : 
        histo.GetXaxis().SetTitle(varInfo['label'] +f" {'N-subjettiness basis' if 'body' in outputName.lower() else ''}")   
        
    
    
    if 'resp' in outputName.lower() or 'prob' in outputName.lower(): 
        histo.GetYaxis().SetTitle('AK8 reco jet '+varInfo['label'].replace('AK8 jet ','') +f" {'basis' if 'body' in outputName.lower() else ''}")
    elif 'cov' in outputName.lower() or 'corr' in outputName.lower() : 
        histo.GetYaxis().SetTitle(varInfo['label'] +f" {'N-subjettiness basis' if 'body' in outputName.lower() else ''}")
    
        
    
    histo.GetYaxis().SetRangeUser(y_min, y_max)
    histo.GetXaxis().SetRangeUser(x_min, x_max)
    
    #histo.SetLabelSize(0.047 * H_ref / Hdw, "X")
    histo.SetLabelOffset(0.012 , "X")#* H_ref / Hdw
    histo.SetLabelOffset(0.012 , "Y")#* H_ref / Hdw
    histo.SetLabelOffset(0.008 , "Z")#* H_ref / Hdw

    histo.SetTickLength(0.03, "Y")  # ?? ok if 1/3
    histo.SetTickLength(0.03, "X")
    histo.SetTickLength(0.03, "Z")

    
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    #ROOT.gStyle.SetNumberContours(255)
    """
    if 'body' in outputName.lower() and ('corr' in histo.GetName().lower() or 'cov' in histo.GetName().lower()) and not('closure' in outputDir.lower()):


        ncont = ROOT.gStyle.GetNumberContours()         # should be 255
        pal  = ROOT.TColor.GetPalette()            # returns a Python list of length ncont

        # pick the palette‐bin whose z‐value is exactly zero
        mid_i = ncont//2
        mid = copy.deepcopy(pal[mid_i])
        print(mid, mid_i)
        # swap in white
        pal[mid_i] = ROOT.TColor.GetColor('#89ADA4')#ROOT.kWhite
        #print( pal[mid_i], mid_i)
    """


    if 'cov' in outputName.lower():
        histo.GetZaxis().SetTitle('Covariance')
        histo.GetZaxis().SetMaxDigits(3)
        
    elif 'corr' in outputName.lower() or ('corr' in outputName.lower() and 'cov' in outputName.lower()):
        histo.GetZaxis().SetTitle('Correlation')
        histo.SetAxisRange(-1., 1.,"Z")
        
    elif 'proba' in outputName.lower():
        histo.GetZaxis().SetTitle('Probability')
        histo.SetAxisRange(0., 1.,"Z")
        
    elif 'resp' in outputName.lower():
        histo.GetZaxis().SetTitle('Events')
        histo.GetZaxis().SetMaxDigits(3)
        
    
    histo.Draw("colz")
    
    ROOT.gPad.Update()
    
    palette = histo.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(0.86)
    palette.SetX2NDC(0.90)
    palette.SetY1NDC(0.+B/H+0.02)
    palette.SetY2NDC(1.-T/H)
    can2D.RedrawAxis()
    can2D.Modified()
    can2D.Update()
    """
    if 'body' in outputName.lower() and ('correl' in histo.GetName().lower() or 'cov6' in histo.GetName().lower()) and not('closure' in outputDir.lower()) and ('normed' in outputLabel.lower()):
        print(outputLabel,outputName)
        
        nX = histo.GetNbinsX()
        
        nY = histo.GetNbinsY()
        
        print(nX,nY,mid,ncont) 
        boxes = []
        for ix in range(1, nX+1):
            x1 = histo.GetXaxis().GetBinLowEdge(ix)
            x2 = histo.GetXaxis().GetBinUpEdge(ix)
            for iy in range(1, nY+1):
                if histo.GetBinContent(ix, iy) == 0:
                    y1 = histo.GetYaxis().GetBinLowEdge(iy)
                    y2 = histo.GetYaxis().GetBinUpEdge(iy)
                    
                    box = ROOT.TBox(x1, y1, x2, y2)
                    #box.SetName(f'{x1}{y1}')
                    box.SetLineColorAlpha(ROOT.kWhite,0.6)
                    box.SetFillColorAlpha(ROOT.kWhite,0.6)
                    box.SetLineWidth(1)
                    box.SetFillStyle(1001)     
                    box.Draw('same')
                    boxes.append(box)
                    
                    #box = ROOT.TBox(x1, y1, x2, y2)
                    #box.SetName(f'{x1}{y1}')
                    #box.SetLineColorAlpha(mid,0.7)
                    #box.SetFillColorAlpha(mid,0.7)
                    #box.SetLineWidth(0)
                    #box.SetFillStyle(3154)     
                    #box.Draw('same')
                    #boxes.append(box)
        ROOT.gPad.Update()
        can2D.Modified()
        can2D.Update()
    """
    CMS_lumi.extraText = ("Simulation " if 'resp' in histo.GetName().lower() or 'prob' in histo.GetName().lower() else "")+"Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = f"{year} fb^{{-1}} (13 TeV)"
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(can2D, 4, 0)
    can2D.Update()
    
    
    if addCorrelation:
        print('|--> Correlation: ', histo.GetCorrelationFactor())
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        textBox.SetTextFont(42)
        textBoxCorr = textBox.Clone()
        textBoxCorr.DrawLatex( 0.15,# if not('body' in ivar) else 0.15, 
                               0.88,
                               #varInfo['bins'][-1]-( .05*(varInfo['bins'][-1]-varInfo['bins'][0]),
                               '#color[8]{Corr. Factor = '+str(np.round(histo.GetCorrelationFactor(),3))+'}' )

    if addCondition:   
                
        conditionNumber = get_condition_number(histo.Clone())
        if conditionNumber<=10.:
            print('|--> Condition Number: ', conditionNumber)
        else:
            print('##############################################')
            print (f' WARNING: Condition Number>10: {conditionNumber} ')
            print('##############################################')
            #print('|--> Condition Number: ', conditionNumber)
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        textBox.SetTextFont(42)

        textBoxCond = textBox.Clone()
        #textBoxCond
        textBoxCond.DrawLatex( 0.15, 
                               0.85,
                               #varInfo['bins'][-1]-( .1*(varInfo['bins'][-1]-varInfo['bins'][0]),
                               '#color[8]{Cond. Number = '+str(np.round(conditionNumber,3))+'}' )

    can2D.SaveAs(outputName)
    if ext.startswith('pdf') and pngToo:
        can2D.SaveAs( outputName.replace('pdf', 'png') )
        if 'cov' in outputName.lower() or 'corr' in outputName.lower() or 'proba' in outputName.lower() or 'rho' in outputName.lower():
            root_macro = outputName.split('.pdf')[0]+'.C'
            can2D.SaveAs(root_macro)
            root_file = outputName.split('.pdf')[0]+'.root'
            outputFile = ROOT.TFile.Open(root_file, "RECREATE")
            histo.SetName(histo.GetName()+'_forHEPData')
            histo.Write()
            can2D.Write()
            outputFile.Close()
    
    can2D.Close()        
    #del(can2D)
    gc.collect()
    
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetNumberContours(255)
    

    #h = histo.DrawCopy("colz")
    
    #if 'body' in outputName.lower() and ('corr' in histo.GetName().lower() or 'cov' in histo.GetName().lower()) and not('closure' in outputDir.lower()):
    #    
    #            
    #    ncont = ROOT.gStyle.GetNumberContours()         # should be 255
    #    pal  = ROOT.TColor.GetPalette()            # returns a Python list of length ncont

        
    #    pal[mid_i] = mid
        
    
        

def drawUncertainties_from_err_shifts_theoryVariations_unitNorm(
                                                                ivar, 
                                                                unfoldHistowoUnc,
                                                                uncerUnfoldHisto, 
                                                                cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, 
                                                                tlegendAlignment, outputName, year, selection, lumi, norming=True 
                                                                ):
    
    unftot = unfoldHistowoUnc.Integral()

    print (f'|> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    colors_cr = [ROOT.kMagenta+2,  ROOT.kBlue-4, 433]  
    colors_syst = get_colour_palette_as_list('vf_8')[1:]#list(reversed())#[ROOT.kBlue+1, ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kGreen+2]
    
    
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
        
    upstyles =   [20,21,34,29,22,23,29,47,33,43,39,41,39,45,117,114,48]
    downstyles = [24,25,28,30,26,32,30,46,27,42,37,40,37,44,38 ,60 , 5 ]
    
    modelVariations = True if 'VariationUNC' in outputName else False
    
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    
    W_ref = 800 #if square else 800
    H_ref = 600 #if square else 600
    extraSpace = 0.02
    W = W_ref
    H = H_ref
    T = 0.07 * H_ref
    B = 0.11 * H_ref
    L = 0.13 * H_ref
    R = 0.03 * H_ref

    canUnc = ROOT.TCanvas('canUnc_relUnc'+ivar, 'canUnc_relUnc'+ivar, 50, 50, W, H)
    canUnc.SetFillColor(0)
    canUnc.SetBorderMode(0)
    canUnc.SetFrameFillStyle(0)
    canUnc.SetFrameBorderMode(0)
    canUnc.SetLeftMargin(L / W + extraSpace)
    canUnc.SetRightMargin(R / W)
    #if with_z_axis:
    #    c.SetRightMargin(B / W + 0.03)
    canUnc.SetTopMargin(T / H)
    canUnc.SetBottomMargin(B / H + 0.02)
    
    #canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  50, 50, 800, 600 )
    #canUnc.SetTopMargin(0.08)
    
    legend=ROOT.TLegend(0.18,0.64,0.88,0.9)

    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.038)
    legend.SetTextFont(42)

    legend.SetBorderSize(0)
    
    
    unfoldHistoTotUnc = unfoldHistowoUnc.Clone('unfoldHistoTotUnc'+ivar)
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc = unfoldHistowoUnc.Clone('unfoldHistoDataStatUnc'+ivar)
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc = unfoldHistowoUnc.Clone('unfoldHistoRMStatUnc'+ivar)
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc = unfoldHistowoUnc.Clone('unfoldHistoBkgSubUnc'+ivar)
    unfoldHistoBkgSubUnc.Sumw2()
    
    
    cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHistoTotUnc, 
                                                       cov_tot.Clone())
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatUnc, 
                                                                   cov_datastat_tot.Clone())
    
    cov_norm_RMStat_np, normed_cov_RMStat = get_normalised_cov(unfoldHistoRMStatUnc, 
                                                               cov_rmstat_tot.Clone())
    cov_norm_BkgSub_np, normed_cov_BkgSub = get_normalised_cov(unfoldHistoBkgSubUnc, 
                                                               cov_bkg_tot.Clone())
    
    
    
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    
    
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoTotUnc, cov_normTot_np)
    unfoldHistoTotUnc.Scale(1.,'width')
    
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoDataStatUnc, cov_norm_dataStat_np)
    unfoldHistoDataStatUnc.Scale(1.,'width')
    
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoRMStatUnc, cov_norm_RMStat_np)
    unfoldHistoRMStatUnc.Scale(1.,'width')
    
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoBkgSubUnc, cov_norm_BkgSub_np)
    unfoldHistoBkgSubUnc.Scale(1.,'width')
    
    
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    

    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    

    dataStatErrHist = unfoldHistoDataStatUnc.Clone('dataStatErrHist'+ivar)
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone('rmStatErrHist'+ivar)
    rmStatErrHist.Sumw2()
    #bkgSubErrHist = unfoldHistoBkgSubUnc.Clone('bkgSubErrHist'+ivar)
    #bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone('totalErrHist'+ivar)
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.056)
    totalErrHist.GetYaxis().SetLabelSize(0.047)    
    totalErrHist.GetYaxis().SetLabelOffset(0.012)    
    totalErrHist.GetYaxis().SetLabelFont(42)    
    totalErrHist.GetYaxis().SetTitleFont(42)
    totalErrHist.GetYaxis().SetTitleOffset(1.1)
    
    set_dynamic_y_range_errRatioHist(totalErrHist,1.35,0.95)
    
    if 'tau' in ivar: totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    else: totalErrHist.GetXaxis().SetTitle(labelX)
    
    totalErrHist.GetXaxis().SetTitleSize(0.056)
    totalErrHist.GetXaxis().SetLabelSize(0.047)    
    totalErrHist.GetXaxis().SetLabelOffset(0.011)    
    totalErrHist.GetXaxis().SetLabelFont(42)    
    totalErrHist.GetXaxis().SetTitleOffset(0.935)
    totalErrHist.GetXaxis().SetTitleFont(42)
    totalErrHist.GetXaxis().SetNdivisions(505)

    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    #totalErrHist.GetYaxis().SetNdivisions(510)

    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3345)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')

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
    
    
    #h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    #bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    #bkgSubErrHist.SetLineWidth(2)
    #h2.SetLineWidth(2)
    #bkgSubErrHist.SetLineStyle(7)
    #h2.SetLineStyle(7)
    #h2.SetLineColor(50)
    #bkgSubErrHist.SetLineColor(50)
    #h2.SetMarkerSize(0)
    #bkgSubErrHist.SetMarkerSize(0)   

    cr_histos = OrderedDict()

    CR1_key=None
    CR2_key=None
    erdOn_key=None

    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            #print (k, col_counter)
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            cr_histos[k].Scale(1./(cr_histos[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
            if 'cr1' in k.lower():
                CR1_key=k
            elif 'cr2' in k.lower():
                CR2_key=k
            elif 'erd' in k.lower():
                erdOn_key=k

            cr_histos[k].SetLineStyle(1)
            cr_histos[k].SetLineWidth(2)
            cr_histos[k].SetMarkerSize(0)
            cr_histos[k].SetFillColor(0)

    cr_keys = [CR1_key, CR2_key, erdOn_key]
    for i, key in enumerate(cr_keys):
        if key:
            color = colors_cr[i % len(colors_cr)]
            cr_histos[key].SetLineColor(color)
            cr_histos[key].SetMarkerColor(color)
            cr_histos[key].Draw('L same')

    if CR1_key: legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')
    if CR2_key: legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')
    if erdOn_key: legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')
    
   
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(normeduncerUnfoldHistoshiftsUp[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                
                up_counter=up_counter+1
    
    col_counter=0
    col_counter_jes=0
    
    syst_sources_up = list(normeduncerUnfoldHistoshiftsUp.keys())
    for i, k in enumerate(syst_sources_up):
        
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsUp[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(color)
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
            elif 'CP5' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'UE tune (CP5)', 'p' )
            elif 'mtop' in k:
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'm_{top}', 'p' )
            

            else: 
                print('else in th. syst comp maker', k)
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
                
    
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(normeduncerUnfoldHistoshiftsDown[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                 
                down_counter=down_counter+1
    
    syst_sources_down = list(normeduncerUnfoldHistoshiftsDown.keys())
    col_counter=0
    for i, k in enumerate(syst_sources_down):
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsDown[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(color)

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            print('else in th. syst comp maker', k)
            normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1} (13 TeV)"
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    canUnc.SaveAs(root_macro)
    
def drawUncertainties_from_err_shifts_unitNorm(ivar, 
                                               unfoldHistowoUnc,
                                               uncerUnfoldHisto, 
                                               cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, 
                                               labelX, tlegendAlignment, 
                                               outputName, year,
                                               selection, 
                                               lumi, with_modelUnc=True, norming=True ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    unftot = unfoldHistowoUnc.Integral()

    print (f'|> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('mod_vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
    colors = colors_list#[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
        
    upstyles =   [20,21,22,29,23,34,47,33,43, 117,114,48]  #39,41, 45,
    downstyles = [24,25,26,30,32,28,46,27,42, 38 ,60 , 5 ]  #37,40, 44,
    
        
    modelkey=None
    JES_key=None
    JER_key=None
    btag_key=None
    btagUncIncluded=False
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    
    
    W_ref = 800 #if square else 800
    H_ref = 600 #if square else 600
    extraSpace = 0.02
    W = W_ref
    H = H_ref
    T = 0.07 * H_ref
    B = 0.11 * H_ref
    L = 0.13 * H_ref
    R = 0.03 * H_ref

    canUnc = ROOT.TCanvas('canUnc_relUnc'+ivar, 'canUnc_relUnc'+ivar, 50, 50, W, H)
    canUnc.SetFillColor(0)
    canUnc.SetBorderMode(0)
    canUnc.SetFrameFillStyle(0)
    canUnc.SetFrameBorderMode(0)
    canUnc.SetLeftMargin(L / W + extraSpace)
    canUnc.SetRightMargin(R / W)
    #if with_z_axis:
    #    c.SetRightMargin(B / W + 0.03)
    canUnc.SetTopMargin(T / H)
    canUnc.SetBottomMargin(B / H + 0.02)

    legend=ROOT.TLegend(0.18,0.64,0.88,0.9)

    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.038)
    legend.SetTextFont(42)

    legend.SetBorderSize(0)
    
    
    unfoldHistoTotUnc = unfoldHistowoUnc.Clone('unfoldHistoTotUnc'+ivar)
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc = unfoldHistowoUnc.Clone('unfoldHistoDataStatUnc'+ivar)
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc = unfoldHistowoUnc.Clone('unfoldHistoRMStatUnc'+ivar)
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc = unfoldHistowoUnc.Clone('unfoldHistoBkgSubUnc'+ivar)
    unfoldHistoBkgSubUnc.Sumw2()
    
    
    cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHistoTotUnc, 
                                                       cov_tot.Clone())
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatUnc, 
                                                                   cov_datastat_tot.Clone())
    
    cov_norm_RMStat_np, normed_cov_RMStat = get_normalised_cov(unfoldHistoRMStatUnc, 
                                                               cov_rmstat_tot.Clone())
    cov_norm_BkgSub_np, normed_cov_BkgSub = get_normalised_cov(unfoldHistoBkgSubUnc, 
                                                               cov_bkg_tot.Clone())
    
    
    
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    
    
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoTotUnc, cov_normTot_np)
    unfoldHistoTotUnc.Scale(1.,'width')
    
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoDataStatUnc, cov_norm_dataStat_np)
    unfoldHistoDataStatUnc.Scale(1.,'width')
    
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoRMStatUnc, cov_norm_RMStat_np)
    unfoldHistoRMStatUnc.Scale(1.,'width')
    
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoBkgSubUnc, cov_norm_BkgSub_np)
    unfoldHistoBkgSubUnc.Scale(1.,'width')
    
    
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone('jesHistoUpMax')
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone('jesHistoDownMax')
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    if 'all' in year:
        jerHistoUpMax = unfoldHistoTotUnc.Clone('jerHistoUpMax')
        jerHistoUpMax.Reset()
        jerHistoDownMax = unfoldHistoTotUnc.Clone('jerHistoDownMax')
        jerHistoDownMax.Reset()

        jerHistoUpMax.Sumw2()
        jerHistoDownMax.Sumw2()
    
    
    #print(uncerUnfoldHisto.keys())
    for k in uncerUnfoldHisto:
        if 'modeltotal'in k.lower() and 'shifthist' in k.lower() and (modelkey==None) and with_modelUnc: 
            modelkey=k
        elif 'jes' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower() and not('const' in k.lower()) and (JES_key==None):
            JES_key=k
            #print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            jesHistoUpMax.Scale(1./(jesHistoUpMax.Integral() if norming else 1.),'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            jesHistoDownMax.Scale(1./(jesHistoDownMax.Integral() if norming else 1.),'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            #print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            jerHistoUpMax.Scale(1./(jerHistoUpMax.Integral() if norming else 1.),'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            jerHistoDownMax.Scale(1./(jerHistoDownMax.Integral() if norming else 1.),'width')
            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            if not('dijet' in selection):
                btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
                btagHistoUpMax.Reset()
                btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
                btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            btagHistoUpMax.Scale(1./(btagHistoUpMax.Integral() if norming else 1.),'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            btagHistoDownMax.Scale(1./(btagHistoDownMax.Integral() if norming else 1.),'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()) and not ('bkg' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #    print(k)
            #    btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(normeduncerUnfoldHistoshiftsUp[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColorAlpha(colors[col_counter], 1 if not('uncl' in text.lower()) else 0.7)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(1 if not('uncl' in text.lower()) else 0.8)
                if 'FSR' in text or 'ISR' in text:
                    normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(downstyles[up_counter])
                else:
                    normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                if 'tau_2_2' in k: print (k,text, up_counter, col_counter,upstyles[up_counter],colors[col_counter])
                col_counter=col_counter+1    
                up_counter=up_counter+1
            
    #up_counter=1
    #down_counter=1
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()) and not ('bkg' in k.lower()):
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
                
            if 'btag' in k.lower(): continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()  if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(normeduncerUnfoldHistoshiftsDown[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColorAlpha(colors[col_counter], 1 if not('uncl' in text.lower()) else 0.7)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(1 if not('uncl' in text.lower()) else 0.8)
                if 'FSR' in text or 'ISR' in text:
                    normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(upstyles[down_counter])
                else:
                    normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                    
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        #modelUnc = normalise_hist(modelUnc.Clone())
        modelUnc.Scale(1./(modelUnc.Integral() if norming else 1.),'width')#
        modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
        modelUnc.SetLineStyle(1)
        modelUnc.SetLineWidth(2)
        modelUnc.SetMarkerSize(0)
        modelUnc.SetLineColor(28)
        #col_counter+=1
        modelUnc.SetFillColor(0)
    
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    #bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    #bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.056)
    totalErrHist.GetYaxis().SetLabelSize(0.047)    
    totalErrHist.GetYaxis().SetLabelOffset(0.012)    
    totalErrHist.GetYaxis().SetLabelFont(42)    
    totalErrHist.GetYaxis().SetTitleFont(42)
    totalErrHist.GetYaxis().SetTitleOffset(1.1)
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.35,0.95)
    
    if 'tau' in ivar: totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    else: totalErrHist.GetXaxis().SetTitle(labelX)
        
    totalErrHist.GetXaxis().SetTitleSize(0.056)
    totalErrHist.GetXaxis().SetLabelSize(0.047)    
    totalErrHist.GetXaxis().SetLabelOffset(0.011)    
    totalErrHist.GetXaxis().SetLabelFont(42)    
    totalErrHist.GetXaxis().SetTitleOffset(0.935)
    totalErrHist.GetXaxis().SetTitleFont(42)
    totalErrHist.GetXaxis().SetNdivisions(505)
   
        
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3345)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')
    
    
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
    
    #h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    #bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    #bkgSubErrHist.SetLineWidth(2)
    #h2.SetLineWidth(2)
    #bkgSubErrHist.SetLineStyle(7)
    #h2.SetLineStyle(7)
    #h2.SetLineColor(50)
    #bkgSubErrHist.SetLineColor(50)
    #h2.SetMarkerSize(0)
    #bkgSubErrHist.SetMarkerSize(0)
    #bkgSubErrHist.Draw('L same ')
    #h2.Draw("L same")
    if with_modelUnc: modelUnc.Draw('L same')

    
    
    for k in otherUncs:
        if ('cr' in k.lower() or 'erd' in k.lower()): continue
        #print(k)
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print ("OtherUncs loop", text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        
        upmax_ibin = 1.
        downmax_ibin = 0.
        diff = 0.
        
        upmax_ibin = jesHistoUpMax.GetBinContent(ibin)
        diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
        
        if (diff>=1. or diff<0.):
            print(f'WARNING: JES total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
            upmax_ibin = 0.
            downmax_ibin = 0.
            diff = 0.
        else:
            downmax_ibin = 1.-diff
        #print ("JES total",ibin, upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if 'all' in year:

        for ibin in range(1,jerHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = jerHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: JER total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("JER total",ibin, upmax_ibin,downmax_ibin)
            jerHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            jerHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if not('dijet' in selection) and btagUncIncluded:

        for ibin in range(1,btagHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = btagHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: b-tagging total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            #print ("b-tagging total",ibin, upmax_ibin,downmax_ibin)
            btagHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            btagHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)
       
    col_counter=9
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(colors[col_counter])
    jesHistoDownMax.SetLineColor(colors[col_counter])
    jesHistoUpMax.SetMarkerColor(colors[col_counter])
    jesHistoDownMax.SetMarkerColor(colors[col_counter])
    jesHistoUpMax.SetMarkerStyle(39)#upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(37)#downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(1)
    jesHistoDownMax.SetMarkerSize(1)
    up_counter+=1
    down_counter+=1
    col_counter+=1

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(1)
        jerHistoDownMax.SetMarkerSize(1)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:# and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(1)
            btagHistoDownMax.SetMarkerSize(1)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #    legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if ('jes' in k.lower() and not('const' in k.lower())) or 'model'in k.lower() or 'tag' in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower() or ('all' in year and 'jer' in k.lower()): 
            continue
        else:
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'l1' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'L1 prefiring', 'p' )
            elif 'unclus' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'MET uncl. energy', 'p' )

            elif 'const' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'neut' in k.lower():
                    text="Neutral ES"
                elif 'charg' in k.lower():
                    text="Charged ES"
                elif 'photon' in k.lower():
                    text="Photon ES"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
            elif 'isr' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], "#alpha_{S}^{ISR}", 'p' )
                
            elif 'fsr' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], "#alpha_{S}^{FSR}", 'p' )

            else: 
                print('else in non-th. syst comp maker', k)
                
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'asandpdf' in k.lower():
                    text = "Scale & PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'PS and Hadr.', 'l' )    
    #legend.AddEntry( bkgSubErrHist, 'Bkg. stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'MC stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1} (13 TeV)"
        
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    canUnc.SaveAs(root_macro)

# deprecated, scaling a bin of the cov by bin width is applied uniformly to all covs, so cancels out in such ratios
def doRelUncPlot_noBW(ivar,
                 year,
                 lumi,
                 sel,
                 variables,
                 allHistos, uncerUnfoldSystCov,
                 outputDir,
                 version,
                 ext='pdf',
                 process='data'
                ):
    
    #Book histos for various summary covariances 
    
    covTot = allHistos['cov'+ivar].Clone()
    covTot.Reset()
    covTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    cov_expsystTot = covTot.Clone()
    cov_expsystTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))

    cov_modelsystTot = covTot.Clone()
    cov_modelsystTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))



    normed_uncerUnfold_systCov = OrderedDict()
    normed_uncerUnfold_systCov_np = OrderedDict()
    normed_unfHisto_withNormedCovErr = OrderedDict()
    

    for key in uncerUnfoldSystCov.keys():
        #print('normed_'+key)
        cov_norm, normed_uncerUnfold_systCov['normed_'+key.replace('173p5','Up').replace('171p5','Down')] = get_normalised_cov(allHistos['unfoldHisto'+ivar].Clone(),uncerUnfoldSystCov[key])

        normed_uncerUnfold_systCov_np['normed_'+key.replace('173p5','Up').replace('171p5','Down')] = cov_norm

        if 'jes' in key.lower() or 'jer' in key.lower() or 'puweight' in key.lower() or 'l1' in key.lower() or 'const' in key.lower() and not('model' in key.lower()):
            cov_expsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
            cov_expsystTot_np+=cov_norm

        elif 'isr' in key.lower() or 'fsr' in key.lower() or 'model' in key.lower():# or 'l1' in key.lower() or 'const' in key.lower():
            #print(key)
            cov_modelsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
            cov_modelsystTot_np+=cov_norm
            
        if not('dijet' in key.lower()):
            
            if 'unclust' in key.lower() or 'btag' in key.lower():
                cov_expsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
                cov_expsystTot_np+=cov_norm

            elif 'cr' in key.lower() or 'erd' in key.lower() or 'hdamp' in key.lower() or 'mtop' in key.lower() or ('tune' in key.lower() and ('up' in key.lower() or 'down' in key.lower())):# or 'l1' in key.lower() or 'const' in key.lower():
                #print(key)
                cov_modelsystTot.Add(normed_uncerUnfold_systCov['normed_'+key.replace('173p5','Up').replace('171p5','Down')])
                cov_modelsystTot_np+=cov_norm
            
            

        covTot.Add(normed_uncerUnfold_systCov['normed_'+key.replace("173p5","Up").replace("171p5","Down")])
        covTot_np+=cov_norm

        unfHisto_temp = allHistos['unfoldHistowoUnc'+ivar].Clone(f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err')

        unfHisto_temp.Scale(1./unfHisto_temp.Integral())
        
        #print(f'unfHisto_with_{key}Err')
        normed_unfHisto_withNormedCovErr[f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err'] = unfHisto_temp.Clone()

        get_th1_normedCovErrors(normed_unfHisto_withNormedCovErr[f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err'], cov_norm)
    
    for ih in  allHistos.keys():
        if ('cov_' in ih and not ('norm' in ih.lower())  and not ('sys' in ih.lower())) :
            #print(ih)
            cov_norm, normed_uncerUnfold_systCov['normed_'+ih] = get_normalised_cov(allHistos['unfoldHisto'+ivar].Clone(),allHistos[ih].Clone())
            covTot.Add(normed_uncerUnfold_systCov['normed_'+ih])
            normed_uncerUnfold_systCov_np['normed_'+ih] = cov_norm
            covTot_np+=cov_norm
            unfHisto_temp = allHistos['unfoldHistowoUnc'+ivar].Clone(f'unfHisto_with_{ih}Err')


            unfHisto_temp.Scale(1./unfHisto_temp.Integral())
            #print(f'unfHisto_with_{ih}Err')
            normed_unfHisto_withNormedCovErr[f'unfHisto_with_{ih}Err'] = unfHisto_temp.Clone()

            get_th1_normedCovErrors(normed_unfHisto_withNormedCovErr[f'unfHisto_with_{ih}Err'], cov_norm)
            
    relativeUncsFromCovsUp = OrderedDict()
    relativeUncsFromCovsDown = OrderedDict()
    relativeUncsFromTUnfCovs = OrderedDict()
    relativeUncsTotals = OrderedDict()


    for key in normed_uncerUnfold_systCov_np.keys():
        if not('up' in key.lower()): continue
        relativeUncsFromCovsUp['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  
        
    for key in normed_uncerUnfold_systCov_np.keys():
        if not('down' in key.lower()): continue
        relativeUncsFromCovsDown['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  

    for key in normed_uncerUnfold_systCov_np.keys():
        if ('up' in key.lower()) or ('down' in key.lower()): continue
        print(key)
        relativeUncsFromTUnfCovs['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  

    relativeUncsTotals = copy.deepcopy(relativeUncsFromTUnfCovs)
    
    unfs_done = []
    for key in normed_uncerUnfold_systCov_np.keys():
        if ('up' in key.lower()):

            up_varn = normed_uncerUnfold_systCov_np[key]
            down_varn = normed_uncerUnfold_systCov_np[key.replace('Up','Down').replace('up','down').replace('UP','DOWN')]
            #print(key, key.replace('Up','Down').replace('up','down').replace('UP','DOWN')) # ,'relTot'+key.replace('Up','total').replace('up','total').replace('UP','total'))

            total = np.sqrt(up_varn**2+down_varn**2)
            relativeUncsTotals['relTot'+key.replace('Up','total').replace('up','total').replace('UP','total')] = compute_relative_uncertainty(total,covTot_np)
            
    data_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_data_'+ivar]
    rel_data_stat_err = compute_relative_uncertainty(data_stat_err,covTot_np)

    mc_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_'+ivar]
    rel_mc_stat_err = compute_relative_uncertainty(mc_stat_err,covTot_np)

    bkg_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_bkg_'+ivar]
    rel_bkg_stat_err = compute_relative_uncertainty(bkg_stat_err,covTot_np)

    stat_plus_expsyst_err =  cov_expsystTot_np+data_stat_err+mc_stat_err+bkg_stat_err
    rel_stat_plus_expsyst_err =  compute_relative_uncertainty(stat_plus_expsyst_err, covTot_np)

    stat_plus_modelsyst_err =  cov_modelsystTot_np+data_stat_err
    rel_stat_plus_modelsyst_err =  compute_relative_uncertainty(stat_plus_modelsyst_err, covTot_np)
    
    modelsyst_err =  cov_modelsystTot_np
    rel_modelsyst_err =  compute_relative_uncertainty(modelsyst_err, covTot_np)
    
    # handling systematics with multiple sub-sources
    jesCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    jerCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    if not('dijet' in sel): 
        btagCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
        CRCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))

    for key in normed_uncerUnfold_systCov_np.keys():
        if 'jes' in key.lower() and not( 'const' in key.lower()):
            jesCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'jer' in key.lower():
            jerCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'btag' in key.lower():
            btagCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'cr1' in key.lower() or 'cr2' in key.lower() or 'erd' in key.lower():
            CRCov_np+=normed_uncerUnfold_systCov_np[key]
            
    rel_jesTot = compute_relative_uncertainty(jesCov_np,covTot_np)
    rel_jerTot = compute_relative_uncertainty(jerCov_np,covTot_np)
    
    rel_unc_dict = OrderedDict()

    #rel_unc_dict['Total'] = rel_jesTot
    rel_unc_dict['JES'] = rel_jesTot
    rel_unc_dict['JER'] = rel_jerTot
    rel_unc_dict['Data stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_data_'+ivar]
    rel_unc_dict['MC stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_'+ivar]
    rel_unc_dict['Bkg stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_bkg_'+ivar]
    rel_unc_dict['Pileup'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_puWeighttotal']
    if 'tau' in ivar: 
        rel_unc_dict['Neutral ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_neutraltotal']
        rel_unc_dict['Charged ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_chargedtotal']
        rel_unc_dict['Photon ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_photontotal']
        
    rel_unc_dict['L1 prefire'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_l1prefiringWeighttotal']
    if not ('dijet' in sel):
        rel_btagTot = compute_relative_uncertainty(btagCov_np, covTot_np)
        
        rel_unc_dict['MET uncl. en.'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_unclustEntotal']
        rel_unc_dict['b-tagging'] = rel_btagTot
    
    make_rel_uncertainty_plot(
                                ivar,
                                year,
                                lumi,
                                selection=sel,
                                dummy_unf_histo=allHistos['unfoldHistowoUnc'+ivar].Clone('dummyUnfHisto'),
                                rel_unc_dict=rel_unc_dict,
                                total_unc=rel_stat_plus_expsyst_err,  
                                total_unc_label="Stat. #oplus Exp.",
                                outfilename=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_RelUNC_exp_'+version+'.'+ext,
                                canvas_title="",
                                x_axis_title=variables[''+ivar]['label'],
                                y_axis_title="Relative uncertainty [%]",
                                y_max=100.0  
    )
    
    rel_unc_dict = OrderedDict()

    rel_unc_dict['ISR'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_isrWeighttotal']
    rel_unc_dict['FSR'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_fsrWeighttotal']
    rel_unc_dict['Data stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_data_'+ivar]
    rel_unc_dict['Shower and had.'] = relativeUncsFromTUnfCovs[f'relnormed_systcov_{ivar}_Physics ModelTotal']
    rel_unc_dict['PDF and #alpha_{S}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_aSandPDFWeighttotal']
    if not('dijet' in sel):
        rel_CRTot = compute_relative_uncertainty(CRCov_np, covTot_np)
        rel_unc_dict['UE tune (CP5)'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_TuneCP5total']
        rel_unc_dict['CR model'] = rel_CRTot#relativeUncsFromTUnfCovs[f'relnormed_systcov_{ivar}_Physics ModelTotal']
        rel_unc_dict['Choice of #m_{top}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_mtoptotal']
        rel_unc_dict['Choice of h_{damp}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_hdamptotal']
        

    
    make_rel_uncertainty_plot(
                                ivar,
                                year,
                                lumi,
                                selection=sel,
                                dummy_unf_histo=allHistos['unfoldHistowoUnc'+ivar].Clone('dummyUnfHisto'),
                                rel_unc_dict=rel_unc_dict,
                                total_unc=rel_stat_plus_modelsyst_err,  
                                total_unc_label="Stat. #oplus Model",
                                outfilename=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_RelUNC_model_'+version+'.'+ext,
                                canvas_title="",
                                x_axis_title=variables[''+ivar]['label'],
                                y_axis_title="Relative uncertainty [%]",
                                y_max=100.
    )
    
    
def doRelUncPlot(ivar,
                 year,
                 lumi,
                 sel,
                 variables,
                 allHistos, uncerUnfoldSystCov,
                 outputDir,
                 version,
                 ext='pdf',
                 process='data'
                ):
    
    #Book histos for various summary covariances 
    
    covTot = allHistos['cov'+ivar].Clone()
    covTot.Reset()
    covTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    cov_expsystTot = covTot.Clone()
    cov_expsystTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))

    cov_modelsystTot = covTot.Clone()
    cov_modelsystTot_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))



    normed_uncerUnfold_systCov = OrderedDict()
    normed_uncerUnfold_systCov_np = OrderedDict()
    normed_unfHisto_withNormedCovErr = OrderedDict()
    

    for key in uncerUnfoldSystCov.keys():
        #print('normed_'+key)
        cov_norm, normed_uncerUnfold_systCov['normed_'+key.replace('173p5','Up').replace('171p5','Down')] = get_normalised_cov_with_binwidth(allHistos['unfoldHisto'+ivar].Clone(),uncerUnfoldSystCov[key])

        normed_uncerUnfold_systCov_np['normed_'+key.replace('173p5','Up').replace('171p5','Down')] = cov_norm

        if 'jes' in key.lower() or 'jer' in key.lower() or 'puweight' in key.lower() or 'l1' in key.lower() or 'const' in key.lower() and not('model' in key.lower()):
            cov_expsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
            cov_expsystTot_np+=cov_norm

        elif 'isr' in key.lower() or 'fsr' in key.lower() or 'model' in key.lower():# or 'l1' in key.lower() or 'const' in key.lower():
            #print(key)
            cov_modelsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
            cov_modelsystTot_np+=cov_norm
            
        if not('dijet' in key.lower()):
            
            if 'unclust' in key.lower() or 'btag' in key.lower():
                cov_expsystTot.Add(normed_uncerUnfold_systCov['normed_'+key])
                cov_expsystTot_np+=cov_norm

            elif 'cr' in key.lower() or 'erd' in key.lower() or 'hdamp' in key.lower() or 'mtop' in key.lower() or ('tune' in key.lower() and ('up' in key.lower() or 'down' in key.lower())):# or 'l1' in key.lower() or 'const' in key.lower():
                print(key)
                cov_modelsystTot.Add(normed_uncerUnfold_systCov['normed_'+key.replace('173p5','Up').replace('171p5','Down')])
                cov_modelsystTot_np+=cov_norm
            
            

        covTot.Add(normed_uncerUnfold_systCov['normed_'+key.replace("173p5","Up").replace("171p5","Down")])
        covTot_np+=cov_norm

        unfHisto_temp = allHistos['unfoldHistowoUnc'+ivar].Clone(f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err')

        unfHisto_temp.Scale(1./unfHisto_temp.Integral(), 'width')
        
        #print(f'unfHisto_with_{key}Err')
        normed_unfHisto_withNormedCovErr[f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err'] = unfHisto_temp.Clone()

        get_th1_normedCovErrors(normed_unfHisto_withNormedCovErr[f'unfHisto_with_{key.replace("173p5","Up").replace("171p5","Down")}Err'], cov_norm)
    
    for ih in  allHistos.keys():
        if ('cov_' in ih and not ('norm' in ih.lower())  and not ('sys' in ih.lower())) :
            #print(ih)
            cov_norm, normed_uncerUnfold_systCov['normed_'+ih] = get_normalised_cov_with_binwidth(allHistos['unfoldHisto'+ivar].Clone(),allHistos[ih].Clone())
            covTot.Add(normed_uncerUnfold_systCov['normed_'+ih])
            normed_uncerUnfold_systCov_np['normed_'+ih] = cov_norm
            covTot_np+=cov_norm
            unfHisto_temp = allHistos['unfoldHistowoUnc'+ivar].Clone(f'unfHisto_with_{ih}Err')


            unfHisto_temp.Scale(1./unfHisto_temp.Integral(), 'width')
            #print(f'unfHisto_with_{ih}Err')
            normed_unfHisto_withNormedCovErr[f'unfHisto_with_{ih}Err'] = unfHisto_temp.Clone()

            get_th1_normedCovErrors(normed_unfHisto_withNormedCovErr[f'unfHisto_with_{ih}Err'], cov_norm)
            
    relativeUncsFromCovsUp = OrderedDict()
    relativeUncsFromCovsDown = OrderedDict()
    relativeUncsFromTUnfCovs = OrderedDict()
    relativeUncsTotals = OrderedDict()


    for key in normed_uncerUnfold_systCov_np.keys():
        if not('up' in key.lower()): continue
        relativeUncsFromCovsUp['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  
        
    for key in normed_uncerUnfold_systCov_np.keys():
        if not('down' in key.lower()): continue
        relativeUncsFromCovsDown['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  

    for key in normed_uncerUnfold_systCov_np.keys():
        if ('up' in key.lower()) or ('down' in key.lower()): continue
        #print(key)
        relativeUncsFromTUnfCovs['rel'+key] = compute_relative_uncertainty(normed_uncerUnfold_systCov_np[key],covTot_np)  

    relativeUncsTotals = copy.deepcopy(relativeUncsFromTUnfCovs)
    
    unfs_done = []
    for key in normed_uncerUnfold_systCov_np.keys():
        if ('up' in key.lower()):

            up_varn = normed_uncerUnfold_systCov_np[key]
            down_varn = normed_uncerUnfold_systCov_np[key.replace('Up','Down').replace('up','down').replace('UP','DOWN')]
            #print(key, key.replace('Up','Down').replace('up','down').replace('UP','DOWN')) # ,'relTot'+key.replace('Up','total').replace('up','total').replace('UP','total'))

            total = np.sqrt(up_varn**2+down_varn**2)
            relativeUncsTotals['relTot'+key.replace('Up','total').replace('up','total').replace('UP','total')] = compute_relative_uncertainty(total,covTot_np)
            
    data_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_data_'+ivar]
    rel_data_stat_err = compute_relative_uncertainty(data_stat_err,covTot_np)

    mc_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_'+ivar]
    rel_mc_stat_err = compute_relative_uncertainty(mc_stat_err,covTot_np)

    bkg_stat_err = normed_uncerUnfold_systCov_np['normed_cov_uncorr_bkg_'+ivar]
    rel_bkg_stat_err = compute_relative_uncertainty(bkg_stat_err,covTot_np)

    stat_plus_expsyst_err =  cov_expsystTot_np+data_stat_err+mc_stat_err+bkg_stat_err
    rel_stat_plus_expsyst_err =  compute_relative_uncertainty(stat_plus_expsyst_err, covTot_np)

    stat_plus_modelsyst_err =  cov_modelsystTot_np+data_stat_err
    rel_stat_plus_modelsyst_err =  compute_relative_uncertainty(stat_plus_modelsyst_err, covTot_np)
    
    modelsyst_err =  cov_modelsystTot_np
    rel_modelsyst_err =  compute_relative_uncertainty(modelsyst_err, covTot_np)
    
    # handling systematics with multiple sub-sources
    jesCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    jerCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
    if not('dijet' in sel): 
        btagCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))
        CRCov_np = np.zeros((covTot.GetNbinsX(),covTot.GetNbinsX()))

    for key in normed_uncerUnfold_systCov_np.keys():
        if 'jes' in key.lower() and not( 'const' in key.lower()):
            jesCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'jer' in key.lower():
            jerCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'btag' in key.lower():
            btagCov_np+=normed_uncerUnfold_systCov_np[key]
        elif 'cr1' in key.lower() or 'cr2' in key.lower() or 'erd' in key.lower():
            CRCov_np+=normed_uncerUnfold_systCov_np[key]
            
    rel_jesTot = compute_relative_uncertainty(jesCov_np,covTot_np)
    rel_jerTot = compute_relative_uncertainty(jerCov_np,covTot_np)
    
    rel_unc_dict = OrderedDict()

    #rel_unc_dict['Total'] = rel_jesTot
    rel_unc_dict['JES'] = rel_jesTot
    rel_unc_dict['JER'] = rel_jerTot
    rel_unc_dict['Statistical'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_data_'+ivar]
    rel_unc_dict['MC stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_'+ivar]
    rel_unc_dict['Bkg stat.'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_bkg_'+ivar]
    rel_unc_dict['Pileup'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_puWeighttotal']
    if 'tau' in ivar: 
        rel_unc_dict['Neutral ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_neutraltotal']
        rel_unc_dict['Charged ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_chargedtotal']
        rel_unc_dict['Photon ES'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_constituentJES_photontotal']
        
    rel_unc_dict['L1 prefire'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_l1prefiringWeighttotal']
    if not ('dijet' in sel):
        rel_btagTot = compute_relative_uncertainty(btagCov_np, covTot_np)
        
        rel_unc_dict['MET uncl. en.'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_unclustEntotal']
        rel_unc_dict['b-tagging'] = rel_btagTot
    
    make_rel_uncertainty_plot(
                                ivar,
                                year,
                                lumi,
                                selection=sel,
                                dummy_unf_histo=allHistos['unfoldHistowoUnc'+ivar].Clone('dummyUnfHisto'),
                                rel_unc_dict=rel_unc_dict,
                                total_unc=rel_stat_plus_expsyst_err,  
                                total_unc_label="Stat #oplus Exp.",
                                outfilename=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_RelUNC_wBW_exp_'+version+'.'+ext,
                                canvas_title="",
                                x_axis_title=variables[''+ivar]['label'],
                                y_axis_title="Relative uncertainty [%]",
                                y_max=100.0  
    )
    
    rel_unc_dict = OrderedDict()

    rel_unc_dict['ISR'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_isrWeighttotal']
    rel_unc_dict['FSR'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_fsrWeighttotal']
    rel_unc_dict['Statistical'] = relativeUncsFromTUnfCovs['relnormed_cov_uncorr_data_'+ivar]
    rel_unc_dict['Shower and hadronisation'] = relativeUncsFromTUnfCovs[f'relnormed_systcov_{ivar}_Physics ModelTotal']
    rel_unc_dict['PDF and #alpha_{S}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_aSandPDFWeighttotal']
    if not('dijet' in sel):
        rel_CRTot = compute_relative_uncertainty(CRCov_np, covTot_np)
        rel_unc_dict['UE tune (CP5)'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_TuneCP5total']
        rel_unc_dict['CR model'] = rel_CRTot#relativeUncsFromTUnfCovs[f'relnormed_systcov_{ivar}_Physics ModelTotal']
        rel_unc_dict['Choice of #m_{top}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_mtoptotal']
        rel_unc_dict['Choice of h_{damp}'] = relativeUncsTotals[f'relTotnormed_systcov_{ivar}_hdamptotal']
        

    
    make_rel_uncertainty_plot(
                                ivar,
                                year,
                                lumi,
                                selection=sel,
                                dummy_unf_histo=allHistos['unfoldHistowoUnc'+ivar].Clone('dummyUnfHisto'),
                                rel_unc_dict=rel_unc_dict,
                                total_unc=rel_stat_plus_modelsyst_err,  
                                total_unc_label="Stat #oplus Model",
                                outfilename=outputDir+ivar+sel+'_from'+('Data' if process.startswith('data') else 'MC')+'_Tunfold_RelUNC_wBW_model_'+version+'.'+ext,
                                canvas_title="",
                                x_axis_title=variables[''+ivar]['label'],
                                y_axis_title="Relative uncertainty [%]",
                                y_max=100.
    )
    
    
def make_rel_uncertainty_plot(ivar,
                              year,
                              lumi,
                              selection,
                              dummy_unf_histo,
                              rel_unc_dict,
                              total_unc=None,
                              total_unc_label="Total",
                              outfilename="relative_uncertainties.pdf",
                              canvas_title="Relative Uncertainties",
                              x_axis_title="Jet Observable",#nSubVariables_dijetSel['Jet_tau_0p5_2']['label'],
                              y_axis_title="Relative uncertainty [%]",
                              y_max=100.0):
    """
    Produce a plot of relative uncertainties vs. bin center, 
    optionally include a band for the total unc. and lines 
    for individual sources.

    Parameters
    -
    
    rel_unc_dict : dict
        Dictionary of { "label": rel_unc_array }, where rel_unc_array has length Nbins
        and each element is the relative uncertainty for that bin. 
        These will be drawn as lines on the plot.
    total_unc : 1D array_like, optional
        If provided, length = Nbins, this will be drawn as a gray band 
        representing the total uncertainty. 
    outfilename : str, optional
        PDF file to save the plot into.
    canvas_title : str, optional
        Title for the top of the canvas.
    x_axis_title : str, optional
        X-axis label.
    y_axis_title : str, optional
        Y-axis label.
    y_max : float, optional
        Maximum on the y-axis for plotting. Adjust as needed.
    """
    #bin_edges = np.array(bin_edges, dtype=float)
    
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
    
    n_bins = dummy_unf_histo.GetNbinsX()
    
    n_bins = dummy_unf_histo.GetNbinsX()
    if total_unc is not None and len(total_unc) != n_bins:
        raise ValueError("Length of total_unc array does not match number of bins in hist_binning!")
    
    for label, arr in rel_unc_dict.items():
        if len(arr) != n_bins:
            raise ValueError(f"Length of array for '{label}' does not match number of bins in hist_binning!")

    
    
    
    W_ref = 600 #if square else 800
    H_ref = 600 #if square else 600
    extraSpace = 0.02
    W = W_ref
    H = H_ref
    T = 0.07 * H_ref
    B = 0.11 * H_ref
    L = 0.13 * H_ref
    R = 0.03 * H_ref

    c = ROOT.TCanvas('c_relUnc'+ivar, 'c_relUnc'+ivar, 50, 50, W, H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(L / W + extraSpace)
    c.SetRightMargin(R / W)
    #if with_z_axis:
    #    c.SetRightMargin(B / W + 0.03)
    c.SetTopMargin(T / H)
    c.SetBottomMargin(B / H + 0.02)
    
    
    # Draw frame and set axis labels
    #h = canv.DrawFrame(x_min, y_min, x_max, y_max)

    #if yTitOffset is None:
    #    y_offset = 1.0 if square else 0.78
    #else:
    #    y_offset = yTitOffset
    
    #c.SetMargin(0.13,0.03,0.12,0.07)  # left, right, bottom, top
    #CMS.SetCmsTextFont(52)
    #CMS.SetCmsTextSize(0.75*0.76)
    
        
    x_min = dummy_unf_histo.GetBinLowEdge(1)
    x_max = dummy_unf_histo.GetBinLowEdge(dummy_unf_histo.GetNbinsX()+1) 
    y_min = 0. 
    y_max = max([y_max]+[10+i*100. for i in total_unc])
    if 'tau' in ivar: 
        x_axis_title = '#'+x_axis_title.split('#')[1] 
        
        
    
    
    frame_histo = dummy_unf_histo.Clone("frameHisto"+ivar)#ROOT.TH1D("frame_histo", "", n_bins, bin_edges)
    frame_histo.Reset("ICE")
    frame_histo.SetTitle("")
    frame_histo.GetXaxis().SetTitle(x_axis_title)   
    frame_histo.GetXaxis().SetTitleFont(42)
    frame_histo.GetXaxis().SetTitleOffset(0.95)
    
    frame_histo.GetYaxis().SetTitle(y_axis_title)
    frame_histo.GetYaxis().SetTitleFont(42)
    frame_histo.GetYaxis().SetTitleOffset( 1.25 )
    frame_histo.GetYaxis().SetTitleSize(0.06-0.004)
    frame_histo.GetXaxis().SetTitleSize(0.06-0.004)
    
    frame_histo.GetYaxis().SetLabelSize(0.05-0.002)
    frame_histo.GetXaxis().SetLabelSize(0.05-0.002)
    
    frame_histo.GetYaxis().SetRangeUser(0.0, y_max)
    frame_histo.GetXaxis().SetRangeUser(x_min, x_max)
    frame_histo.Draw("AXIS")


    legend = ROOT.TLegend(0.65, 0.50, 0.90, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.036 if not('Model' in total_unc_label) else 0.030)

    total_band = None
    if total_unc is not None:
        x_vals = []
        y_vals = []
        ex_vals = []
        ey_vals = []
        
        total_band = frame_histo.Clone("total_unc_band")
        total_band.Reset("ICE")
        
        for i in range(1, n_bins+1):
            val_percent = total_unc[i-1]*100.0
            total_band.SetBinContent(i, val_percent)
            

        '''
                    ROOT.TGraphErrors(n_bins, 
                                       np.array(x_vals, dtype=float), 
                                       np.array(y_vals, dtype=float),
                                       np.array(ex_vals, dtype=float),
                                       np.array(ey_vals, dtype=float))
        '''
        total_band.SetFillColor(17)
        
        total_band.SetFillColorAlpha(17, 0.6)
        total_band.SetLineColor(17)
        total_band.SetMarkerColor(17)
        
        total_band.SetLineWidth(1)
        #total_band.SetFillStyle(3254)  
        total_band.Draw("hist same")      

        legend.AddEntry(total_band, total_unc_label, "f")

    #colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta+1,
    #          ROOT.kOrange+1, ROOT.kAzure+2, ROOT.kTeal+1, ROOT.kViolet+1, ROOT.kGray+2,
    #          ROOT.kMagenta              
    #         ]
    
    colors = get_colour_palette_as_list()
    colors.insert(2,ROOT.kBlack)
    styles = [1, 2, 3, 4, 5, 6, 7, 9, 10, 2, 3]  

    graphs = []
    color_index = 0
    style_index = 0

    for label, unc_array in rel_unc_dict.items():
        h = frame_histo.Clone(label.replace(' ',''))
        
        for i in range(1, n_bins+1):
            h.SetBinContent(i,unc_array[i-1]* 100.0)
            h.SetBinError(i,0.000001)
            

        h.SetMarkerStyle(1)
        h.SetMarkerSize(0)
        h.SetMarkerColor(colors[color_index % len(colors)])
        h.SetLineColor(colors[color_index % len(colors)])
        h.SetLineStyle(styles[style_index % len(styles)])
        h.SetLineWidth(2)

        color_index += 1
        style_index += 1

        h.Draw("E1 same")
        legend.AddEntry(h, label, "l")
        graphs.append(h)

    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(c, 4, 0)
    
    c.Update()
    
    legend.Draw()
    
    c.SaveAs(outfilename)
    
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12) 
    

def drawUncertainties_from_err_shifts_theoryVariations(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, 
                                                       cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, norming=True ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    colors_cr = [ROOT.kMagenta+2,  ROOT.kBlue-4, 433]  
    colors_syst = get_colour_palette_as_list('vf_8')[1:]#list(reversed())#[ROOT.kBlue+1, ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kGreen+2]
    
    #colors = get_colour_palette_as_list('vf_10')
    #[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    #ROOT.gStyle.SetPalette(len(colors),array('i', colors))
        
    upstyles =   [20,21,34,29,22,23,29,47,33,43,39,41,39,45,117,114,48]
    downstyles = [24,25,28,30,26,32,30,46,27,42,37,40,37,44,38 ,60 , 5 ]
    
    modelVariations = True if 'VariationUNC' in outputName else False
    
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    #otherUncs = OrderedDict()
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)
    
    #if tlegendAlignment.startswith('right'): 
    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoNoNorm.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    

    up_counter=0#1
    down_counter=0#1
    col_counter=0#1
    col_counter_jes=0
    

    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    #bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    #bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    '''
    if not('dijet' in selection): 
        totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    else:
        if 'all' in year:
            totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        
        if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
    '''
    set_dynamic_y_range_errRatioHist(totalErrHist,1.3)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')

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
    #rmStatErrHist.Draw('L same ')
    #h1.Draw("L same")
    #h.Delete()
    
    #h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    #bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    #bkgSubErrHist.SetLineWidth(2)
    #h2.SetLineWidth(2)
    #bkgSubErrHist.SetLineStyle(7)
    #h2.SetLineStyle(7)
    #h2.SetLineColor(50)
    #bkgSubErrHist.SetLineColor(50)
    #h2.SetMarkerSize(0)
    #bkgSubErrHist.SetMarkerSize(0)
    #bkgSubErrHist.Draw('L same ')
    #h2.Draw("L same")
    

    cr_histos = OrderedDict()

    CR1_key=None
    CR2_key=None
    erdOn_key=None
    #set_palette_from_list(colors_cr)

    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            #print (k, col_counter)
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            #cr_histos[k] = normalise_hist(cr_histos[k].Clone())
            cr_histos[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
            if 'cr1' in k.lower():
                CR1_key=k
            elif 'cr2' in k.lower():
                CR2_key=k
            elif 'erd' in k.lower():
                erdOn_key=k

            cr_histos[k].SetLineStyle(1)
            cr_histos[k].SetLineWidth(2)
            cr_histos[k].SetMarkerSize(0)
            #cr_histos[k].SetLineColor(colors[col_counter])
            #col_counter+=1
            cr_histos[k].SetFillColor(0)

    cr_keys = [CR1_key, CR2_key, erdOn_key]
    for i, key in enumerate(cr_keys):
        if key:
            color = colors_cr[i % len(colors_cr)]
            cr_histos[key].SetLineColor(color)
            cr_histos[key].SetMarkerColor(color)
            cr_histos[key].Draw('L same')

    if CR1_key: legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')
    if CR2_key: legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')
    if erdOn_key: legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')
    
    #cr_histos[CR1_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')

    #cr_histos[CR2_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')

    #cr_histos[erdOn_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')

   
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                up_counter=up_counter+1
    
    col_counter=0
    col_counter_jes=0
    
    syst_sources_up = list(normeduncerUnfoldHistoshiftsUp.keys())
    for i, k in enumerate(syst_sources_up):
        
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsUp[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(color)
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
            elif 'CP5' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'UE tune (CP5)', 'p' )
            elif 'mtop' in k:
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'm_{top}', 'p' )
            

            else: 
                print('else in th. syst comp maker', k)
                if 'asandpdf' in k.lower():
                    text = "#alpha_{S} & PDF wt."
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
                
    
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                down_counter=down_counter+1
    
    #set_palette_from_list(colors_syst)
    syst_sources_down = list(normeduncerUnfoldHistoshiftsDown.keys())
    col_counter=0
    for i, k in enumerate(syst_sources_down):
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsDown[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(color)

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k or 'CP5' in k or 'mtop' in k:
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")

            else: 
                print('else in th. syst comp maker', k)
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    #legend.AddEntry( bkgSubErrHist, 'Bkg. stat.', 'l' )    
    #legend.AddEntry( rmStatErrHist, 'MC stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
def drawUncertainties_from_err_shifts(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, with_modelUnc=True, norming=False ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('mod_vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
    colors = colors_list#[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
        
    upstyles =   [20,21,22,29,23,34,47,33,43, 117,114,48]  #39,41, 45,
    downstyles = [24,25,26,30,32,28,46,27,42, 38 ,60 , 5 ]  #37,40, 44,
    
        
    modelkey=None
    JES_key=None
    JER_key=None
    btag_key=None
    btagUncIncluded=False
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)

    #if tlegendAlignment.startswith('right'): 
    #    legend=ROOT.TLegend(0.2,0.65,0.8,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone('jesHistoUpMax')
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone('jesHistoDownMax')
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    if 'all' in year:
        jerHistoUpMax = unfoldHistoTotUnc.Clone('jerHistoUpMax')
        jerHistoUpMax.Reset()
        jerHistoDownMax = unfoldHistoTotUnc.Clone('jerHistoDownMax')
        jerHistoDownMax.Reset()

        jerHistoUpMax.Sumw2()
        jerHistoDownMax.Sumw2()
    
    
    #print(uncerUnfoldHisto.keys())
    for k in uncerUnfoldHisto:
        if 'modeltotal'in k.lower() and 'shifthist' in k.lower() and (modelkey==None) and with_modelUnc: 
            modelkey=k
        elif 'jes' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower() and not('const' in k.lower()) and (JES_key==None):
            JES_key=k
            print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            jesHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            jesHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            jerHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            jerHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            if not('dijet' in selection):
                btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
                btagHistoUpMax.Reset()
                btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
                btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            btagHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            btagHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #    print(k)
            #    btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(1)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                if 'tau_2_2' in k: print (k,text, up_counter, col_counter,upstyles[up_counter],colors[col_counter])
                col_counter=col_counter+1    
                up_counter=up_counter+1
            
    #up_counter=1
    #down_counter=1
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
                
            if 'btag' in k.lower(): continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()  if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(1)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        #modelUnc = normalise_hist(modelUnc.Clone())
        modelUnc.Scale(1./(unftot if norming else 1.),'width')#
        modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
        modelUnc.SetLineStyle(1)
        modelUnc.SetLineWidth(2)
        modelUnc.SetMarkerSize(0)
        modelUnc.SetLineColor(28)
        #col_counter+=1
        modelUnc.SetFillColor(0)
    
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    #bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    #bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    #if not('dijet' in selection): 
    #    totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    #else:
    #    if 'all' in year:
    #        totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    
    #    if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
   
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.4,0.95)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')
    
    
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
    
    #h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    #bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    #bkgSubErrHist.SetLineWidth(2)
    #h2.SetLineWidth(2)
    #bkgSubErrHist.SetLineStyle(7)
    #h2.SetLineStyle(7)
    #h2.SetLineColor(50)
    #bkgSubErrHist.SetLineColor(50)
    #h2.SetMarkerSize(0)
    #bkgSubErrHist.SetMarkerSize(0)
    #bkgSubErrHist.Draw('L same ')
    #h2.Draw("L same")
    if with_modelUnc: modelUnc.Draw('L same')

    
    
    for k in otherUncs:
        if ('cr' in k.lower() or 'erd' in k.lower()): continue
        #print(k)
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print ("OtherUncs loop", text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        
        upmax_ibin = 1.
        downmax_ibin = 0.
        diff = 0.
        
        upmax_ibin = jesHistoUpMax.GetBinContent(ibin)
        diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
        
        if (diff>=1. or diff<0.):
            print(f'WARNING: JES total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
            upmax_ibin = 0.
            downmax_ibin = 0.
            diff = 0.
        else:
            downmax_ibin = 1.-diff
        print ("JES total",ibin, upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if 'all' in year:

        for ibin in range(1,jerHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = jerHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: JER total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("JER total",ibin, upmax_ibin,downmax_ibin)
            jerHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            jerHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if not('dijet' in selection) and btagUncIncluded:

        for ibin in range(1,btagHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = btagHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: b-tagging total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("b-tagging total",ibin, upmax_ibin,downmax_ibin)
            btagHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            btagHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)
       
    col_counter=9
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(colors[col_counter])
    jesHistoDownMax.SetLineColor(colors[col_counter])
    jesHistoUpMax.SetMarkerColor(colors[col_counter])
    jesHistoDownMax.SetMarkerColor(colors[col_counter])
    jesHistoUpMax.SetMarkerStyle(39)#upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(37)#downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(1)
    jesHistoDownMax.SetMarkerSize(1)
    up_counter+=1
    down_counter+=1
    col_counter+=1

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(1)
        jerHistoDownMax.SetMarkerSize(1)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:# and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(1)
            btagHistoDownMax.SetMarkerSize(1)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #    legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if ('jes' in k.lower() and not('const' in k.lower())) or 'model'in k.lower() or 'tag' in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower() or ('all' in year and 'jer' in k.lower()): 
            continue
        else:
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'l1' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'L1 prefiring', 'p' )
            elif 'unclus' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'MET uncl. en.', 'p' )

            elif 'const' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'neut' in k.lower():
                    text="Neutral ES"
                elif 'charg' in k.lower():
                    text="Charged ES"
                elif 'photon' in k.lower():
                    text="Photon ES"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )


            else: 
                print('else in non-th. syst comp maker', k)
                
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'asandpdf' in k.lower():
                    text = "#alpha_{S} & PDF wt."
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'PS and Hadr.', 'l' )    
    #legend.AddEntry( bkgSubErrHist, 'Bkg. stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'MC stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
def plotSysComparison2( nomHisto, dictUncHistos, outputName, labelX='', 
                       log=False, version='', ext='pdf', year='2017', 
                       outputDir='Plots/',
                       #sys_pref = 'jes',
                       mode = 'onlyJES', sysList=[],pngToo=False,
                      ):
    """Draw comparison of input systematic uncertainties/variations vs. nominal"""
    colors = [ 95, 38, 6, 7, 8, 42, 50, 218, 225, 30, 16,51, 61, 67, 89, 133, 142, 208, 36, 2, 144, 225, 198, 190, 83, 167, 207, 209, 212, 216,  227, 150, 93, 40]

    if 'jes' in mode.lower():
        jes_uncorr_list = [
                                '_jesAbsolute_2016', '_jesBBEC1_2016', '_jesEC2_2016', '_jesHF_2016', '_jesRelativeSample_2016',
                                '_jesAbsolute_2017', '_jesBBEC1_2017', '_jesEC2_2017', '_jesHF_2017', '_jesRelativeSample_2017',
                                '_jesAbsolute_2018', '_jesBBEC1_2018', '_jesEC2_2018', '_jesHF_2018', '_jesRelativeSample_2018'
                              ] 
    elif 'jer' in mode.lower() and year=='all':
        jes_uncorr_list = [
                                '_jer_2016', '_jer_2017', '_jer_2018'
                            ] 
        print(sysList)
    jes_corr_list = ['_jesAbsolute', '_jesBBEC1', '_jesEC2', '_jesFlavorQCD', '_jesHF', '_jesRelativeBal']
    
    if '18' in year and '_jesHEMIssue' in sysList: 
        jes_corr_list.append('_jesHEMIssue')
        
    elif 'all' in year and '_jesHEMIssue_2018' in sysList and 'jes' in mode.lower(): 
        jes_uncorr_list.append('_jesHEMIssue_2018')
        #print(sysList)
        #print(jes_uncorr_list)
        
    
    
    outputFileName = outputName+'_'+version+'.'+ext
    print ('Processing plots for sys comparisons......', outputFileName)
    
    legend=ROOT.TLegend(0.15,0.7,0.9,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.022)
    legend.SetBorderSize(0)

    multiGraph = ROOT.TMultiGraph()
    nomHisto.Sumw2()
    
    
    nomHisto=nomHisto.Clone()
    #print("Before scaling, integral:", nomHisto.Integral())

    nomHisto.Scale(1.,'width')
    
    #print("After normalising and scaling to b.w., nom integral:", nomHisto.Integral())
    
    gnom = ROOT.TGraphAsymmErrors() #nomHisto.Clone()
    gnom.Divide(nomHisto,nomHisto,'pois')
    
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetLineWidth(3)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    multiGraph.Add( gnom )
    
    
    dictShifts = {}
    col_counter=0
    colUp_counter=0
    colDown_counter=0
    dictGraphs = {}
    #print(dictUncHistos.keys())
    for ih in dictUncHistos.keys():
        
        dictUncHistos[ih].Sumw2()
        if (('jes' in mode.lower()) and not('jes' in ih.lower())) or (('jer' in mode.lower()) and not('jer' in ih.lower())):
            continue
        elif not('jes' in mode.lower()) and ('jes' in ih.lower()) and not('const' in ih.lower()):
            continue
        #print(ih)
        #print("Before scaling, integral:", ih,dictUncHistos[ih].Integral())
        dictShifts[ih] = dictUncHistos[ih].Clone()

        dictShifts[ih].Scale(1.,'width')
        dictGraphs[ih] = ROOT.TGraphAsymmErrors()
        dictGraphs[ih].Divide(dictShifts[ih],nomHisto,'pois')
        
        if 'model' in 'ih' or 'erd' in ih or 'mtop' in ih or 'CR' in ih:
            col_counter = colUp_counter
            colUp_counter+=1
            colDown_counter+=1
            dictGraphs[ih].SetLineStyle( 2 )#dictShifts[ih]
            
        else:
            if 'up' in ih.lower():
                col_counter = colUp_counter
                dictGraphs[ih].SetLineStyle( 3 if not(('isr' in ih) or ('fsr' in ih)) else 2)#dictShifts[ih]
                #colors= colors_up
                colUp_counter+=1
                
            elif 'down' in ih.lower():
                col_counter = colDown_counter
                dictGraphs[ih].SetLineStyle( 2 if not(('isr' in ih) or ('fsr' in ih)) else 3)#dictShifts[ih]
                #colors = colors_down
                colDown_counter+=1

        if not col_counter==len(colors)-1:
            dictGraphs[ih].SetLineColor( colors[col_counter] )#dictShifts[ih]
        else:
            col_counter = col_counter-len(colors)+2
            if 'up' in ih.lower(): 
                colUp_counter = col_counter-len(colors)+2
            if 'down' in ih.lower(): 
                colDown_counter = col_counter-len(colors)+2
            dictGraphs[ih].SetLineColor( colors[col_counter] )#dictShifts[ih]


        dictGraphs[ih].SetMarkerStyle(0)
        dictGraphs[ih].SetMarkerSize(0)
        dictGraphs[ih].SetLineWidth( 1 )
        
        y=year if not('all' in year or '+' in year) else ''
        #print(ih)
        
        stringtocheck=ih.split('_' )[1] if not('const' in ih) else (ih.split('_' )[1]+ih.split('_' )[2]).replace('JES', ' ES ')#if (year=='all') else ih#+'_'+(year if not('+' in year) else '_fullRunII')
            
        if ('jes' in ih and not('const' in ih)) or ( 'jer' in ih and year=='all'):
            flagUncorr=False
            for j in jes_uncorr_list:
                #print(j)
                if j in ih:
                    stringtocheck=ih[1:]
                    #print("Sys. comp plot for uncorr sources, working on:", j,ih,stringtocheck)
                    flagUncorr=True
                    break
            
        if not('isr' in stringtocheck or 'fsr' in stringtocheck): 
            legend.AddEntry( dictGraphs[ih], stringtocheck,'l' )#+'_'+y
        else: 
            stringtocheck = stringtocheck.replace('Down','Up')  if 'Down' in stringtocheck else stringtocheck.replace('Up','Down')
            legend.AddEntry( dictGraphs[ih], stringtocheck,'l' )#+'_'+y

        multiGraph.Add(dictGraphs[ih])
        
        
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    
    multiGraph.GetYaxis().SetTitle( 'Variation/Nominal' )
    multiGraph.GetXaxis().SetTitle( labelX )
    if not('jes' in mode.lower() or 'jer' in mode.lower() or 'const' in outputFileName.lower()):
        multiGraph.SetMaximum( 1.7 )
        multiGraph.SetMinimum( 0.5 )
    
    else:
        if not('const' in outputFileName.lower() or 'jes' in mode.lower()):
            multiGraph.SetMaximum( 1.10 )
            multiGraph.SetMinimum( 0.96 )
        elif 'const' in outputFileName.lower():
            multiGraph.SetMaximum( 1.1 )
            multiGraph.SetMinimum( 0.95 )
        else:
            
            multiGraph.SetMaximum( 1.06 )
            multiGraph.SetMinimum( 0.97 )
        
        
    
    multiGraph.Draw('ALP')

    
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    canUnc.Update()
    #canUnc.cd()
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    if 'pdf' in outputFileName and pngToo: canUnc.SaveAs( outputDir + outputFileName.replace('png', 'pdf') )
    del canUnc

    
    


def makePSplot_simple(purity,stability,
                      accepGen,#=None,
                      fakeReco,#=None,
                      variables,
                      lumi,
                      var,outputDir, ext,
                      dictHistos=OrderedDict(),
                      
                      bins=[0.,1.],
                      year='2017',
                      #ext,#='pdf',
                      sel='_WSel',
                      signalLabelBegin='TTToSemiLeptonic'):

    if not os.path.exists(outputDir): os.makedirs(outputDir)

    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc")]
    extraSpace = 0.01
    W_ref = 600 #if square else 800
    H_ref = 600 #if square else 600
    extraSpace = 0.02
    W = W_ref
    H = H_ref
    T = 0.07 * H_ref
    B = 0.11 * H_ref
    L = 0.13 * H_ref
    R = 0.03 * H_ref

    canvas = ROOT.TCanvas('canvasPurity_Stability'+var+year, 'canvasPurity_Stability'+var+year,  50, 50, W, H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(L / W + extraSpace)
    canvas.SetRightMargin(R / W)
    #if with_z_axis:
    #    c.SetRightMargin(B / W + 0.03)
    canvas.SetTopMargin(T / H)
    canvas.SetBottomMargin(B / H + 0.02)
    
    
        
    #canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)
    p=ROOT.TH1D("Purity"+var+year,";;",len(bins)-1,array('d',bins))
    s=ROOT.TH1D("Stability"+var+year,";;",len(bins)-1,array('d',bins))
    a=accepGen.Clone()
    f=fakeReco.Clone()
    for i in range(len(purity)):
        p.SetBinContent(i+1,purity[i])
        s.SetBinContent(i+1,stability[i])
    
    x_min = p.GetBinLowEdge(1)
    x_max = p.GetBinLowEdge(p.GetNbinsX()+1) 
    y_min = 0. 
    y_max = 1.
    
    if 'tau' in var: 
        x_axis_title = '#'+variables[var]['label'].split('#')[1] 
    else:
        x_axis_title = variables[var]['label']
    print(x_axis_title)
    
    #legend.SetNColumns( 4 )
    frame_histo = p.Clone("purityFrame"+var+year)#ROOT.TH1D("frame_histo", "", n_bins, bin_edges)
    frame_histo.Reset("ICE")
    frame_histo.SetTitle("")
    frame_histo.GetXaxis().SetTitle(x_axis_title)   
    frame_histo.GetXaxis().SetTitleFont(42)
    frame_histo.GetXaxis().SetTitleOffset(0.90)
    
    frame_histo.GetYaxis().SetTitle("")
    frame_histo.GetYaxis().SetTitleFont(42)
    frame_histo.GetYaxis().SetTitleOffset( 1.25 )
    frame_histo.GetYaxis().SetTitleSize(0.06-0.006)
    frame_histo.GetXaxis().SetTitleSize(0.06-0.006)
    
    frame_histo.GetYaxis().SetLabelSize(0.05-0.004)
    frame_histo.GetXaxis().SetLabelSize(0.05-0.004)
    
    frame_histo.GetYaxis().SetRangeUser(0.0, 1.1)
    frame_histo.GetXaxis().SetRangeUser(x_min, x_max)
    frame_histo.Draw("AXIS")


    legend = ROOT.TLegend(0.55, 0.75, 0.90, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.032)# if not('Model' in total_unc_label) else 0.030)
    #legend.SetNColumns(3)



    p.SetMaximum( 1.2 )
    p.SetMinimum( 0. )
    
    p.SetTitle("")
    p.GetXaxis().SetTitle(x_axis_title)   
    p.GetXaxis().SetTitleFont(42)
    p.GetXaxis().SetTitleOffset(0.93)
    
    p.GetYaxis().SetTitle("")
    p.GetYaxis().SetTitleFont(42)
    #p.GetYaxis().SetTitleOffset( 1.25 )
    p.GetYaxis().SetTitleSize(0.06-0.006)
    p.GetXaxis().SetTitleSize(0.06-0.006)
    
    p.GetYaxis().SetLabelSize(0.05-0.004)
    p.GetXaxis().SetLabelSize(0.05-0.004)
    
    p.GetYaxis().SetLabelOffset(0.012)
    p.GetXaxis().SetLabelOffset(0.011)

    
    p.GetYaxis().SetRangeUser(0.0, 1.2)
    p.GetXaxis().SetRangeUser(x_min, x_max)
    p.GetXaxis().SetNdivisions(505)
    p.GetXaxis().SetNdivisions(505)
    
    p.SetLineColor( ROOT.kBlack )
    p.SetMarkerColor( ROOT.kBlack )
    p.SetMarkerSize( 0.5 )
    p.SetLineWidth( 2 )

    s.SetLineColor( ROOT.kMagenta )
    s.SetMarkerColor( ROOT.kMagenta )
    s.SetMarkerSize( 0.5 )
    s.SetLineWidth( 2 )
    #print(variables[var]['label'] )
    #p.GetXaxis().SetTitle( variables[var]['label'] )
    legend.AddEntry( p, 'Purity', 'l' )
    legend.AddEntry( s,  'Stability' , 'l' )
    p.Draw('hist')
    s.Draw('hist same')
    
    #if not('tau' in var):
    if not(a==None or f==None):
        a.SetLineStyle(2)
        f.SetLineStyle(2)
        a.SetLineWidth(1)
        f.SetLineWidth(1)
        a.SetLineColor(colors[0])
        f.SetLineColor(colors[1])

        legend.AddEntry( a,  'Acceptance (gen)' , 'l' )
        legend.AddEntry( f, 'Fake rate (reco)', 'l' )


        a.Draw('hist same')
        f.Draw('hist same')

    dictHistos[ 'purityGraph_'+var ] = p.Clone()
    dictHistos[ 'stabilityGraph_'+var ] = s.Clone()
    dictHistos[ 'acceptanceRateGraph_'+var ] = a.Clone()
    dictHistos[ 'fakeRateGraph_'+var ] = f.Clone()


    legend.Draw()
    CMS_lumi.extraText = "Simulation Preliminary"
    #CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('135' if 'dijet' in sel else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in sel else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(canvas, 4, 0)

    canvas.Update()
    canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+year+'.'+ext)
    if ext.startswith('pdf'):
        canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.png')
    return 1#numBins


def createCanvasPads():
    c = ROOT.TCanvas("c", "canvas", 800, 600)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()
 
    return c, pad1, pad2



##########################################################################
#################### Helpers for plot agglomeration ######################
##########################################################################


def png_to_pdf(png_files, output_file):
    images = []
    for file in png_files:
        img = Image.open(file)
        if img.mode == 'RGBA':
            img = img.convert('RGB')
        images.append(img)

    images[0].save(output_file, save_all=True, append_images=images[1:])

def merge_pdfs(pdf_files, output_file):
    merger = PdfMerger()
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(output_file)
    merger.close()

def combine_images_and_pdfs(input_files, output_file):
    png_files = [file for file in input_files if file.lower().endswith('.png')]
    pdf_files = [file for file in input_files if file.lower().endswith('.pdf')]

    temp_pdf_files = []
    if png_files:
        temp_pdf = "temp_images.pdf"
        png_to_pdf(png_files, temp_pdf)
        temp_pdf_files.append(temp_pdf)
    
    all_pdf_files = pdf_files + temp_pdf_files
    if all_pdf_files:
        merge_pdfs(all_pdf_files, output_file)

    for temp_file in temp_pdf_files:
        os.remove(temp_file)

##########################################################################
####################### Helpers for plots, etc. ##########################
##########################################################################

def set_dynamic_y_range(histList, extra_margin=1.3):
    """
    Given a list of histograms (or a single histogram), determine the min and max values
    to adjust the Y-axis range and require adding some buffer space above the max for legends.

    Parameters:
        histList (list or ROOT.TH1): List of histograms or a single histogram.
        extra_margin (float): Factor by which to multiply the max for spacing.

    Returns:
        None (the first histogram in histList will have its Y range updated)
    """
    if not isinstance(histList, list):
        histList = [histList]

    minVal = 1e9
    maxVal = -1e9

    for h in histList:
        for ibin in range(1, h.GetNbinsX()+1):
            val = h.GetBinContent(ibin)
            if val < minVal: minVal = val
            if val > maxVal: maxVal = val

    hMain = histList[0]
    hMain.GetYaxis().SetRangeUser(minVal * 0.95, maxVal * extra_margin)    
    
    
def set_dynamic_y_range_errRatioHist(histList, extra_margin=1.4, bottom_margin=0.9):
    """
    Given a list of error bar histograms (ie, hist w. y-errors divide by same hist with no y-errors), 
    determine the min and max values thereof to adjust the Y-axis range and require
    adding some buffer space above the max for legends.

    Parameters:
        histList (list or ROOT.TH1): List of histograms or a single histogram.
        extra_margin (float): Factor by which to multiply the max for blank y-spacing in 
        top part of plot.

    Returns:
        None ( first histogram in histList will have Y-range updated)
    """
    if not isinstance(histList, list):
        histList = [histList]

    minVal = 1e9
    maxVal = -1e9

    for h in histList:
        for ibin in range(1, h.GetNbinsX()+1):
            val = 1.-h.GetBinError(ibin)
            if val < minVal: 
                minVal = val
            val = 1.+h.GetBinError(ibin)
            if val > maxVal: 
                maxVal = val
            #print(val)

    hMain = histList[0]
    hMain.GetYaxis().SetRangeUser(minVal * bottom_margin, maxVal * extra_margin)        
    
def get_colour_palette_as_list(palette_requested=None):
    if palette_requested=='vf_10' or palette_requested==None:
        hex_list = ["#3f90da", "#ffa90e",  "#b9ac70", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#bd1f01", "#717581", "#92dadd"]
    elif palette_requested=='mod_vf_10':# or palette_requested==None:
        hex_list = ["#94a4a2", "#3f90da",  "#bd1f01", "#ffa90e", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#964a8b", "#92dadd"]
    elif palette_requested=='vf_8':
        hex_list = ["#1845fb", "#ff5e02", "#c91f16", "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"]
    elif palette_requested=='vf_6':
        hex_list = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
    else:
        hex_list = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
    colour_list = []
    for i in hex_list:
        colour_list.append(ROOT.TColor.GetColor(i))
    return colour_list


def set_palette_from_list(color_list):
    arr = array('i', color_list)
    ROOT.gStyle.SetPalette(len(color_list), arr)
    

def createRecoBins(genBins):
    recoBins=[genBins[0]]
    for i in range(1,len(genBins)):
        recoBins.append(np.round(genBins[i-1]+(genBins[i]-genBins[i-1])/2.,3))
        recoBins.append(genBins[i])
        #recoBins.append(b)
    return recoBins



#from https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/print_bottom_line_test.py#L112C1-L135C15
def get_null_bins(h):
    null_bins = []
    if isinstance(h, (ROOT.TH1, ROOT.TH2)):
        if isinstance(h, ROOT.TH2):
            h_proj = h.ProjectionX()
        else:
            h_proj = h
        for ix in range(1, h_proj.GetNbinsX()+1):
            if h_proj.GetBinContent(ix) == 0:
                null_bins.append(ix)
        return null_bins
    else:
        proj = h.sum(axis=0)
        return np.where(proj == 0)[0]


def remove_null_bins(arr, null_bins):
    if len(arr.shape) > 1:
        for ax in range(len(arr.shape)):
            if arr.shape[ax] > 1:
                arr = np.delete(arr, null_bins, axis=ax)
    else:
        arr = np.delete(arr, null_bins)
    return arr


def bottomLineTest( ivar, dataHisto, dataHistoLabel, MCHisto, covMatrix, varInfo, outputLabel, outputDir, rebin=1,
                    ext='png', version='Vn', selection='_dijet', process='data', no_null_bins=True):#, ignore_UF=True):
    #based on https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/my_unfolder.py#L2501
    
    print(f"Data nbins: {dataHisto.GetNbinsX()}",f"MC nbins: {MCHisto.GetNbinsX()}",f"cov nbinsX: {covMatrix.GetNbinsX()}",f"cov nbinsY: {covMatrix.GetNbinsY()}")
    if (dataHisto.GetNbinsX()!=MCHisto.GetNbinsX()):
        print("!!! ERROR:(dataHisto.GetNbinsX()!=MCHisto.GetNbinsX()) !!!" )
        dataBins=[]
        MCbins=[]
        for i in range(1,dataHisto.GetNbinsX()+2):#!=MCHisto.GetNbinsX()):
            dataBins.append(dataHisto.GetBinLowEdge(i))
        for i in range(1,MCHisto.GetNbinsX()+2):
            MCbins.append(MCHisto.GetBinLowEdge(i))
        
        #print("Data bins", dataBins )
        #print("MC bins", MCbins )
            
    #print(dataHisto.GetNbinsX(),MCHisto.GetNbinsX(),covMatrix.GetNbinsX(),covMatrix.GetNbinsY())
    
    
    if not(rebin==1):
        MCHisto.Rebin( rebin )  
        dataHisto.Rebin( rebin )
        covMatrix.Rebin2D( rebin, rebin )
    
    if isinstance(covMatrix,ROOT.TH2): 
        cov_arr,_ = th2_to_ndarray(covMatrix.Clone())
    else:
        cov_arr = covMatrix
    if isinstance(dataHisto,ROOT.TH1): 
        data_arr,_ = th1_to_ndarray(dataHisto.Clone())
    else:
        data_arr = dataHisto
    if isinstance(MCHisto,ROOT.TH1): 
        mc_arr,_ = th1_to_ndarray(MCHisto.Clone())
    else:
        mc_arr = MCHisto
    #print(cov_arr.shape,data_arr.shape,mc_arr.shape)
    # check symmetry
    #print("max asymmetry:", np.max(np.abs(cov_arr - cov_arr.T)))

    # compute eigenvalues
    eigvals = np.linalg.eigvalsh(cov_arr)  
    #print("cov eigenvalues:", eigvals)

    
    
    if no_null_bins:
        
        null_bins = get_null_bins(cov_arr)
        print(f'null bins {dataHistoLabel} space:', null_bins)
        
        cov_arr = remove_null_bins(cov_arr, null_bins)
        mc_arr = remove_null_bins(mc_arr, null_bins)
        data_arr = remove_null_bins(data_arr, null_bins)    
    
    
    try:
        delta = data_arr - mc_arr
        #print(f'Data b.c.: {data_arr}, MC b.c.: {mc_arr}, Residuals:{delta}')        
    except ValueError:
        print(f'Data b.c.: {data_arr}, MC b.c.: {mc_arr}')
        print("chi2 cannot be calculated since something is wrong with the data mc event th2->ndarrays, please check what's going on")
        return 1,1,1
    
    # check condition
    print("cond number:", np.linalg.cond(cov_arr))
    print("cov_arr shape:", cov_arr.shape)
    print("  data shape:", data_arr.shape)
    print("   MC  shape:",  mc_arr.shape)
    print("delta shape:", delta.shape)
    #print("first 5 δ entries:", delta.ravel()[:])
    #print("cov:", cov_arr[:,:])
    
    try: 
        v_inv = np.linalg.inv(cov_arr)
    except np.linalg.LinAlgError:
        print("Using pseudo-inverse instead since true inv. operation via np.linalg.inv() failed")
        v_inv = np.linalg.pinv(cov_arr)#, rcond=1E-30)
        
    inter = v_inv.dot(delta.T)
    chi2 = delta.dot(inter)[0][0]
    
    data_nonzero = len([i for i in range(data_arr[0].shape[0]) if data_arr[0][i]>0])
    mc_nonzero = len([i for i in range(mc_arr[0].shape[0]) if mc_arr[0][i]>0])
    
    
    ndof = max(data_nonzero,mc_nonzero)#delta.shape[1]# # only consider n bins where at least one has data - if both 0, don't count it
    print(1.-scipy.stats.chi2.cdf(chi2, int(ndof)))
    
    p = 1.-scipy.special.gammainc(chi2/2.,ndof/2.)#1.-scipy.stats.chi2.cdf(chi2, int(ndof))
        
    print(f'chi2 for {dataHistoLabel}, ndf, p, chi2/ndf = ', np.round(chi2,5), ndof, p, np.round(chi2/ndof,5))
    
    return chi2, ndof, p

def calc_chi2_stats(one_hist, other_hist, cov_matrix):
    one_vec = one_hist 
    other_vec = other_hist
    delta = one_vec - other_vec
    if isinstance(cov_matrix, ROOT.TH2):
        v, _ = th2_to_ndarray(cov_matrix)
    else:
        v = cov_matrix
    # print("delta:", delta)
    # v = np.diag(np.diag(v))  
    # print("v:", v)
    try:
        v_inv = np.linalg.inv(v)
    except np.linalg.LinAlgError:
        print("Trying pseudo-inverse instead")
        v_inv = np.linalg.pinv(v, rcond=1E-30)
    inter = v_inv.dot(delta.T)
    # print("parts:", delta * inter.T)
    chi2 = delta.dot(inter)[0][0]
    ndof = delta.shape[1]
    p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
    return chi2, ndof, p


def fold_generator_level(hist_truth, probability_matrix, bins_reco, oflow=False):
    
    # Convert ROOT TH1 to vector
    gen_vec, gen_vec_err = th1_to_ndarray(hist_truth, oflow_x=oflow)
    
    #convert TH2 to ndarray for PM
    probaM, _ = th2_to_ndarray(probability_matrix, oflow)

    # Multiply array (vector) with PM (transpose from row vec to column vec as necessary)
    folded_vec = probaM.dot(gen_vec.T)

    # Convert array to TH1
    folded_mc_truth = ndarray_to_th1(folded_vec.T, has_oflow_x=oflow, offset=0., bins=bins_reco)

    # Err. prop.: if y = Ax, with covariance matrices Vyy and Vxx,respectively, then Vyy = (A*Vxx)*A^T
    vxx, _ = th2_to_ndarray(make_diag_cov_hist_from_errors(hist_truth, inverse=False), oflow)
    result = probaM.dot(vxx)
    
    folded_covariance = result.dot(probaM.T)
    print(folded_covariance.shape)
    folded_errors = make_hist_from_diagonal_errors(folded_covariance,bins=bins_reco)
    update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_mc_truth)
    return folded_mc_truth


def DoUnfolding(Response,Reco):
    tunfolder = ROOT.TUnfoldDensity(Response,
                                    ROOT.TUnfold.kHistMapOutputHoriz,
                                    ROOT.TUnfold.kRegModeCurvature, 
                                    ROOT.TUnfold.kEConstraintNone, 
                                    ROOT.TUnfoldDensity.kDensityModeBinWidth)
    tunfolder.SetInput(Reco)
    tunfolder.DoUnfold(0.)
    return tunfolder.GetOutput("MC_unfolded")

def get_folded_unfolded(folded, unfolded, cov_tot, probaM, oflow=True):
    # don't use getfoldedoutput, because it doesn't have the updated errors from the total error matrix
    # so we'll have to do it ourselves
    # 1. Make unfolded hist into TVector/TMatrix

    # 2. Make response 2d hist into matrix

    # 3. Multiply the two, convert to TH1

    # Get the TUnfold one for reference, although its errors will be wrong
    
    #print(folded.GetNbinsX(),unfolded.GetNbinsX())
    
    
    #probability matrix is simply the normalised RM, 
    #so need to consider that there is an UF bin on the y-axis for miss(gen)-corrections
    probaM, _ = th2_to_ndarray(probaM, oflow_x=False, oflow_y=False, uflow_x=False, uflow_y=False)
    #print(probaM.shape, "with UF on y axis")
    
    folded_unfolded_tunfold = folded.Clone()

    # Get unfolded results as array
    unfolded_vector, _ = th1_to_ndarray(unfolded, oflow_x=False)
    #print(unfolded_vector.shape)
    
    # Multiply
    # first, transpose from row vec to column vec
    folded_vec = probaM.dot(unfolded_vector.T)
    #print(folded_vec.shape)
    
    
    bins = array('d',[folded.GetBinLowEdge(i) for i in range(1,folded.GetNbinsX()+2)])

    
    # Convert vector to TH1
    folded_unfolded = ndarray_to_th1(folded_vec.T, has_oflow_x=False, offset=0., bins=bins)

    # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
    # respectively, then Vyy = (A*Vxx)*A^T
    unfolded_covariance_matrix, _ = th2_to_ndarray(cov_tot , oflow_x=False, oflow_y=False, uflow_x=False, uflow_y=False )
    result = probaM.dot(unfolded_covariance_matrix)
    #print(result.shape)
    
    folded_covariance = result.dot(probaM.T)
    #print(folded_covariance.shape)
    
    
    folded_errors = make_hist_from_diagonal_errors(folded_covariance,bins=bins)
    
    #print(folded_errors.GetBinLowEdge(4))
    #print(folded_unfolded.GetBinLowEdge(4))
    
    update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_unfolded)

    return folded_unfolded.Clone(folded.GetName()+'errorFixed')

def CrossClosure(response1,reco1,response2,reco2):
    unf11=DoUnfolding(response1,reco1)
    unf12=DoUnfolding(response2,reco1)
    unf21=DoUnfolding(response1,reco2)
    unf22=DoUnfolding(response2,reco2)
    return unf11,unf12,unf21,unf22

def SelfClosure(response1,reco2,response2,reco1):
    unf21=DoUnfolding(response1,reco2)
    unf12=DoUnfolding(response2,reco1)
    return unf21,unf12


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





def generate_latex_table(observables, data, column_titles, filename=None):
    """
    Generates a LaTeX table based on the input data.

    Args:
    observables (list of str): List of row identifiers (first column, like "Observable").
    data (list of lists): List of lists containing the data for each row (each list represents a row).
    column_titles (list of str): List of column titles (like 'Observable', '$\chi^2$', 'ndf', '$\frac{\chi^2}{ndf}$').
    filename (str): If provided, saves the LaTeX table to a file with the given name. Otherwise, it prints it.

    Returns:
    str: LaTeX table code.
    """
    
    latex_table = "\\begin{table}[ht!]\n\\centering\n\\begin{tabular}{|" + "c|" * len(column_titles) + "}\n\\hline\n"
    
    latex_table += " & ".join(column_titles) + " \\\\\n\\hline\n"
    
    for i, observable in enumerate(observables):
        row_data = " & ".join(map(str, data[i]))
        latex_table += f"{observable} & {row_data} \\\\\n"
    
    latex_table += "\\hline\n\\end{tabular}\n\\caption{The table caption here.}\n\\label{tab:the_label}\n\\end{table}"

    if filename:
        with open(filename, 'w') as file:
            file.write(latex_table)
        print(f"LaTeX table saved to {filename}.")
    else:
        print(latex_table)
    
    return latex_table


def extendTH1(h, extendUF=True, extendOF=True):
    """
    Given a TH1 (e.g. TH1D), return a new TH1 with extra bins
    that include the underflow (UF) and/or overflow (OF) entries.
    
    The extra bin(s) will have the same width as the nominal first
    (for underflow) and last (for overflow) bins.
    
    Parameters:
      h         : The original TH1 histogram.
      extendUF  : If True, add an extra (leftmost) bin for underflow.
      extendOF  : If True, add an extra (rightmost) bin for overflow.
    
    Returns:
      A new TH1 histogram with the “extended” x‐axis.
    """
    # Number of nominal bins
    n = h.GetNbinsX()
    axis = h.GetXaxis()
    x_min = axis.GetXmin()  # lower edge of first nominal bin
    x_max = axis.GetXmax()  # upper edge of last nominal bin
    first_bin_width = axis.GetBinWidth(1)
    last_bin_width = axis.GetBinWidth(n)
    
    # Get the original bin edges.
    # If the histogram was created with nonuniform binning, GetXbins() returns
    # a TArrayD of bin edges (of size n+1). Otherwise it is empty.
    bins_array = axis.GetXbins()
    if bins_array.GetSize() > 0:
        # Non-uniform binning: get edges for bins 1..n+1.
        orig_edges = [axis.GetBinLowEdge(i) for i in range(1, n+2)]
    else:
        # Uniform binning.
        orig_edges = [x_min + i*(x_max - x_min)/n for i in range(0, n+1)]
    
    # Build the new edge array.
    new_edges = []
    if extendUF:
        new_edges.append(x_min - first_bin_width)
    new_edges.extend(orig_edges)
    if extendOF:
        new_edges.append(x_max + last_bin_width)
    new_edges_arr = array('d', new_edges)
    
    h_ext = ROOT.TH1D(h.GetName() + "_ext", h.GetTitle() + " (extended)", len(new_edges_arr)-1, new_edges_arr)
    h_ext.Sumw2()  # preserve Sumw2
    
    # Offset in bin numbering: if extended UF then new bin 1 is the UF bin.
    offset = 1 if extendUF else 0
    
    # Copy underflow if extended.
    if extendUF:
        h_ext.SetBinContent(1, h.GetBinContent(0))
        h_ext.SetBinError(1, h.GetBinError(0))
    else:
        h_ext.SetBinContent(0, h.GetBinContent(0))
        h_ext.SetBinError(0, h.GetBinError(0))
    
    # Copy the nominal bins.
    for i in range(1, n+1):
        new_bin = i + offset
        h_ext.SetBinContent(new_bin, h.GetBinContent(i))
        h_ext.SetBinError(new_bin, h.GetBinError(i))
    
    # Copy overflow if extended.
    if extendOF:
        new_bin = n + 1 + offset
        h_ext.SetBinContent(new_bin, h.GetBinContent(n+1))
        h_ext.SetBinError(new_bin, h.GetBinError(n+1))
    else:
        h_ext.SetBinContent(h_ext.GetNbinsX()+1, h.GetBinContent(n+1))
        h_ext.SetBinError(h_ext.GetNbinsX()+1, h.GetBinError(n+1))
    
    return h_ext
        


def extendTH2(h, extendUF_x=True, extendOF_x=True, extendUF_y=True, extendOF_y=True):
    """
    Create a new TH2D whose binning in x and y optionally extends underflow/overflow
    into the visible range. Corner bins (x-flow, y-flow) follow the chosen setting on
    each axis independently:
      - If an axis is extended for underflow/overflow, that flow bin becomes visible
        along that axis.
      - If an axis is NOT extended, that flow bin remains in underflow/overflow for
        that axis.
    Thus, a partially extended corner (e.g. x=overflow is extended, y=underflow is not)
    ends up in the new histogram at (x=lastVisibleBin, y=underflowBin=0), without merging
    into the first visible y bin.

    Parameters:
      h          : The original TH2 histogram (TH2D assumed).
      extendUF_x : If True, x underflow becomes the first visible bin in x.
      extendOF_x : If True, x overflow becomes the last visible bin in x.
      extendUF_y : If True, y underflow becomes the first visible bin in y.
      extendOF_y : If True, y overflow becomes the last visible bin in y.

    Returns:
      A new TH2D with the extended axes and contents/errors correctly placed, including
      partial-flow corners.
    """
    import math
    from array import array

    n_x = h.GetNbinsX()
    n_y = h.GetNbinsY()

    axisX = h.GetXaxis()
    axisY = h.GetYaxis()

    # Original bin edges (x)
    x_min = axisX.GetXmin()
    x_max = axisX.GetXmax()
    binsX = axisX.GetXbins()
    if binsX.GetSize() > 0:
        orig_edges_x = [axisX.GetBinLowEdge(i) for i in range(1, n_x+2)]
    else:
        # uniform binning
        orig_edges_x = [x_min + i*(x_max - x_min)/n_x for i in range(0, n_x+1)]
    first_bin_width_x = axisX.GetBinWidth(1)
    last_bin_width_x  = axisX.GetBinWidth(n_x)

    # Original bin edges (y)
    y_min = axisY.GetXmin()
    y_max = axisY.GetXmax()
    binsY = axisY.GetXbins()
    if binsY.GetSize() > 0:
        orig_edges_y = [axisY.GetBinLowEdge(i) for i in range(1, n_y+2)]
    else:
        # uniform binning
        orig_edges_y = [y_min + i*(y_max - y_min)/n_y for i in range(0, n_y+1)]
    first_bin_width_y = axisY.GetBinWidth(1)
    last_bin_width_y  = axisY.GetBinWidth(n_y)

    # Build the new bin-edge arrays
    new_edges_x = []
    if extendUF_x:
        new_edges_x.append(orig_edges_x[0] - first_bin_width_x)
    new_edges_x.extend(orig_edges_x)
    if extendOF_x:
        new_edges_x.append(orig_edges_x[-1] + last_bin_width_x)

    new_edges_y = []
    if extendUF_y:
        new_edges_y.append(orig_edges_y[0] - first_bin_width_y)
    new_edges_y.extend(orig_edges_y)
    if extendOF_y:
        new_edges_y.append(orig_edges_y[-1] + last_bin_width_y)

    new_edges_x_arr = array('d', new_edges_x)
    new_edges_y_arr = array('d', new_edges_y)

    # Number of visible bins in new histogram
    n_new_x = len(new_edges_x_arr) - 1
    n_new_y = len(new_edges_y_arr) - 1

    # Create new histogram
    h_ext = ROOT.TH2D(h.GetName()+"_ext", h.GetTitle()+" (extended)",
                      n_new_x, new_edges_x_arr,
                      n_new_y, new_edges_y_arr)
    h_ext.Sumw2()

    #
    # 1) Helper to map old bin index (0..n+1) -> new bin index
    #    (0..n_new+1). If that axis is extended for underflow,
    #    old underflow(0) -> new bin 1, else -> 0. Similarly for
    #    overflow (n+1). Nominal bins map to [1..n] or [2..n+1].
    #
    def map_axis_bin(old_bin, n, extendUF, extendOF):
        # old_bin can be 0..(n+1)
        # new_n = # of visible bins = n + (1 if UF extended) + (1 if OF extended)
        n_visible = n + (1 if extendUF else 0) + (1 if extendOF else 0)

        if old_bin == 0:   # underflow
            return 1 if extendUF else 0
        elif old_bin == n+1:  # overflow
            return n_visible if extendOF else (n_visible + 1)
        else:
            # nominal bin => shift by +1 if one extended the UF
            offset = 1 if extendUF else 0
            newb = old_bin + offset
            return newb

    #
    # 2) Fill all bins (including corners) in ONE pass.
    #    Loop over `old' bins [0..n_x+1, 0..n_y+1] and add them
    #    into the new histogram bin that corresponds.
    #    This automatically handles partial corners without merging.
    #
    for old_i in range(0, n_x+2):  # x in [0..n_x+1]
        for old_j in range(0, n_y+2):  # y in [0..n_y+1]
            c  = h.GetBinContent(old_i, old_j)
            ce = h.GetBinError(old_i, old_j)
            if c == 0 and ce == 0:
                continue  # skip empty to save a bit of time

            # map to new bin indices
            i_new = map_axis_bin(old_i, n_x, extendUF_x, extendOF_x)
            j_new = map_axis_bin(old_j, n_y, extendUF_y, extendOF_y)

            # accumulate (in case multiple old bins map to the same new bin)
            old_c  = h_ext.GetBinContent(i_new, j_new)
            old_ce = h_ext.GetBinError(i_new, j_new)
            new_c  = old_c + c
            new_ce = math.sqrt(old_ce*old_ce + ce*ce)
            h_ext.SetBinContent(i_new, j_new, new_c)
            h_ext.SetBinError(i_new, j_new, new_ce)

    return h_ext





def extendTH2_old(h, extendUF_x=True, extendOF_x=True, extendUF_y=True, extendOF_y=True):
    """
    Given a TH2 (e.g. TH2D), return a new TH2 with extra bins along the x- and/or y-axes
    that include the underflow and/or overflow entries.
    
    For each axis the extra bin has the same width as the first (for underflow)
    or last (for overflow) nominal bin.
    
    Parameters:
      h          : The original TH2 histogram.
      extendUF_x : If True, add an extra x-bin at the left for underflow.
      extendOF_x : If True, add an extra x-bin at the right for overflow.
      extendUF_y : If True, add an extra y-bin at the bottom for underflow.
      extendOF_y : If True, add an extra y-bin at the top for overflow.
    
    Returns:
      A new TH2 histogram with extended x- and y-axes.
    """
    # Process the X axis 
    n_x = h.GetNbinsX()
    axisX = h.GetXaxis()
    x_min = axisX.GetXmin()
    x_max = axisX.GetXmax()
    first_bin_width_x = axisX.GetBinWidth(1)
    last_bin_width_x = axisX.GetBinWidth(n_x)
    binsX = axisX.GetXbins()
    if binsX.GetSize() > 0:
        orig_edges_x = [axisX.GetBinLowEdge(i) for i in range(1, n_x+2)]
    else:
        orig_edges_x = [x_min + i*(x_max-x_min)/n_x for i in range(0, n_x+1)]
    new_edges_x = []
    if extendUF_x:
        new_edges_x.append(x_min - first_bin_width_x)
    new_edges_x.extend(orig_edges_x)
    if extendOF_x:
        new_edges_x.append(x_max + last_bin_width_x)
    new_edges_x_arr = array('d', new_edges_x)
    new_n_x = len(new_edges_x_arr) - 1

    #  Process the Y axis 
    n_y = h.GetNbinsY()
    axisY = h.GetYaxis()
    y_min = axisY.GetXmin()  # for Y axis, use GetXmin()/GetXmax() as well
    y_max = axisY.GetXmax()
    first_bin_width_y = axisY.GetBinWidth(1)
    last_bin_width_y = axisY.GetBinWidth(n_y)
    binsY = axisY.GetXbins()
    if binsY.GetSize() > 0:
        orig_edges_y = [axisY.GetBinLowEdge(i) for i in range(1, n_y+2)]
    else:
        orig_edges_y = [y_min + i*(y_max-y_min)/n_y for i in range(0, n_y+1)]
    new_edges_y = []
    if extendUF_y:
        new_edges_y.append(y_min - first_bin_width_y)
    new_edges_y.extend(orig_edges_y)
    if extendOF_y:
        new_edges_y.append(y_max + last_bin_width_y)
    new_edges_y_arr = array('d', new_edges_y)
    new_n_y = len(new_edges_y_arr) - 1

    # Create new TH2 histogram (here TH2D is assumed)
    h_ext = ROOT.TH2D(h.GetName() + "_ext", h.GetTitle() + " (extended)",
                      new_n_x, new_edges_x_arr,
                      new_n_y, new_edges_y_arr)
    h_ext.Sumw2()

    # Offsets in bin numbering: if underflow is extended then nominal bins start at new bin index 2.
    offset_x = 1 if extendUF_x else 0
    offset_y = 1 if extendUF_y else 0

    # Loop over the bins of the new histogram (which now includes the flow bins as actual bins)
    for i_new in range(1, new_n_x+1):
        # Map new bin i_new to original x bin index:
        if extendUF_x and i_new == 1:
            orig_i = 0      # original underflow bin in x
        elif extendOF_x and i_new == new_n_x:
            orig_i = n_x + 1  # original overflow bin in x
        else:
            orig_i = i_new - offset_x  # nominal bin (1..n_x)
        for j_new in range(1, new_n_y+1):
            if extendUF_y and j_new == 1:
                orig_j = 0      # original underflow bin in y
            elif extendOF_y and j_new == new_n_y:
                orig_j = n_y + 1  # original overflow bin in y
            else:
                orig_j = j_new - offset_y
            # Get the content and error from the original histogram.
            content = h.GetBinContent(orig_i, orig_j)
            err = h.GetBinError(orig_i, orig_j)
            new_bin = h_ext.GetBin(i_new, j_new)
            h_ext.SetBinContent(new_bin, content)
            h_ext.SetBinError(new_bin, err)
    return h_ext


def computeOOA(resp2D):
    """
    Computes per-bin out of acceptance correction factors of the form:
       (SB-only gen counts) / (nom+SB gen counts)
    by summing over all reco bins for each gen bin.
    
   
    """
    effCorr = []
    for i in range(0,resp2D.GetNbinsX()+1):
        # Numerator: integral over Y= SB region only
        # Denominator: integral over Y= nominal + SB region
        
        num = resp2D.Integral(i, i+1, resp2D.GetNbinsY(), resp2D.GetNbinsY()+1) # "SB" portion in Y 
        den = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY() + 1) # "nom + SB" portion

        if den != 0:
            effCorr.append(num / den)
        else:
            effCorr.append(1.0)  
    return effCorr




def computeAcceptance(resp2D):
    """
    Computes per-bin acceptance correction factors of the form:
       (nom-only gen counts) / (nom+SB gen counts)
    by summing over all reco bins for each gen bin.
    
   
    """
    effCorr = []
    for i in range(0,resp2D.GetNbinsX()+1):
        # Numerator: integral over Y= nominal region only
        # Denominator: integral over Y= nominal + SB region
        
        num = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY())     # "nom" portion
        den = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY() + 1) # "nom + SB" portion

        if den != 0:
            effCorr.append(num / den)
        else:
            effCorr.append(1.0)  
    return effCorr


def applyAcceptanceCorrection(hist1D, effCorr):
    """
    Scales each bin of 'hist1D' by the corresponding acceptance factor
    stored in 'effCorr'.  
    
    
    """
    
    for i in range(0, hist1D.GetNbinsX() + 2):
        scaleFactor = effCorr[i - 1]  
        val  = hist1D.GetBinContent(i)
        err  = hist1D.GetBinError(i)
        hist1D.SetBinContent(i, val * scaleFactor)
        hist1D.SetBinError(i,  err * scaleFactor)
        
        
def applyInvAcceptanceCorrection(hist1D, effCorr):
    """
    Scales each bin of 'hist1D' by the corresponding acceptance factor
    stored in 'effCorr'.  
    
    
    """
    
    for i in range(0, hist1D.GetNbinsX() + 2):
        scaleFactor = effCorr[i - 1]  
        val  = hist1D.GetBinContent(i)
        err  = hist1D.GetBinError(i)
        hist1D.SetBinContent(i, val * (1.-scaleFactor))
        hist1D.SetBinError(i,  err * (1.-scaleFactor))
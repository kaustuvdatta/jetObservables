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
from root_numpy import array2hist, hist2array
#import histoHelpers
from histoHelpers import *
import os
import glob
import sys
import math
#sys.path.insert(0,'../test/')
#from Drawhistogram_dijetSel import *
#from datasets_dijetSel_RunIISummer20UL_nomWts import dictSamples, checkDict
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
#import pdb
lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)

##########################################################################

#################### Helpers/plotters for unfoldings #####################

##########################################################################


###################### Unfolding plotting functions ######################

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

    
def drawDataMCReco( ivar, selection, year, lumi, process,
                    dataJetHisto, nominal_recoJetHisto, 
                    alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto,
                    labelX, maxX, tlegendAlignment, outputName ):
    """docstring for drawUnfold"""
    print ("Drawing Data/MC")
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    if 'pt' in ivar.lower():
        can.SetLogy()
        outputName.replace('DataMC','DataMC_logPlot')
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
    print(dataJetHisto,nominal_recoJetHisto,alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto)
    
    dataHisto = normalise_hist(dataJetHisto.Clone())
    recoHisto = normalise_hist(nominal_recoJetHisto.Clone())
    alt0recoHisto = normalise_hist(alt0_recoJetHisto.Clone())
    alt1recoHisto = normalise_hist(alt1_recoJetHisto.Clone())
    alt2recoHisto = normalise_hist(alt2_recoJetHisto.Clone())
    
    dataHisto.Scale(1, 'width')  ### divide by bin width
    #dataHisto.Scale(1/dataHisto.Integral(), 'width')  ### divide by bin width
    dataHisto.SetMarkerStyle(8)
    dataHisto.SetMarkerSize(2)
    dataHisto.SetMarkerColor(ROOT.kBlack)
    dataHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( dataHisto, 'Data', 'pe' )
    
    recoHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    #genJetHisto.Scale(1/genJetHisto.Integral(), 'width')  ### divide by bin width
    recoHisto.SetLineWidth(2)
    recoHisto.SetLineColor(ROOT.kRed)
    recoHisto.SetMarkerColor(ROOT.kRed)
    recoHisto.SetMarkerStyle(25)
    legend.AddEntry( recoHisto, 'MG5-MLM+Pythia8', 'lp' )

    print(labelX)

    try: 
        dataHisto.GetYaxis().SetTitle( '#frac{1}{dN} #frac{dN}{d#'+labelX.split('#')[1]+'}' )
    except IndexError: 
        try: 
            #print('#frac{1}{dN} #frac{dN}{d'+labelX.split('AK8 jet ')[1]+'}') 
            dataHisto.GetYaxis().SetTitle( '#frac{1}{dN} #frac{dN}{d'+labelX.split('AK8 jet ')[1]+'}' )
        except IndexError:
            #print('#frac{1}{dN} #frac{dN}{d'+labelX.split('good_')[1]+'}') 
            dataHisto.GetYaxis().SetTitle('#frac{1}{dN} #frac{dN}{d'+labelX.split('_')[1]+'}' )
    print(dataHisto.GetYaxis().GetTitle())
    #dataHisto.GetYaxis().SetTitleOffset(0.95)
    dataHisto.GetYaxis().SetTitleSize(0.05)
    dataHisto.SetMaximum( 1.6*max([ recoHisto.GetMaximum(), dataHisto.GetMaximum()] )  )
    dataHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), dataHisto.GetMaximum()] ) )

    dataHisto.Draw( "E")
    recoHisto.Draw( "histe same")

    alt0recoHisto.Scale(1, 'width')  ### divide by bin width
    alt0recoHisto.SetLineWidth(2)
    alt0recoHisto.SetLineColor(ROOT.kBlue)
    alt0recoHisto.SetMarkerColor(ROOT.kBlue)
    alt0recoHisto.SetMarkerStyle(25)
    legend.AddEntry( alt0recoHisto, 'MG5-MLM+Herwig7', 'lp' )
    alt0recoHisto.Draw("histe same")
    
    alt1recoHisto.Scale(1, 'width')  ### divide by bin width
    alt1recoHisto.SetLineWidth(2)
    alt1recoHisto.SetLineColor(ROOT.kViolet)
    alt1recoHisto.SetMarkerColor(ROOT.kViolet)
    alt1recoHisto.SetMarkerStyle(25)
    legend.AddEntry( alt1recoHisto, 'MG5+Pythia8', 'lp' )
    alt1recoHisto.Draw("histe same")
    
    
    alt2recoHisto.Scale(1, 'width')  ### divide by bin width
    alt2recoHisto.SetLineWidth(2)
    alt2recoHisto.SetLineColor(ROOT.kGray+4)
    alt2recoHisto.SetMarkerColor(ROOT.kGray+4)
    alt2recoHisto.SetMarkerStyle(25)
    legend.AddEntry( alt2recoHisto, 'Pythia8', 'lp' )
    alt2recoHisto.Draw("histe same")
    
    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.045)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
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
    
    
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    #print (labelX)
    tmpPad2.GetYaxis().SetTitle( "Sim./Data" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(-0.5, 2.5 )
    tmpPad2.GetYaxis().CenterTitle()
    try: tmpPad2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    except IndexError: tmpPad2.GetXaxis().SetTitle(labelX)
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
    
    
    ratio_nominal = ROOT.TGraphAsymmErrors()
    ratio_nominal.Divide( recoHisto, dataHisto, 'pois' )
    ratio_nominal.SetLineColor(ROOT.kRed)
    ratio_nominal.SetMarkerColor(ROOT.kRed)
    try:  ratio_nominal.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    except IndexError: ratio_nominal.GetXaxis().SetTitle( labelX )
    #ratio_nominal.GetYaxis().SetTitle( "Sim./Data" )
    #ratio_nominal.GetYaxis().SetTitleOffset( 0.5 )
    #ratio_nominal.GetYaxis().SetRangeUser(0., 2. )
    #ratio_nominal.GetYaxis().CenterTitle()
    ratio_nominal.GetXaxis().SetLabelSize(0.12)
    ratio_nominal.GetXaxis().SetTitleSize(0.12)
    #ratio_nominal.GetYaxis().SetLabelSize(0.12)
    #ratio_nominal.GetYaxis().SetTitleSize(0.12)
    ratio_nominal.GetXaxis().SetNdivisions(505)
    #ratio_nominal.GetYaxis().SetNdivisions(505)
    ratio_nominal.SetMarkerStyle(25)
    ratio_nominal.Draw('P0 same')
    
    ratio_alt0MC = ROOT.TGraphAsymmErrors()
    ratio_alt0MC.Divide( alt0recoHisto, dataHisto, 'pois' )
    ratio_alt0MC.SetLineColor(ROOT.kBlue)
    ratio_alt0MC.SetMarkerColor(ROOT.kBlue)
    ratio_alt0MC.SetMarkerStyle(25)
    
    ratio_alt0MC.Draw('P0 same')
    
    ratio_alt1MC = ROOT.TGraphAsymmErrors()
    ratio_alt1MC.Divide( alt1recoHisto, dataHisto, 'pois' )
    ratio_alt1MC.SetLineColor(ROOT.kViolet)
    ratio_alt1MC.SetMarkerColor(ROOT.kViolet)
    ratio_alt1MC.SetMarkerStyle(25)
    
    ratio_alt1MC.Draw('P0 same')
    
    ratio_alt2MC = ROOT.TGraphAsymmErrors()
    ratio_alt2MC.Divide( alt2recoHisto, dataHisto, 'pois' )
    ratio_alt2MC.SetLineColor(ROOT.kGray+4)
    ratio_alt2MC.SetMarkerColor(ROOT.kGray+4)
    ratio_alt2MC.SetMarkerStyle(25)
    
    ratio_alt2MC.Draw('P0 same')
    
    ratioLegend=ROOT.TLegend(0.20,0.85,0.8,0.95)
    ratioLegend.SetTextSize(0.06)
    ratioLegend.SetNColumns(4)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_nominal, 'MG5-MLM+P8', 'lp' )
    ratioLegend.AddEntry( ratio_alt0MC, 'MG5-MLM+H7', 'lp' )
    ratioLegend.AddEntry( ratio_alt1MC, 'MG5+P8', 'lp' )
    ratioLegend.AddEntry( ratio_alt2MC, 'P8+P8', 'lp' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)
                    

def drawUnfold( ivar, selection, process, year, lumi, dataJetHisto, genJetHisto, unfoldHisto, unfoldHistowoUnc, altMCHisto,
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
    
    #bins = variables[ivar]['bins']

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
    legend.AddEntry( genJetHisto, 'MG5-MLM+Pythia8', 'lp' )

   

    unfoldHisto.GetYaxis().SetTitle( '#frac{1}{dN} #frac{dN}{d#'+labelX.split('#')[1]+'}' )
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
    legend.AddEntry( altMCHisto, 'MG5-MLM+Herwig7', 'lp' )
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

    if selection.startswith("_dijet"): seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
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
    tmpPad2.GetYaxis().SetRangeUser(0., 2. )
    #tmpPad2.GetXaxis().SetRangeUser(unfoldHisto.GetBinLowEdge(1),unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2) )   
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
    ratio_datastatUnc.GetYaxis().SetRangeUser(0.25, 2.25 )
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
    
def drawClosures(ivar, selection, process, year, lumi, genJetHisto, genJetHistoCross, unfoldHisto, unfoldHistoCross,
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
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.55,0.67,0.85,0.9)

    else: legend=ROOT.TLegend(0.20,0.67,0.40,0.9)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
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
    legend.AddEntry( unfoldHisto, ('MG5-MLM+Pythia8 (self-closure)' if process.startswith('MCSelfClosure') else 'MG5-MLM+P8 unf. w/ MG5-MLM+P8'), 'pe' )
    
    #genJetHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    genJetHisto.Scale(1/genJetHisto.Integral(), 'width')  ### divide by bin width
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(ROOT.kBlue)
    genJetHisto.SetMarkerStyle(0)
    genJetHisto.SetLineStyle(2)
    legend.AddEntry( genJetHisto, 'MG5-MLM+Pythia8 (gen)', 'lp' )
    unfoldHisto.GetYaxis().SetTitle( '#frac{1}{dN} #frac{dN}{d#'+labelX.split('#')[1]+'}' )
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
        legend.AddEntry( unfoldHistoCross, 'MG5-MLM+P8 unf. w/ MG5-MLM+H7', 'pe' )

        #genJetHisto.Scale(1, 'width')
        #genJetHisto.Scale(scaleFactor)
        genJetHistoCross.Scale(1/genJetHistoCross.Integral(), 'width')  ### divide by bin width
        genJetHistoCross.SetLineWidth(2)
        genJetHistoCross.SetLineColor(ROOT.kMagenta)
        genJetHistoCross.SetMarkerStyle(0)
        genJetHistoCross.SetLineStyle(2)
        legend.AddEntry( genJetHistoCross, 'MG5-MLM+H7 (gen)', 'lp' )
        
        unfoldHistoCross.Draw( "E same")
        genJetHistoCross.Draw( "histe same")
        

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.045)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
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
    tmpPad2.GetYaxis().SetRangeUser(0.5,1.5 )
    tmpPad2.GetXaxis().SetRangeUser(unfoldHisto.GetBinLowEdge(1),unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2) )
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
        #hRatioUp.GetXaxis().SetLimits(0.,unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        hRatioUp.Draw('P0')

    else:
        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfoldHistoCross, unfoldHisto, 'pois' )
        hRatioUp2.SetLineColor(ROOT.kBlack)
        hRatioUp2.SetMarkerColor(ROOT.kBlack)
        hRatioUp2.SetLineWidth(2)
        hRatioUp2.SetMarkerStyle(25)
        #hRatioUp2.GetXaxis().SetLimits(0.,unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
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
        dictHistos['combData'].GetYaxis().SetTitle( '#frac{1}{dN} #frac{dN}{d#tau_{X}}' )
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
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or 'JER' in text or ('PU' in text and 'DAMP' not in text) or 'PDF' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2 if not('L1' in text) else 1)
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
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or 'JER' in text or ('PU' in text and 'DAMP' not in text) or 'PDF' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[dummy])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2 if not('L1' in text) else 1)
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
    modelUnc.SetLineWidth(2)
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
    totalErrHist.GetYaxis().SetRangeUser(0.,2.0)
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(ROOT.kGray+3,0.9)
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3154)
    totalErrHist.SetLineColor(1)
    totalErrHist.Draw(' E5')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(3)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.8)
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
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
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
    print ('Processing plots for sys comparisons......', outputFileName)

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
    for ih in dictUncHistos.keys():
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

        
        if 'jes' in ih: 
            dictShifts[ih].SetLineStyle( 1 )

        dictShifts[ih].SetMarkerStyle(0)
        dictShifts[ih].SetLineWidth( 1 )
        
        y=year if not('all' in year or '+' in year) else 'all'
        stringtocheck=ih.split(y+'_')[1]#+'_'+(year if not('+' in year) else '_fullRunII')
        #print(stringtowrite)
        #if 'jes' in stringtocheck and ('2016' in stringtocheck or'2017' in stringtocheck or '2018' in stringtocheck): 
            #print (ih, ih.split('_')[1])
        legend.AddEntry( dictShifts[ih], stringtocheck+'_'+y, 'l' )
        #else: legend.AddEntry( dictShifts[ih], ih.split('_')[1] , 'l' )
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

def createRecoBins(genBins):
    recoBins=[genBins[0]]
    for i in range(1,len(genBins)):
        recoBins.append(np.round(genBins[i-1]+(genBins[i]-genBins[i-1])/2.,3))
        recoBins.append(genBins[i])
        #recoBins.append(b)
    return recoBins

def makePSplot_simple(purity,stability,variables,var,outputDir,ext='pdf',dictHistos=OrderedDict(),bins=[0.,1.],sel='_dijetSel',year='2017',signalLabelBegin='QCD_HT_MG5+P8'):
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    colorPallete = [ 0, 2, 4, 8, 12, 28, 30 ]

    ROOT.gStyle.SetPadRightMargin(0.05)
    canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)
    #canvas.cd()
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
    p.GetXaxis().SetTitle( variables[var]['label'] )
    p.Draw('hist')
    s.Draw('hist same')
    
    dictHistos[ 'purityGraph_'+var ] = p.Clone()
    dictHistos[ 'stabilityGraph_'+var ] = s.Clone()

    '''
    numBins = numBins + dictHistos[ 'purityGraph_'+var ].GetNbinsX()
    numBinsList.append(numBins)
    multigraph = ROOT.THStack()
    for i in dictHistos:
        if i.endswith(var): multigraph.Add( dictHistos[i]  )

    ROOT.gStyle.SetPadRightMargin(0.05)
    canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)

    multigraph.Draw("hist nostack")
    multigraph.GetXaxis().SetTitle( variables[var]['label'] )
    multigraph.GetYaxis().SetTitle( 'Percentage' )
    multigraph.SetMaximum( 1.1 )
    multigraph.GetYaxis().SetTitleOffset(0.8)
    '''
    #canvas.PlaceLegend()
    legend.Draw()
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canvas, 4, 0)
    '''
    textBox=ROOT.TLatex()
    textBox.SetNDC()
    textBox.SetTextSize(0.04)
    textBox.SetTextFont(62) ### 62 is bold, 42 is normal
    '''
    #textBox.DrawLatex(0.65, 0.75, sel.split('Sel')[0].split('_')[1]+' Selection' )
    canvas.Update()
    canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.'+ext)
    if ext.startswith('pdf'):
        canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.png')
    return 1#numBins
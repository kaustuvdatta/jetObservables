import copy, pprint, array, bisect, scipy,os, sys, glob, math, random, gc
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
from unfoldingPlottersAndHelpers import *

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
######################Helpers for combined unfoldings ######################
#############################################################################
def plot_combined_and_1D_MCClosures( combined_altMC_truth_hist,
                                     unfolded_combined_nom,
                                     unfolded_combined_alt,
                                     combined_nomMC_truth_hist,
                                     varDict,
                                     cov_tot, cov_datastat_tot,
                                     cov_tot_cross, cov_datastat_tot_cross,genBinMap,
                                     n_obs=25,
                                     selection='_dijetSel',
                                     process='MCCrossClosure', 
                                     labelX='',
                                     outputDir = '../Plots_April25_dijetSel_ApprovalChecks/dijetSel/',
                                     outputFilename = 'combined_MCCrossClosure',
                                     tlegendAlignment='right',                                 
                                     year='all',
                                     ext='.pdf',
                                     maxYFactor=1.4,
                                     nomMCHisto_label = None, 
                                     altMCHisto_label = None, 
                                     noNorm = False,version='April25',
                                     signalLabel = 'MLMQCD_HT2000toInf',
                                     alt0SignalLabel = 'H7MLMQCD_HT2000toInf',
                                     return1DHistDict = False
                                    ):

    if return1DHistDict:
        dict_1DHists = OrderedDict()
        for var in varDict.keys():
            dict_1DHists[var] = OrderedDict()
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 800
    H_ref = 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw
    #can = ROOT.TCanvas('can'+'CrossClosure', 'can'+'CrossClosure',  10, 10, 800, 500 )
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1) 
    
    can = ROOT.TCanvas('canUnfolding', 'canUnfolding',  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
    
    pad1 = ROOT.TPad("pad1", "Main",0, Hdw / H, 1, 1, -1)
    
    
    #pad1.SetPad(0, Hdw / H, 1, 1)
    pad1.SetLeftMargin(L)
    pad1.SetRightMargin(R)
    pad1.SetTopMargin(Tup)
    pad1.SetBottomMargin(Bup)
    
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.60,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.60,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.040)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    

    if not(noNorm): 
        cov_normTot_np, normed_covTot = get_normalised_cov_combined(unfolded_combined_nom, cov_tot.Clone(), genBinMap)
        cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov_combined(unfolded_combined_nom, cov_datastat_tot.Clone(), genBinMap)
    
        unfolded_combined_nom = normalize_combined_TH1_by_blocks(unfolded_combined_nom, genBinMap, noNorm=False)
        get_th1_normedCovErrors(unfolded_combined_nom, cov_normTot_np)
    
    physical_unfoldHisto = physical_histograms_from_combined_Ndim(combined_hist=unfolded_combined_nom.Clone(unfolded_combined_nom.GetName() + '_cloned_physical'),
                                                                  bin_map=genBinMap,
                                                                  varDict=varDict,
                                                                  new_hist_prefix='unfoldHisto',
                                                                ) 
    
    if process.startswith('MCCrossClosure'): 
        
        
        if not(noNorm): 
            cov_normTot_np_cross, normed_covTot_cross = get_normalised_cov_combined(unfolded_combined_alt.Clone(), 
                                                                       cov_tot_cross.Clone(), genBinMap)
            cov_norm_dataStat_np_cross, normed_cov_dataStat_cross = get_normalised_cov_combined(unfolded_combined_alt.Clone(), cov_datastat_tot_cross.Clone(), genBinMap)
            unfolded_combined_alt = normalize_combined_TH1_by_blocks(unfolded_combined_alt, genBinMap, noNorm=False)
            get_th1_normedCovErrors(unfolded_combined_alt, cov_normTot_np_cross)

        physical_unfoldHistoCross = physical_histograms_from_combined_Ndim(combined_hist=unfolded_combined_alt.Clone(unfolded_combined_alt.GetName() + '_cloned_physical'),
                                                                             bin_map=genBinMap,
                                                                             varDict=varDict,
                                                                             new_hist_prefix='unfoldHistoCross',
                                                                            ) 
        print(combined_nomMC_truth_hist.Integral(), unfolded_combined_nom.Integral(),unfolded_combined_alt.Integral(),combined_altMC_truth_hist.Integral())
        

    unfolded_combined_nom.SetMarkerStyle(4)
    unfolded_combined_nom.SetMarkerSize(0.5)
    unfolded_combined_nom.SetMarkerColor(ROOT.kRed)
    unfolded_combined_nom.SetLineColor(ROOT.kRed)
    unfolded_combined_nom.SetLineWidth(1)
    
    legend.AddEntry( unfolded_combined_nom, (f'{nomMCHisto_label} (closure)' if process.startswith('MCSelfClosure') else f'#splitline{{{nomMCHisto_label} unf. with }}{{{nomMCHisto_label.replace("-FXFX","")} }}'), 'pe' )
    
    
    if not(noNorm): 
        combined_nomMC_truth_hist = normalize_combined_TH1_by_blocks(combined_nomMC_truth_hist, genBinMap, noNorm=False)
    
    physical_genJetHisto = physical_histograms_from_combined_Ndim(combined_hist=combined_nomMC_truth_hist.Clone(),
                                                                  bin_map=genBinMap,
                                                                  varDict=varDict,
                                                                  new_hist_prefix = signalLabel+'_gen',#'MLMQCD_HT2000toInf',
                                                                  withSuff=True, sys='_nom'
                                                                 ) 
    combined_nomMC_truth_hist.SetLineWidth(1)
    combined_nomMC_truth_hist.SetLineColor(ROOT.kBlue)
    combined_nomMC_truth_hist.SetMarkerStyle(0)
    #combined_nomMC_truth_hist.SetMarkerSize(0.5)
    combined_nomMC_truth_hist.SetLineStyle(2)
    legend.AddEntry( combined_nomMC_truth_hist, f'{nomMCHisto_label} (gen)', 'lp' )
    
    unfolded_combined_nom.GetXaxis().SetTitleOffset(999)    
    unfolded_combined_nom.GetXaxis().SetLabelOffset(999)    

    
    #if 'tau' in labelX :
    unfolded_combined_nom.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#tau_{N}^{(#beta)}}' if not noNorm else 'N_{events}')
    
            
    unfolded_combined_nom.GetYaxis().SetTitleOffset(extraSpace+0.9*Hup/H_ref)    
    unfolded_combined_nom.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    unfolded_combined_nom.GetYaxis().SetLabelSize(0.048* H_ref / Hup)
    unfolded_combined_nom.GetYaxis().SetTitleFont(42)

    #unfolded_combined_nom.SetTickLength(0.03, "XY")  #?? ok if 1/3
    unfolded_combined_nom.GetXaxis().SetRangeUser(-0.5, unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1))#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfolded_combined_nom.GetYaxis().SetTickLength(0.03)#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfolded_combined_nom.GetXaxis().SetTickLength(0.03)#* H_ref / Hdw)#, "X")

    unfolded_combined_nom.SetMaximum( maxYFactor*max([ combined_nomMC_truth_hist.GetMaximum(), unfolded_combined_nom.GetMaximum(), unfolded_combined_alt.GetMaximum() if 'Cross' in process else unfolded_combined_nom.GetMaximum() ] )  )
    
    unfolded_combined_nom.Draw( "AXIS")
    can.Update()
    can.Modified()
    unfolded_combined_nom.Draw( "PE1 SAME")    
    combined_nomMC_truth_hist.Draw( "histE1 same")
    
    if not process.startswith('MCSelfClosure'):
        
        #unfolded_combined_alt.Scale(1./( (1./n_obs) * unfolded_combined_alt.Integral() if not noNorm else 1.),'width')
        
        unfolded_combined_alt.SetMarkerStyle(26)
        unfolded_combined_alt.SetMarkerSize(0.5)
        unfolded_combined_alt.SetMarkerColor(ROOT.kRed+4)
        unfolded_combined_alt.SetLineColor(ROOT.kRed+4)
        unfolded_combined_alt.SetLineWidth(1)
        legend.AddEntry( unfolded_combined_alt, f'#splitline{{{nomMCHisto_label} unf. with }}{{{altMCHisto_label.replace("-FXFX","")} }}', 'pe')

        if not(noNorm): 
            combined_altMC_truth_hist = normalize_combined_TH1_by_blocks(combined_altMC_truth_hist, genBinMap, noNorm=False)
        
        physical_genJetHistoCross = physical_histograms_from_combined_Ndim(combined_hist=combined_altMC_truth_hist.Clone(),
                                                                      bin_map=genBinMap,
                                                                      varDict=varDict,
                                                                      new_hist_prefix = alt0SignalLabel+'_gen',#'MLMQCD_HT2000toInf',
                                                                      withSuff=True, sys='_nom'
                                                                     ) 
        combined_altMC_truth_hist.SetLineWidth(1)
        combined_altMC_truth_hist.SetLineColor(ROOT.kMagenta)
        combined_altMC_truth_hist.SetMarkerStyle(0)
        combined_altMC_truth_hist.SetLineStyle(2)
        legend.AddEntry( combined_altMC_truth_hist, f'{altMCHisto_label.replace("-FXFX","")} (gen)', 'lp')
        
        unfolded_combined_alt.Draw( "PE1 same")
        combined_altMC_truth_hist.Draw( "histE1 same")

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if "dijet" in selection: 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif "W" in selection: seltext = 'Boosted W-enriched'
    elif "top" in selection: seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if "dijet" in selection: seltext = 'p_{T}>200 GeV' 
    elif "W" in selection: seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif "top" in selection: seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.80, seltext )

    
    legend.Draw()
    for var in varDict.keys():
        outName = outputDir+var+selection+f"_fromMBody_{'MCSelfClosure' if process.startswith('MCSelfClosure') else 'MCCrossClosure'}"+signalLabel+f'_TUnfold_{"NO_NORM_" if noNorm else ""}'+version+'.pdf'
    


        drawClosures1DfromNDim( var, selection, 
                                process, year, lumi, 
                                physical_genJetHisto[var].Clone(), 
                                physical_genJetHistoCross[var].Clone() if 'cross' in process.lower() else None, 
                                physical_unfoldHisto[var].Clone(), 
                                physical_unfoldHistoCross[var].Clone() if 'cross' in process.lower() else None,
                                #ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, 
                                labelX=varDict[var]['label'], 
                                maxX=varDict[var]['bins'][-1],
                                tlegendAlignment=varDict[var]['alignLeg'], 
                                outputName=outName,version=version,
                                nomMCHisto_label = nomMCHisto_label, 
                                altMCHisto_label = altMCHisto_label if 'cross' in process.lower() else None ,
                                noNorm = True )
        if return1DHistDict:
            dict_1DHists[var]['selection'] = selection 
            dict_1DHists[var]['process'] = process 
            dict_1DHists[var]['year'] = year 
            dict_1DHists[var]['lumi'] = lumi
            dict_1DHists[var]['genJetHisto'] = physical_genJetHisto[var].Clone() 
            dict_1DHists[var]['genJetHisto'].Sumw2()
            dict_1DHists[var]['genJetHisto'].SetDirectory(0)
            
            dict_1DHists[var]['unfoldHisto'] = physical_unfoldHisto[var].Clone()
            dict_1DHists[var]['unfoldHisto'].Sumw2()
            dict_1DHists[var]['unfoldHisto'].SetDirectory(0)
            
            #dict_1DHists[var]['ratioUncHisto'] = ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, 
            dict_1DHists[var]['labelX'] = varDict[var]['label'] 
            dict_1DHists[var]['maxX'] = varDict[var]['bins'][-1]
            dict_1DHists[var]['tlegendAlignment'] = varDict[var]['alignLeg'][-1] 
            dict_1DHists[var]['nomMCHisto_label'] = nomMCHisto_label 
            dict_1DHists[var]['noNorm'] = True 
            
            if 'cross' in process.lower():
                dict_1DHists[var]['genJetHistoCross'] = physical_genJetHistoCross[var].Clone() if 'cross' in process.lower() else None
                dict_1DHists[var]['genJetHistoCross'].Sumw2()
                dict_1DHists[var]['genJetHistoCross'].SetDirectory(0)
                dict_1DHists[var]['unfoldHistoCross'] = physical_unfoldHistoCross[var].Clone() if 'cross' in process.lower() else None
                dict_1DHists[var]['unfoldHistoCross'].Sumw2()
                dict_1DHists[var]['unfoldHistoCross'].SetDirectory(0)
                dict_1DHists[var]['altMCHisto_label'] = altMCHisto_label if 'cross' in process.lower() else None 
                
        
    
    can.cd()
    pad1.cd()
    #if process.startswith('data'):
    CMS_lumi.extraText = "Simulation Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #   CMS_lumi.extraText = "Simulation Preliminary"
    #   CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2", "Ratio",0, 0, 1, Hdw / H,-1);
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
    
    tmpPad2= pad2.DrawFrame( -0.5, 0., unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1), 1.9 )
    print (labelX)
    
    tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    tmpPad2.GetXaxis().SetRangeUser(-0.5, unfolded_combined_nom.GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1) )
    
    tmpPad2.GetYaxis().SetTitleOffset(extraSpace + (0.8 ) * Hdw / H_ref)
    tmpPad2.GetXaxis().SetTitleOffset(0.92)
    tmpPad2.SetTitleSize(0.054 * H_ref / Hdw, "Y")
    tmpPad2.SetLabelSize(0.048 * H_ref / Hdw, "Y")
    tmpPad2.SetTitleSize(0.054 * H_ref / Hdw, "X")
    tmpPad2.SetLabelSize(0.048 * H_ref / Hdw, "X")
    tmpPad2.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    tmpPad2.GetXaxis().SetTitle( f'{labelX} N-subjettiness basis' )
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Unf.}" if 'Self' in process else "#frac{Unf. alt.}{Unf. nom.}" )
    tmpPad2.GetYaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    tmpPad2.GetYaxis().SetNdivisions(505)
    tmpPad2.GetYaxis().CenterTitle()
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    pad2.cd()
    
    if 'Self' in process:

        hRatioUp = ROOT.TGraphAsymmErrors()
        hRatioUp.Divide( combined_nomMC_truth_hist, unfolded_combined_nom, 'pois' )
        hRatioUp.SetLineColor(ROOT.kBlack)
        hRatioUp.SetMarkerColor(ROOT.kBlack)
        hRatioUp.SetLineWidth(1)
        hRatioUp.SetMarkerStyle(25)
        hRatioUp.SetMarkerSize(0.5)
        
    
        #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
        hRatioUp.GetYaxis().SetTitleOffset(extraSpace + (0.8 ) * Hdw / H_ref)
        hRatioUp.GetXaxis().SetTitleOffset(0.92)
        hRatioUp.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
        hRatioUp.GetYaxis().SetLabelSize(0.048 * H_ref / Hdw)#, "Y")
        hRatioUp.GetXaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelSize(0.048 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelOffset(0.012 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetTitle( f'{labelX} N-subjettiness basis' )
        hRatioUp.GetYaxis().SetTitle( "#frac{Sim.}{Unf.}" )
        hRatioUp.GetYaxis().SetTitleFont(42)
        hRatioUp.GetXaxis().SetTitleFont(42)

        #tmpPad2.GetXaxis().SetTitle(nameXaxis)
        #tmpPad2.GetYaxis().SetTitle(nameRatio)

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp.GetYaxis().SetNdivisions(505)
        hRatioUp.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp,1.5,0.5)
        
        hRatioUp.GetXaxis().SetLimits( -0.05,unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1))#unfolded_combined_nom.GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+2))
        
        hRatioUp.Draw('AP')
        hRatioUp.GetYaxis().SetRangeUser(0.2,1.8 )
        hRatioUp.GetXaxis().SetRangeUser(-0.5, unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1) )
        pad2.Update()
        can.Update()
        
    else:
        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfolded_combined_alt, unfolded_combined_nom, 'pois' )
        hRatioUp2.SetLineColor(ROOT.kBlack)
        hRatioUp2.SetMarkerColor(ROOT.kBlack)
        hRatioUp2.SetLineWidth(1)
        hRatioUp2.SetMarkerStyle(25)
        hRatioUp2.SetMarkerSize(0.5)
        
        #hRatioUp2.GetXaxis().SetLimits(0.,maxX)#unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        
        #hRatioUp2.GetYaxis().SetRangeUser(0.3,1.9 )
    
        #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
        hRatioUp2.GetYaxis().SetTitleOffset(extraSpace + (0.9 ) * Hdw / H_ref)
        hRatioUp2.GetXaxis().SetTitleOffset(0.92)
        hRatioUp2.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
        hRatioUp2.GetYaxis().SetLabelSize(0.048 * H_ref / Hdw)#, "Y")
        hRatioUp2.GetXaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "X")
        hRatioUp2.GetXaxis().SetLabelSize(0.048 * H_ref / Hdw)#, "X")
        hRatioUp2.GetXaxis().SetLabelOffset(0.012 * H_ref / Hdw)#, "X")
        hRatioUp2.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis' )
        hRatioUp2.GetYaxis().SetTitle( "#frac{Unf. alt.}{Unf. nom.}" )
        hRatioUp2.GetYaxis().SetTitleFont(42)
        hRatioUp2.GetXaxis().SetTitleFont(42)

        #tmpPad2.GetXaxis().SetTitle(nameXaxis)
        #tmpPad2.GetYaxis().SetTitle(nameRatio)

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp2.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp2.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp2.GetYaxis().SetNdivisions(505)
        hRatioUp2.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp2,1.5,0.5)
        
        hRatioUp2.GetXaxis().SetLimits( -0.05,unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1))#unfolded_combined_nom.GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+2))
        
        hRatioUp2.Draw('AP')
        hRatioUp2.GetYaxis().SetRangeUser(0.2,1.8 )
        hRatioUp2.GetXaxis().SetRangeUser(-0.5, unfolded_combined_nom.GetXaxis().GetBinLowEdge(unfolded_combined_nom.GetNbinsX()+1) )
        pad2.Update()
        can.Update()
    
    png = (outputDir+outputFilename+ext).split('.pdf')[0]+'.png'
    
    can.SaveAs(outputDir+outputFilename+ext)
    can.SaveAs(png)
    
    if return1DHistDict:
        return dict_1DHists


def drawClosures1DfromNDim(ivar, selection, process, year, lumi, genJetHisto, genJetHistoCross, unfoldHisto, unfoldHistoCross,
                             #ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, 
                             labelX, maxX, tlegendAlignment, 
                             outputName, version, nomMCHisto_label = None, altMCHisto_label = None, noNorm = True ):
    
    if process.startswith('MCCrossClosure'):
    
        genJetHistoCross.SetTitle("") 
        unfoldHistoCross.SetTitle("")
        
    else:

        genJetHisto.SetTitle("") 
        unfoldHisto.SetTitle("")

    
    """docstring for drawClosures"""
    print ("Drawing unfolding closure")
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
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
    #   else: legend=ROOT.TLegend(0.20,0.61,0.52,0.89)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.60,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.60,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    
        
    
    
    #unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.))
    #if not(noNorm): get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    
    
    if process.startswith('MCCrossClosure'): 
        print(genJetHisto.Integral(), unfoldHisto.Integral(),unfoldHistoCross.Integral(),genJetHistoCross.Integral())
        
                
        #unfoldHistoCross.Scale(1./(unfoldHistoCross.Integral() if not(noNorm) else 1.))
        #if not(noNorm): get_th1_normedCovErrors(unfoldHistoCross, cov_normTot_np_cross)
        unfoldHistoCross.Scale(1.,'width')

    unfoldHisto.Scale(1.,'width')

    
    
    #unfoldHisto.Scale(1./(unfoldHisto.Integral() if not noNorm else 1.),'width')
    unfoldHisto.SetMarkerStyle(4)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    unfoldHisto.SetLineWidth(2)
    
    legend.AddEntry( unfoldHisto, (f'{nomMCHisto_label} (closure)' if process.startswith('MCSelfClosure') else f'#splitline{{{nomMCHisto_label} unf. with }}{{{nomMCHisto_label.replace("-FXFX","")} }}'), 'pe' )
    
    
    genJetHisto.Scale(1.,'width')
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

    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3
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

        genJetHistoCross.Scale(1.,'width')
        
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
    
    if "dijet" in selection: 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.20
    elif "W" in selection: seltext = 'Boosted W-enriched'
    elif "top" in selection: seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if "dijet" in selection: seltext = 'p_{T}>200 GeV' 
    elif "W" in selection: seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif "top" in selection: seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.80, seltext )

    
    legend.Draw()
    #if process.startswith('data'):
    CMS_lumi.extraText = "Simulation Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #   CMS_lumi.extraText = "Simulation Preliminary"
    #   CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
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

    #Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp2.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp2.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp2.GetYaxis().SetNdivisions(505)
        hRatioUp2.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp2,1.5,0.5)
        
        
        hRatioUp2.Draw('P0')
    
    
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)

def drawClosures1DfromNDim_vs_singleObs(
                             ivar, 
                             selection, 
                             process, 
                             year, 
                             lumi, 
                             genJetHisto, 
                             genJetHistoCross, 
                             unfoldHisto, 
                             unfoldHistoCross,
                             unfoldHisto_old, 
                             unfoldHistoCross_old,
                             labelX, maxX, tlegendAlignment, 
                             outputName, version, 
                             nomMCHisto_label = None, altMCHisto_label = None, noNorm = True ):
    
    if process.startswith('MCCrossClosure'):
    
        genJetHistoCross.SetTitle("") 
        unfoldHistoCross.SetTitle("")
        
    else:

        genJetHisto.SetTitle("") 
        unfoldHisto.SetTitle("")

    
    """docstring for drawClosures"""
    print ("Drawing unfolding closure")
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
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
    #   else: legend=ROOT.TLegend(0.20,0.61,0.52,0.89)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.50,0.91-0.02,0.88)

    else: legend=ROOT.TLegend(0.19,0.50,0.45-0.02,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.028)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    
    
    
    if process.startswith('MCCrossClosure'): 
        print(genJetHisto.Integral(), unfoldHisto.Integral(),unfoldHistoCross.Integral(),genJetHistoCross.Integral())
      
        unfoldHistoCross.Scale(1.,'width')

        #old histos for model dep

        unfoldHistoCross_old.Scale(1.,'width')
        unfoldHistoCross_old.SetMarkerStyle(27)
        unfoldHistoCross_old.SetLineStyle(2)
        unfoldHistoCross_old.SetMarkerColor(38)
        unfoldHistoCross_old.SetLineColor(38)
        unfoldHistoCross_old.SetLineWidth(2)
        

    unfoldHisto.Scale(1.,'width')
    
    unfoldHisto.SetMarkerStyle(4)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    unfoldHisto.SetLineWidth(2)
    
    legend.AddEntry( unfoldHisto, (f'{nomMCHisto_label} (closure)' if process.startswith('MCSelfClosure') else f'#splitline{{{nomMCHisto_label} unf. with }}{{{nomMCHisto_label.replace("-FXFX","")} }}'), 'lpe' )
    

    unfoldHisto_old.Scale(1.,'width')
    unfoldHisto_old.SetMarkerStyle(4)
    unfoldHisto_old.SetLineStyle(2)
    unfoldHisto_old.SetMarkerColor(ROOT.kGray+1)
    unfoldHisto_old.SetLineColor(ROOT.kGray+1)
    unfoldHisto_old.SetLineWidth(2)
    
    


    genJetHisto.Scale(1.,'width')
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

    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3
    unfoldHisto.SetMaximum( 1.7*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum(), unfoldHistoCross.GetMaximum() if 'Cross' in process else 1.7*unfoldHisto.GetMaximum() ] )  )
    
    unfoldHisto.Draw( "AXIS")
    can.Update()
    can.Modified()
    genJetHisto.Draw( "histE same")
    unfoldHisto_old.Draw( "histE same")    
    unfoldHisto.Draw( "E same")    

    if not process.startswith('MCSelfClosure'):
        
        #unfoldHistoCross.Scale(1./(unfoldHistoCross.Integral() if not noNorm else 1.),'width')
        
        unfoldHistoCross.SetMarkerStyle(26)
        #unfoldHistoCross.SetMarkerSize(2)
        unfoldHistoCross.SetMarkerColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineWidth(2)
        legend.AddEntry( unfoldHistoCross, f'#splitline{{{nomMCHisto_label} unf. with }}{{{altMCHisto_label.replace("-FXFX","")} }}', 'lpe')

        genJetHistoCross.Scale(1.,'width')
        
        genJetHistoCross.SetLineWidth(2)
        genJetHistoCross.SetLineColor(ROOT.kMagenta)
        genJetHistoCross.SetMarkerStyle(0)
        genJetHistoCross.SetLineStyle(2)
        legend.AddEntry( genJetHistoCross, f'{altMCHisto_label.replace("-FXFX","")} (gen)', 'lp')
        
        legend.AddEntry(unfoldHistoCross_old, 
                        f'#splitline{{{nomMCHisto_label} unf. with }}{{{altMCHisto_label.replace("-FXFX","")} (1D unf.) }}',
                        'lpe')
        legend.AddEntry(unfoldHisto_old, 
                    f'{nomMCHisto_label} (closure)' if process.startswith('MCSelfClosure') else f'#splitline{{{nomMCHisto_label} unf. with }}{{{nomMCHisto_label.replace("-FXFX","")} (1D unf.)}}', 
                    'lpe')
        
        unfoldHistoCross.Draw( "E same")
        genJetHistoCross.Draw( "histE same")
        unfoldHistoCross_old.Draw( "histE same")

        

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if "dijet" in selection: 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.20
    elif "W" in selection: seltext = 'Boosted W-enriched'
    elif "top" in selection: seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if "dijet" in selection: seltext = 'p_{T}>200 GeV' 
    elif "W" in selection: seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif "top" in selection: seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.51+dijetOffset ), 0.80, seltext )

    
    legend.Draw()
    #if process.startswith('data'):
    CMS_lumi.extraText = "Simulation Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #   CMS_lumi.extraText = "Simulation Preliminary"
    #   CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
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
    
    tmpPad2.GetYaxis().SetRangeUser(0.4,1.8 )
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

    #Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp.GetYaxis().SetNdivisions(505)
        hRatioUp.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp,1.5,0.5)
        
        
        hRatioUp.Draw('P0')

    else:
        hRatioUp = ROOT.TGraphAsymmErrors()
        hRatioUp.Divide( unfoldHistoCross, unfoldHisto, 'pois' )
        hRatioUp.SetLineColor(ROOT.kBlack)
        hRatioUp.SetMarkerColor(ROOT.kBlack)
        hRatioUp.SetLineWidth(2)
        hRatioUp.SetMarkerStyle(25)
        
        #hRatioUp2.GetXaxis().SetLimits(0.,maxX)#unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        
        #hRatioUp2.GetYaxis().SetRangeUser(0.3,1.9 )
    
        #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
        hRatioUp.GetYaxis().SetTitleOffset(extraSpace + (1.15 ) * Hdw / H_ref)
        hRatioUp.GetXaxis().SetTitleOffset(0.92)
        hRatioUp.GetYaxis().SetTitleSize(0.054 * H_ref / Hdw)#, "Y")
        hRatioUp.GetYaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "Y")
        hRatioUp.GetXaxis().SetTitleSize(0.056 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelSize(0.047 * H_ref / Hdw)#, "X")
        hRatioUp.GetXaxis().SetLabelOffset(0.012 * H_ref / Hdw)#, "X")
        #hRatioUp2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        if 'tau' in labelX: 
            hRatioUp.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
        else:
            hRatioUp.GetXaxis().SetTitle( labelX)
        hRatioUp.GetYaxis().SetTitle( "#frac{Unf. alt.}{Unf. nom.}" )
        hRatioUp.GetYaxis().SetTitleFont(42)
        hRatioUp.GetXaxis().SetTitleFont(42)

        #tmpPad2.GetXaxis().SetTitle(nameXaxis)
        #tmpPad2.GetYaxis().SetTitle(nameRatio)

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp.GetYaxis().SetNdivisions(505)
        hRatioUp.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp2,1.5,0.5)
        
        

        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfoldHistoCross_old, unfoldHisto_old, 'pois' )
        hRatioUp2.SetLineColor(38)
        hRatioUp2.SetMarkerColor(38)
        hRatioUp2.SetLineWidth(2)
        hRatioUp2.SetMarkerStyle(8)
        hRatioUp2.SetLineStyle(3)
        
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

        #Set tick lengths to match original (these are fractions of axis length)
        hRatioUp2.GetYaxis().SetTickLength(0.03 * H_ref / Hup)#, "Y")  #?? ok if 1/3
        hRatioUp2.GetXaxis().SetTickLength(0.03 * H_ref / Hdw)#, "X")

        #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
        hRatioUp2.GetYaxis().SetNdivisions(505)
        hRatioUp2.GetYaxis().CenterTitle()
        
        #set_dynamic_y_range_errRatioHist(hRatioUp2,1.5,0.5)
        
        
        hRatioUp.Draw('P0')
        hRatioUp2.Draw('P0 same')

        ratioLegend=ROOT.TLegend(0.19,0.70,0.85,0.88)
        ratioLegend.SetTextSize(0.077)
        ratioLegend.SetTextFont(42)
        ratioLegend.SetNColumns(2)
        ratioLegend.SetFillStyle(0)#ColorAlpha(10,0.6)
        ratioLegend.SetBorderSize(0)
        ratioLegend.AddEntry( hRatioUp, '#frac{Unf. alt.}{Unf. nom.}', 'lpe' )
        ratioLegend.AddEntry( hRatioUp2, '#frac{Unf. alt.}{Unf. nom.} 1D unf. (old)', 'lpe' )
        ratioLegend.Draw()

    
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)

def drawUnfold_Ndim(ivar, process, 
                    lumi, 
                    varDict,
                    genBinMap,
                    dataJetHisto, 
                    genJetHisto, 
                    unfoldHisto, 
                    unfoldHistoStatUnc, 
                    unfoldHistowoUnc,
                    altMCHisto, 
                    foldHisto, 
                    recoJetHisto, 
                    cov_tot, 
                    cov_datastat_tot, 
                    labelX, 
                    maxX, 
                    tlegendAlignment, 
                    outputDir,
                    outputName,
                    year='all',
                    selection='_dijetSel',
                    altMC1Histo = None, 
                    altMC2Histo = None, 
                    altMC1Histo_label = None, 
                    altMC2Histo_label = None, 
                    nomMCHisto_label = None, altMCHisto_label = None,
                    extraMC=False,
                    includeFSR = False, fsrUpHisto = None, fsrDownHisto=False, 
                    noNorm=False,
                    signalLabel = 'MLMQCD_HT2000toInf',
                    alt0SignalLabel = 'H7MLMQCD_HT2000toInf',
                    alt1SignalLabel = 'QCD_HT2000toInf',
                    alt2SignalLabel = 'QCD_Pt_3200toInf',
                    fsrLabel = 'sysMLMQCD_fsrWeight_HT2000toInf',
                    n_obs=25
                   ):
    
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 800
    H_ref = 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
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
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.70,0.60,0.90,0.90)

    else: legend=ROOT.TLegend(0.19,0.60,0.45-0.02,0.88)
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
                                                       cov_tot.Clone(),ndim=n_obs)
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatErr, 
                                                                   cov_datastat_tot.Clone(),ndim=n_obs)
    
    
    
    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()

    unfoldHistowoUnc.Scale(1./((1./n_obs) * unfoldHistowoUnc.Integral() if not(noNorm) else 1.))#,'width')

    unfoldHistoDataStatErr.Scale(1./((1./n_obs) * unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.))
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_norm_dataStat_np)
    else:
        cov_abs_dataStat_np, _ = th2_to_ndarray(cov_datastat_tot.Clone())
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_abs_dataStat_np)
        
    #unfoldHistoDataStatErr.Scale(1.,'width')

    unfoldHisto.Scale(1./((1./n_obs) * unfoldHisto.Integral() if not(noNorm) else 1.))
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    else:
        cov_absTot_np, _ = th2_to_ndarray(cov_tot.Clone())
        get_th1_normedCovErrors(unfoldHisto, cov_absTot_np)
        
    #unfoldHisto.Scale(1.,'width')

    
    dataJetHisto.Scale(1./((1./n_obs) * dataJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    genJetHisto.Scale(1./((1./n_obs) * genJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    altMCHisto.Scale(1./((1./n_obs) * altMCHisto.Integral() if not(noNorm) else 1.))#,'width')
    foldHisto.Scale(1./((1./n_obs) * foldHisto.Integral() if not(noNorm) else 1.))#,'width')
    recoJetHisto.Scale(1./((1./n_obs) * recoJetHisto.Integral() if not(noNorm) else 1.))#,'width')
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        fsrUpHisto.Scale(1./((1./n_obs) * fsrUpHisto.Integral() if not(noNorm) else 1.))#,'width')
        fsrDownHisto.Sumw2()
        fsrDownHisto.Scale(1./((1./n_obs) * fsrDownHisto.Integral() if not(noNorm) else 1.))#,'width')
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        altMC1Histo.Scale(1./((1./n_obs) * altMC1Histo.Integral() if not(noNorm) else 1.))#,'width')
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo.Scale(1./((1./n_obs) * altMC2Histo.Integral() if not(noNorm) else 1.))#,'width')
            
            altMC2Histo.SetTitle("")

    
    

    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetLineWidth(1)
    unfoldHisto.SetMarkerSize(0.5)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(1)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(0.5)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(0.55)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(23)


        fsrDownHisto.SetMarkerSize(0.55)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(22)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )
    
    unfoldHisto.GetXaxis().SetTitleOffset(999)    
    unfoldHisto.GetXaxis().SetLabelOffset(999)    

    #unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+1.15*Hup/H_ref)    
    #unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    #unfoldHisto.GetYaxis().SetLabelSize(0.047* H_ref / Hup)
    
    unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+0.92*Hup/H_ref)    
    unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelSize(0.048* H_ref / Hup)
    unfoldHisto.GetYaxis().SetTitleFont(42)

    #unfolded_combined_nom.SetTickLength(0.03, "XY")  #?? ok if 1/3
    unfoldHisto.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetXaxis().GetBinLowEdge(unfoldHisto.GetNbinsX()+1))#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfoldHisto.GetYaxis().SetTickLength(0.03)#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfoldHisto.GetXaxis().SetTickLength(0.03)#* H_ref / Hdw)#, "X")

   
    if 'body' in labelX: 
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#tau_{N}^{(#beta)}}' if not noNorm else 'N_{events}')
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
    #unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( (1.6)*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3
    
    if noNorm:
        unfoldHisto.GetYaxis().SetMaxDigits(3)  #?? ok if 1/3
        ROOT.TGaxis.SetExponentOffset(-0.06, 0.005, "y")
    unfoldHisto.Draw( "AXIS")
    
    
    can.Update()
    
    
    
    #altMCHisto.Scale(1, 'width')  ###divide by bin width
    altMCHisto.SetLineWidth(1)
    altMCHisto.SetMarkerSize(0.5)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 
        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')

        
        
    legend.AddEntry( altMCHisto, 'MG5-MLM+H7' if 'dijet' in selection else'PWHG+H7','lpe')#'aMC@NLO+Pythia8', 'lp' )
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            #altMC2Histo.Scale(1, 'width')  ###divide by bin width
            altMC2Histo.SetLineWidth(1)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(0.5)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1X0 same")
        else:
            #altMC1Histo.Scale(1, 'width')  ###divide by bin width
            altMC1Histo.SetLineWidth(1)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(0.5)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            legend.AddEntry( altMC1Histo, altMC1Histo_label.replace('-FXFX',''),'lpe')#'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1X0 same")

        
    genJetHisto.Draw( "histE1X0 same")
    altMCHisto.Draw("histE1X0 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1X0 same")
        fsrDownHisto.Draw("PE1X0 same")

    unfoldHisto.Draw( "E1X0 same")
    
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
    
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.80, seltext )

    
    legend.Draw()
    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #   CMS_lumi.extraText = "Simulation Preliminary"
    #   CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
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
    
    can.cd()
    pad2 = ROOT.TPad("pad2", "Ratio",0, 0, 1, Hdw / H,-1);
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
    
    tmpPad2= pad2.DrawFrame( -0.5, 0., unfoldHisto.GetXaxis().GetBinLowEdge(unfoldHisto.GetNbinsX()+1), 1.9 )
    print (labelX)
    
    tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    tmpPad2.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )

    #tmpPad2.GetYaxis().CenterTitle()
    #tmpPad2.SetLabelSize(0.13, 'x')
    #tmpPad2.SetTitleSize(0.12, 'x')
    #tmpPad2.SetLabelSize(0.12, 'y')
    #tmpPad2.SetTitleSize(0.12, 'y')
    #mpPad2.SetNdivisions(505, 'x')
    #mpPad2.SetNdivisions(505, 'y')
    
    tmpPad2.GetYaxis().SetTitleOffset(extraSpace + (0.8 ) * Hdw / H_ref)
    tmpPad2.GetXaxis().SetTitleOffset(0.92)
    tmpPad2.SetTitleSize((0.056 if noNorm else 0.054) * H_ref / Hdw, "Y")
    tmpPad2.SetLabelSize(0.048 * H_ref / Hdw, "Y")
    tmpPad2.SetTitleSize(0.054 * H_ref / Hdw, "X")
    tmpPad2.SetLabelSize(0.048 * H_ref / Hdw, "X")
    tmpPad2.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    tmpPad2.GetXaxis().SetTitle( f'{labelX} N-subjettiness basis' )
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}"  )

    tmpPad2.GetYaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    tmpPad2.GetYaxis().SetNdivisions(505)
    tmpPad2.GetYaxis().CenterTitle()
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    pad2.cd()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.1,2.1 )
    
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (0.92 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.95)
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.048 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.048 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.012 * H_ref / Hdw, "X")
    ratio_totalUnc.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')#'#'+labelX.split('#')[1] )
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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
    #set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.7,0.4)
    
    ratio_totalUnc.GetXaxis().SetRangeUser(0., 2.)#-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )

    ratio_totalUnc.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerSize(0.5)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerSize(0.5)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerSize(0.55)
        hRatio3.SetMarkerStyle(23)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerSize(0.55)
        hRatio4.SetMarkerStyle(22)
        #hRatio3.SetMarkerSize(1)
        

        #hRatio4.SetMarkerSize(1)
        
        
    #hRatio.SetMarkerSize(1)
    hRatio.Draw('P0X0 same')
    
    #hRatio2.SetMarkX0erSize(1)
    hRatio2.Draw('P0 same')
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerSize(0.5)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
        #hRatio5.SetMarkerSize(1)
        hRatio5.Draw('P0X0 same')

    if includeFSR:
        hRatio3.Draw('P0X0 same')
        hRatio4.Draw('P0X0 same')
    
        
    
    
    ratioLegend=ROOT.TLegend(0.15,0.78,0.60,0.88)
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
    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)     
    
    
def drawUnfold_Ndim_Plot1D(ivar, process, 
                    lumi, 
                    varDict,
                    genBinMap,
                    dataJetHisto, 
                    genJetHisto, 
                    unfoldHisto, 
                    unfoldHistoStatUnc, 
                    unfoldHistowoUnc,
                    altMCHisto, 
                    foldHisto, 
                    recoJetHisto, 
                    cov_tot, 
                    cov_datastat_tot, 
                    labelX, 
                    maxX, 
                    tlegendAlignment, 
                    outputDir,
                    outputName,
                    year='all',
                    selection='_dijetSel',
                    altMC1Histo = None, 
                    altMC2Histo = None, 
                    altMC1Histo_label = None, 
                    altMC2Histo_label = None, 
                    nomMCHisto_label = None, altMCHisto_label = None,
                    extraMC=False,
                    includeFSR = False, fsrUpHisto = None, fsrDownHisto=False, 
                    noNorm=False,version='April25',
                    signalLabel = 'MLMQCD_HT2000toInf',
                    alt0SignalLabel = 'H7MLMQCD_HT2000toInf',
                    alt1SignalLabel = 'QCD_HT2000toInf',
                    alt2SignalLabel = 'QCD_Pt_3200toInf',
                    fsrLabel = 'sysMLMQCD_fsrWeight_HT2000toInf',
                    n_obs=25, 
                    return1DHistDict = False
                                    ):

    if return1DHistDict:
        dict_1DHists = OrderedDict()
        for var in varDict.keys():
            dict_1DHists[var] = OrderedDict()
            
    
    
    
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 800
    H_ref = 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
    Tup = T_ref * H_ref / Hup
    Tdw = M_ref * H_ref / Hdw
    Bup = 0.022
    Bdw = B_ref * H_ref / Hdw

    can = ROOT.TCanvas('canUnfoldingCombined'+ivar, 'canUnfolding'+ivar,  50, 50, W, H)
    can.SetFillColor(0)
    can.SetBorderMode(0)
    can.SetFrameFillStyle(0)
    can.SetFrameBorderMode(0)
    can.SetFrameLineColor(0)
    can.SetFrameLineWidth(0)
    

    pad1 = ROOT.TPad("pad1Combined"+ivar, "Main",0, Hdw / H, 1, 1, -1)
    
    
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
    
    legend=ROOT.TLegend(0.72,0.60,0.92,0.90)

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
    
    #cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHisto, 
    #                                                  cov_tot.Clone(),ndim=n_obs)
    #cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoDataStatErr, 
    #                                                              cov_datastat_tot.Clone(),ndim=n_obs)
    
    cov_normTot_np, normed_covTot = get_normalised_cov_combined(unfoldHisto, cov_tot.Clone(), genBinMap)
    #cov_datastat_tot.Clone#get_normalised_cov(unfoldHisto, 
    #cov_tot.Clone(),ndim=n_obs)
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov_combined(unfoldHisto, cov_datastat_tot.Clone(), genBinMap)
    #print(cov_normTot_np)
    #print(cov_norm_dataStat_np)
    

    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()
    
    unfoldHistowoUnc = normalize_combined_TH1_by_blocks(unfoldHistowoUnc, genBinMap, noNorm=noNorm)
    unfoldHistoDataStatErr = normalize_combined_TH1_by_blocks(unfoldHistoDataStatErr, genBinMap, noNorm=noNorm)
    
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_norm_dataStat_np)
    else:
        cov_abs_dataStat_np, _ = th2_to_ndarray(cov_datastat_tot.Clone())
        get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_abs_dataStat_np)
        

    unfoldHisto = normalize_combined_TH1_by_blocks(unfoldHisto, genBinMap, noNorm=noNorm)
    if not(noNorm): 
        get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    else:
        cov_absTot_np, _ = th2_to_ndarray(cov_tot.Clone())
        get_th1_normedCovErrors(unfoldHisto, cov_absTot_np)
        
    #unfoldHisto.Scale(1.,'width')

    
    dataJetHisto = normalize_combined_TH1_by_blocks(dataJetHisto, genBinMap, noNorm=False)
    genJetHisto = normalize_combined_TH1_by_blocks(genJetHisto, genBinMap, noNorm=False)
    altMCHisto = normalize_combined_TH1_by_blocks(altMCHisto, genBinMap, noNorm=False)
    foldHisto = normalize_combined_TH1_by_blocks(foldHisto, genBinMap, noNorm=False)
    recoJetHisto = normalize_combined_TH1_by_blocks(recoJetHisto, genBinMap, noNorm=False)

    if includeFSR:
        fsrUpHisto.Sumw2()
        fsrUpHisto = normalize_combined_TH1_by_blocks(fsrUpHisto, genBinMap, noNorm=False)
        fsrDownHisto.Sumw2()
        fsrDownHisto = normalize_combined_TH1_by_blocks(fsrDownHisto, genBinMap, noNorm=False)

    if extraMC:
        altMC1Histo.Sumw2()
        altMC1Histo = normalize_combined_TH1_by_blocks(altMC1Histo, genBinMap, noNorm=False)
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo = normalize_combined_TH1_by_blocks(altMC2Histo, genBinMap, noNorm=False)
            altMC2Histo.SetTitle("")
    
    
    

    
    #for var in varDict.keys(): 
    physical_unfoldHistoDataStatErr = physical_histograms_from_combined_Ndim(combined_hist=unfoldHistoDataStatErr.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix='unfoldHistoDataStatUnc',
                                                                ) 
    physical_unfoldHisto = physical_histograms_from_combined_Ndim(combined_hist=unfoldHisto.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix='unfoldHisto',
                                                                ) 

    physical_genJetHisto = physical_histograms_from_combined_Ndim(combined_hist=genJetHisto.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix = signalLabel+'_gen',#'MLMQCD_HT2000toInf',
                                                                     withSuff=True, sys='_nom'
                                                                ) 
    physical_altMCHisto = physical_histograms_from_combined_Ndim(combined_hist=altMCHisto.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix = alt0SignalLabel+'_gen',#'H7MLMQCD_HT2000toInf',
                                                                     withSuff=True, sys='_nom'
                                                                ) 
    physical_fsrUpHisto = physical_histograms_from_combined_Ndim(combined_hist=fsrUpHisto.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix=fsrLabel+'_gen',
                                                                     withSuff=True, sys='_fsrWeightUp'
                                                                ) 
    physical_fsrDownHisto = physical_histograms_from_combined_Ndim(combined_hist=fsrDownHisto.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix=fsrLabel+'_gen',
                                                                     withSuff=True, sys='_fsrWeightDown'
                                                                ) 
    if extraMC:
        physical_altMC1Histo = physical_histograms_from_combined_Ndim(combined_hist=altMC1Histo.Clone(),
                                                                     bin_map=genBinMap,
                                                                     varDict=varDict,
                                                                     new_hist_prefix = alt1SignalLabel+'_gen',#'QCD_HT2000toInf',
                                                                     withSuff=True, sys='_nom'
                                                                    ) 
        if 'dijet' in selection and altMC2Histo:

            physical_altMC2Histo = physical_histograms_from_combined_Ndim(combined_hist=altMC2Histo.Clone(),
                                                                 bin_map=genBinMap,
                                                                 varDict=varDict,
                                                                 new_hist_prefix = alt2SignalLabel+'_gen',#'QCD_Pt_3200toInf',
                                                                     withSuff=True, sys='_nom'
                                                                ) 
            
    for var in varDict.keys():
        outName = outputDir+var+selection+'_fromMBody_Data'+signalLabel+'_TUnfold_'+version+'.pdf'
        #.split('6bodyOC_dijetSel_fromDataMLMQCD_HT2000toInf_TUnfold__March25_Ndim_dataWithBkgCorr.pdf')
        print(outName,varDict[var]['alignLeg'])
        drawUnfolded1DfromNDim(var, selection, process, year, lumi,
                               physical_genJetHisto[var].Clone(), 
                               physical_unfoldHistoDataStatErr[var].Clone(),
                               physical_unfoldHisto[var].Clone(),
                               physical_altMCHisto[var].Clone(), 
                               labelX=varDict[var]['label'], 
                               maxX=varDict[var]['bins'][-1],
                               tlegendAlignment=varDict[var]['alignLeg'], 
                               outputName=outName,
                               altMC1Histo = physical_altMC1Histo[var].Clone(), 
                               altMC2Histo = physical_altMC2Histo[var].Clone() if 'dijet' in selection and altMC2Histo else None, 
                               altMC1Histo_label=altMC1Histo_label,
                               altMC2Histo_label=altMC2Histo_label if 'dijet' in selection and altMC2Histo else None, 
                               nomMCHisto_label=nomMCHisto_label,
                               altMCHisto_label=altMCHisto_label,
                               extraMC=True, includeFSR = True, 
                               fsrUpHisto = physical_fsrUpHisto[var].Clone(),
                               fsrDownHisto=physical_fsrDownHisto[var].Clone(), noNorm=True, version=version
                              )
        if return1DHistDict:
            dict_1DHists[var]['selection'] = selection 
            dict_1DHists[var]['process'] = process 
            dict_1DHists[var]['year'] = year 
            dict_1DHists[var]['lumi'] = lumi
            dict_1DHists[var]['genJetHisto'] = physical_genJetHisto[var].Clone() 
            dict_1DHists[var]['genJetHisto'].Sumw2()
            dict_1DHists[var]['genJetHisto'].SetDirectory(0)
            
            dict_1DHists[var]['unfoldHistoDataStatErr'] = physical_unfoldHistoDataStatErr[var].Clone()
            dict_1DHists[var]['unfoldHistoDataStatErr'].Sumw2()
            dict_1DHists[var]['unfoldHistoDataStatErr'].SetDirectory(0)
            
            dict_1DHists[var]['unfoldHisto'] = physical_unfoldHisto[var].Clone()
            dict_1DHists[var]['unfoldHisto'].Sumw2()
            dict_1DHists[var]['unfoldHisto'].SetDirectory(0)
            
            dict_1DHists[var]['unfoldHistowoUnc'] = physical_unfoldHisto[var].Clone()
            for ibin in range(1,dict_1DHists[var]['unfoldHistowoUnc'].GetNbinsX()+1):
                dict_1DHists[var]['unfoldHistowoUnc'].SetBinError(ibin,0)
                
            dict_1DHists[var]['unfoldHistowoUnc'].Sumw2()
            dict_1DHists[var]['unfoldHistowoUnc'].SetDirectory(0)
            
            dict_1DHists[var]['altMCHisto'] = physical_altMCHisto[var].Clone()
            dict_1DHists[var]['altMCHisto'].Sumw2()
            dict_1DHists[var]['altMCHisto'].SetDirectory(0)
            
            dict_1DHists[var]['altMC1Histo'] = physical_altMC1Histo[var].Clone()
            dict_1DHists[var]['altMC1Histo'].Sumw2()
            dict_1DHists[var]['altMC1Histo'].SetDirectory(0)
            
            dict_1DHists[var]['fsrUpHisto'] = physical_fsrUpHisto[var].Clone()
            dict_1DHists[var]['fsrUpHisto'].Sumw2()
            dict_1DHists[var]['fsrUpHisto'].SetDirectory(0)
            
            dict_1DHists[var]['fsrDownHisto'] = physical_fsrDownHisto[var].Clone()
            dict_1DHists[var]['fsrDownHisto'].Sumw2()
            dict_1DHists[var]['fsrDownHisto'].SetDirectory(0)
            
            
            
            dict_1DHists[var]['nomMCHisto_label'] = nomMCHisto_label 
            dict_1DHists[var]['altMCHisto_label'] = altMCHisto_label 
            dict_1DHists[var]['altMC1Histo_label'] = altMC1Histo_label 
            
            
            if 'dijet' in selection:
                dict_1DHists[var]['altMC2Histo'] = physical_altMC2Histo[var].Clone() if altMC2Histo!=None else None
                dict_1DHists[var]['altMC2Histo'].Sumw2()
                dict_1DHists[var]['altMC2Histo'].SetDirectory(0)
                dict_1DHists[var]['altMC2Histo_label'] = altMC2Histo_label if altMC2Histo!=None else None
            

            #dict_1DHists[var]['ratioUncHisto'] = ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, 
            dict_1DHists[var]['labelX'] = varDict[var]['label'] 
            dict_1DHists[var]['maxX'] = varDict[var]['bins'][-1]
            dict_1DHists[var]['tlegendAlignment'] = varDict[var]['alignLeg'][-1] 
            
            dict_1DHists[var]['noNorm'] = True 
            

    can.cd()
    pad1.cd()
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetLineWidth(1)
    unfoldHisto.SetMarkerSize(0.5)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(1)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(0.5)
    #genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(0.55)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(23)


        fsrDownHisto.SetMarkerSize(0.55)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(22)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )
    
    unfoldHisto.GetXaxis().SetTitleOffset(999)    
    unfoldHisto.GetXaxis().SetLabelOffset(999)    
    #unfoldHisto.GetYaxis().SetLabelOffset(0.011 * H_ref / Hup)    

    #unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+1.15*Hup/H_ref)    
    #unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    #unfoldHisto.GetYaxis().SetLabelSize(0.047* H_ref / Hup)
    
    unfoldHisto.GetYaxis().SetTitleOffset(extraSpace+0.92*Hup/H_ref)    
    unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
    unfoldHisto.GetYaxis().SetLabelSize(0.046* H_ref / Hup)
    unfoldHisto.GetYaxis().SetTitleFont(42)


    #unfolded_combined_nom.SetTickLength(0.03, "XY")  #?? ok if 1/3
    unfoldHisto.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetXaxis().GetBinLowEdge(unfoldHisto.GetNbinsX()+1))#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfoldHisto.GetYaxis().SetTickLength(0.03)#* H_ref / Hup)#, "Y")  #?? ok if 1/3
    unfoldHisto.GetXaxis().SetTickLength(0.03)#* H_ref / Hdw)#, "X")

   
    if 'body' in labelX: 
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#tau_{N}^{(#beta)}}' if not noNorm else 'N_{events}')
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
    #unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( (1.75 if not('dijet' in selection) else 1.55 )*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3
    
    if noNorm:
        unfoldHisto.GetYaxis().SetMaxDigits(3)  #?? ok if 1/3
        ROOT.TGaxis.SetExponentOffset(-0.06, 0.005, "y")
    unfoldHisto.Draw( "AXIS")
    
    
    can.Update()
    
    
    
    #altMCHisto.Scale(1, 'width')  ###divide by bin width
    altMCHisto.SetLineWidth(1)
    altMCHisto.SetMarkerSize(0.5)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    #altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 
        #due to historial reasons, confusing fsrDown/Up naming corresponds to energy scale down/up, ie, \alpha_S^FSR up/down!
        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')

        
        
        
    legend.AddEntry( altMCHisto, 'MG5-MLM+H7' if 'dijet' in selection else'PWHG+H7','lpe')#'aMC@NLO+Pythia8', 'lp' )
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            #altMC2Histo.Scale(1, 'width')  ###divide by bin width
            altMC2Histo.SetLineWidth(1)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            #altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(0.5)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            #altMC1Histo.Scale(1, 'width')  ###divide by bin width
            altMC1Histo.SetLineWidth(1)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            #altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(0.5)
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
    
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.15 if tlegendAlignment.startswith('right') else 0.53+dijetOffset ), 0.80, seltext )
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    
    legend.Draw()
    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
    else:
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+ year
    #else:
    #   CMS_lumi.extraText = "Simulation Preliminary"
    #   CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
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
    
    
    
    tmpPad2= pad2.DrawFrame( -0.5, 0., unfoldHisto.GetXaxis().GetBinLowEdge(unfoldHisto.GetNbinsX()+1), 2.1 )
    print (labelX)
    
    if not 'dijet' in selection:
        tmpPad2.GetYaxis().SetRangeUser(0.2,1.99 )
    else:
        tmpPad2.GetYaxis().SetRangeUser(0.4,1.9 )

    tmpPad2.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )

    #tmpPad2.GetYaxis().CenterTitle()
    #tmpPad2.SetLabelSize(0.13, 'x')
    #tmpPad2.SetTitleSize(0.12, 'x')
    #tmpPad2.SetLabelSize(0.12, 'y')
    #tmpPad2.SetTitleSize(0.12, 'y')
    #mpPad2.SetNdivisions(505, 'x')
    #mpPad2.SetNdivisions(505, 'y')
    
    tmpPad2.GetYaxis().SetTitleOffset(extraSpace + (0.8 ) * Hdw / H_ref)
    tmpPad2.GetXaxis().SetTitleOffset(0.92)
    tmpPad2.SetTitleSize((0.054 if noNorm else 0.053) * H_ref / Hdw, "Y")
    tmpPad2.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    tmpPad2.SetTitleSize(0.054 * H_ref / Hdw, "X")
    tmpPad2.SetLabelSize(0.046 * H_ref / Hdw, "X")
    tmpPad2.SetLabelOffset(0.011 * H_ref / Hdw, "X")
    #tmpPad2.SetLabelOffset(0.011 * H_ref / Hdw, "Y")
    tmpPad2.GetXaxis().SetTitle( f'{labelX} N-subjettiness basis' )
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}"  )

    tmpPad2.GetYaxis().SetTitleFont(42)
    tmpPad2.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    tmpPad2.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    tmpPad2.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
    tmpPad2.GetYaxis().SetNdivisions(505)
    tmpPad2.GetYaxis().CenterTitle()
    
    
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    pad2.cd()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.8)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )
    
    if not 'dijet' in selection:
        ratio_totalUnc.GetYaxis().SetRangeUser(0.2,1.99 )
    else:
        ratio_totalUnc.GetYaxis().SetRangeUser(0.4,1.9 )
    
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (0.92 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.95)
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.011 * H_ref / Hdw, "X")
    #ratio_totalUnc.SetLabelOffset(0.011 * H_ref / Hdw, "Y")
    ratio_totalUnc.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')#'#'+labelX.split('#')[1] )
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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
    #set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.7,0.4)
    
    #ratio_totalUnc.GetXaxis().SetRangeUser(0., 2.)#-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )

    ratio_totalUnc.GetXaxis().SetRangeUser(-0.5, unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+1) )
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerSize(0.5)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerSize(0.5)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerSize(0.5)
        hRatio3.SetMarkerStyle(23)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerSize(0.5)
        hRatio4.SetMarkerStyle(22)
        
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerSize(0.5)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
        #hRatio5.SetMarkerSize(1)
        
    hRatio5.Draw('PE1 same')
    #hRatio.SetMarkerSize(1)
    hRatio.Draw('PE1 same')
    
    #hRatio2.SetMarkerSize(1)
    hRatio2.Draw('PE1 same')
    
    
    if includeFSR:
        #hRatio3.SetMarkerSize(1)
        hRatio3.Draw('PE1 same')

        #hRatio4.SetMarkerSize(1)
        hRatio4.Draw('PE1 same')
    
    
    
    ratioLegend=ROOT.TLegend(0.15,0.78,0.60,0.88)
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
    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12) 
    
    if return1DHistDict:
        return dict_1DHists

def drawUnfolded1DfromNDim(ivar, selection, process, year, lumi,
                           genJetHisto, 
                           unfoldHistoDataStatUnc, unfoldHistoTotUnc, 
                           altMCHisto, 
                           labelX, maxX, tlegendAlignment, outputName, version,
                           altMC1Histo = None, altMC2Histo = None, 
                           altMC1Histo_label = None, altMC2Histo_label = None, 
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
    #Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
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
    
    unfoldHistowoUnc = unfoldHistoTotUnc.Clone('unfoldHistowoUnc'+ivar)
    unfoldHistowoUnc.Sumw2()
    for i in range(1,unfoldHistowoUnc.GetNbinsX()+1):
        unfoldHistowoUnc.SetBinError(i, 0)
    
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
    
    
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())
    
    
    

    #genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    #altMCHisto.Sumw2()

    unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')

    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.), 'width')
    
    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.),'width')

    
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.),'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.),'width')   
    
    
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
    unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
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
    unfoldHisto.SetMaximum(  1.7*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )
    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3

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
        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')

        
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'PWHG+H7','lpe')#
    
    
    
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
            CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
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
    ratio_datastatUnc.SetFillStyle(3245)#if not('dijet' in selection) else 3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.94)
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.011 * H_ref / Hdw, "X")
    if 'tau' in ivar: ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: ratio_totalUnc.GetXaxis().SetTitle(labelX)
        
    
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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
    
    
    
    hRatio.SetMarkerSize(1)
    hRatio.Draw('PE1 same')
    
    hRatio2.SetMarkerSize(1)
    hRatio2.Draw('PE1 same')
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
        hRatio5.SetMarkerSize(1)
        hRatio5.Draw('PE1 same')

    if includeFSR:
        hRatio3.SetMarkerSize(1)
        hRatio3.Draw('PE1 same')

        hRatio4.SetMarkerSize(1)
        hRatio4.Draw('PE1 same')
    
    
    ratioLegend=ROOT.TLegend(0.19,0.78,0.69,0.88)
    ratioLegend.SetTextSize(0.09)
    ratioLegend.SetTextFont(42)
    ratioLegend.SetNColumns(2)
    ratioLegend.SetFillStyle(0)#ColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    ratioLegend.AddEntry( ratio_totalUnc, 'Total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    ratioLegend.Draw()
    print(outputName)
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    can.SaveAs(root_macro)

    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)    


def drawUnfolded1DfromNDim_vs_singleObs(ivar, selection, process, year, lumi,
                                        genJetHisto, 
                                        unfoldHistoDataStatUnc, unfoldHistoTotUnc, 
                                        unfoldHisto_old,
                                        altMCHisto, 
                                        labelX, maxX, tlegendAlignment, outputName, version,
                                        altMC1Histo = None, altMC2Histo = None, 
                                        altMC1Histo_label = None, altMC2Histo_label = None, 
                                        nomMCHisto_label = None, altMCHisto_label = None,
                                        extraMC=False, includeFSR = False, 
                                        fsrUpHisto = None, fsrDownHisto=None, noNorm=False
                                      ):
    

    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    #ROOT.gStyle.SetPadRightMargin(0.04)
    #ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    extraSpace = 0.02
    #Set canvas dimensions and margins
    W_ref = 700 #if square else 800
    H_ref = 600 #if square else 500
    #Set bottom pad relative height and relative margin
    F_ref = 1.0 / 3.0
    M_ref = 0.03
    #Set reference margins
    T_ref = 0.07
    B_ref = 0.13
    L = 0.15 #if square else 0.12
    R = 0.05
    #Calculate total canvas size and pad heights
    W = W_ref
    H = int(H_ref * (1 + (1 - T_ref - B_ref) * F_ref + M_ref))
    Hup = H_ref * (1 - B_ref)
    Hdw = H - Hup
    #references for T, B, L, R
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
    
    unfoldHistowoUnc = unfoldHistoTotUnc.Clone('unfoldHistowoUnc'+ivar)
    unfoldHistowoUnc.Sumw2()
    for i in range(1,unfoldHistowoUnc.GetNbinsX()+1):
        unfoldHistowoUnc.SetBinError(i, 0)
    
    unfoldHisto = unfoldHistoTotUnc.Clone('unfoldHisto'+ivar)
    unfoldHisto.Sumw2()
 
    unfoldHistoDataStatErr=unfoldHistoDataStatUnc.Clone('unfoldHistoStatUnc'+ivar)
    unfoldHistoDataStatErr.Sumw2()
    
    
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    
    unfoldHisto_old.SetTitle("")
    print("unfoldHisto_old.Integral()",unfoldHisto_old.Integral())
    
    unfoldHistoDataStatErr.SetTitle("")
    
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    
    
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())
    
    
    

    #genJetHisto.Sumw2()
    #unfoldHistowoUnc.Sumw2()
    #altMCHisto.Sumw2()

    unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')
    unfoldHisto_old.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')

    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.), 'width')
    
    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.),'width')

    
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.),'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.),'width')   
    
    
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
    unfoldHisto.GetYaxis().SetTitleSize(0.054* H_ref / Hup)
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
    unfoldHisto.SetMaximum(  1.7*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )
    unfoldHisto.SetTickLength(0.03, "XY")  #?? ok if 1/3

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
        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} up", 'pe')

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "#alpha_{S}^{FSR} down", 'pe')

        
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'PWHG+H7','lpe')#
    
    
    
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
    
    unfoldHisto_old.SetMarkerStyle(27)
    unfoldHisto_old.SetMarkerSize(1)
    unfoldHisto_old.SetMarkerColor(ROOT.kGray+1)
    unfoldHisto_old.SetLineColor(ROOT.kGray+1)
    unfoldHisto_old.SetLineStyle(2)
    unfoldHisto_old.SetLineWidth(2)
    legend.AddEntry( unfoldHisto_old, 'Data (1D unf.)', 'lpe' )
    unfoldHisto_old.Draw( "histE1 same")

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
            CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
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
    ratio_datastatUnc.SetFillStyle(3245)#if not('dijet' in selection) else 3245)
    #ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    #ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    #ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )
    
    #ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset(extraSpace + (1.13 ) * Hdw / H_ref)
    ratio_totalUnc.GetXaxis().SetTitleOffset(0.94)
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "Y")
    ratio_totalUnc.SetTitleSize(0.054 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelSize(0.046 * H_ref / Hdw, "X")
    ratio_totalUnc.SetLabelOffset(0.011 * H_ref / Hdw, "X")
    if 'tau' in ivar: ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: ratio_totalUnc.GetXaxis().SetTitle(labelX)
        
    
    ratio_totalUnc.GetYaxis().SetTitle( "Ratio")#frac{Old.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleFont(42)
    ratio_totalUnc.GetXaxis().SetTitleFont(42)

    #tmpPad2.GetXaxis().SetTitle(nameXaxis)
    #tmpPad2.GetYaxis().SetTitle(nameRatio)

    #Set tick lengths to match original (these are fractions of axis length)
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hup, "Y")  #?? ok if 1/3
    ratio_totalUnc.SetTickLength(0.03 * H_ref / Hdw, "X")

    #Reduce divisions to match smaller height (default n=510, optim=kTRUE)
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
    if not('dijet' in selection):
        set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.42,0.6)
    else:
        set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.37,0.75)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   
    hRatioVsOld = ROOT.TGraphAsymmErrors()
    hRatioVsOld.Divide( unfoldHisto_old, unfoldHisto, 'pois' )
    hRatioVsOld.SetLineColor(ROOT.kGray+2)#ROOT.kRed)
    hRatioVsOld.SetMarkerColor(ROOT.kGray+2)#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatioVsOld.SetMarkerStyle(27)
    hRatioVsOld.SetMarkerSize(1)
    hRatioVsOld.Draw('P0 same')
    
    ratioLegend=ROOT.TLegend(0.19,0.70,0.85,0.88)
    ratioLegend.SetTextSize(0.077)
    ratioLegend.SetTextFont(42)
    ratioLegend.SetNColumns(3)
    ratioLegend.SetFillStyle(0)#ColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    ratioLegend.AddEntry( ratio_totalUnc, 'Total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    ratioLegend.AddEntry( hRatioVsOld, '#frac{1D unf. (old) }{1D from comb. unf.}', 'pe' )
    ratioLegend.Draw()
    print(outputName)
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    root_macro = outputName.split('.pdf')[0]+'.C'
    can.SaveAs(root_macro)

    #ROOT.gStyle.SetPadRightMargin(0.09)     ##reseating
    #ROOT.gStyle.SetPadLeftMargin(0.12)        

def drawUncertainties_from_err_shifts_1DfromNDim(
                                                    ivar,    
                                                    unfoldHistoTotUnc,
                                                    unfoldHistoDataStatUnc,
                                                    unfoldHistowoUnc,
                                                    unfoldHistoRMStatUnc,
                                                    unfoldHistoBkgSubUnc,
                                                    uncerUnfoldHisto, 
                                                    labelX, 
                                                    tlegendAlignment, 
                                                    outputName, year, 
                                                    selection, lumi, with_modelUnc=True 
                                                ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    unftot = unfoldHistowoUnc.Integral()

    print (f'|------> Procesing uncertainty plot for {ivar} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
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
    #   c.SetRightMargin(B / W + 0.03)
    canUnc.SetTopMargin(T / H)
    canUnc.SetBottomMargin(B / H + 0.02)

    legend=ROOT.TLegend(0.18,0.64,0.88,0.9)

    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.038)
    legend.SetTextFont(42)

    legend.SetBorderSize(0)
    
    
    
    unfoldHistoTotUnc.Scale(1.,'width')
    unfoldHistowoUnc.Scale(1.,'width')
    unfoldHistoDataStatUnc.Scale(1.,'width')
    unfoldHistoRMStatUnc.Scale(1.,'width')
    unfoldHistoBkgSubUnc.Scale(1.,'width')

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
            print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            jesHistoUpMax.Scale(1.,'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            jesHistoDownMax.Scale(1.,'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            jerHistoUpMax.Scale(1.,'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            jerHistoDownMax.Scale(1.,'width')
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
            btagHistoUpMax.Scale(1.,'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            btagHistoDownMax.Scale(1.,'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()) and not ('bkg' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #   print(k)
            #   btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1.,'width')#./(unftot if norming else 1.)
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
            normeduncerUnfoldHistoshiftsDown[k].Scale(1.,'width')#./(unftot if norming else 1.)
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
        modelUnc.Scale(1.,'width')#
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
    totalErrHist.GetYaxis().SetTitleOffset(1.08)
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.35,0.95 if ('dijet' in selection) else 0.88)
    
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
        
        
        if btagUncIncluded:#and not(btag_key!=None):
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
        #   legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
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
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], "ISR", 'p' )
                
            elif 'fsr' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], "FSR", 'p' )

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
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
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

def drawUncertainties_from_err_shifts_theoryVariations_1DfromNDim(
                                                                ivar,    
                                                                unfoldHistoTotUnc,
                                                                unfoldHistoDataStatUnc,
                                                                unfoldHistowoUnc,
                                                                unfoldHistoRMStatUnc,
                                                                unfoldHistoBkgSubUnc,
                                                                uncerUnfoldHisto, 
                                                                labelX, 
                                                                tlegendAlignment, 
                                                                outputName, year, 
                                                                selection, lumi 
                                                                ):
    
    unftot = unfoldHistowoUnc.Integral()

    print (f'|------> Procesing theory/model variation uncertainty plot for {ivar}, 1D slice from Ndim  ')
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
    #   c.SetRightMargin(B / W + 0.03)
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
    
    
    up_counter=0
    down_counter=0
    col_counter=0
    #col_counter_jes=0
    
    
    unfoldHistoTotUnc.Scale(1.,'width')
    unfoldHistowoUnc.Scale(1.,'width')
    unfoldHistoDataStatUnc.Scale(1.,'width')
    unfoldHistoRMStatUnc.Scale(1.,'width')
    unfoldHistoBkgSubUnc.Scale(1.,'width')

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
    totalErrHist.GetYaxis().SetTitleSize(0.056)
    totalErrHist.GetYaxis().SetLabelSize(0.047)    
    totalErrHist.GetYaxis().SetLabelOffset(0.012)    
    totalErrHist.GetYaxis().SetLabelFont(42)    
    totalErrHist.GetYaxis().SetTitleFont(42)
    totalErrHist.GetYaxis().SetTitleOffset(1.08)
    
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.35,0.95 if ('dijet' in selection) else 0.88)
    
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
    dataStatErrHist.SetFillStyle(3245)
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
    ##h2.SetLineWidth(2)
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
            cr_histos[k].Scale(1.,'width')
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
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            normeduncerUnfoldHistoshiftsUp[k].Scale(1.,'width')
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                
                up_counter=up_counter+1
    
    
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
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            normeduncerUnfoldHistoshiftsDown[k].Scale(1.,'width')
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

            normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    if year=='all': 
        #if 'dijet' in selection:
        CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"
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

def drawUncertainties_from_err_shifts_theoryVariations_Ndim_unitNorm(ivar, 
                                                                     unfoldHistoTotUnc, 
                                                                     unfoldHistowoUnc,
                                                                     unfoldHistoDataStatUnc, 
                                                                     unfoldHistoRMStatUnc, 
                                                                     unfoldHistoBkgSubUnc, 
                                                                     uncerUnfoldHisto, 
                                                                     cov_tot, 
                                                                     cov_datastat_tot, 
                                                                     cov_rmstat_tot, 
                                                                     cov_bkg_tot, 
                                                                     labelX, 
                                                                     tlegendAlignment, 
                                                                     outputName, 
                                                                     unftot, 
                                                                     selection, 
                                                                     genBinMap, 
                                                                     varDict, 
                                                                     outputDir,
                                                                     norming=True,
                                                                     year='all',
                                                                     n_obs=25,
                                                                     
                                                                    ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    
    print (f'|------> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
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
    #   legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistoTotUnc = unfoldHistowoUnc.Clone('unfoldHistoTotUnc'+ivar)
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc = unfoldHistowoUnc.Clone('unfoldHistoDataStatUnc'+ivar)
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc = unfoldHistowoUnc.Clone('unfoldHistoRMStatUnc'+ivar)
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc = unfoldHistowoUnc.Clone('unfoldHistoBkgSubUnc'+ivar)
    unfoldHistoBkgSubUnc.Sumw2()
    
    #print(cov_normTot_np)
    #print(cov_norm_dataStat_np)
    
    cov_normTot_np, normed_covTot = get_normalised_cov_combined(unfoldHistoTotUnc, cov_tot.Clone(), genBinMap)
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov_combined(unfoldHistoDataStatUnc, cov_datastat_tot.Clone(), genBinMap)
    
    
    cov_norm_RMStat_np, normed_cov_RMStat = get_normalised_cov_combined(unfoldHistoRMStatUnc, 
                                                               cov_rmstat_tot.Clone(), genBinMap)
    cov_norm_BkgSub_np, normed_cov_BkgSub = get_normalised_cov_combined(unfoldHistoBkgSubUnc, 
                                                               cov_bkg_tot.Clone(), genBinMap)
    
    
    
    unfoldHistowoUnc = normalize_combined_TH1_by_blocks(unfoldHistowoUnc, genBinMap, noNorm=False)
    unfoldHistoDataStatUnc = normalize_combined_TH1_by_blocks(unfoldHistoDataStatUnc, genBinMap, noNorm=False)
   
    unfoldHistoTotUnc = normalize_combined_TH1_by_blocks(unfoldHistoTotUnc, genBinMap, noNorm=False)
    unfoldHistoRMStatUnc = normalize_combined_TH1_by_blocks(unfoldHistoRMStatUnc, genBinMap, noNorm=False)
    unfoldHistoBkgSubUnc = normalize_combined_TH1_by_blocks(unfoldHistoBkgSubUnc, genBinMap, noNorm=False)
    
    #unfoldHistoTotUnc.Scale(1./((1./n_obs)*unfoldHistoTotUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoTotUnc, cov_normTot_np)
    #unfoldHistoTotUnc.Scale(1.,'width')
    
    #unfoldHistoDataStatUnc.Scale(1./(1./n_obs)*(unfoldHistoDataStatUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoDataStatUnc, cov_norm_dataStat_np)
    #unfoldHistoDataStatUnc.Scale(1.,'width')
    
    #unfoldHistoRMStatUnc.Scale(1./((1./n_obs)*unfoldHistoRMStatUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoRMStatUnc, cov_norm_RMStat_np)
    #unfoldHistoRMStatUnc.Scale(1.,'width')
    
    #unfoldHistoBkgSubUnc.Scale(1./((1./n_obs)*unfoldHistoBkgSubUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoBkgSubUnc, cov_norm_BkgSub_np)
    #unfoldHistoBkgSubUnc.Scale(1.,'width')
    
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

    physical_unfoldHistoTotUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistoTotUnc.Clone(),
                                                                        bin_map=genBinMap,
                                                                        varDict=varDict,
                                                                        new_hist_prefix='unfoldHistoTotUnc',
                                                                        )  
    
    physical_unfoldHistoDataStatUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistoDataStatUnc.Clone(),
                                                                        bin_map=genBinMap,
                                                                        varDict=varDict,
                                                                        new_hist_prefix='unfoldHistoDataStatUnc',
                                                                        )     
    physical_unfoldHistowoUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistowoUnc.Clone(),
                                                                     bin_map=genBinMap,
                                                                     varDict=varDict,
                                                                     new_hist_prefix='unfoldHistowoUnc',
                                                                    )     
    physical_unfoldHistoRMStatUnc = physical_histograms_from_combined_Ndim(combined_hist=unfoldHistoRMStatUnc.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix='unfoldHistoRMStatUnc',
                                                                    )     
    physical_unfoldHistoBkgSubUnc = physical_histograms_from_combined_Ndim(combined_hist=unfoldHistoBkgSubUnc.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix='unfoldHistoBkgSubUnc',
                                                                    )     

    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    
    set_dynamic_y_range_errRatioHist(totalErrHist,1.35,0.90)
    
    totalErrHist.GetXaxis().SetTitle(f'{labelX}')#6-body_OC')
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

    physical_uncerUnfoldHisto = OrderedDict()


    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            #cr_histos[k].Scale(1./((1./n_obs)*cr_histos[k].Integral() if norming else 1.))#,'width')#./(unftot if norming else 1.)         
            
            cr_histos[k] = normalize_combined_TH1_by_blocks(cr_histos[k], genBinMap, noNorm=False)


            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=cr_histos[k].Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'cr_histos_gen_{k}',
                                                                    )     

            
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
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k].Scale(1./((1./n_obs)*normeduncerUnfoldHistoshiftsUp[k].Integral() if norming else 1.))#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = normalize_combined_TH1_by_blocks(normeduncerUnfoldHistoshiftsUp[k], genBinMap, noNorm=False)

            if k in physical_uncerUnfoldHisto.keys():
                print(f"WARNING: sys key: {k} in sys slicing dict being overwritten!!!!!!!!!")

            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'normeduncerUnfoldHistoshiftsUp_gen_{k}',
                                                                    )     
            
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)
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
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k].Scale(1./((1./n_obs)*normeduncerUnfoldHistoshiftsDown[k].Integral() if norming else 1.))#,'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = normalize_combined_TH1_by_blocks(normeduncerUnfoldHistoshiftsDown[k], genBinMap, noNorm=False)
            
            if k in physical_uncerUnfoldHisto.keys():
                print(f"WARNING: sys key: {k} in sys slicing dict being overwritten!!!!!!!!!")

            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'normeduncerUnfoldHistoshiftsDown_gen_{k}',
                                                                    )     

            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
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

    for var in varDict.keys():
        outName = outputName.replace('6bodyOC_',var+'_').replace('_dataWithBkgCorr','')
        print(outName)
        physical_uncerUnfoldHistos = OrderedDict()
        for k in physical_uncerUnfoldHisto:
            key = k.replace('6bodyOC_',var+'_')
            physical_uncerUnfoldHistos[key] = physical_uncerUnfoldHisto[k][var].Clone(physical_uncerUnfoldHisto[k][var].GetName()+'_1DSlice')


        drawUncertainties_from_err_shifts_theoryVariations_1DfromNDim(
                                                                        ivar=var,    
                                                                        unfoldHistoTotUnc=physical_unfoldHistoTotUnc[var].Clone(physical_unfoldHistoTotUnc[var].GetName()+'_1DSlice'),
                                                                        unfoldHistoDataStatUnc=physical_unfoldHistoDataStatUnc[var].Clone(physical_unfoldHistoDataStatUnc[var].GetName()+'_1DSlice'),
                                                                        unfoldHistowoUnc=physical_unfoldHistowoUnc[var].Clone(physical_unfoldHistowoUnc[var].GetName()+'_1DSlice'),
                                                                        unfoldHistoRMStatUnc=physical_unfoldHistoRMStatUnc[var].Clone(physical_unfoldHistoRMStatUnc[var].GetName()+'_1DSlice'),
                                                                        unfoldHistoBkgSubUnc=physical_unfoldHistoBkgSubUnc[var].Clone(physical_unfoldHistoBkgSubUnc[var].GetName()+'_1DSlice'),
            
                                                                        uncerUnfoldHisto=physical_uncerUnfoldHistos, 
                                                                        labelX=varDict[var]['label'], 
                                                                        tlegendAlignment=varDict[var]['alignLeg'], 
                                                                        outputName=outName, 
                                                                        year=year, 
                                                                        selection=selection, 
                                                                        lumi=lumi 
                                                                    )
        
    #legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    #legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
        

    
        
def drawUncertainties_from_err_shifts_Ndim_unitNorm(ivar, 
                                                    unfoldHistoTotUnc, 
                                                    unfoldHistowoUnc, 
                                                    unfoldHistoDataStatUnc, 
                                                    unfoldHistoRMStatUnc, 
                                                    unfoldHistoBkgSubUnc, 
                                                    uncerUnfoldHisto, 
                                                    cov_tot, 
                                                    cov_datastat_tot, 
                                                    cov_rmstat_tot, 
                                                    cov_bkg_tot, 
                                                    labelX, 
                                                    tlegendAlignment, 
                                                    outputName, 
                                                    unftot, 
                                                    selection,
                                                    genBinMap,varDict,
                                                    outputDir, 
                                                    with_modelUnc=True, 
                                                    norming=False,
                                                    year='all',
                                                    n_obs=25,
                                                    
                                                ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
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
    #   legend=ROOT.TLegend(0.2,0.65,0.8,0.9)
    #else: 
    #   legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistoTotUnc = unfoldHistowoUnc.Clone('unfoldHistoTotUnc'+ivar)
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc = unfoldHistowoUnc.Clone('unfoldHistoDataStatUnc'+ivar)
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc = unfoldHistowoUnc.Clone('unfoldHistoRMStatUnc'+ivar)
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc = unfoldHistowoUnc.Clone('unfoldHistoBkgSubUnc'+ivar)
    unfoldHistoBkgSubUnc.Sumw2()
    
    
    cov_normTot_np, normed_covTot = get_normalised_cov_combined(unfoldHistoTotUnc, cov_tot.Clone(), genBinMap)
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov_combined(unfoldHistoDataStatUnc, cov_datastat_tot.Clone(), genBinMap)
    
    
    cov_norm_RMStat_np, normed_cov_RMStat = get_normalised_cov_combined(unfoldHistoRMStatUnc, 
                                                               cov_rmstat_tot.Clone(), genBinMap)
    cov_norm_BkgSub_np, normed_cov_BkgSub = get_normalised_cov_combined(unfoldHistoBkgSubUnc, 
                                                               cov_bkg_tot.Clone(), genBinMap)
    
    
    
    unfoldHistowoUnc = normalize_combined_TH1_by_blocks(unfoldHistowoUnc, genBinMap, noNorm=False)
    unfoldHistoDataStatUnc = normalize_combined_TH1_by_blocks(unfoldHistoDataStatUnc, genBinMap, noNorm=False)
   
    unfoldHistoTotUnc = normalize_combined_TH1_by_blocks(unfoldHistoTotUnc, genBinMap, noNorm=False)
    unfoldHistoRMStatUnc = normalize_combined_TH1_by_blocks(unfoldHistoRMStatUnc, genBinMap, noNorm=False)
    unfoldHistoBkgSubUnc = normalize_combined_TH1_by_blocks(unfoldHistoBkgSubUnc, genBinMap, noNorm=False)
    
    #unfoldHistoTotUnc.Scale(1./((1./n_obs)*unfoldHistoTotUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoTotUnc, cov_normTot_np)
    #unfoldHistoTotUnc.Scale(1.,'width')
    
    #unfoldHistoDataStatUnc.Scale(1./(1./n_obs)*(unfoldHistoDataStatUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoDataStatUnc, cov_norm_dataStat_np)
    #unfoldHistoDataStatUnc.Scale(1.,'width')
    
    #unfoldHistoRMStatUnc.Scale(1./((1./n_obs)*unfoldHistoRMStatUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoRMStatUnc, cov_norm_RMStat_np)
    #unfoldHistoRMStatUnc.Scale(1.,'width')
    
    #unfoldHistoBkgSubUnc.Scale(1./((1./n_obs)*unfoldHistoBkgSubUnc.Integral() if norming else 1.))
    if norming: get_th1_normedCovErrors(unfoldHistoBkgSubUnc, cov_norm_BkgSub_np)
    #unfoldHistoBkgSubUnc.Scale(1.,'width')
    
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

    physical_uncerUnfoldHisto = OrderedDict()

    
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
            #jesHistoUpMax.Scale(1./((1./n_obs)*jesHistoUpMax.Integral() if norming else 1.))#,'width')
            jesHistoUpMax = normalize_combined_TH1_by_blocks(jesHistoUpMax, genBinMap, noNorm=False)
            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=jesHistoUpMax.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'jesHistoMax_gen_{k}',
                                                                    )     
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            #jesHistoDownMax.Scale(1./((1./n_obs)*jesHistoDownMax.Integral() if norming else 1.))#,'width')
            jesHistoDownMax = normalize_combined_TH1_by_blocks(jesHistoDownMax, genBinMap, noNorm=False)
            
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            #jerHistoUpMax.Scale(1./((1./n_obs)*jerHistoUpMax.Integral() if norming else 1.))#,'width')
            jerHistoUpMax = normalize_combined_TH1_by_blocks(jerHistoUpMax, genBinMap, noNorm=False)

            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=jerHistoUpMax.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'jerHistoMax_gen_{k}',
                                                                    )     

            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            #jerHistoDownMax.Scale(1./((1./n_obs)*jerHistoDownMax.Integral() if norming else 1.))#,'width')

            jerHistoDownMax = normalize_combined_TH1_by_blocks(jerHistoDownMax, genBinMap, noNorm=False)

            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            #btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
            #btagHistoUpMax.Reset()
            #btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
            #btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            #btagHistoUpMax.Scale(1./((1./n_obs)*btagHistoUpMax.Integral() if norming else 1.))#,'width')
            btagHistoUpMax = normalize_combined_TH1_by_blocks(btagHistoUpMax, genBinMap, noNorm=False)
            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=btagHistoUpMax.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'btagHistoMax_gen_{k}',
                                                                    )     
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            #btagHistoDownMax.Scale(1./((1./n_obs)*btagHistoDownMax.Integral() if norming else 1.))#,'width')
            btagHistoDownMax = normalize_combined_TH1_by_blocks(btagHistoDownMax, genBinMap, noNorm=False)
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()) and not ('bkg' in k.lower()):#and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #   print(k)
            #   btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k] = normalize_combined_TH1_by_blocks(normeduncerUnfoldHistoshiftsUp[k], genBinMap, noNorm=False)
            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'normeduncerUnfoldHistoshiftsUp_gen_{k}',
                                                                    )     
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)#if not('L1' in text) else 1)
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
            normeduncerUnfoldHistoshiftsDown[k] = normalize_combined_TH1_by_blocks(normeduncerUnfoldHistoshiftsDown[k], genBinMap, noNorm=False)

            physical_uncerUnfoldHisto[k] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'normeduncerUnfoldHistoshiftsDown_gen_{k}',
                                                                    )  

            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)#if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        modelUnc = normalize_combined_TH1_by_blocks(modelUnc, genBinMap, noNorm=False)
        physical_uncerUnfoldHisto[modelkey] = physical_histograms_from_combined_Ndim(
                                                                    combined_hist=modelUnc.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix=f'modelUnc_gen_{k}',
                                                                    )  
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

    physical_unfoldHistoTotUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistoTotUnc.Clone(),
                                                                        bin_map=genBinMap,
                                                                        varDict=varDict,
                                                                        new_hist_prefix='unfoldHistoTotUnc',
                                                                        )  
    
    physical_unfoldHistoDataStatUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistoDataStatUnc.Clone(),
                                                                        bin_map=genBinMap,
                                                                        varDict=varDict,
                                                                        new_hist_prefix='unfoldHistoDataStatUnc',
                                                                        )     
    physical_unfoldHistowoUnc = physical_histograms_from_combined_Ndim(  combined_hist=unfoldHistowoUnc.Clone(),
                                                                     bin_map=genBinMap,
                                                                     varDict=varDict,
                                                                     new_hist_prefix='unfoldHistowoUnc',
                                                                    )     
    physical_unfoldHistoRMStatUnc = physical_histograms_from_combined_Ndim(combined_hist=unfoldHistoRMStatUnc.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix='unfoldHistoRMStatUnc',
                                                                    )     
    physical_unfoldHistoBkgSubUnc = physical_histograms_from_combined_Ndim(combined_hist=unfoldHistoBkgSubUnc.Clone(),
                                                                    bin_map=genBinMap,
                                                                    varDict=varDict,
                                                                    new_hist_prefix='unfoldHistoBkgSubUnc',
                                                                    )     
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    #if not('dijet' in selection): 
    #   totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    #else:
    #   if 'all' in year:
    #       totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
    #   else:
    #       totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #   
    #   if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
    #       totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #   else:
    #       totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
   
    
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.3,0.95 if ('dijet' in selection) else 0.9)
    
    totalErrHist.GetXaxis().SetTitle(f'{labelX}')#'6-body_OC')
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
    #
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
            #print ("JER total",ibin, upmax_ibin,downmax_ibin)
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
    jesHistoUpMax.SetMarkerSize(2.0)
    jesHistoDownMax.SetMarkerSize(2.0)
    up_counter+=1
    down_counter+=1
    col_counter+=1


    for var in varDict.keys():
        outName = outputName.replace('6bodyOC_',var+'_').replace('_dataWithBkgCorr','')
        print(outName)
        physical_uncerUnfoldHistos = OrderedDict()
        for k in physical_uncerUnfoldHisto:
            key = k.replace('6bodyOC_',var+'_')
            physical_uncerUnfoldHistos[key] = physical_uncerUnfoldHisto[k][var].Clone(physical_uncerUnfoldHisto[k][var].GetName()+'_1DSlice')


        drawUncertainties_from_err_shifts_1DfromNDim(
                                                    ivar=var,    
                                                    unfoldHistoTotUnc=physical_unfoldHistoTotUnc[var].Clone(physical_unfoldHistoTotUnc[var].GetName()+'_1DSlice'),
                                                    unfoldHistoDataStatUnc=physical_unfoldHistoDataStatUnc[var].Clone(physical_unfoldHistoDataStatUnc[var].GetName()+'_1DSlice'),
                                                    unfoldHistowoUnc=physical_unfoldHistowoUnc[var].Clone(physical_unfoldHistowoUnc[var].GetName()+'_1DSlice'),
                                                    unfoldHistoRMStatUnc=physical_unfoldHistoRMStatUnc[var].Clone(physical_unfoldHistoRMStatUnc[var].GetName()+'_1DSlice'),
                                                    unfoldHistoBkgSubUnc=physical_unfoldHistoBkgSubUnc[var].Clone(physical_unfoldHistoBkgSubUnc[var].GetName()+'_1DSlice'),

                                                    uncerUnfoldHisto=physical_uncerUnfoldHistos, 
                                                    labelX=varDict[var]['label'], 
                                                    tlegendAlignment=varDict[var]['alignLeg'], 
                                                    outputName=outName, 
                                                    year=year, 
                                                    selection=selection, 
                                                    lumi=lumi 
                                                )

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(2.0)
        jerHistoDownMax.SetMarkerSize(2.0)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:#and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(2.0)
            btagHistoDownMax.SetMarkerSize(2.0)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #   legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
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
                    text = f"Scale and PDF"
                
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'Shower & hadronization', 'l' )    
    #legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'MC stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    


def combine_hist_2D_withUF(
                            aBlankTH2 = None,
                            RM_lists_dict=[],
                            gen_lists_dict=[],
                            truereco_lists_dict=[],
                
                            ):
    combined_response_matrix = aBlankTH2.Clone()
    combined_response_matrix.Reset()
    combined_response_matrix.Sumw2()
    gen_bin_offset = 0
    reco_bin_offset = 0
    
    if len(gen_lists_dict)==0:
        for i in range(len(RM_lists_dict)):
            gen_lists_dict.append(RM_lists_dict[i].ProjectionX(f'{i}_projX')) #just need for bin counting

    if len(truereco_lists_dict)==0:
        for i in range(len(RM_lists_dict)):
            truereco_lists_dict.append(RM_lists_dict[i].ProjectionY(f'{i}_projY')) 
    
    for i in range(len(RM_lists_dict)):
        rm = RM_lists_dict[i]
        gen_bins = gen_lists_dict[i].GetNbinsX()+1
        reco_bins = truereco_lists_dict[i].GetNbinsX()+2

        for g in range(1,gen_bins):
            #print(g)#,gb,rb)
            for r in range(1,reco_bins):
                combined_response_matrix.SetBinContent(gen_bin_offset + g, reco_bin_offset + r, rm.GetBinContent(g, r))
                combined_response_matrix.SetBinError(gen_bin_offset + g, reco_bin_offset + r, rm.GetBinError(g, r))
            #print(r)    
            #include reco underflows for each observable (corrections for reconstruction inefficiency (missgen))
            combined_response_matrix.SetBinContent(gen_bin_offset + g, 0, rm.GetBinContent(g, 0))
            combined_response_matrix.SetBinError(gen_bin_offset + g, 0, rm.GetBinError(g, 0))

            #rb+=1
            #gb+=1
        #print(g)#,gb,rb)
        gen_bin_offset += gen_bins
        reco_bin_offset += reco_bins
    gc.collect()    
    return combined_response_matrix

def combine_all_histogram_types(allVarsDict, varList, 
                                histTypes = [
                                                "reco", "truereco", "fakereco",
                                                "reco_genBin", "truereco_genBin", "fakereco_genBin",
                                                "gen", "accepgen", "missgen",
                                                "respWithMiss"
                                            ],
                                samplePrefLabel="MLMQCD_HT2000toInf",  
                                sel="_dijetSel",
                                sysName = '_nom',
                                extraRecoGap=2,     
                                extraGenGap=1,      
                                verbose=False,
                                combined_var_name = "6bodyOC"
                                                            ):
    """
    Combines all individual observables' histograms (one dict entry per obs) into a set of combined histograms plus record a binMap for each category.

    Inputs
    ----------
    allVarsDict : dict[str -> dict[str->TH1 or TH2]]
        Outer dict: keys are obs names. 
        Value: a dict mapping histogramName->histogramObject for that obs.
        Example: allVars_signalHistos or allVars_dataHistos, etc.
    varList : list of str
        Which obs to combine (keys in allVarsDict).
    histTypes : list of str
        The histogram types we look foruse (provided in allVars dicts): reco, gen, respWithMiss, etc.
    samplePrefLabel : str
        The prefix used for sample+hist naming, e.g. "MLMQCD_HT2000toInf" or "data" or "sysMLMQCD_..." 
        so for each obs hists in dict have keys like "MLMQCD_HT2000toInf_recoJet_tau_0p25_1_nom_dijetSel", etc.
    sel : str ---  "_dijetSel, _WSel, _topSel"
    extraRecoGap : int
        Number of empty bins to insert between sub-vars in the combined 'reco' dimension/axis.
    extraGenGap : int
        Number of empty bins to insert between sub-vars in the combined 'gen' dimension/axis.
    verbose : bool
        Print debug if True.

    Returns
    -------
    A dictionary with:
      {
        "combined_{histType}": TH1D (if found),
        "binMap_{histType}": { ivar: (startBin,endBin) },
        also, e.g.,
        "binMap_respWithMiss": { ivar: (genStart, genEnd, recoStart, recoEnd) }
      }
    """

    outDict = {}
    #Store one combined histogram (if found) per type of histogram in input dict.
    #Keep a bin map for each type, a la binMap_reco[ivar] = (startBin, endBin) for  1D histos
    #or for the 2D (resp. matrix), binMap_respWithMiss[ivar] = (genStart, genEnd, recoStart, recoEnd).

    #Keep track of how many bins each type needs in the reco dimension vs gen dimension for generalisation's sake
    #(currently hard-set to 2x finer reco bins in same global bin range as gen bins which have unit width, ie, reco bin widths uniformly 0.5)
    #For 1D reco-scheme-binned histograms, increment totalRecoBins_{type}.
    #For 1D gen-scheme-binned histograms, increment totalGenBins_{type}.
    #For 2D resp. increment on both axes.
    totalRecoBins = {}
    totalGenBins  = {}

    for t in histTypes:
        totalRecoBins[t] = 0
        totalGenBins[t]  = 0

    #Build dict of which histName for each (ivar, type) to combine them afterwards.
    foundHists = {}
    hasResp = 0 
    for ivar in varList:
        subDict = allVarsDict[ivar]  
        for t in histTypes:
            # helper to build expected hist name strings a la the pattern used for 1D unfoldings

            def build_expected_name(t, ivar, sel):
                """Return the suffix for the histo naming convention."""
                
                if "_genBin" in t:
                    #naming e.g. "MLMQCD_HT2000toInf_recoJet_tau_0p25_1_nom_dijetSel_genBin"
                    mainT = t.replace("_genBin", "")  
                    return f"_{mainT}{ivar}{sysName}{sel}_genBin"
                else:
                    
                    return f"_{t}{ivar}{sysName}{sel}"

            suffix = build_expected_name(t, ivar, sel)
            #print(suffix,samplePrefLabel)
            
            for histKey in subDict.keys():
                if histKey.startswith(samplePrefLabel) and suffix in histKey:
                    #print(histKey)
                    hObj = subDict[histKey]
                    foundHists[(ivar, t)] = hObj
                    if t == "respWithMiss":
                        hasResp+=1
                        #2D
                        nX = hObj.GetNbinsX()
                        nY = hObj.GetNbinsY()
                        
                        totalGenBins[t]  += (nX + extraGenGap)
                        totalRecoBins[t] += (nY + extraRecoGap)
                    else:
                        #1D
                        nBins = hObj.GetNbinsX()
                        if ("gen" in t.lower()) or (t in ["gen", "accepgen", "missgen"]):
                            totalGenBins[t] += (nBins + extraGenGap)
                        else:
                            totalRecoBins[t] += (nBins + extraRecoGap)

                    break  

    if verbose:
        print("[combine_all_histogram_types] Summed bin counts:")
        for t in histTypes:
            if totalRecoBins[t] or totalGenBins[t]:
                print(f"  type={t}, totalRecoBins={totalRecoBins[t]}, totalGenBins={totalGenBins[t]}")

    
    binMaps = {}
    for t in histTypes:
        if not any((ivar, t) in foundHists for ivar in varList):
            continue
        if t == "respWithMiss":
            nx = totalGenBins[t]
            ny = totalRecoBins[t]
            if nx < 1 or ny < 1:
                continue
            h2 = ROOT.TH2D(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}",
                           nx, 0, nx,
                           ny, 0, ny // 2 if ny>2 else ny)  #or just ny
            h2.Sumw2()
            outDict[f"combined_{t}"] = h2
            binMaps[t] = {}
        else:
            
            if ("gen" in t) or (t in ["gen", "accepgen", "missgen"]):
                nb = totalGenBins[t]
                h1 = ROOT.TH1D(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}", nb, 0, nb)
                h1.Sumw2()
                outDict[f"combined_{t}"] = h1
                binMaps[t] = {}
            else:
                nb = totalRecoBins[t]
                h1 = ROOT.TH1D(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}", nb, 0, nb // 2 if nb>2 else nb)
                h1.Sumw2()
                outDict[f"combined_{t}"] = h1
                binMaps[t] = {}

    #Fill all combined histos, for each obs find the hist of a certain type, t, and if it exists, offset bins, copy contents/errors from original histos
    offsets_reco = {t:0 for t in histTypes}
    offsets_gen  = {t:0 for t in histTypes}
    #rm_gen_bin_offset = 0 
    #rm_reco_bin_offset = 0
    
    
    for ivar in varList:
        for t in histTypes:
            if (ivar, t) not in foundHists:
                continue
            hObj = foundHists[(ivar, t)]
            if t == "respWithMiss":
                #h2_comb = outDict[f"combined_{t}"]
                #nX = hObj.GetNbinsX()+1
                #nY = hObj.GetNbinsY()+2
                
                #gxOff = rm_gen_bin_offset#offsets_gen[t]
                #ryOff = rm_reco_bin_offset#offsets_reco[t]
                #binMaps[t][ivar] = (gxOff+1, gxOff+nX, ryOff+1, ryOff+nY)
                
                
                #for gx in range(1, nX):
                #   for ry in range(1, nY):
                #       c = hObj.GetBinContent(gx, ry)
                #       e = hObj.GetBinError(gx, ry)
                #       h2_comb.SetBinContent(gxOff + gx, ryOff + ry, c)
                #       h2_comb.SetBinError(gxOff + gx, ryOff + ry, e)
                #handle misreconstructed gen in reco UF
                #h2_comb.SetBinContent(gxOff + gx, 0, hObj.GetBinContent(g, 0))
                #h2_comb.SetBinError(gxOff + gx, 0, hObj.GetBinError(g, 0))

                #offsets_gen[t]  += (nX + extraGenGap)
                #offsets_reco[t] += (nY + extraRecoGap)
                #rm_gen_bin_offset += nX
                #rm_reco_bin_offset += nY
                continue
            else:
                #1D
                h1_comb = outDict[f"combined_{t}"]
                nBins   = hObj.GetNbinsX()
                if ("gen" in t) or (t in ["gen", "accepgen", "missgen"]):
                    baseOff = offsets_gen[t]
                    binMaps[t][ivar] = (baseOff+1, baseOff + nBins)
                    for iBin in range(1, nBins+1):
                        c = hObj.GetBinContent(iBin)
                        e = hObj.GetBinError(iBin)
                        h1_comb.SetBinContent(baseOff + iBin, c)
                        h1_comb.SetBinError(baseOff + iBin, e)
                    #Insert gap bins between observables
                    for g in range(nBins+1, nBins+1 + extraGenGap):
                        h1_comb.SetBinContent(baseOff + g, 0)
                        h1_comb.SetBinError(baseOff + g, 0)
                    offsets_gen[t] += (nBins + extraGenGap)
                else:
                    baseOff = offsets_reco[t]
                    binMaps[t][ivar] = (baseOff+1, baseOff + nBins)
                    for iBin in range(1, nBins+1):
                        c = hObj.GetBinContent(iBin)
                        e = hObj.GetBinError(iBin)
                        h1_comb.SetBinContent(baseOff + iBin, c)
                        h1_comb.SetBinError(baseOff + iBin, e)
                    for g in range(nBins+1, nBins+1 + extraRecoGap):
                        h1_comb.SetBinContent(baseOff + g, 0)
                        h1_comb.SetBinError(baseOff + g, 0)
                    offsets_reco[t] += (nBins + extraRecoGap)
    #2D
    if hasResp>1:
        RM_lists_dict = []
        gen_lists_dict = []
        truereco_lists_dict = [] 

        for ivar in varList:
            for t in histTypes:
                if (ivar, t) not in foundHists:
                    continue
                if '_gen' in t and not ('genBin' in t):
                    gen_lists_dict.append(foundHists[(ivar, t)].Clone())
                elif '_truereco' in t and not ('genBin' in t):
                    truereco_lists_dict.append(foundHists[(ivar, t)].Clone())
                elif 'respWithMiss' in t:
                    RM_lists_dict.append(foundHists[(ivar, t)].Clone())


        outDict[f"combined_{t}"]  = combine_hist_2D_withUF( outDict[f"combined_{t}"].Clone(outDict[f"combined_{t}"].GetName()+'_blankClone'),
                                                            RM_lists_dict,
                                                            gen_lists_dict,
                                                            truereco_lists_dict,
                                                          
                                                          )

    finalDict = {}
    for t in histTypes:
        if f"combined_{t}" in outDict:
            
            if not('genBin' in t): 
                finalDict[f"combined_{samplePrefLabel}_{t}{combined_var_name}{sysName}{sel}"] = outDict[f"combined_{t}"]
            else:
                finalDict[f"combined_{samplePrefLabel}_{t.split('_genBin')[0]}{combined_var_name}{sysName}{sel}_genBin"] = outDict[f"combined_{t}"]

    for t in binMaps:
        
        if binMaps[t]: 
            if not('genBin' in t):
                finalDict[f"binMap_{samplePrefLabel}_{t}{combined_var_name}{sysName}{sel}"] = binMaps[t]
            else:
                finalDict[f"binMap_{samplePrefLabel}_{t.split('_genBin')[0]}{combined_var_name}{sysName}{sel}_genBin"] = binMaps[t]
    return finalDict

def build_combined_covariance_matrix(hist_list,    
                                     correlation_matrix,
                                     use_off_diag_corr=False,
                                     name_hist='cov_combined'
                                    ):
    """
    Inputs:
      hist_list           : list of TH1s,
      correlation_matrix  : 2D numpy array describing correlation among observables corresponding to input hists in list,
      use_off_diag_corr   : bool; if False, off-diagonal blocks are set to zero instead using correlations from 
                            correlation_matrix[i, j].

    Returns:
      A TH2D (combined_cov_hist) representing the combined covariance matrix
      for all histograms after combining them in one big 1D histo.
      Diagonal (blocks) contains each 1D histogram's bin variances.
      Off-diagonal blocks incorporate correlations between observables or are zeroed on use_off_diag_corr input value (default=False),
      If zero, return just one big diagonal matrix
    """

    #For each histogram, build its diagonal (co)variance array
    #(ie, bin error^2 per bin along diagonal of new combined cov). 
    cov_matrices = []
    n_bins_total = 0

    for reco_hist in hist_list:
        n_bins = reco_hist.GetNbinsX() + 2
        n_bins_total += n_bins
        
        #for blocks on a given diagonal
        cov_matrix = np.zeros((n_bins, n_bins))
        for i in range(1, n_bins):
            error = reco_hist.GetBinError(i)
            cov_matrix[i-1, i-1] = error**2
        cov_matrices.append(cov_matrix)
        
    print(n_bins_total)
    
    #Create square, combined covariance as a TH2 with dimension n_bins_total, 
    #over a range of global bins [0, n_bins_total//2] used for reco axes in other
    #1-/2-D combined hists
    combined_cov_hist = ROOT.TH2D(
        name_hist,
        "Combined Covariance Matrix",
        n_bins_total, 0, n_bins_total/2,
        n_bins_total, 0, n_bins_total/2
    )
    combined_cov_hist.Sumw2()

    #Fill the diagonal blocks from each individual observables' 1-D, bin-wise variances
    bin_offset = 0
    for cov_matrix in cov_matrices:
        n_bins = cov_matrix.shape[0]
        for i in range(n_bins):
            for j in range(n_bins):
                combined_cov_hist.SetBinContent(
                    bin_offset + i,
                    bin_offset + j,
                    cov_matrix[i, j]
                )
        bin_offset += n_bins

    #Fill the off-diagonal blocks, and entries in blocks, using the correlation_matrix for the observables
    #unless 'use_off_diag_corr' is False
    bin_offset_i = 0
    for i in range(len(cov_matrices)):
        cov_matrix_i = cov_matrices[i]
        n_bins_i = cov_matrix_i.shape[0]
        bin_offset_j = 0
        for j in range(len(cov_matrices)):
            cov_matrix_j = cov_matrices[j]
            n_bins_j = cov_matrix_j.shape[0]

            if i != j:
                for k in range(n_bins_i):
                    for l in range(n_bins_j):
                        if use_off_diag_corr:
                            corr = correlation_matrix[i, j]
                        else:
                            corr = 0.0
                        combined_cov_value = corr * np.sqrt(cov_matrix_i[k, k] * cov_matrix_j[l, l])
                        combined_cov_hist.SetBinContent(
                            bin_offset_i + k + 1,
                            bin_offset_j + l + 1,
                            combined_cov_value
                        )
            bin_offset_j += n_bins_j
        bin_offset_i += n_bins_i

    return combined_cov_hist

def rand_int_as_string(nMin=1,nMax=999):
    return str(random.randint(nMin, nMax))

def physical_histograms_from_combined_Ndim( combined_hist, 
                                            bin_map, 
                                            varDict, 
                                            new_hist_prefix="split", 
                                            withSuff=False, 
                                            sys='_nom',
                                            sel='_dijetSel'):
    """
    Input:
    -----------
    combined_hist : ROOT.TH1D
        The combined 1D histogram that contains all the individual histograms merged together.
    bin_map : dict
        A dictionary mapping each observable (ivar) to a tuple (start_bin, end_bin) that indicates
        the bin range (1-indexed) in the combined histogram corresponding to that observable.
    new_hist_prefix : str, optional
        A prefix for naming the split histograms (default is "split").
        
    Return:
    --------
    physical_hists : dict
        A dictionary mapping each observable name to a new TH1D histogram with the corresponding
        bin contents and errors copied from the combined histogram.
    """
    physical_hists = {}
    
    
    #Loop over each observable in the bin map
    for var, (start_bin, end_bin) in bin_map.items():
        #Calculate number of bins in combined histo for this observable
        nBins = end_bin - start_bin + 1
        
        #Create a new TH1D to store the physically-binned distributions
        
        #some redundant code here, FIXME?!
        
        genBins = varDict[var]['bins']
        recoBins = varDict[var]['bins_reco']
        bins=None
        if 'gen' in new_hist_prefix or 'unfold' in new_hist_prefix: 
            bins=genBins
            offset=1
        elif 'data' in new_hist_prefix or 'reco' in new_hist_prefix:
            bins=recoBins
            offset=2
        nBins2 = len(bins)-1
        #print(nBins,nBins2)
        assert(nBins==nBins2)
        if withSuff:
            #print(f"{new_hist_prefix}{var}{sys}{sel}_fromND")
            new_hist = ROOT.TH1D(f"{new_hist_prefix}{var}{sys}{sel}_fromND"+rand_int_as_string(), f"{new_hist_prefix}{var}{sys}{sel}_fromND", len(bins)-1, array( 'd', bins))
        else:
            #print(f"{new_hist_prefix}{var}_fromND")
            new_hist = ROOT.TH1D(f"{new_hist_prefix}{var}_fromND"+rand_int_as_string(), f"{new_hist_prefix}{var}_fromND", len(bins)-1, array( 'd', bins))
        new_hist.Sumw2()  
        bcounter=1
        for i in range(start_bin, end_bin+1):
            comb_bin = i
            content = combined_hist.GetBinContent(comb_bin)
            error = combined_hist.GetBinError(comb_bin)
            new_hist.SetBinContent(bcounter, content)
            new_hist.SetBinError(bcounter, error)
            bcounter+=1
        #print(bcounter)
        physical_hists[var] = new_hist
        
    return physical_hists


def th2_to_numpy_array(th2):
    """
    Convert a ROOT.TH2 (e.g. TH2D) into a NumPy 2D array of bin contents.
    Only the main (in-range) bins are included (i.e. overflow/underflow bins are omitted).
    
    Parameters:
    -----------
    th2 : ROOT.TH2
        The ROOT histogram to convert.
    
    Returns:
    --------
    np.ndarray
        A 2D NumPy array of shape (nbins_x, nbins_y) containing the bin contents.
    """
    nbins_x = th2.GetNbinsX()
    nbins_y = th2.GetNbinsY()
    arr = np.empty((nbins_x, nbins_y), dtype=float)
    #ROOT bins are 1-indexed. Loop over the in-range bins.
    for ix in range(1, nbins_x+1):
        for iy in range(1, nbins_y+1):
            arr[ix-1, iy-1] = th2.GetBinContent(ix, iy)
    return arr

def physical_covariances_from_combined_Ndim(combined_cov, bin_map, varDict,
                                              new_cov_prefix="split_cov", sys='_nom', sel='_WSel'):
    """
    Extract per-observable covariance blocks from a combined covariance matrix.
    
    Parameters:
    -----------
    combined_cov : ROOT.TH2 or numpy.ndarray
        The combined covariance matrix (as a TH2 or already as a NumPy array).
    bin_map : dict
        A dictionary mapping each observable to a tuple (start_bin, end_bin)
        that indicates the bin range (1-indexed) in the combined histogram corresponding to that observable.
    varDict : dict
        A dictionary with binning information for each observable.
        For each observable key, it should have keys like 'bins' (for gen/unfolded) or 'bins_reco' (for data/reco).
    new_cov_prefix : str, optional
        A prefix for naming the split covariance matrices (default is "split_cov").
    sys : str, optional
        Suffix to append to the names (default is '_nom').
    sel : str, optional
        Additional suffix (default is '_dijetSel').
        
    Returns:
    --------
    physical_covs : dict
        A dictionary mapping each observable to its corresponding covariance block as a new ROOT.TH2D.
    """
    physical_covs = {}
    
    
    cov_full = th2_to_numpy_array(combined_cov)
    

    #Loop over each observable using the bin map.
    #N.B.: the bin_map uses 1-indexed bin numbers, so subtract 1 for NumPy slicing.
    for var, (start_bin, end_bin) in bin_map.items():
        nBins = end_bin - start_bin + 1  #number of physical bins for this observable
        
        #Extract the corresponding block from the combined covariance array.
        block = cov_full[start_bin - 1:end_bin, start_bin - 1:end_bin].copy()
        
        #Pick appropriate physical binning from dict of variables.
        if ('gen' in new_cov_prefix) or ('unfold' in new_cov_prefix):
            bins = varDict[var]['bins']
        elif ('data' in new_cov_prefix) or ('reco' in new_cov_prefix):
            bins = varDict[var]['bins_reco']
        else:
            bins = varDict[var]['bins']
            
        nBins2 = len(bins) - 1
        assert(nBins == nBins2), f"Mismatch in bin count for {var}: {nBins} vs {nBins2}"
        
        #Create a new ROOT.TH2D to store this covariance block.
        name = f"{new_cov_prefix}{var}{sys}{sel}"
        new_cov = ROOT.TH2D(name, name, nBins, array('d', bins), nBins, array('d', bins))
        new_cov.Sumw2()
        
        #Fill the new TH2D with the covariance block contents.
        for i in range(nBins):
            for j in range(nBins):
                new_cov.SetBinContent(i + 1, j + 1, block[i, j])
                
        physical_covs[var] = new_cov.Clone()
    
    return physical_covs


def numpy_to_hist2D_Ndim(numpy_array, input_hist):
    """
    Convert a 2D numpy array to a ROOT.TH2D histogram, preserving variable bin widths.
    
    If the input_hist is a TH2 with more than one bin in Y,
    its x- and y-axes are used for the binning.
    Otherwise, if the input_hist is a TH1 (or a TH2 with only one y-bin),
    its x-axis is used as the global binning for both dimensions.
    
    Parameters:
      numpy_array  : 2D numpy.ndarray 
                     The array with shape (nBins, nBins) holding the desired bin contents.
      input_hist: ROOT.TH1 or ROOT.TH2D
                     The histogram from which to extract the binning scheme.
      
    Returns:
      hist_out : ROOT.TH2D
                 A new TH2D with variable bin widths and bin contents set from numpy_array.
    """
    if input_hist.GetNbinsY() > 1:
        nBinsX = input_hist.GetNbinsX()
        nBinsY = input_hist.GetNbinsY()
        xaxis = input_hist.GetXaxis()
        yaxis = input_hist.GetYaxis()
        xedges = [xaxis.GetBinLowEdge(i) for i in range(1, nBinsX+2)]
        yedges = [yaxis.GetBinLowEdge(i) for i in range(1, nBinsY+2)]
    else:
        #If input_hist is a TH1 (or a TH2 with 1 bin along Y), assume its x-axis defines
        #the global binning for both axes.
        nBins = input_hist.GetNbinsX()
        xaxis = input_hist.GetXaxis()
        xedges = [xaxis.GetBinLowEdge(i) for i in range(1, nBins+2)]
        #Use the same edges for Y.
        nBinsX = nBins
        nBinsY = nBins
        yedges = xedges

    if numpy_array.shape != (nBinsX, nBinsY):
        raise ValueError("Shape of numpy_array {} does not match expected shape ({}, {}) from histogram."
                         .format(numpy_array.shape, nBinsX, nBinsY))
    
    hist_out = ROOT.TH2D(input_hist.GetName(), input_hist.GetTitle(),
                         nBinsX, array('d', xedges),
                         nBinsY, array('d', yedges))
    hist_out.Sumw2()
    
    for ix in range(1, nBinsX+1):
        for iy in range(1, nBinsY+1):
            hist_out.SetBinContent(ix, iy, numpy_array[ix-1, iy-1])
    
    return hist_out


def compute_block_jacobian(hist, start_bin, end_bin, withOF=False):
    """
    Compute the Jacobian for a block corresponding to one observable.
    
    
    For bins i in this block, assume
         v_i = c_i / N_block, where N_block = sum of bin contents in (start_bin, end_bin).
         
    
    """
    #Compute summed yield for one block in combined unfolding.
    N_block = 0.
    for i in range(start_bin, end_bin+1):
        c = hist.GetBinContent(i)
        if c < 0:
            c = 0.
        N_block += c

    nBins = end_bin - start_bin + 1
    jac_block = np.zeros((nBins, nBins))
    for i in range(nBins):
        c_i = hist.GetBinContent(start_bin + i)
        if c_i < 0:
            c_i = 0.
        for j in range(nBins):
            if i == j:
                jac_block[i, j] = (N_block - c_i) / (N_block**2)
            else:
                jac_block[i, j] = - c_i / (N_block**2)
    return jac_block

def get_normalised_cov_combined(unfHisto, ematrix, bin_map):
    """
    Compute the normalized covariance matrix for a combined 1D unfolded histogram.
    
    1) compute a block-diagonal Jacobian by computing a separate Jacobian for each observable
       (using its own yield, summed over (start_bin,end_bin) in bin_map), 
    2) transform the relevant covariance matrix with the block-diagonal global Jacobian, for the combined unfolding,
       effectively applying a simple normalising transformation to the corresponding block for a given obs. 
       and off-diagonals contain correlated entries between bins of various observables in the normalised
       combined covariance matrix.
    
    """
    cov_abs = th2_to_numpy_array(ematrix)
    
    nbins_total = unfHisto.GetNbinsX()
    cov_norm = np.zeros((nbins_total, nbins_total))
    globalJac_norm = np.zeros((nbins_total, nbins_total))

    
    for var, (start_bin, end_bin) in bin_map.items():
        nBins = end_bin - start_bin + 1
        jac = compute_block_jacobian(unfHisto, start_bin, end_bin)
        
        #Extract block from the combined covariance.
        #block = cov_abs[start_bin - 1 : end_bin, start_bin - 1 : end_bin].copy()
        
        #block_norm = jac @ block @ jac.T
        
        globalJac_norm[start_bin - 1 : end_bin, start_bin - 1 : end_bin] = jac#block_norm

    cov_norm = globalJac_norm @ cov_abs @ globalJac_norm.T
    
    ematrix_norm_th2 = numpy_to_hist2D(cov_norm, ematrix.Clone(ematrix.GetName()+'_jacTrafo'))
    
    return cov_norm, ematrix_norm_th2

def normalize_combined_TH1_by_blocks(hist, bin_map, noNorm=False):
    """
    Normalizes a combined TH1 (stitched from all individual observables)
    on an obs.-by-obs./block-by-block basis. For each observable (defined by (start_bin, end_bin)
    in bin_map): 
    1) computes the block integral, 
    2) divides each bin’s content and error in that block by the block integral.
    3) error computation is anyway overwritten when setting the bin errors correctly
       from the normalised covariances
    Parameters:
      hist    : ROOT.TH1D
                The combined TH1 to be normalized.
      bin_map : dict
                A dictionary mapping each observable to a tuple (start_bin, end_bin)
                (1-indexed) that defines the bin range in the combined histogram.
      noNorm  : bool, optional
                If True, no normalization is applied (the histogram is returned as is).

    Returns:
      hist    : ROOT.TH1D
                The same histogram with each block normalized by its own integral.
    """
    if noNorm:
        return hist

    #Loop over each observable block in bin_map.
    for var, (start_bin, end_bin) in bin_map.items():
        #Compute the integral (yield) for this unfolded observable/block.
        block_integral = 0.
        for i in range(start_bin, end_bin + 1):
            block_integral += hist.GetBinContent(i)
        
        if block_integral > 0:
            for i in range(start_bin, end_bin + 1):
                content = hist.GetBinContent(i)
                error   = hist.GetBinError(i)
                new_content = content / block_integral
                new_error   = error / block_integral
                hist.SetBinContent(i, new_content)
                hist.SetBinError(i, new_error)
    return hist

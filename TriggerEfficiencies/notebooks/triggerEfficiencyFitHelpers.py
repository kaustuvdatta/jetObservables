from ROOT import *
import os, math, sys, argparse, gc, copy, json, pprint
from array import array
from random import gauss
import numpy as np
from collections import OrderedDict

import histoHelpers
from histoHelpers import *

sys.path.insert(0,'../../')
from datasets_dijetSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict

sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle


tdrstyle.setTDRStyle()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptStat(0)


def fitTriggerEfficiencies_perYear(inFileSample, sample, year, triggerDict, triggerPassList, triggerDenomList, nameList, lumi, outputDir,
                                   minFit=60., maxFit=1200., ext='.png', version='_Aug23',
                                   xlabel='Leading jet pt [GeV]', unzoomed = False,
                                   xmin=100, xmax=1500, rebin=1, fitRefine=True, verbose=False):
    
    ############# Set up plot style #############
    CMS_lumi.lumi_13TeV = f"{year} (13 TeV)"
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText = "Preliminary"

    CMS_lumi.relPosX = 0.1

    iPos = 0

    H_ref = 720; 
    W_ref = 1024; 
    W = W_ref
    H  = H_ref

    iPeriod = 4

    # references for T, B, L, R
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref

    x1_l = 0.9
    y1_l = 0.9

    dx_l = 0.25
    dy_l = 0.25
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l
    #############################################
    
    
    for triggerPass,triggerDenom,name in zip(triggerPassList,triggerDenomList,nameList):

        outputFileName = outputDir+name+'_'+triggerDenom+"_"+triggerPass+'_'+sample+'_TriggerEfficiency'+version+ext
        if verbose: print (f'Processing.......{outputFileName} for {triggerPass,triggerDenom,name}')


        threshold = int(triggerPass.split('AK8PFJet')[1])

        can = TCanvas('c2'+triggerPass, 'c2'+triggerPass,W,H)
        can.SetTicks(1)
        can.SetFillColor(0)
        can.SetBorderMode(0)
        can.SetBorderSize(0)
        can.SetFrameFillStyle(0)
        can.SetFrameBorderMode(0)
        #can.SetFrameLineWidth(0)
        can.SetLeftMargin( L/W )
        can.SetRightMargin( R/W )
        can.SetTopMargin( T/H )
        can.SetBottomMargin( B/H )
        can.SetTickx(0)
        can.SetTicky(0)
        #can.cd()

        can.cd()


        histos={}
        histos[ 'denomOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_'+triggerDenom+'_only' ).Clone()
        #histos[ 'denomOnly'+triggerPass ].Scale(triggerDict[triggerDenom]['prescale'])
        histos[ 'PassingOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_'+triggerPass+'_simulated' ).Clone()
        #histos[ 'PassingOnly'+triggerPass ].Scale(triggerDict[triggerPass]['prescale'])


        histos[ 'denomOnly'+triggerPass ].Rebin(rebin)
        histos[ 'PassingOnly'+triggerPass ].Rebin(rebin)
        binWidth = histos[ 'denomOnly'+triggerPass ].GetBinWidth(1)

        histos[ 'denomOnly'+triggerPass ].SetLineWidth(2)
        histos[ 'denomOnly'+triggerPass ].SetLineColor(kGray)

        histos[ 'eff'+triggerPass ] = histos[ 'PassingOnly'+triggerPass ].Clone("Efficiency")#ROOT.TEfficiency( histos[ 'PassingOnly'+triggerPass ].Clone(), histos[ 'denomOnly'+triggerPass ].Clone() )
        histos[ 'eff'+triggerPass ].Divide(histos[ 'eff'+triggerPass ],histos[ 'denomOnly'+triggerPass ],1,1,"B")


        histos[ 'eff'+triggerPass ].SetMarkerStyle(8)
        if not(unzoomed):
            histos[ 'eff'+triggerPass ].SetMaximum(1.005)
            histos[ 'eff'+triggerPass ].SetMinimum(0.99)
        else:
            histos[ 'eff'+triggerPass ].SetMaximum(1.5)
            histos[ 'eff'+triggerPass ].SetMinimum(0.)
        
        histos[ 'eff'+triggerPass ].GetXaxis().SetRangeUser(0, min(6*threshold, 2000))
        #histos[ 'eff'+triggerPass ].SetTitleOffset(1.3,'X')

        title = f'{triggerPass} relative to {triggerDenom}'    
        histos[ 'eff'+triggerPass ].SetTitle(title+"; Leading jet p_{T} [GeV]; Efficiency (#epsilon)")
        histos[ 'eff'+triggerPass ].GetYaxis().CenterTitle()
        histos[ 'eff'+triggerPass ].GetYaxis().SetTitleSize(0.06)
        histos[ 'eff'+triggerPass ].GetXaxis().SetTitleSize(0.05)
        histos[ 'eff'+triggerPass ].GetYaxis().SetTitleOffset(0.9)
        histos[ 'eff'+triggerPass ].GetXaxis().SetTitleOffset(1.15)
        histos[ 'eff'+triggerPass ].Draw()

        #histos[ 'eff'+triggerPass ].Draw()
        # Fit to obtain turn-ons

        #lower_threshold = int(triggerPass.split('AK8PFJet')[1])*1.2
        lower_threshold = int(triggerPass.split('AK8PFJet')[1])/1.1
        higher_threshold = int(triggerPass.split('AK8PFJet')[1])*3.5

        if verbose: print(triggerPass,triggerPass.split('AK8PFJet')[1],lower_threshold,higher_threshold)

        fit = ROOT.TF1("eff_%s" % name, '([0] + 0.5 * (1-[0]) * (1 + erf((x-[1])/[2])))', lower_threshold, higher_threshold)
        # set parameter names
        fit.SetParName(0, 'a')
        fit.SetParName(1, 'mu')
        fit.SetParName(2, 'sigma')
        #fit.SetParName(3, 'N')
        # set parameter limits
        fit.SetParLimits(1, lower_threshold, higher_threshold)  # enforce +ve parameter values
        fit.SetParLimits(1, 1, 1000)  # enforce +ve parameter values
        fit.SetParLimits(2, 20, 500)
        # fit.SetParLimits(3, 0.00001, 100)
        #fit.SetParLimits(3, 0.9, 1.1)
        fit.SetLineColor(ROOT.kBlack)
        fit.SetLineWidth(1)
        # Set starting values
        fit.SetParameter('a', 0)
        fit.SetParameter('mu', int(triggerPass.split('AK8PFJet')[1]))
        fit.SetParameter('sigma', int(triggerPass.split('AK8PFJet')[1])/10)
        #fit.SetParameter('N', 1)
        fit.SetNpx(1000)
        fit.SetLineColor(ROOT.kRed)
        fit_result = histos[ 'eff'+triggerPass ].Fit(fit, 'RSEQ')
        histos[ 'eff'+triggerPass ].Draw()

        #cl1=histos[ 'PassingOnly'+triggerPass ].Clone()
        #cl2=histos[ 'denomOnly'+triggerPass ].Clone()
        #cl1.Scale(20./cl1.Integral())
        #cl2.Scale(20./cl2.Integral())
        #cl1.Draw("hist same")
        #cl2.Draw("hist same")

        ROOT.gPad.Modified()
        ROOT.gPad.Update()

        
        good_eff = 0.995# * fit.GetParameter("N")
        good_eff_pt = fit.GetX(good_eff)
        triggerDict[triggerPass]['turn-on_1'] = good_eff_pt
        if verbose: print("turn-on, unrefined fit",good_eff_pt,triggerPass)
        
        if not(fitRefine):
            
            stats_box = histos[ 'eff'+triggerPass ].FindObject("stats")
            stats_box.SetFillColor(ROOT.kWhite)
            stats_box.SetBorderSize(0)
            stats_box.SetFillStyle(0)
            stats_box.SetX1NDC(0.62)
            stats_box.SetX2NDC(0.88)
            stats_box.SetY1NDC(0.75)
            stats_box.SetY2NDC(0.88)
            #line.Draw("same")
            CMS_lumi.CMS_lumi(can, iPeriod, iPos)
            can.cd()
            can.Update()
            can.RedrawAxis()
            frame = can.GetFrame()
            frame.Draw()

            eff_text = ROOT.TPaveText(0.63, 0.65, 0.88, 0.73, "NDC")
            eff_text.AddText("#epsilon = 0.995 @ p_{T} = %3.f GeV" % good_eff_pt)
            eff_text.SetFillStyle(0)
            eff_text.SetBorderSize(0)
            eff_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            eff_text.Draw()
            
            extra = ROOT.TPaveText(0.32,0.8,0.77,1.1, "brNDC")
            extra.SetTextSize(0.037)
            extra.AddText(title)
            extra.SetFillColor(0)
            extra.SetBorderSize(0)
            extra.Draw()

            can.Update()
            can.Draw()
            #input("Press Enter to end")
            can.SaveAs(outputFileName)
            can.Close()
            
            
        else:
            # Update fit by increasing starting point to better capture effects in the high efficiency region
            fit_factor = 0.95 
            fit.SetRange(fit.GetX(fit_factor), higher_threshold*1.)
            fit_result = histos[ 'eff'+triggerPass ].Fit(fit, 'RSEMQ')
            ROOT.gPad.Modified()
            ROOT.gPad.Update()


            #fit_result.Draw()#histos[ 'eff'+triggerPass ].Draw()
            # Draw fit stats
            stats_box = histos[ 'eff'+triggerPass ].FindObject("stats")
            stats_box.SetFillColor(ROOT.kWhite)
            stats_box.SetBorderSize(0)
            stats_box.SetFillStyle(0)
            stats_box.SetX1NDC(0.62)
            stats_box.SetX2NDC(0.88)
            stats_box.SetY1NDC(0.75)
            stats_box.SetY2NDC(0.88)
            #line.Draw("same")
            CMS_lumi.CMS_lumi(can, iPeriod, iPos)
            can.cd()
            can.Update()
            can.RedrawAxis()
            frame = can.GetFrame()
            frame.Draw()


            # Add in info about 99% relative efficiency
            good_eff = 0.995# * fit.GetParameter("N")
            good_eff_pt = fit.GetX(good_eff)
            #info['good_eff_pt'] = good_eff_pt
            #info['good_eff_pt_err'] = 0
            if verbose: print("turn-on, refined fit", good_eff_pt,triggerPass)
            eff_text = ROOT.TPaveText(0.63, 0.65, 0.88, 0.73, "NDC")
            eff_text.AddText("#epsilon = 0.995 @ p_{T} = %3.f GeV" % good_eff_pt)
            eff_text.SetFillStyle(0)
            eff_text.SetBorderSize(0)
            eff_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            eff_text.Draw()


            extra = ROOT.TPaveText(0.32,0.8,0.77,1.1, "brNDC")
            extra.SetTextSize(0.037)
            extra.AddText(title)
            extra.SetFillColor(0)
            extra.SetBorderSize(0)
            extra.Draw()
            triggerDict[triggerPass]['turn-on_2'] = good_eff_pt

            can.Update()
            can.Draw()
            #input("Press Enter to end")
            can.SaveAs(outputFileName)
            can.Close()
    return triggerDict
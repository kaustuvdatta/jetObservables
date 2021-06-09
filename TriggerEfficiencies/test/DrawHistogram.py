#!/usr/bin/env python
'''
File: DrawHistogram.py
Author: Alejandro Gomez Espinosa
Description: My Draw histograms. Check for options at the end.
'''

from ROOT import *
import time, os, math, sys
from array import array
import argparse
sys.path.insert(0,'../../Unfolding/python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
from collections import OrderedDict

#gROOT.Reset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()

gStyle.SetOptStat(0)


xline = array('d', [0,1000])
yline = array('d', [0.99,.99])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

def plotTriggerEfficiency( inFileSample, sample, triggerDenom, denomScale, triggerPass, passScale, name, minFit, maxFit, xlabel, xmin, xmax, rebin, log):
    """docstring for plot"""

    outputFileName = name+'_'+triggerDenom+"_"+triggerPass+'_'+sample+'_TriggerEfficiency'+args.version+'.'+args.ext
    print 'Processing.......', outputFileName

    histos = {}

    if args.type.startswith('simulated'):
        histos[ 'denomOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_'+triggerDenom+'_only' ).Clone()
        histos[ 'PassingOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_'+triggerPass+'_simulated' ).Clone()
    else:
        histos[ 'denomOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_HLT_'+triggerDenom+'_scaled' ).Clone()
        #histos[ 'denomOnly'+triggerPass ].Scale( denomScale )
        histos[ 'PassingOnly'+triggerPass ] = inFileSample.Get( 'TriggerEfficiencies/jet1Pt_HLT_'+triggerPass+'_scaled' ).Clone()
        histos[ 'PassingOnly'+triggerPass ].Scale( passScale )

    histos[ 'denomOnly'+triggerPass ].Rebin(rebin)
    histos[ 'PassingOnly'+triggerPass ].Rebin(rebin)
    if args.type.startswith('simulated'):
        histos[ 'eff'+triggerPass ] = TEfficiency( histos[ 'PassingOnly'+triggerPass ].Clone(), histos[ 'denomOnly'+triggerPass ].Clone() )
    else:
        histos[ 'eff'+triggerPass ] = TGraphAsymmErrors( histos[ 'PassingOnly'+triggerPass ].Clone(), histos[ 'denomOnly'+triggerPass ].Clone(), 'pois'  )

    binWidth = histos[ 'denomOnly'+triggerPass ].GetBinWidth(1)

    legend=TLegend(0.65,0.55,0.90,0.70)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    histos[ 'denomOnly'+triggerPass ].SetLineWidth(2)
    histos[ 'denomOnly'+triggerPass ].SetLineColor(kGray)
    #histos[ 'PassingOnly'+triggerPass ].SetLineWidth(2)
    #histos[ 'PassingOnly'+triggerPass ].SetLineColor(kBlue-4)

    tdrStyle.SetPadTickY(0)
    can = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    pad1 = TPad("pad1", "Histo",0,0.00,1.00,1.00,-1)
    pad2 = TPad("pad2", "",0,0.00,1.00,1.00,-1)
    pad2.SetFillStyle(4000)
    pad1.Draw()
    pad1.cd()
    if log: pad1.SetLogy()

    legend.AddEntry( histos[ 'denomOnly'+triggerPass ], triggerDenom, 'l' )
    legend.AddEntry( histos[ 'PassingOnly'+triggerPass ], triggerPass, 'l' )
    #histos[ 'denomOnly'+triggerPass ].SetMinimum(10)
    histos[ 'denomOnly'+triggerPass ].GetXaxis().SetRangeUser( xmin, xmax )
    histos[ 'denomOnly'+triggerPass ].SetMaximum( histos[ 'denomOnly'+triggerPass ].GetMaximum()*1.5 )
    histos[ 'denomOnly'+triggerPass ].Draw('histe')
    histos[ 'denomOnly'+triggerPass ].GetYaxis().SetTitleSize(0.06)
    histos[ 'denomOnly'+triggerPass ].GetYaxis().SetTitleOffset(0.8)
    histos[ 'denomOnly'+triggerPass ].GetYaxis().SetLabelSize(0.06)
    histos[ 'denomOnly'+triggerPass ].GetXaxis().SetTitleOffset(0.8)
    histos[ 'denomOnly'+triggerPass ].GetXaxis().SetTitleSize(0.06)
    histos[ 'denomOnly'+triggerPass ].GetXaxis().SetLabelSize(0.05)
    histos[ 'PassingOnly'+triggerPass ].Draw('histe same')
    histos[ 'denomOnly'+triggerPass ].GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    histos[ 'denomOnly'+triggerPass ].GetXaxis().SetTitle( 'Leading Jet p_{T} [GeV]' )

    CMS_lumi.CMS_lumi(pad1, 4, 0)
    gPad.Update()
    xMin = pad1.GetUxmin()
    xMax = pad1.GetUxmax()
    dx = (xMax - xMin) / 0.8
    ymin = 0 #histos[ 'eff'+triggerPass ].CreateGraph().GetHistogram().GetMinimum()
    ymax = 1 #histos[ 'eff'+triggerPass ].CreateGraph().GetHistogram().GetMaximum()
    dy = (ymax - ymin) / 0.8
    pad2.Range(xMin-0.1*dx, ymin-0.1*dy, xMax+0.1*dx, ymax+0.1*dy)
    pad2.Draw()
    pad2.cd()
    histos[ 'eff'+triggerPass ].SetMarkerStyle(8)
    #histos[ 'eff'+triggerPass ].SetLineWidth(2)
    #histos[ 'eff'+triggerPass ].SetLineColor(kBlue-4)
    histos[ 'eff'+triggerPass ].SetFillStyle(1001)
    histos[ 'eff'+triggerPass ].Draw()
    legend.AddEntry( histos[ 'eff'+triggerPass ], 'Efficiency', 'p' )
    legend.Draw()
    pad2.Update()
#    histos[ 'eff'+triggerPass ].GetPaintedGraph().GetYaxis().SetLabelSize(0.06)
#    histos[ 'eff'+triggerPass ].GetPaintedGraph().GetXaxis().SetLabelSize(0.06)
#    histos[ 'eff'+triggerPass ].GetPaintedGraph().GetYaxis().SetTitleSize(0.06)
#    histos[ 'eff'+triggerPass ].GetPaintedGraph().GetYaxis().SetTitleOffset(0.8)
    if args.type.startswith('simulated'):
        histos[ 'eff'+triggerPass ].GetPaintedGraph().GetYaxis().SetLabelOffset(9999)
        histos[ 'eff'+triggerPass ].GetPaintedGraph().GetXaxis().SetLimits( xmin, xmax )
        histos[ 'eff'+triggerPass ].GetPaintedGraph().GetYaxis().SetTickLength( 0 )
    else:
        histos[ 'eff'+triggerPass ].GetYaxis().SetLabelOffset(9999)
        histos[ 'eff'+triggerPass ].GetXaxis().SetLimits( xmin, xmax )
        histos[ 'eff'+triggerPass ].GetYaxis().SetTickLength( 0 )
        histos[ 'eff'+triggerPass ].SetMinimum(0)
        histos[ 'eff'+triggerPass ].SetMaximum(1.2)
    gPad.Update()

    newAxis = TGaxis( xmax,ymin,xmax, ymax,ymin,ymax,502,"+L") #,0.03 )
    #newAxis.SetLabelOffset(0.1)
    #newAxis.SetLineColor(kRed)
    #newAxis.SetLabelColor(kRed)
    newAxis.SetTitleSize(0.06)
    newAxis.SetTitleOffset(0.6)
    newAxis.SetLabelSize(0.05)
    #newAxis.SetLabelOffset(0.06)
    #newAxis.SetTickLength(0.01)
    #newAxis.SetLabelSize(0.03)
    newAxis.SetTitle("Efficiency")
    newAxis.Draw()
    line.Draw('same')

##    errF = TF1('errF', '(1+TMath::Erf((x-[0])/[1]))', minFit, maxFit )
#    errF = TF1('errF', '[3]*( [2] + (0.5 * ( 1-[2]) * ( 1 + TMath::Erf( (x-[0])/[1] ) )))', minFit, maxFit )  ## Mass
#    errF.SetParLimits(0, 0, 200 )
#    errF.SetParLimits(1, 0, 100 )
#    errF.SetParLimits(2, -10, 10 )
#    errF.SetParLimits(3, .9, 1.1 )
##    histos[ 'eff'+triggerPass ].SetStatisticOption(TEfficiency.kFWilson)
#    for i in range(5): histos[ 'eff'+triggerPass ].Fit(errF, 'MIR')
#    errF.Draw("same")
##    #for i in range(5): Efficiency.Fit('errF', 'MIR')

    can.SaveAs( 'Plots/'+outputFileName )
    del can

#    #### Fitting
#    #errF = TF1('errF', '0.5*(1+TMath::Erf((x-[0])/[1]))', 500, 1500 )
#    #errF = TF1('errF', '0.5*(1+TMath::Erf(([0]*x-[1])/[2]))', 400, 1000 )  ## HT
#    #errF = TF1('errF', '0.5*(1+TMath::Erf(([0]*x-[1])/[2]))', 0, 100 )  ## Mass
#    #Efficiency.SetStatisticOption(TEfficiency.kFWilson)
#    #for i in range(5): eff.Fit(errF, '+')
#    #for i in range(5): Efficiency.Fit('errF', 'MIR')
#    #print '&'*10, '900', errF.Eval(900)
#    #print '&'*10, '1000', errF.Eval(1000)
#    gStyle.SetOptFit(1)
#    '''
#    errF.SetLineColor(kRed)
#    errF.SetLineWidth(2)
#    errF.Draw('sames')
#    can1.Update()
#    st1 = Efficiency.GetListOfFunctions().FindObject("stats")
#    st1.SetX1NDC(.60);
#    st1.SetX2NDC(.90);
#    st1.SetY1NDC(.20);
#    st1.SetY2NDC(.50);
##	#eff.Draw("same")
#    can1.Modified()
#    '''
#
    return histos[ 'eff'+triggerPass ]


###################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='simple', help='Process to draw, example.' )
    parser.add_argument('-t', '--type', action='store', default='simulated', help='Using prescales or simulated' )
    parser.add_argument('-d', '--dataset', action='store', default='all', help='Dataset: JetHT, SingleMuon, etc.' )
    parser.add_argument('-v', '--version', action='store', default='v01', help='Version of the files' )
    parser.add_argument('-C', '--cut', action='store', default='_cutDEta', help='cut, example: cutDEta' )
    parser.add_argument('-s', '--single', action='store', default='all', help='single histogram, example: massAve_cutDijet.' )
    parser.add_argument('-l', '--lumi', action='store', default='15.5', help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('--year', action='store', default='2017', help='Extension of plots.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if not os.path.exists('Plots/'): os.makedirs('Plots/')

    plotList = [
        [ 'simple', 'AK8Jet1Pt', 'Leading jet pt [GeV]', 100, 1000, 1, True],
#        [ 'simple', 'AK8Jet2Pt', '2nd Leading jet pt [GeV]', 100, 800, 1, True],
            ]

    if 'all' in args.single: Plots = [ x[1:] for x in plotList if x[0] in args.proc ]
    else: Plots = [ y[1:] for y in plotList if ( ( y[0] in args.proc ) and ( y[1].startswith(args.single ) ))  ]


    bkgFiles = {}
    signalFiles = {}
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = '13 TeV, '+args.year

    triggerList = [
#                    ['AK8PFJet80', 16419.911,   'AK8PFJet60', 82183.849,  150, 250 ],
                    ['AK8PFJet140', 0.0085,   'AK8PFJet80', 16419.911, 150, 250 ],
                    ['AK8PFJet200', 0.067,    'AK8PFJet140', 1560.153, 100, 800 ],
                    ['AK8PFJet260', 0.42,     'AK8PFJet200', 219.705, 100, 800 ],
                    ['AK8PFJet320', 1.,     'AK8PFJet260', 88.426, 100, 800 ],
                    ['AK8PFJet400', 1.,      'AK8PFJet320', 33.827, 100, 800  ],
                    ['AK8PFJet450', 1.,      'AK8PFJet400', 5.404, 100, 800  ],
                    ['AK8PFJet500', 1.,          'AK8PFJet450', 4.299, 100, 800  ],
                    ['AK8PFJet550', 1.,          'AK8PFJet500', 1, 100, 800  ]
                    ]

    processingSamples = {}
    if args.year.startswith('2017'): processingSamples[ 'JetHT2017' ] = [ TFile.Open('/eos/home-a/algomez/tmpFiles/jetObservables/triggerEfficiencies/Plots/v05/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017ALL.root'), 0 ]
    else: processingSamples[ 'JetHT2018' ] = [ TFile.Open('/eos/home-a/algomez/tmpFiles/jetObservables/triggerEfficiencies/Plots/v05/triggerEfficiencies_histograms_MiniAOD_JetHTRun2018ALL.root'), 0 ]

    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    for i in Plots:
        for isam, samFile in processingSamples.iteritems():
            for q, it in enumerate(triggerList):
                plotTriggerEfficiency( samFile[0], isam, it[2], it[3], it[0], it[1],  i[0], it[4], it[5], i[1], i[2], i[3], i[4], i[5] )

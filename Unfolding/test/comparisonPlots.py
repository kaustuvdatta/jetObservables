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

##################################################
def plotSignalBkg( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False,
                      addRatioFit=False, ext='png', outputDir='Plots/', legendAlignment='right' ):
    """function to plot s and b histos"""

    outputFileName = name+'_AnalysisPlots_'+args.version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    if legendAlignment.startswith('right'): legend=ROOT.TLegend(0.60,(0.75 if args.selection.startswith('dijet') else 0.65),0.90,0.90)
    else: legend=ROOT.TLegend(0.20,(0.75 if args.selection.startswith('dijet') else 0.65),0.50,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    dictHistos = OrderedDict()
    for ivar in nSubVariables.keys():
        if name.split('_')[0].endswith( ivar.split('_')[0] ):
            dictHistos[ivar] = {
                    'data' : dictFiles[ivar].Get( 'dataMinusBkgsGenBin' ),
                    'unfoldHisto' : dictFiles[ivar].Get( 'unfoldHisto'+ivar )
                    }

    #legend.AddEntry( dataHistos[ 'DATA' ], 'Data', 'lep' )

    pad = {}
    canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 2000, 750 )
    canvas[outputFileName].Divide( len(dictHistos), 1 )
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    dummy = 0
    for ivar, ihistos in dictHistos.items():
        dummy=dummy+1
        canvas[outputFileName].cd(dummy)
        pad[ivar+'1'] = ROOT.TPad(ivar+'1', "Fit",0,0.207,1.00,1.00,-1)
        pad[ivar+'2'] = ROOT.TPad(ivar+'2', "Pull",0,0.00,1.00,0.30,-1);
        pad[ivar+'1'].Draw()
        pad[ivar+'2'].Draw()

        pad[ivar+'1'].cd()
        if log: pad[ivar+'1'].SetLogy()
#        if 'tau' in axisX:
#            dataHistos['DATA'].GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+axisX.split('#')[1]+'}' )
#        else:
#            dataHistos['DATA'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )
#        dataHistos[ 'DATA' ].SetMarkerStyle(8)
#        dataHistos['DATA'].GetYaxis().SetTitleSize( 0.05 )
#        dataHistos['DATA'].GetYaxis().SetTitleOffset( 1.2 )
#        dataHistos['DATA'].GetXaxis().SetTitle( axisX )
#        if xmax: dataHistos['DATA'].GetXaxis().SetRangeUser( xmin, xmax )
#        dataHistos['DATA'].SetMaximum( hBkg.GetMaximum()*1.2 )
#        #dataHistos['DATA'].SetMinimum( 2. )
#
#        #stackHisto.SetMinimum( 0.1 )
#        if not args.selection.startswith('dijet'):
#            hBkg.SetLineColor( bkgFiles['TTToSemiLeptonic'][1]['color'] )
#            hBkg.SetFillColor( 0 )
#            hBkg.SetLineWidth(2)
#        else:
#            hBkg.SetLineWidth(2)

        dictHistos[ivar][ 'data' ].Draw('hist')
        dictHistos[ivar][ 'unfoldHisto' ].Draw('hist same')

        CMS_lumi.extraText = " Preliminary"
        CMS_lumi.CMS_lumi( pad[ivar+'1'], 4, 0)
        legend.Draw()

        pad[ivar+'2'].cd()
        pad[ivar+'2'].SetGrid()
        pad[ivar+'2'].SetTopMargin(0)
        pad[ivar+'2'].SetBottomMargin(0.3)

        tmppad= pad[ivar+'2'].DrawFrame(xmin,0.5,xmax,1.5)
        tmppad.GetYaxis().SetTitle( "Data/Sim." )
        tmppad.GetXaxis().SetTitle(  axisX )
        tmppad.GetYaxis().SetTitleOffset( 0.5 )
        tmppad.GetYaxis().CenterTitle()
        tmppad.SetLabelSize(0.12, 'x')
        tmppad.SetTitleSize(0.12, 'x')
        tmppad.SetLabelSize(0.12, 'y')
        tmppad.SetTitleSize(0.12, 'y')
        tmppad.SetNdivisions(505, 'x')
        tmppad.SetNdivisions(505, 'y')
        pad[ivar+'2'].Modified()
        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dictHistos[ivar][ 'data' ], dictHistos[ivar][ 'unfoldHisto' ], 'pois' )
        hRatio.Draw('P')

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    del canvas[outputFileName]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='bkgData', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
    parser.add_argument('-s', '--selection', action='store', default='dijetSel', help='Selection: dijetSel, WSel, topSel.' )
    parser.add_argument('-v', '--version', action='store', default='v0', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    #parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    VER = args.version.split('_')[1] if '_' in args.version else args.version
    args.lumi = checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[args.year]['lumi']

    tmpFile = '/afs/cern.ch/user/a/algomez/cernbox/JetObservables/Updates/20210510/Plots/Unfold/'
    dictFiles = OrderedDict()
    for ivar, varInfo in nSubVariables.items():
        if args.selection.startswith('dijet') and ivar.startswith(('Jet1', 'Jet2')):
            dictFiles[ivar] = ROOT.TFile.Open( tmpFile+'/'+args.year+'/'+ivar+'/MCClosure/outputHistograms_QCD_HT2000toInf.root' )

    outputDir = 'Plots/'#+args.selection+'/'+('Resolution' if args.process.startswith('reso') else 'Basic')+'/'+args.year
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+args.year

    plotList = []
    plotList.append( [ 'recoJet1', nSubVariables['Jet1_tau21']['label'], nSubVariables['Jet1_tau21']['bins'][0], nSubVariables['Jet1_tau21']['bins'][-1], nSubVariables['Jet1_tau21']['bins'], nSubVariables['Jet1_tau21']['alignLeg'] ] )

    for i in plotList:
        if ( 'bkgData' in args.process ):
            plotSignalBkg( i[0]+'_'+args.selection, i[2], i[3], i[4],
                            log=args.log, axisX=i[1], legendAlignment=i[5], outputDir=outputDir)

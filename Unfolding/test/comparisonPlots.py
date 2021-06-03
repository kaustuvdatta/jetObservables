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
def makePlots( name, ext='png', outputDir='Plots/' ):
    """function to plot s and b histos"""

    numBins = 0
    dictHistos = OrderedDict()
    for ivar in nSubVariables.keys():
        if name.endswith( ivar.split('_')[0] ):
            dictHistos[ivar] = {
                    'data' : dictFiles[ivar].Get( 'dataMinusBkgsGenBin' ),
                    'unfold' : dictFiles[ivar].Get( 'unfoldHisto'+ivar ),
                    'unfoldwoUnc' : dictFiles[ivar].Get( 'unfoldHistowoUnc'+ivar ),
                    'fold' : dictFiles[ivar].Get( 'foldHisto'+ivar ),
                    'gen' : dictFiles[ivar].Get( mainQCD+'_accepgen'+ivar+'_'+args.selection ),
                    'reco' : dictFiles[ivar].Get( mainQCD+'_reco'+ivar+'_nom_'+args.selection ),
                    'data'+secondQCD : dictFiles[ivar+secondQCD].Get( 'dataMinusBkgsGenBin' ),
                    'unfold'+secondQCD : dictFiles[ivar+secondQCD].Get( 'unfoldHisto'+ivar ),
                    'unfoldwoUnc'+secondQCD : dictFiles[ivar+secondQCD].Get( 'unfoldHistowoUnc'+ivar ),
                    'fold'+secondQCD : dictFiles[ivar+secondQCD].Get( 'foldHisto'+ivar ),
                    'gen'+secondQCD : dictFiles[ivar+secondQCD].Get( secondQCD+'_accepgen'+ivar+'_'+args.selection ),
                    'reco'+secondQCD : dictFiles[ivar+secondQCD].Get( secondQCD+'_reco'+ivar+'_nom_'+args.selection )
                    }
            print(dictHistos)
            numBins = numBins + dictHistos[ivar]['data'].GetNbinsX()

            ############################################
            ## Self Clsoure
            otherHistos = OrderedDict()
            mainHisto = dictHistos[ivar]['unfold']
            if args.process.startswith('MCClosure'):
                mainHistoLabel = 'Unfold'
                otherHistos['Accepted Gen'] = dictHistos[ivar]['gen']
                otherHistos['Folded'] = dictHistos[ivar]['fold']
            elif args.process.startswith('MCSelfClosure'):
                mainHistoLabel = 'MC reco'
                otherHistos['Accepted Gen'] = dictHistos[ivar]['gen']
                otherHistos['Folded'] = dictHistos[ivar]['fold']
                otherHistos['True Reco'] = dictHistos[ivar]['reco']
            elif args.process.startswith('data'):
                mainHistoLabel = 'Data'
                otherHistos['Accepted Gen'] = dictHistos[ivar]['gen']
            drawUnfold( ivar, mainHisto, mainHistoLabel, otherHistos, nSubVariables[ivar], outputDir=outputDir )

    combinePlots( name+'_'+args.selection, dictHistos, numBins, outputDir=outputDir, axisX='' )

##########################################################################
def drawUnfold( ivar, mainHisto, mainHistoLabel, otherHistos, varInfo, outputDir ):
    """docstring for drawUnfold"""

    outputDir=outputDir+'/'+ivar+'/'+args.process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+args.selection+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_Tunfold_'+args.version+'.'+args.ext

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    if varInfo['alignLeg'].startswith('right'): legend=ROOT.TLegend(0.65,0.65,0.90,0.88)
    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    dataIntegral = mainHisto.Integral()
    mainHisto.Scale(1, 'width')  ### divide by bin width
    #######mainHisto.Scale(1/mainHisto.Integral(), 'width')  ### divide by bin width
    mainHisto.SetLineWidth(2)
    mainHisto.SetLineStyle(2)
    mainHisto.SetLineColor(ROOT.kMagenta)
    legend.AddEntry( mainHisto, mainHistoLabel, 'l' )

    tmpMax = [ mainHisto.GetMaximum() ]
    for ilabel, ih in otherHistos.items():
        #ih.Scale(1, 'width')
        #ih.Scale(scaleFactor)
        ih.Scale(dataIntegral/ih.Integral(), 'width')  ### divide by bin width
        ih.SetLineWidth(2)
        ih.SetLineColor(1)
        tmpMax.append( ih.GetMaximum() )
        legend.AddEntry( ih, ilabel, 'l' )


    mainHisto.GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+varInfo['label'].split('#')[1]+'}' )
    #mainHisto.GetYaxis().SetTitleOffset(0.95)
    mainHisto.GetYaxis().SetTitleSize(0.05)
    mainHisto.SetMaximum( 1.2*max( tmpMax )  )

    mainHisto.Draw( "histe")
    for ih in otherHistos: otherHistos[ih].Draw("histe same")
    mainHisto.Draw( "histe same")

    legend.Draw()
    CMS_lumi.relPosX = 0.12
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

    can.SaveAs(outputName)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)


##################################################
def combinePlots( name, dictHistos, numBins, axisX='', outputDir='Plots/'):
    """docstring for combinePlots"""

    outputFileName = name+'_combinePlots_'+args.version+'.'+args.ext
    if args.log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.50,0.80,0.90,0.90)
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.07)

    dictHistos['combData'] = ROOT.TH1F('combData', 'combData', numBins, 0, numBins)
    legend.AddEntry( dictHistos[ 'combData' ], 'Data', 'lep' )
    dictHistos['combUnfold'] = ROOT.TH1F('combUnfold', 'combUnfold', numBins, 0, numBins)
    legend.AddEntry( dictHistos[ 'combUnfold' ], 'Unfold', 'lep' )

    tmpNbin = 0
    for ivar,ih in dictHistos.items():
        if ivar.startswith('Jet'):
            for ibin in range(1, ih['data'].GetNbinsX()+1):
                tmpNbin = tmpNbin+1
                dictHistos['combData'].SetBinContent( tmpNbin, ih['data'].GetBinContent(ibin) )
                dictHistos['combData'].SetBinError( tmpNbin, ih['data'].GetBinError(ibin) )
                dictHistos['combUnfold'].SetBinContent( tmpNbin, ih['unfold'].GetBinContent(ibin) )
                dictHistos['combUnfold'].SetBinError( tmpNbin, ih['unfold'].GetBinError(ibin) )


    canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 750 )
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.10)
    pad1 = ROOT.TPad(ivar+'1', "Fit",0.,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad(ivar+'2', "Pull",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    if args.log: pad1.SetLogy()
    dictHistos['combData'].GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#tau_{X}}' )
    dictHistos['combData'].GetYaxis().SetTitleSize( 0.05 )
    dictHistos['combData'].GetYaxis().SetTitleOffset( 0.8 )
    dictHistos['combData'].SetMaximum( dictHistos['combUnfold'].GetMaximum()*1.2 )
    dictHistos['combData'].SetMinimum( 0.001 )
    dictHistos['combData'].SetLineColor( ROOT.kBlack )
    dictHistos['combData'].SetLineWidth( 2 )

    dictHistos['combUnfold'].SetLineColor( ROOT.kMagenta )
    dictHistos['combUnfold'].SetLineWidth( 2 )
    dictHistos['combData'].Draw('E')
    dictHistos[ 'combUnfold' ].Draw('E same')

    CMS_lumi.extraText = " Preliminary"
    CMS_lumi.relPosX = 0.07
    CMS_lumi.CMS_lumi( pad1, 4, 0)
    legend.Draw()

    pad2.cd()
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)

    tmppad= pad2.DrawFrame(0,0.,numBins,1.95)
    tmppad.GetYaxis().SetTitle( "Data/Unfold" )
    tmppad.GetXaxis().SetTitle(  axisX )
    tmppad.GetYaxis().SetTitleOffset( 0.3 )
    tmppad.GetYaxis().CenterTitle()
    tmppad.SetLabelSize(0.12, 'x')
    tmppad.SetTitleSize(0.12, 'x')
    tmppad.SetLabelSize(0.12, 'y')
    tmppad.SetTitleSize(0.12, 'y')
    tmppad.SetNdivisions(505, 'x')
    tmppad.SetNdivisions(505, 'y')
    pad2.Modified()
    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( dictHistos['combData'], dictHistos[ 'combUnfold' ], 'pois' )
    hRatio.Draw('P')

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    del canvas[outputFileName]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='MCClosure', dest='process', help='Process to draw, example: data, MCClosure, MCSelfClosure.' )
    parser.add_argument('-s', '--selection', action='store', default='dijetSel', help='Selection: dijetSel, WSel, topSel.' )
    parser.add_argument('-v', '--version', action='store', default='v02', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    #parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )
    parser.add_argument("--QCDHT", action='store_true', dest="QCDHT", default=False, help="For dijet sel, signal QCD is HT or not" )
    parser.add_argument("--inputFolder", action='store', dest="inputFolder", default="", help="input folder" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    VER = args.version.split('_')[1] if '_' in args.version else args.version
    args.lumi = checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[args.year]['lumi']

    mainQCD = 'QCD_HT2000toInf' if args.QCDHT else 'QCD_Pt_3200toInf'
    secondQCD =  'QCD_Pt_3200toInf' if args.QCDHT else 'QCD_HT2000toInf'
    dictFiles = OrderedDict()
    for ivar, varInfo in nSubVariables.items():
        if args.selection.startswith('dijet') and ivar.startswith(('Jet1', 'Jet2')):
            dictFiles[ivar] = ROOT.TFile.Open( args.inputFolder+'Plots/'+args.selection+'/Unfold/'+args.year+'/'+ivar+'/'+args.process+'/outputHistograms_'+mainQCD+'.root' )
            dictFiles[ivar+secondQCD] = ROOT.TFile.Open( args.inputFolder+'Plots/'+args.selection+'/Unfold/'+args.year+'/'+ivar+'/'+args.process+'/outputHistograms_'+secondQCD+'.root' )

    outputDir = 'Plots/'+args.selection+'/Unfold/'+args.year
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    CMS_lumi.extraText = 'Preliminary' if args.process.startswith('data') else "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+args.year

    plotList = [ 'recoJet1', 'recoJet2' ] if args.selection.startswith('dijet') else [ 'recoJet' ]

    for i in plotList:
        makePlots( i, outputDir=outputDir)

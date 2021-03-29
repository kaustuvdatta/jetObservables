#!/usr/bin/env python

import ROOT
import time, os, math, sys, copy
from array import array
import argparse
from collections import OrderedDict
import subprocess
#from makeDatacards import dataFiles, bkgFiles, sigFiles
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrStyle
from datasets import checkDict, dictSamples

####gReset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()
tdrStyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)

selection = {}
#selection['SL_presel'] = [ 'nlep > 0', 'nJets > 3', 'nDeepCSVM > 1' ]

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
        histo.Draw('hist')
        histo2.Draw('hist same')
    if not axisY: histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: histo.GetXaxis().SetTitle( axisX )

    legend.Draw()

    canvas[name].SaveAs( outputDir+outName )
    #del can

######################################################
def plotSysComparison( nomHisto, upHisto, downHisto, outputName, syst, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ):
    """docstring for plot"""

    outputFileName = outputName+'_'+syst+'SystPlots_'+version+'.'+ext
    print 'Processing.......', outputFileName

    binWidth = nomHisto.GetBinWidth(1)

    legend=ROOT.TLegend(0.70,0.75,0.90,0.87)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.AddEntry( nomHisto, 'Nominal' , 'l' )
    legend.AddEntry( upHisto, syst+'Up', 'l' )
    legend.AddEntry( downHisto, syst+'Down', 'l' )

    nomHisto.SetLineColor(ROOT.kBlack)
    nomHisto.SetLineWidth(2)
    upHisto.SetLineColor(ROOT.kRed-4)
    upHisto.SetLineWidth(2)
    upHisto.SetLineStyle(2)
    downHisto.SetLineColor(ROOT.kBlue)
    downHisto.SetLineWidth(2)
    downHisto.SetLineStyle(2)

    hRatioUp = ROOT.TGraphAsymmErrors()
    hRatioUp.Divide( nomHisto, upHisto, 'pois' )
    hRatioUp.SetLineColor(ROOT.kRed-4)
    hRatioUp.SetLineWidth(2)
    hRatioDown = ROOT.TGraphAsymmErrors()
    hRatioDown.Divide( nomHisto, downHisto, 'pois' )
    hRatioDown.SetLineColor(ROOT.kBlue)
    hRatioDown.SetLineWidth(2)

    #tdrStyle.SetPadRightMargin(0.05)
    #tdrStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('c1'+outputName, 'c1'+outputName,  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1"+outputName, "Fit",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+outputName, "Pull",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    if log: pad1.SetLogy()
    nomHisto.Draw("E")
    upHisto.Draw('hist same')
    downHisto.Draw("hist same ")
    #hData.SetMaximum( 1.2* max( hData.GetMaximum(), hBkg.GetMaximum() )  )
    #if 'pt' in label: hData.SetMinimum( 1 )
    #hData.GetYaxis().SetTitleOffset(1.2)
    #if xmax: hData.GetXaxis().SetRangeUser( xmin, xmax )
    nomHisto.GetXaxis().SetTitle( labelX )
    nomHisto.GetYaxis().SetTitle( 'Normalized / '+str(int(binWidth)) )

    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+year
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    legend.Draw()

    pad2.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame( 0, 0.5, 1, 1.5 )
    tmpPad2.GetXaxis().SetTitle( nomHisto.GetXaxis().GetTitle()  )
    tmpPad2.GetYaxis().SetTitle( "Ratio to Nom" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    #hRatioUp.SetMarkerStyle(8)
    hRatioUp.Draw('P')
    hRatioDown.Draw('P same')

    can.SaveAs( outputDir + outputFileName )
    del can

##################################################

######################################################
def plotQuality( nameInRoot, label, xmin, xmax, rebinX, labX, labY, log, moveCMSlogo=False, fitRatio=False ):
    """docstring for plot"""

    outputFileName = nameInRoot+'_'+args.ttbarDecay+'_dataQualityPlots_'+args.version+'_'+args.year+'.'+args.ext
    print 'Processing.......', outputFileName

    histos = {}

    for idataLabel, idata in dataFile.iteritems():
        try:
            histos[ 'Data' ].Add( idata.Get( 'jetObservables/'+nameInRoot ) )
        except (KeyError, AttributeError) as e:
            histos[ 'Data' ] = idata.Get( 'jetObservables/'+nameInRoot.split('Total')[0] )

    histos[ 'Bkg' ] = histos[ 'Data' ].Clone()
    histos[ 'Bkg' ].Reset()
    for isamLabel, isam in bkgFiles.iteritems():
        histos[ isamLabel ] = isam[0].Get( 'jetObservables/'+nameInRoot )
        histos[ isamLabel ].Scale( isam[1] )
        histos[ 'Bkg' ].Add( histos[ isamLabel ] )

    if rebinX != 1:
        histos[ 'Data' ].Rebin( rebinX )
        histos[ 'Bkg' ].Rebin( rebinX )
    hData = histos[ 'Data' ].Clone()
    hBkg = histos[ 'Bkg' ].Clone()

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( hData, hBkg, 'pois' )
    hRatioStatErr = hBkg.Clone()
    hRatioStatErr.Divide( hBkg )
    hRatioStatErr.SetFillColor(kBlack)
    hRatioStatErr.SetFillStyle(3004)

    binWidth = histos['Data'].GetBinWidth(1)

    if (labY < 0.5) and ( labX < 0.5 ): legend=ROOT.TLegend(0.20,0.50,0.50,0.62)
    elif (labX < 0.5): legend=ROOT.TLegend(0.20,0.75,0.50,0.87)
    else: legend=ROOT.TLegend(0.70,0.75,0.90,0.87)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.AddEntry( hData, 'DATA' , 'ep' )
    legend.AddEntry( hBkg, 'All Bkg', 'lp' )

    hBkg.SetLineColor(kRed-4)
    hBkg.SetLineWidth(2)
    #hBkg.SetFillColor(kBlack)
    hBkg.SetFillStyle(3004)
    hData.SetMarkerStyle(8)

    tdrStyle.SetPadRightMargin(0.05)
    tdrStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('c1', 'c1',  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    if log: pad1.SetLogy()
    hData.Draw("E")
    hBkg.Draw('hist same E1')
    hData.Draw("same E")
    hData.SetMaximum( 1.2* max( hData.GetMaximum(), hBkg.GetMaximum() )  )
    if 'pt' in label: hData.SetMinimum( 1 )
    #hData.GetYaxis().SetTitleOffset(1.2)
    if xmax: hData.GetXaxis().SetRangeUser( xmin, xmax )
    #hData.GetYaxis().SetTitle( 'Normalized' )
    #hData.GetYaxis().SetTitle( 'Normalized / '+str(int(binWidth))+' GeV' )
    hData.GetYaxis().SetTitle( ( 'Events / '+str(int(binWidth))+' GeV' if nameInRoot in [ 'massAve', 'HT', 'jet1Pt', 'jet2Pt', 'MET' ] else 'Events' ) )

    #CMS_lumi.relPosX = 0.13
    if moveCMSlogo:
        CMS_lumi.cmsTextOffset = 0.1
        CMS_lumi.relPosX = 0.15
    else:
        CMS_lumi.cmsTextOffset = 0.0
        CMS_lumi.relPosX = 0.13
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    #labelAxis( name, hData, '' )
    legend.Draw()
    #setSelection( selection[ args.ttbarDecay+'_'+args.cut ], labX, labY )

    pad2.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame(xmin, ( 0.5 if fitRatio else 0.5), xmax,1.5)
    #labelAxis( name.replace( args.cut, ''), tmpPad2, ( 'softDrop' if 'Puppi' in args.grooming else Groom ) )
    tmpPad2.GetXaxis().SetTitle( label )
    tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    hRatio.SetMarkerStyle(8)
    hRatio.Draw('P')
    hRatioStatErr.Draw('same e2')
    if fitRatio:
        fitLine = TF1( 'fitLine', 'pol1', 0, 2 ) #800, 5000)
        hRatio.Fit( 'fitLine', 'MIR')
        fitLine.Draw("same")
        pad2.Update()
        st1 = hRatio.GetListOfFunctions().FindObject("stats")
        st1.SetX1NDC(.65)
        st1.SetX2NDC(.95)
        st1.SetY1NDC(.75)
        st1.SetY2NDC(.95)
        #st1.SetTextColor(kRed)
        pad2.Modified()

    can.SaveAs( 'Plots/Basic/'+ outputFileName.replace('Plots', ( 'Fit' if fitRatio else '') ) )
    del can

##################################################
def plotSignalBkg( name, xmin, xmax, rebinX, axisX='', axisY='', labX=0.92, labY=0.50, log=False,
                      addRatioFit=False, Norm=False, ext='png' ):
    """function to plot s and b histos"""

    outputFileName = name+'_PlusBkg_AnalysisPlots_'+args.version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    if Norm: outputFileName = outputFileName.replace('Plots','Plots_Normalized')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.60,0.80,0.90,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    dataHistos = {}
    if args.selection.startswith('dijet'):
        tmpdataHistos = {}
        for it in checkDict( 'JetHT', dictSamples )[args.year]['triggerList']:
            tmpdataHistos[ it ] = dataFile['data'].Get( 'jetObservables/'+name.replace( args.selection, it+'_'+args.selection ) )
            tmpdataHistos[ it ].Scale( checkDict( 'JetHT', dictSamples )[args.year]['triggerList'][it] )
        dataHistos[ 'DATA' ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
        dataHistos[ 'DATA' ].Reset()
        for i in tmpdataHistos: dataHistos['DATA'].Add( tmpdataHistos[i] )
    else:
        dataHistos[ 'DATA' ] = dataFile['data'].Get( 'jetObservables/'+name )
    if rebinX > 1: dataHistos[ 'DATA' ] = dataHistos[ 'DATA' ].Rebin( rebinX )
    legend.AddEntry( dataHistos[ 'DATA' ], 'Data', 'lep' )
    if Norm: dataHistos[ 'DATA' ].Scale( 1 /dataHistos['DATA'].Integral() )

    bkgHistos = OrderedDict()
    maxList = []
    if len(bkgFiles) > 0:
        for bkgSamples in bkgFiles:
            bkgHistos[ bkgSamples ] = bkgFiles[ bkgSamples ][0].Get( 'jetObservables/'+name )
            bkgHistos[ bkgSamples ].SetTitle(bkgSamples)
            bkgHistos[ bkgSamples ].Scale( args.lumi*bkgFiles[ bkgSamples ][1]['XS'] / bkgFiles[ bkgSamples ][1][args.year]['nGenWeights'] )
            print(bkgSamples, round(bkgHistos[ bkgSamples ].Integral(), 2) )
            if rebinX > 1: bkgHistos[ bkgSamples ] = bkgHistos[ bkgSamples ].Rebin( rebinX )

            if Norm:
                bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][1]['color'] )
                bkgHistos[ bkgSamples ].SetLineWidth( 2 )
                try: bkgHistos[ bkgSamples ].Scale( 1 / bkgHistos[ bkgSamples ].Integral() )
                except ZeroDivisionError: pass
                maxList.append( bkgHistos[ bkgSamples ].GetMaximum() )
            else:
                #bkgHistos[ bkgSamples ].SetFillStyle( 1001 )
                #bkgHistos[ bkgSamples ].SetFillColor( int(bkgFiles[ bkgSamples ][1]['color']) )
                bkgHistos[ bkgSamples ].SetLineColor( bkgFiles[ bkgSamples ][1]['color'] )


#    #### Merging samples
#    for bkg in bkgFiles:
#        if bkg.endswith(('WZ','ZZ')):
#            bkgHistos['WW'].Add( bkgHistos[bkg] )
#            bkgHistos.pop(bkg, None)
#        elif bkg.startswith('ST_t'):
#            bkgHistos['ST_s-channel'].Add( bkgHistos[bkg] )
#            bkgHistos.pop(bkg, None)
#        elif bkg.startswith('QCD') and not bkg.endswith('Inf'):
#            bkgHistos['QCD_Pt_3200toInf'].Add( bkgHistos[bkg] )
#            bkgHistos.pop(bkg, None)
#        else:
#            legend.AddEntry( bkgHistos[ bkg ], bkg, 'l' if Norm else 'f' )

    hBkg = bkgHistos[next(iter(bkgHistos))].Clone()
    hBkg.Reset()
    binWidth = dataHistos['DATA'].GetBinWidth(1)

    if not Norm:

        stackHisto = ROOT.THStack('stackHisto'+name, 'stack'+name)
        for samples in bkgHistos:
            stackHisto.Add( bkgHistos[ samples ].Clone() )
            hBkg.Add( bkgHistos[ samples ].Clone() )
            legend.AddEntry( bkgHistos[ samples ], bkgFiles[ samples ][1]['label'], 'le' )

        canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name,  10, 10, 750, 750 )
        #tdrStyle.SetPadRightMargin(0.05)
        #tdrStyle.SetPadLeftMargin(0.15)
        pad1 = ROOT.TPad("pad1", "Fit",0,0.207,1.00,1.00,-1)
        pad2 = ROOT.TPad("pad2", "Pull",0,0.00,1.00,0.30,-1);
        pad1.Draw()
        pad2.Draw()

        pad1.cd()
        #if log and not args.final: pad1.SetLogy()
        if log: pad1.SetLogy()
        dataHistos[ 'DATA' ].SetMarkerStyle(8)
        dataHistos['DATA'].GetYaxis().SetTitle( 'Events / '+str(binWidth) )
        dataHistos['DATA'].GetYaxis().SetTitleSize( 0.05 )
        dataHistos['DATA'].GetYaxis().SetTitleOffset( 1.2 )
        dataHistos['DATA'].GetXaxis().SetTitle( axisX )
        if xmax: dataHistos['DATA'].GetXaxis().SetRangeUser( xmin, xmax )
        dataHistos['DATA'].SetMaximum( hBkg.GetMaximum()*1.2 )
        dataHistos['DATA'].SetMinimum( 2. )

        dataHistos[ 'DATA' ].Draw('E')
        CMS_lumi.extraText = "Preliminary"
        stackHisto.Draw('hist same')

        #stackHisto.SetMinimum( 0.1 )
        #hBkg.SetFillStyle(0)
        hBkg.SetLineColor(ROOT.kBlack)
        hBkg.SetLineStyle(1)
        hBkg.SetLineWidth(1)
        #hBkg.SetFillStyle(3004)
        #hBkg.SetFillColor( kRed )
        #hBkg.Draw("same")


        #CMS_lumi.relPosX = 0.14
        CMS_lumi.CMS_lumi( pad1, 4, 0)
        legend.Draw()

        pad2.cd()
        pad2.SetGrid()
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.3)

        tmpPad2= pad2.DrawFrame(xmin,0.5,xmax,1.5)
        tmpPad2.GetYaxis().SetTitle( "Data/Bkg" )
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
        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dataHistos[ 'DATA' ], hBkg, 'pois' )
        hRatio.SetMarkerStyle(8)
        hRatio.Draw('P')
        canvas[outputFileName].SaveAs( 'Plots/Basic/'+outputFileName )

    else:

        tdrStyle.SetPadRightMargin(0.05)
        canvas[outputFileName]= ROOT.TCanvas('c1', 'c1', 750, 500 )
        if log: canvas[outputFileName].SetLogy()
        signalHistos[next(iter(signalHistos))].GetYaxis().SetTitleOffset(1.0)
        signalHistos[next(iter(signalHistos))].GetYaxis().SetTitle( ( 'Normalized / '+str(int(binWidth))+' GeV' ) )
        if xmax: signalHistos[next(iter(signalHistos))].GetXaxis().SetRangeUser( xmin, xmax )
        signalHistos[next(iter(signalHistos))].Draw('hist')
        for signalSamples in signalHistos: signalHistos[ signalSamples ].Draw('hist same')
        for bkgSamples in bkgHistos: bkgHistos[ bkgSamples ].Draw('hist same')
        if 'DATA' in args.process:
                dataHistos[ 'DATA' ].SetMarkerStyle(8)
                dataHistos[ 'DATA' ].Draw('same')
                CMS_lumi.extraText = ""#"Preliminary"
        signalHistos[next(iter(signalHistos))].SetMaximum( 1.1 * max( maxList ) )

        if not 'DATA' in args.process: CMS_lumi.lumi_13TeV = ''
        CMS_lumi.relPosX = 0.11
        CMS_lumi.CMS_lumi(canvas[outputFileName], 4, 0)
        legend.Draw()

        canvas[outputFileName].SaveAs( 'Plots/Basic/'+outputFileName )
    del canvas[outputFileName]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proc', action='store', default='signalBkg', dest='process', help='Process to draw, example: 1D, 2D, MC.' )
    parser.add_argument('-s', '--selection', action='store', default='dijetSel', help='Selection: dijetSel, WSel, topSel.' )
    parser.add_argument('-v', '--version', action='store', default='v0', help='Version: v01, v02.' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    #parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('-u', '--unc', action='store', default='JES', dest='unc',  help='Type of uncertainty' )
    parser.add_argument('-L', '--log', action='store_true', default=False, dest='log',  help='Plot in log scale (true) or not (false)' )
    parser.add_argument('-n', '--norm', action='store_true', default=False, dest='norm',  help='Normalized plot (true) or not (false)' )
    parser.add_argument('-F', '--addFit', action='store_true', default=False, dest='addFit',  help='Plot fit in ratio plot.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)


    VER = args.version.split('_')[1] if '_' in args.version else args.version
    dataFile = {}
    dataFile['data'] = ROOT.TFile.Open(checkDict( 'JetHT', dictSamples )[args.year]['skimmerHisto'])
    args.lumi = checkDict( 'JetHT', dictSamples )[args.year]['lumi']

    bkgFiles = {}
    bkgFiles['QCDHerwig'] = [
                            ROOT.TFile.Open( checkDict( 'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7', dictSamples )[args.year]['skimmerHisto'] ),
                            checkDict( 'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7', dictSamples )
            ]
#    #bkgFiles['TTJets'] = [ ROOT.TFile.Open('Rootfiles/TTJets_amcatnloFXFX-pythia8_2017UL.root'), 722.8, 298034800561.50, kWhite ] #from xsdb
#    bkgFiles['TTToSemiLeptonic'] = [ ROOT.TFile.Open('Rootfiles/TTToSemileptonic_powheg_pythia8_2017UL_v1recheck.root'), 365.34, 30155734578.508015, kWhite]
#    bkgFiles['TTTo2L2Nu'] = [ ROOT.TFile.Open('Rootfiles/TTTo2L2Nu_powheg_pythia8_2017UL.root'), 88.29 , 4778168320.403400, kGreen ]
#
#    bkgFiles['WJetsToLNu'] = [ ROOT.TFile.Open('Rootfiles/WJetsToLNu_madgraphMLM_pythia8_2017UL.root'), 5.368e+04, 10501241224509.000000, kMagenta]
#
#    bkgFiles['ST_s-channel'] = [ ROOT.TFile.Open('Rootfiles/ST_s-channel_amcatnlo_pythia8_2017UL.root'), 3.549e+00, 70465639.296516, kBlue ]
#    bkgFiles['ST_t-channel_top'] = [ ROOT.TFile.Open('Rootfiles/ST_t-channel_top_powheg_pythia8_2017UL.root'), 1.197e+02, 655197632.514000, kBlue]
#    bkgFiles['ST_t-channel_antitop'] = [ ROOT.TFile.Open('Rootfiles/ST_t-channel_antitop_powheg_pythia8_2017UL.root'), 7.174e+01, 265124234.844000, kBlue]
#    #bkgFiles['ST_t-channel_antitop'] = [ ROOT.TFile.Open('Rootfiles/ST_t-channel_antitop_powheg_pythia8_2017UL.root'), 7.174e+01, 461187168.978000, kBlue]
#    bkgFiles['ST_tW_top'] = [ ROOT.TFile.Open('Rootfiles/ST_tW_top_powheg_pythia8_2017UL.root'), 3.245e+01, 325819206.029700, kBlue]
#    #bkgFiles['ST_tW_top'] = [ ROOT.TFile.Open('Rootfiles/ST_tW_top_powheg_pythia8_2017UL.root'), 3.245e+01, 651638412.059400, kBlue]
#    bkgFiles['ST_tW_antitop'] = [ ROOT.TFile.Open('Rootfiles/ST_tW_antitop_powheg_pythia8_2017UL.root'), 3.251e+01, 298789163.881200, kBlue]
#
#    bkgFiles['QCD_Pt_170to300'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_170to300_pythia8_2017UL.root'), 1.025e+05, 29522100.000000, kBlack ]
#    bkgFiles['QCD_Pt_300to470'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_300to470_pythia8_2017UL.root'), 6.762e+03, 57365509.961543, kBlack]
#    bkgFiles['QCD_Pt_470to600'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_470to600_pythia8_2017UL.root'), 5.461e+02, 27559688.974355, kBlack]
#    bkgFiles['QCD_Pt_600to800'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_600to800_pythia8_2017UL.root'), 1.549e+02, 65128701.570204, kBlack]
#    bkgFiles['QCD_Pt_800to1000'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_800to1000_pythia8_2017UL.root'), 2.597e+01, 39261600.000000, kBlack]
#    bkgFiles['QCD_Pt_1000to1400'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_1000to1400_pythia8_2017UL.root'), 7.398e+00, 19967700.000000, kBlack]
#    bkgFiles['QCD_Pt_1400to1800'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_1400to1800_pythia8_2017UL.root'), 6.423e-01, 5434800.000000, kBlack]
#    bkgFiles['QCD_Pt_1800to2400'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_1800to2400_pythia8_2017UL.root'), 8.671e-02, 2999700.000000, kBlack]
#    bkgFiles['QCD_Pt_2400to3200'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_2400to3200_pythia8_2017UL.root'), 5.193e-03, 1919400.000000, kBlack]
#    bkgFiles['QCD_Pt_3200toInf'] = [ ROOT.TFile.Open('Rootfiles/QCD_Pt_3200toInf_pythia8_2017UL.root'), 1.340e-04, 800000.000000, kBlack]
#    #bkgFiles['QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8'] = [ ROOT.TFile.Open('Rootfiles/'), 1.347e+09, ]
#    #bkgFiles['QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7'] = [ ROOT.TFile.Open('Rootfiles/'), 7.149e+08, ]
#
#    bkgFiles['WW'] = [ ROOT.TFile.Open('Rootfiles/WW_pythia8_2017UL.root'), 7.577e+01, 7876265.259367, kYellow ]
#    bkgFiles['WZ'] = [ ROOT.TFile.Open('Rootfiles/WZ_pythia8_2017UL.root'), 1.21e+00, 3970000.000000, kYellow ]
#    bkgFiles['ZZ'] = [ ROOT.TFile.Open('Rootfiles/ZZ_pythia8_2017UL.root'), 2.748e+00,1981800.000000, kYellow  ]

    if not os.path.exists('Plots/Basic'): os.makedirs('Plots/Basic')
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, "+args.year

    plotList = [
            [ 'qual', 'leadAK8JetTau21', 'Leading AK8 jet #tau_{21}', 0, 1, 2, 0.85, 0.70, True, False ],
            [ 'qual', 'recoJetPt', 'Leading AK8 jet pt [GeV]', 100, 1500, 2, 0.85, 0.70, True, False ],
            [ 'signalBkg', 'recoJet1_sortedPt_pt', 'Leading AK8 jet pt [GeV]', 100, 1500, 5 ] ,
            [ 'signalBkg', 'recoJet2_sortedPt_pt', '2nd Leading AK8 jet pt [GeV]', 100, 1500, 5 ],
            ]

    for ijet in [ ('Jet1', 'Outer'), ('Jet2', 'Central') ]:
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_pt_nom', ijet[1]+' AK8 jet pt [GeV]', 100, 1500, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_eta_nom', ijet[1]+' AK8 jet eta', -3, 3, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau21_nom', ijet[1]+' AK8 jet #tau_{21}', 0, 1, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau32_nom', ijet[1]+' AK8 jet #tau_{32}', 0, 1, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_0p5_1_nom', ijet[1]+' AK8 jet #tau_{1}^{0.5}', 0, 1, 10 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_0p5_2_nom', ijet[1]+' AK8 jet #tau_{2}^{0.5} ', 0, 0.8, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_0p5_3_nom', ijet[1]+' AK8 jet #tau_{3}^{0.5} ', 0, 0.6, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_0p5_4_nom', ijet[1]+' AK8 jet #tau_{4}^{0.5} ', 0, 0.6, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_1_1_nom', ijet[1]+' AK8 jet #tau_{1}^{1}', 0, 1, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_1_2_nom', ijet[1]+' AK8 jet #tau_{2}^{1} ', 0, 0.6, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_1_3_nom', ijet[1]+' AK8 jet #tau_{3}^{1} ', 0, 0.4, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_1_4_nom', ijet[1]+' AK8 jet #tau_{4}^{1} ', 0, 0.4, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_2_1_nom', ijet[1]+' AK8 jet #tau_{1}^{2}', 0, 1, 5 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_2_2_nom', ijet[1]+' AK8 jet #tau_{2}^{2} ', 0, 0.4, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_2_3_nom', ijet[1]+' AK8 jet #tau_{3}^{2} ', 0, 0.4, 20 ] )
        plotList.append( [ 'signalBkg', 'reco'+ijet[0]+'_tau_2_4_nom', ijet[1]+' AK8 jet #tau_{4}^{2} ', 0, 0.4, 20 ] )

    Plots = [ x[1:] for x in plotList if ( ( args.process in x[0] ) )  ]
    #else: Plots = [ y[1:] for y in plotList if ( ( args.process in y[0] ) and ( y[1] in args.single ) )  ]

    for i in Plots:
        if ( 'signalBkg' in args.process ):
            plotSignalBkg( i[0]+'_'+args.selection, i[2], i[3], i[4], log=args.log, axisX=i[1], Norm=args.norm)
        elif ( 'qual' in args.process ):
            plotQuality(
                i[0]+'_'+args.selection, i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8],
                fitRatio=args.addFit )
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

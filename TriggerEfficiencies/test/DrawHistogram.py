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
yline = array('d', [0.98,.98])
line = TGraph(2, xline, yline)
line.SetLineColor(kRed)

def plotTriggerEfficiency( inFileSample, sample, triggerDenom, denomScale, triggerPass, passScale, name, minFit, maxFit, xlabel, xmin, xmax, rebin, log):
    """docstring for plot"""

    outputFileName = name+'_'+triggerDenom+"_"+triggerPass+'_'+sample+'_TriggerEfficiency'+args.version+'.'+args.extension
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

    legend=TLegend(0.50,0.75,0.90,0.90)
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

    legend.AddEntry( histos[ 'denomOnly'+triggerPass ], 'SingleMu trigger', 'l' )
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

    CMS_lumi.CMS_lumi(pad1, 4, 0)
    #legend.Draw()
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

def plotDiffEff( listOfEff, name ):
	"""docstring for plotDiffEff"""

	can = TCanvas('c1', 'c1',  10, 10, 750, 500 )

	legend=TLegend(0.60,0.25,0.90,0.40)
	legend.SetFillStyle(0)
	legend.SetTextSize(0.04)

	dummy = 1
	for sample in listOfEff:
		legend.AddEntry( listOfEff[ sample ], sample, 'l' )

		listOfEff[ sample ].SetMarkerStyle(8)
		listOfEff[ sample ].SetLineWidth(2)
		listOfEff[ sample ].SetLineColor(dummy)
		listOfEff[ sample ].GetYaxis().SetTitle("Efficiency")
		listOfEff[ sample ].GetYaxis().SetLabelSize(0.06)
		listOfEff[ sample ].GetXaxis().SetLabelSize(0.06)
		listOfEff[ sample ].GetYaxis().SetTitleSize(0.06)
		listOfEff[ sample ].GetYaxis().SetTitleOffset(0.8)
		listOfEff[ sample ].SetMinimum(0.8)
		listOfEff[ sample ].SetMaximum(1.05)
		#listOfEff[ sample ].GetXaxis().SetLimits( 850, 950 )
		listOfEff[ sample ].GetXaxis().SetLimits( 0, 200 )
		if dummy == 1:
			labelAxis( name, listOfEff[ sample ], 'Pruned')
			listOfEff[ sample ].Draw()
		else:
			listOfEff[ sample ].Draw('same')
		dummy+=1

	legend.Draw('same')
	CMS_lumi.lumi_13TeV = ""
	CMS_lumi.relPosX = 0.11
	CMS_lumi.cmsTextSize = 0.7
	CMS_lumi.extraOverCmsTextSize = 0.6
	CMS_lumi.CMS_lumi(can, 4, 0)
	can.SaveAs( 'Plots/'+name+'_DiffEfficiencies.'+args.extension )
	del can

###### Rebin 2D plots
def Rebin2D( h1, rebinx, rebiny ):
	"""docstring for Rebin2D"""

	tmph1 = h1.Clone()
	nbinsx = h1.GetXaxis().GetNbins()
	nbinsy = h1.GetYaxis().GetNbins()
	xmin  = h1.GetXaxis().GetXmin()
	xmax  = h1.GetXaxis().GetXmax()
	ymin  = h1.GetYaxis().GetXmin()
	ymax  = h1.GetYaxis().GetXmax()
	nx = nbinsx/rebinx
	ny = nbinsy/rebiny
	#.SetLineColor( i )

	for biny in range( 1, nbinsy):
		for binx in range(1, nbinsx):
			ibin1 = h1.GetBin(binx,biny)
			h1.SetBinContent( ibin1, 0 )

	for biny in range( 1, nbinsy):
		by = tmph1.GetYaxis().GetBinCenter( biny )
		iy = h1.GetYaxis().FindBin(by)
		for binx in range(1, nbinsx):
			bx = tmph1.GetXaxis().GetBinCenter(binx)
			ix  = h1.GetXaxis().FindBin(bx)
			bin = tmph1.GetBin(binx,biny)
			ibin= h1.GetBin(ix,iy)
			cu  = tmph1.GetBinContent(bin)
			h1.AddBinContent(ibin,cu)
	return h1


def plot2DTriggerEfficiency( inFileSample, dataset, triggerDenom, name, xlabel, ylabel, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, labX, labY ):
	"""docstring for plot"""

	outputFileName = name+'_'+args.cut+'_'+triggerDenom+"_"+args.trigger+'_'+dataset+'_'+'TriggerEfficiency'+args.version+'.'+args.extension
	print 'Processing.......', outputFileName

	print args.trigger+'TriggerEfficiency/'+name+'Denom'+args.cut
	rawDenom = inFileSample.Get( args.trigger+'TriggerEfficiency/'+name+'Denom'+args.cut )
	Denom = Rebin2D( rawDenom, rebinx, rebiny )

	rawPassing = inFileSample.Get( args.trigger+'TriggerEfficiency/'+name+'Passing'+args.cut )
	Passing = Rebin2D( rawPassing, rebinx, rebiny )


	'''
	if ( TEfficiency.CheckConsistency( Passing, Denom ) ): Efficiency = TEfficiency( Passing, Denom )
	else:
		print '--- Passing and Denom are inconsistent.'
		#sys.exit(0)
	'''

	Efficiency = Denom.Clone()
	Efficiency.Reset()
	Efficiency.Divide( Passing, Denom, 1, 1, 'B' )


	tdrStyle.SetPadRightMargin(0.12)
	can = TCanvas('c1', 'c1',  10, 10, 1000, 750 )
	gStyle.SetPaintTextFormat("4.2f")
	Efficiency.SetMarkerSize(0.01)
	Efficiency.SetMaximum(1)
	Efficiency.SetMinimum(0)
	Efficiency.Draw('colz')
	Efficiency.Draw('same text')
	#gPad.Update()
	Efficiency.GetYaxis().SetTitleOffset(1.0)
	Efficiency.SetMarkerSize(2)
	Efficiency.GetXaxis().SetRange( int(Xmin/(rebinx)), int(Xmax/(rebinx)) )
	Efficiency.GetXaxis().SetTitle( xlabel )
	Efficiency.GetYaxis().SetTitle( ylabel )
	Efficiency.GetYaxis().SetRange( int(Ymin/(rebiny)), int(Ymax/(rebiny)) )
	#gPad.Update()

	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)

	can.SaveAs( 'Plots/'+outputFileName )
	del can

#	##### 1D Efficiency
#	newEfficiency = TGraphAsymmErrors( newPassing, newDenom, 'cp'  )
#	binWidth = rebinx
#
#	#legend=TLegend(0.50,0.75,0.90,0.90)
#	#legend.SetFillStyle(0)
#	#legend.SetTextSize(0.04)
#
#	#gStyle.SetOptFit(1)
#	can1 = TCanvas('can1', 'can1',  10, 10, 750, 500 )
#	newEfficiency.SetMarkerStyle(8)
#	newEfficiency.SetMarkerColor(kGray)
#	newEfficiency.SetMinimum(-0.15)
#	#newEfficiency.SetMinimum(0.7)
#	newEfficiency.SetMaximum(1.15)
#	newEfficiency.GetYaxis().SetLabelSize(0.05)
#	newEfficiency.GetXaxis().SetLabelSize(0.05)
#	newEfficiency.GetYaxis().SetTitleSize(0.06)
#	newEfficiency.GetYaxis().SetTitleOffset(0.8)
#	newEfficiency.GetXaxis().SetTitleOffset(0.8)
#	#newEfficiency.GetXaxis().SetLimits( 400, 1200 )
#	#newEfficiency.GetXaxis().SetLimits( 700, 1050 )
#	newEfficiency.GetXaxis().SetLimits( 0, 200 ) #xmin, xmax )
#	newEfficiency.Draw('AP')
#	labelAxis( name, newEfficiency, 'Pruned')
#	CMS_lumi.relPosX = 0.11
#	CMS_lumi.cmsTextSize = 0.7
#	CMS_lumi.extraOverCmsTextSize = 0.6
#	CMS_lumi.CMS_lumi(can1, 4, 0)
#
#	outputFileName = outputFileName.replace( 'jet1Pt', ''  )
#	can1.SaveAs( 'Plots/'+outputFileName )
#	del can1
#
#	return Efficiency


def diffplotTriggerEfficiency( name, diffTriggers, labelX, xmin, xmax, rebin, log ):
    """docstring for plot"""

    outputFileName = name+'_'+'_'.join([ q.split('_')[0] for q in diffTriggers])+'_'+(next(iter(diffTriggers)).split('_')[1])+'_TriggerEfficiency_'+args.version+'.'+args.extension
    print 'Processing.......', outputFileName

    diffEffDenom = {}
    diffEffPassing = {}
    diffEff = OrderedDict()
    for trigger in diffTriggers:
        print trigger+'/'+name+'Denom'
        diffEffDenom[ trigger ] = Samples[next(iter(Samples))][0].Get( trigger+'/'+name+'Denom' )
        diffEffDenom[ trigger ].Rebin(rebin)
        diffEffPassing[ trigger ] = Samples[next(iter(Samples))][0].Get( trigger+'/'+name+'Passing' )
        diffEffPassing[ trigger ].Rebin(rebin)
        diffEff[ trigger ] = TGraphAsymmErrors( diffEffPassing[trigger], diffEffDenom[trigger], 'cp'  )
        #Efficiency = TEfficiency( Passing, Denom )

    #binWidth = DenomOnly.GetBinWidth(1)

    legend=TLegend(0.50,0.25,0.90,0.55)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    can1 = TCanvas('c1', 'c1',  10, 10, 750, 500 )
    dummy=1
    for q in diffEff:
        diffEff[q].SetMarkerStyle(8)
        diffEff[q].SetMarkerColor(dummy)
        legend.AddEntry( diffEff[ q ], q.split('PF')[1], 'pl' )
        dummy+=1

    #diffEff[diffEff.iterkeys().next()].SetMinimum(0.8)
    diffEff[diffEff.iterkeys().next()].SetMaximum(1.15)
    diffEff[diffEff.iterkeys().next()].GetYaxis().SetLabelSize(0.05)
    diffEff[diffEff.iterkeys().next()].GetXaxis().SetLabelSize(0.05)
    diffEff[diffEff.iterkeys().next()].GetYaxis().SetTitleSize(0.06)
    diffEff[diffEff.iterkeys().next()].GetYaxis().SetTitleOffset(0.8)
    diffEff[diffEff.iterkeys().next()].GetXaxis().SetTitleOffset(0.8)
    diffEff[diffEff.iterkeys().next()].GetXaxis().SetLimits( xmin, xmax )
    diffEff[diffEff.iterkeys().next()].GetXaxis().SetTitle( labelX )
    diffEff[diffEff.iterkeys().next()].GetYaxis().SetTitle( 'Efficiency' )
    diffEff[diffEff.iterkeys().next()].Draw('AP')
    for q in diffEff: diffEff[q].Draw("P same")
    CMS_lumi.relPosX = 0.11
    CMS_lumi.cmsTextSize = 0.7
    CMS_lumi.extraOverCmsTextSize = 0.6
    CMS_lumi.CMS_lumi(can1, 4, 0)
    legend.Draw()

    can1.SaveAs( 'Plots/'+outputFileName )
    del can1

def plot2D( inFile, sample, name, titleXAxis, titleXAxis2, Xmin, Xmax, rebinx, Ymin, Ymax, rebiny, legX, legY ):
	"""docstring for plot"""

	outputFileName = name+args.cut+'_'+args.trigger+'_'+sample+'_'+'TriggerEfficiency'+args.version+'.'+args.extension
	print 'Processing.......', outputFileName
	h1 = inFile.Get( args.trigger+'TriggerEfficiency/'+name+args.cut )
	tmph1 = h1.Clone()

	### Rebinning
	nbinsx = h1.GetXaxis().GetNbins()
	nbinsy = h1.GetYaxis().GetNbins()
	xmin  = h1.GetXaxis().GetXmin()
	xmax  = h1.GetXaxis().GetXmax()
	ymin  = h1.GetYaxis().GetXmin()
	ymax  = h1.GetYaxis().GetXmax()
	nx = nbinsx/rebinx
	ny = nbinsy/rebiny
	h1.SetBins( nx, xmin, xmax, ny, ymin, ymax )

	for biny in range( 1, nbinsy):
		for binx in range(1, nbinsx):
			ibin1 = h1.GetBin(binx,biny)
			h1.SetBinContent( ibin1, 0 )

	for biny in range( 1, nbinsy):
		by = tmph1.GetYaxis().GetBinCenter( biny )
		iy = h1.GetYaxis().FindBin(by)
		for binx in range(1, nbinsx):
			bx = tmph1.GetXaxis().GetBinCenter(binx)
			ix  = h1.GetXaxis().FindBin(bx)
			bin = tmph1.GetBin(binx,biny)
			ibin= h1.GetBin(ix,iy)
			cu  = tmph1.GetBinContent(bin)
			h1.AddBinContent(ibin,cu)

	h1.GetXaxis().SetTitle( titleXAxis )
	h1.GetYaxis().SetTitleOffset( 1.0 )
	h1.GetYaxis().SetTitle( titleXAxis2 )

	if (Xmax or Ymax):
		h1.GetXaxis().SetRangeUser( Xmin, Xmax )
		h1.GetYaxis().SetRangeUser( Ymin, Ymax )

	tdrStyle.SetPadRightMargin(0.12)
	can = TCanvas('c1', 'c1',  750, 500 )
	can.SetLogz()
	#h1.SetMaximum(100)
	#h1.SetMinimum(0.1)
	h1.Draw('colz')

	CMS_lumi.extraText = "Simulation Preliminary"
	CMS_lumi.relPosX = 0.13
	CMS_lumi.CMS_lumi(can, 4, 0)

	can.SaveAs( 'Plots/'+outputFileName )
	del can

###################################################################
def simpleComparison( name, listComparison, labelX, xmin, xmax, rebin, log, PU ):

    for ifile in Samples:

        dictHistos = {}
        legend=TLegend( (0.60 if len(listComparison)>1 else 0.5 ),0.65,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)

        for i,k in enumerate(listComparison):
            dictHistos[ k+ifile ] = Samples[ifile][0].Get(k+'/'+name)
            dictHistos[ k+ifile ].SetLineColor( i+i+2 )
            dictHistos[ k+ifile ].SetLineWidth( 2 )
            dictHistos[ k+ifile ].Rebin( rebin  )
            try:
                #dictHistos[ k+ifile ].Scale( 1/dictHistos[ next(iter(dictHistos)) ].Integral() )
                if not name.startswith(('num','hltnum')): dictHistos[ k+ifile ].Scale( 1/dictHistos[ k+ifile ].Integral() )
                if name.startswith(('reco', 'numReco')): label = 'slimmedJetsPuppi '+k.split('HT')[1]
                elif name.startswith('gen'): label = 'slimmedGenJets '+k.split('HT')[1]
                else: label = k.split('PF')[1].replace('HT_', ' ')
                legend.AddEntry( dictHistos[ k+ifile ], label, 'l' )
            except ZeroDivisionError: dictHistos.pop(k+ifile)

        outputFileName = name+'_'+ifile+'_'.join([ x.split('PF')[1] for x in listComparison])+'_simpleComparison_'+args.version+'.'+args.extension
        print 'Processing.......', outputFileName

        can1 = TCanvas(outputFileName, outputFileName,  10, 10, 750, 500 )
        if log: can1.SetLogy()
        dictHistos[next(iter(dictHistos))].GetYaxis().SetLabelSize(0.05)
        dictHistos[next(iter(dictHistos))].GetXaxis().SetLabelSize(0.05)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitleSize(0.06)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitleOffset(0.9)
        dictHistos[next(iter(dictHistos))].GetXaxis().SetTitleOffset(0.8)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitle( ('Normalized' if not name.startswith(('num','hltnum')) else 'Events' )+' / '+str(dictHistos[next(iter(dictHistos))].GetBinWidth(1)) )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetTitle( labelX )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetRangeUser( xmin, xmax )
        dictHistos[next(iter(dictHistos))].SetMaximum( 1.2*max([ dictHistos[x].GetMaximum() for x in dictHistos  ]) )
        #dictHistos[next(iter(dictHistos))].SetMaximum( 10e4 )
        #dictHistos[next(iter(dictHistos))].SetMinimum( 1 )
        dictHistos[next(iter(dictHistos))].Draw('e')
        for ih in dictHistos:
            dictHistos[ih].Draw('histe same')
        CMS_lumi.lumi_13TeV = ifile+' ' #.split('_')[1]+' '
        CMS_lumi.relPosX = 0.11
        CMS_lumi.cmsTextSize = 0.7
        CMS_lumi.extraOverCmsTextSize = 0.7
        CMS_lumi.CMS_lumi(can1, 4, 0)
        legend.Draw()

        can1.SaveAs( 'Plots/'+outputFileName )
        del can1
        del dictHistos

###################################################################
def meanResponse( name, listComparison, labelY, labelX, xmin, xmax, rebin, log ):

    for ifile in Samples:
        dictHistos = {}
        outputFileName = name+'_'+ifile+'_'.join([ x.split('PF')[1] for x in listComparison])+'_meanResponse_'+args.version+'.'+args.extension
        print 'Processing.......', outputFileName

        legend=TLegend(0.50,0.65,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)

        for i,k in enumerate(listComparison):
            print k+'/'+name
            tmp = Samples[ifile][0].Get(k+'/'+name)
            dictHistos[ k+ifile ] = tmp.ProfileX()
            dictHistos[ k+ifile ].SetName(k+ifile)
            dictHistos[ k+ifile ].SetLineColor( i+1 )
            if not isinstance( rebin, list ): dictHistos[ k+ifile ].Rebin( rebin )
            else: dictHistos[ k+ifile ].Rebin(len(rebin)-1, dictHistos[ k+ifile ].GetName(), array('d', rebin) )
	    dictHistos[ k+ifile ].SetMarkerStyle(8)
	    dictHistos[ k+ifile ].SetMarkerColor( i+1 )
            legend.AddEntry( dictHistos[ k+ifile ], k.replace('Response', ''), 'pl' )


        can1 = TCanvas(outputFileName, outputFileName,  10, 10, 750, 500 )
        dictHistos[dictHistos.iterkeys().next()].SetMinimum(0)
        dictHistos[dictHistos.iterkeys().next()].SetMaximum(2)
        dictHistos[dictHistos.iterkeys().next()].GetYaxis().SetLabelSize(0.05)
        dictHistos[dictHistos.iterkeys().next()].GetXaxis().SetLabelSize(0.05)
        dictHistos[dictHistos.iterkeys().next()].GetYaxis().SetTitleSize(0.06)
        dictHistos[dictHistos.iterkeys().next()].GetYaxis().SetTitleOffset(0.8)
        dictHistos[dictHistos.iterkeys().next()].GetXaxis().SetTitleOffset(0.8)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitle( 'Response '+labelY )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetTitle( labelX )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetRangeUser( xmin, xmax )
        dictHistos[next(iter(dictHistos))].Draw('e')
        for ih in dictHistos:
            dictHistos[ih].Draw('e same')
        #labelAxis( name, diffEff[diffEff.iterkeys().next()], 'SoftDrop')
        CMS_lumi.lumi_13TeV = ifile #.split('_')[1]+' '
        CMS_lumi.relPosX = 0.11
        CMS_lumi.cmsTextSize = 0.7
        CMS_lumi.extraOverCmsTextSize = 0.6
        CMS_lumi.CMS_lumi(can1, 4, 0)
        legend.Draw()

        can1.SaveAs( 'Plots/'+outputFileName )
        del can1



###################################################################
def simpleComparisonTwoFiles( files, histName, listComparison, labelX, xmin, xmax, rebin, log, PU ):

    for k in listComparison:

        dictHistos = {}
        legend=TLegend( 0.5,0.75,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)
        labels = []
        for i, ifile in enumerate(files):
            if histName.startswith('resp'):
                tmp = ifile[0].Get(k+'/'+histName)
                dictHistos[ k+ifile[1] ] = tmp.ProfileX(k+ifile[1])
                dictHistos[ k+ifile[1] ].SetMarkerStyle(8)
                dictHistos[ k+ifile[1] ].SetMarkerColor( i+i+2 )
            else: dictHistos[ k+ifile[1] ] = ifile[0].Get(k+'/'+histName)
            dictHistos[ k+ifile[1] ].SetLineColor( i+i+2 )
            dictHistos[ k+ifile[1] ].SetLineWidth( 2 )
            dictHistos[ k+ifile[1] ].Rebin( rebin  )
            labels.append(ifile[1])
            try:
                if not histName.startswith('resp'): dictHistos[ k+ifile[1] ].Scale( 1/dictHistos[ k+ifile[1] ].Integral() )
                label = k.split('PF')[1].replace('HT_', ' ')+' '+ifile[1]
                legend.AddEntry( dictHistos[ k+ifile[1] ], label, 'pl' )
            except ZeroDivisionError: dictHistos.pop(k+ifile[1])

        outputFileName = histName+'_'+k+'_'.join(labels)+'_simpleComparison_'+args.version+'.'+args.extension
        print 'Processing.......', outputFileName

        can1 = TCanvas(outputFileName, outputFileName,  10, 10, 750, 500 )
        if log: can1.SetLogy()
        dictHistos[next(iter(dictHistos))].GetYaxis().SetLabelSize(0.05)
        dictHistos[next(iter(dictHistos))].GetXaxis().SetLabelSize(0.05)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitleSize(0.06)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitleOffset(0.9)
        dictHistos[next(iter(dictHistos))].GetXaxis().SetTitleOffset(0.8)
        dictHistos[next(iter(dictHistos))].GetYaxis().SetTitle( ('Response ' if histName.startswith('resp') else 'Normalized / ')+str(dictHistos[next(iter(dictHistos))].GetBinWidth(1)) )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetTitle( labelX )
        dictHistos[next(iter(dictHistos))].GetXaxis().SetRangeUser( xmin, xmax )
        dictHistos[next(iter(dictHistos))].SetMaximum( 1.2*max([ dictHistos[x].GetMaximum() for x in dictHistos  ]) )
        dictHistos[next(iter(dictHistos))].Draw('e')
        for ih in dictHistos:
            dictHistos[ih].Draw('e same')
        CMS_lumi.lumi_13TeV = ' ' #.split('_')[1]+' '
        CMS_lumi.relPosX = 0.11
        CMS_lumi.cmsTextSize = 0.7
        CMS_lumi.extraOverCmsTextSize = 0.7
        CMS_lumi.CMS_lumi(can1, 4, 0)
        legend.Draw()

        can1.SaveAs( 'Plots/'+outputFileName )
        del can1
        del dictHistos

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
    parser.add_argument('-e', '--extension', action='store', default='png', help='Extension of plots.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    if not os.path.exists('Plots/'): os.makedirs('Plots/')

    plotList = [

        [ '2D', 'jet1SDMassvsPt', 'Leading SD Jet Mass [GeV]', 'Leading Jet Pt [GeV]', 20, 200, 20, 400, 800, 50, 0.85, 0.2],

        [ 'simple', 'AK8Jet1Pt', 'Leading jet pt [GeV]', 100, 800, 1, True],
#        [ 'simple', 'AK8Jet2Pt', '2nd Leading jet pt [GeV]', 100, 800, 1, True],

        [ 'eff', 'recoLeadJetPt', 'Leading jet pt [GeV]', 100, 1000, 10, True],
        [ 'eff', 'recoHT', 'HT [GeV]', 500, 1500, 10, True],
            ]

    if 'all' in args.single: Plots = [ x[1:] for x in plotList if x[0] in args.proc ]
    else: Plots = [ y[1:] for y in plotList if ( ( y[0] in args.proc ) and ( y[1].startswith(args.single ) ))  ]


    bkgFiles = {}
    signalFiles = {}
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = '13 TeV'

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

    Samples = {}
    #Samples[ 'SingleMuon2017B' ] = [ TFile.Open('Rootfiles/triggerEfficiencies_histograms.root'), 0 ]
    Samples[ 'JetHT2017' ] = [ TFile.Open('/eos/home-a/algomez/tmpFiles/jetObservables/triggerEfficiencies/Plots/v05/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017ALL.root'), 0 ]
    #Samples[ 'JetHT2018' ] = [ TFile.Open('/eos/home-a/algomez/tmpFiles/jetObservables/triggerEfficiencies/Plots/v05/triggerEfficiencies_histograms_MiniAOD_JetHTRun2018ALL.root'), 0 ]


    processingSamples = {}
    if 'all' in args.dataset:
        for sam in Samples: processingSamples[ sam ] = Samples[ sam ]
    else:
        for sam in Samples:
            if sam.startswith( args.dataset ): processingSamples[ sam ] = Samples[ sam ]

    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    for i in Plots:
        for isam, samFile in processingSamples.iteritems():
            for q, it in enumerate(triggerList):
                plotTriggerEfficiency( samFile[0], isam, it[2], it[3], it[0], it[1],  i[0], it[4], it[5], i[1], i[2], i[3], i[4], i[5] )
#        if args.proc.startswith('simple'):
#            #cuts = [ 'Pt10', 'Pt20','Pt30','Pt40','Pt50', 'Pt10Eta5', 'Pt20Eta5', 'Pt30Eta5', 'Pt40Eta5', 'Pt50Eta5' ]
#            cuts = [ 'Pt10', 'Pt30','Pt50', 'Pt10Eta5', 'Pt30Eta5', 'Pt50Eta5' ]
#            for q in cuts:
#                PUlist = [''] if i[0].startswith(('reco', 'gen', 'numReco')) else [ '', 'CHS', 'PUPPI', 'SK' ]
#                simpleComparison( i[0], ['ResponseHLTPF'+pu+'HT_'+q for pu in PUlist  ], i[1], i[2], i[3], i[4], i[5], i[5] )
#        elif args.proc.startswith('meanRes'):
#            cuts = [ 'Pt10', 'Pt30','Pt50', 'Pt10Eta5', 'Pt30Eta5', 'Pt50Eta5' ]
#            for q in cuts:
#                meanResponse( i[0], ['ResponseHLTPF'+pu+'HT_'+q for pu in [ '', 'CHS', 'PUPPI', 'SK' ] ], i[1], i[2], i[3], i[4], i[5], i[6] )
#        elif args.proc.startswith('compa'):
#            cuts = [ 'Pt10','Pt30','Pt50', 'Pt10Eta5', 'Pt30Eta5', 'Pt50Eta5' ]
#            for q in cuts:
#                PUlist = [ '', 'CHS', 'PUPPI' ]
#                simpleComparisonTwoFiles(
#                        [ [TFile('reRunHLTwithAnalyzer_MC_PixelVertices.root'), 'PixelVertices'],
#                        [ TFile('reRunHLTwithAnalyzer_MC_VerticesPF.root'), 'VerticesPF'] ],
#                        i[0], ['ResponseHLTPF'+pu+'HT_'+q for pu in PUlist  ], i[1], i[2], i[3], i[4], i[5], i[5] )
#
#        elif args.proc.startswith('eff'):
#            cuts = [ 'Pt10', 'Pt30','Pt50', 'Pt10Eta5', 'Pt30Eta5', 'Pt50Eta5' ]
#            HT = [ 'HT850', 'HT950', 'HT1050' ]
#            PT = [ 'Pt350', 'Pt400', 'Pt450', 'Pt500' ]
#            for q in cuts:
#                for ht in HT:
#                    diffplotTriggerEfficiency( i[0], ['EffHLTPF'+pu+'HT_'+q+ht for pu in [ '', 'CHS', 'PUPPI', 'SK' ] ], i[1], i[2], i[3], i[4], i[5] )
#            for pt in PT:
#                diffplotTriggerEfficiency( i[0], ['EffHLTPF'+pu+'Jet_'+pt for pu in [ '', 'CHS', 'PUPPI', 'SK' ] ], i[1], i[2], i[3], i[4], i[5] )
#
#
#    ''' MAYBE THIS IS NEEDED LATER< DONT ERASE IT
#        elif '2D' in args.proc:
#                plot2DTriggerEfficiency( TFile.Open('Rootfiles/'+processingSamples[sam][0]),
#                                        sam,
#                                        BASEDTrigger,
#                                        i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
#
#        elif '2d' in args.proc:
#                plot2D( TFile.Open('Rootfiles/'+processingSamples[sam][0]),
#                                sam,
#                                i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], i[10] )
#    '''

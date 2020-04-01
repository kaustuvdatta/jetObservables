#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from ROOT import *
import ROOT
from datasets import *
from array import array
import numpy as np
from DrawHistogram import plotSimpleComparison
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
####gReset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)


def runTUnfold( dataFile, sigFiles, bkgFiles, variables, sel, sysUncert ):
    """docstring for createDataCards"""

    ### Getting input histos
    dataHistos = loadHistograms( dataFile, variables, sel, isMC=False )
    signalHistos = loadHistograms( sigFiles, variables, sel, sysUnc=sysUncert )
    bkgHistos = loadHistograms( bkgFiles, variables, sel, sysUnc=sysUncert ) if args.process.startswith('data') else {}

    for ivar in variables:

        print '|------> Unfolding '+ivar

        ######## Cross check: plotting data vs all MC Scaled
        print '|------> Cross check: plotting data vs all MC'
        allBkgHisto = dataHistos['data_recoJet'+ivar+'_nom'+sel].Clone()
        allBkgHisto.Reset()
        allBkgHistoNorm = allBkgHisto.Clone()
        for ibkg in bkgHistos:
            if ibkg.endswith('_recoJet'+ivar+'_nom'+sel): allBkgHisto.Add( bkgHistos[ibkg].Clone() )
            if ibkg.endswith('_recoJet'+ivar+'_nom'+sel+'_Normalized'): allBkgHistoNorm.Add( bkgHistos[ibkg].Clone() )
        allMCHisto = allBkgHisto.Clone()
        allMCHisto.Add( signalHistos[ next(iter(sigFiles))+'_recoJet'+ivar+'_nom'+sel ].Clone() )
        plotSimpleComparison( dataHistos['data_recoJet'+ivar+'_nom'+sel].Clone(), 'data', allMCHisto, 'allBkgs', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))+"_nom", rebinX=variables[ivar][0], version=sel+'_'+args.version  )
        allMCHistoNorm = allBkgHistoNorm.Clone()
        allMCHistoNorm.Add( signalHistos[ next(iter(sigFiles))+'_recoJet'+ivar+'_nom'+sel+"_Normalized" ].Clone() )
        plotSimpleComparison( dataHistos['data_recoJet'+ivar+'_nom'+sel+"_Normalized"].Clone(), 'data', allMCHistoNorm, 'allBkgs', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))+"_nom_Normalized", rebinX=variables[ivar][0], version=sel+'_'+args.version  )

        ######## Cross check: plotting response matrix
        tdrStyle.SetPadRightMargin(0.12)
        print '|------> Cross check: plotting response matrix for signal'
        can2D = TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[next(iter(sigFiles))+'_respJet'+ivar+'_nom'+sel].Draw("colz")
        can2D.SaveAs('Plots/'+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))+sel+'_responseMatrix'+args.version+'.png')
        can2DNorm = TCanvas(ivar+'can2DNorm', ivar+'can2DNorm', 750, 500 )
        signalHistos[next(iter(sigFiles))+'_respJet'+ivar+'_nom'+sel+'_Normalized'].Draw("colz")
        can2DNorm.SaveAs('Plots/'+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))+sel+'Normalized_responseMatrix'+args.version+'.png')


        ######## TUnfold part
        print '|------> TUnfolding starts:'

        ##### Defining options for TUnfold
        tunfolder = ROOT.TUnfoldDensity(
                                            #signalHistos[next(iter(sigFiles))+'_respJet'+ivar+'_nom'+sel+'_Normalized'], ### response matrix
                                            signalHistos[next(iter(sigFiles))+'_respJet'+ivar+'_nom'+sel], ### response matrix
                                            ROOT.TUnfold.kHistMapOutputHoriz,  #### kHistMapOutputVert if x->reco and y->gen, kHistMapOutputHoriz if x->gen and y->reco
                                            #ROOT.TUnfold.kRegModeCurvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                            #ROOT.TUnfold.kEConstraintNone,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. (Need to look into this!!)
                                            #ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                            )

        ##### Defining input (data recoJet )
        tunfolder.SetInput( dataHistos[ 'data_recoJet'+ivar+'_nom'+sel ].Clone() )
        #tunfolder.SetInput( dataHistos[ 'data_recoJet'+ivar+'_nom'+sel+'_Normalized' ].Clone() )

        ###### Removing bkgs from data
        for ibkg in bkgHistos:
            if ibkg.endswith('_recoJet'+ivar+'_nom'+sel):
                print '|--------> Removing this bkg: ', ibkg
                tunfolder.SubtractBackground( bkgHistos[ibkg], ibkg )

        ###### Adding SYS unc
        for sys in [ 'jesTotal', 'jer', 'pu' ]:
            for upDown in [ 'Up', 'Down' ]:
                tunfolder.AddSysError(
                                    signalHistos[next(iter(sigFiles))+'_respJet'+ivar+'_'+sys+upDown+sel],
                                    sys+upDown,
                                    ROOT.TUnfold.kHistMapOutputHoriz,
                                    ROOT.TUnfoldSys.kSysErrModeShift  #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                    )

        ###### Running the unfolding
        nScan=50
        tauMin=0.0
        tauMax=0.0
        iBest=0

        logTauX = ROOT.MakeNullPointer(ROOT.TSpline)
        logTauY = ROOT.MakeNullPointer(ROOT.TSpline)
        lCurve = ROOT.MakeNullPointer(ROOT.TGraph)
        ## this method scans the parameter tau and finds the kink in the L curve finally, the unfolding is done for the best choice of tau
        tunfolder.ScanLcurve(nScan,tauMin,tauMax,lCurve,logTauX,logTauY)
        #########################


        ###### Plot unfolding results
        canvasStack = ROOT.TCanvas('canvasStack', 'canvasStack', 750, 500)

        legend=TLegend(0.70,0.70,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)

        genJetHisto = signalHistos[ next(iter(sigFiles))+'_genJet'+ivar+sel ].Clone()
        #genJetHisto = signalHistos[ next(iter(sigFiles))+'_genJet'+ivar+sel+'_Normalized' ].Clone()
        #genJetHisto.Scale(1, 'width')  ### divide by bin width
        genJetHisto.SetLineWidth(2)
        genJetHisto.SetLineColor(2)
        legend.AddEntry( genJetHisto, 'Gen level', 'l' )

        unfoldHisto = tunfolder.GetOutput("unfolded")
        #unfoldHisto.Scale(1, 'width')  ### divide by bin width
        unfoldHisto.SetMarkerStyle(22)
        unfoldHisto.SetMarkerColor(4)
        legend.AddEntry( unfoldHisto, 'Unfolding', 'pe' )

        recoJetHisto = signalHistos[ next(iter(sigFiles))+'_recoJet'+ivar+'_nom'+sel+'_genBin' ].Clone()
        #recoJetHisto.Scale(1/recoJetHisto.Integral(), 'width')  ### divide by bin width
        recoJetHisto.SetLineWidth(2)
        recoJetHisto.SetLineColor(kBlue)
        legend.AddEntry( recoJetHisto, 'reco level', 'l' )

        hs = ROOT.THStack("hs", "hs")
        hs.Add( genJetHisto, "hist")
        hs.Add( recoJetHisto, "hist same")
        hs.Add( unfoldHisto, "e")
        hs.Draw("nostack")
        hs.GetXaxis().SetTitle( genJetHisto.GetXaxis().GetTitle() )
        hs.GetYaxis().SetTitle( 'dN/d#sigma' )
        hs.GetYaxis().SetTitleOffset(0.8)

        legend.Draw()
        CMS_lumi.extraText = "Simulation Preliminary"
        #CMS_lumi.lumi_13TeV = str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV, 2016"
        CMS_lumi.lumi_13TeV = "13 TeV, 2016"
        CMS_lumi.relPosX = 0.11
        CMS_lumi.CMS_lumi(canvasStack, 4, 0)
        canvasStack.SaveAs('Plots/'+ivar+sel+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_Tunfold_'+args.version+'.png')

##########################################################################
def loadHistograms( samples, variables, sel, sysUnc=[], isMC=True ):
    """docstring for loadHistograms"""

    SYSUNC = [ '_nom' ] + [ s+u for u in ['Up', 'Down'] for s in sysUnc ]

    allHistos = {}
    for var in variables:
        for isam in samples:
            tmpList = [ 'recoJet'+var+syst+sel for syst in SYSUNC ]
            if isMC: tmpList = tmpList + [ 'genJet'+var+sel ] + [ 'respJet'+var+syst+sel for syst in SYSUNC ]
            for ih in tmpList:
                allHistos[isam+'_'+ih] = samples[isam][0].Get( 'jetObservables/'+ih )
                if isMC:
                    tmpIsam = 'TT' if isam.startswith('data') else isam
                    MCScale = checkDict( tmpIsam, dictSamples )['XS'] * args.lumi / checkDict( tmpIsam, dictSamples )['2016']['nGenWeights']
                    allHistos[isam+'_'+ih].Scale( MCScale )

                if not ih.startswith('resp'):
                    if len(variables[var])==1:
                        if ih.startswith('reco'):
                            allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih].Clone()
                            allHistos[isam+'_'+ih+'_genBin'].Rebin( variables[var][0][1] )
                            allHistos[isam+'_'+ih].Rebin( variables[var][0][0] )
                        else: allHistos[isam+'_'+ih].Rebin( variables[var][0][1] )
                    else: allHistos[isam+'_'+ih] = allHistos[isam+'_'+ih].Rebin( len(variables[var])-1, allHistos[isam+'_'+ih].GetName()+"_Rebin", array( 'd', variables[var] ) )
                else:
                    if len(variables[var])==1: allHistos[isam+'_'+ih].Rebin2D( variables[var][0][1], variables[var][0][0] )
                    else:
                        #### fancy way to create variable binning TH2D
                        tmpHisto = TH2F( allHistos[isam+'_'+ih].GetName()+isam+"_Rebin", allHistos[isam+'_'+ih].GetName()+isam+"_Rebin", len(variables[var])-1, array( 'd', variables[var]), len(variables[var])-1, array( 'd', variables[var]) )

                        tmpArrayContent = np.zeros((len(variables[var]), len(variables[var])))
                        tmpArrayError = np.zeros((len(variables[var]), len(variables[var])))

                        for biny in range( 1, allHistos[isam+'_'+ih].GetNbinsY()+1 ):
                            by = allHistos[isam+'_'+ih].GetYaxis().GetBinCenter( biny )
                            for binx in range( 1, allHistos[isam+'_'+ih].GetNbinsX()+1 ):
                                bx = allHistos[isam+'_'+ih].GetXaxis().GetBinCenter(binx)
                                for iY in range( len(variables[var])-1 ):
                                    for iX in range( len(variables[var])-1 ):
                                        if (by<variables[var][iY+1] and by>variables[var][iY]) and (bx<variables[var][iX+1] and bx>variables[var][iX]):
                                            jbin = allHistos[isam+'_'+ih].GetBin(binx,biny)
                                            tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + allHistos[isam+'_'+ih].GetBinContent( jbin )
                                            tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + TMath.Power( allHistos[isam+'_'+ih].GetBinError( jbin ), 2 )

                        for biny in range( 1, tmpHisto.GetNbinsY()+1 ):
                            for binx in range( 1, tmpHisto.GetNbinsX()+1 ):
                                tmpHisto.SetBinContent( tmpHisto.GetBin(binx,biny), tmpArrayContent[binx-1][biny-1] )

                        allHistos[isam+'_'+ih] = tmpHisto

                allHistos[isam+'_'+ih+'_Normalized'] = allHistos[isam+'_'+ih].Clone()
                try: allHistos[isam+'_'+ih+'_Normalized'].Scale( 1/allHistos[isam+'_'+ih+'_Normalized'].Integral() )
                except ZeroDivisionError: continue

    return allHistos



###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="data", help="Process to unfold: data or MC." )
    parser.add_argument("-s", "--selection", action='store', dest="selection", default="_topSel", help="Selection to unfold: _dijetSel, _WSel, _topSel" )
    ##parser.add_argument("-r", "--runCombine", action='store_true', dest="runCombine", help="Run combine (true)" )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=35920., help='Luminosity, example: 1.' )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    inputFolder='Rootfiles/'+args.version
    dataFile = {}
    if args.process.startswith('MC'): dataFile['data'] = [ TFile( inputFolder+'/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root' ), 'Data', 'kBlack' ]
    else: dataFile['data'] = [ TFile( inputFolder+'/jetObservables_histograms_SingleMuonRun2016ALL.root' ), 'Data', 'kBlack' ]

    sigFiles = {}
    sigFiles['TTJets'] = [ TFile( inputFolder+'/jetObservables_histograms_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]
    #sigFiles['TT'] = [ TFile( inputFolder+'/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]

    bkgFiles = {}
    bkgFiles['ST_s-channel_4f_InclusiveDecays'] = [ TFile( inputFolder+'/jetObservables_histograms_ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t-channel_antitop'] = [ TFile( inputFolder+'/jetObservables_histograms_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t-channel_top'] = [ TFile( inputFolder+'/jetObservables_histograms_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_antitop'] = [ TFile( inputFolder+'/jetObservables_histograms_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_top'] = [ TFile( inputFolder+'/jetObservables_histograms_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['WJets'] = [ TFile( inputFolder+'/jetObservables_histograms_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root' ), 'WJets', 'kCyan' ]
    ##bkgFiles[] = [ '', TFile( inputFolder+'/jetObservables_histograms_'+ibkg+'.root' ), '', 'kMagenta' ]

    variables = {}
    variables[ 'Tau21' ] = [ (10,20) ]  ### (reco,gen)
    #variables[ 'Tau21' ] = [ 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1. ]

    sysUncert = [ '_jesTotal', '_jer', '_pu' ]

    runTUnfold( dataFile, sigFiles, bkgFiles, variables, args.selection, sysUncert )


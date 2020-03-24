#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from ROOT import *
import ROOT
from datasets import *
from array import array
import numpy as np
from DrawHistogram import plotSimpleComparison, plotUnfoldCombined
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
####gReset()
gROOT.SetBatch()
gROOT.ForceStyle()
tdrstyle.setTDRStyle()
gStyle.SetOptStat(0)


def createDatacards( dataFile, sigFiles, bkgFiles, variables ):
    """docstring for createDataCards"""

    ### Getting input histos
    dataHistos = loadHistograms( dataFile, variables, isMC=False )
    signalHistos = loadHistograms( sigFiles, variables )
    if args.process.startswith('data'): bkgHistos = loadHistograms( bkgFiles, variables )
    else:
        for ih in signalHistos: signalHistos[ih].Scale( dataHistos['Tau21_data_recoJet'].Integral()/signalHistos[ih].Integral() )

    for ivar in variables:

        print '|------> Unfolding '+ivar

        ######## Cross check: plotting data vs all MC
        print '|------> Cross check: plotting data vs all MC'
        allBkgHisto = dataHistos[ivar+'_data_recoJet'].Clone()
        allBkgHisto.Reset()
        if args.process.startswith('data'):
            for ibkg in bkgHistos:
                if ibkg.startswith(ivar) and ibkg.endswith('recoJet'):
                    allBkgHisto.Add( bkgHistos[ibkg].Clone() )
        allMCHisto = allBkgHisto.Clone()
        allMCHisto.Add( signalHistos[ ivar+'_'+next(iter(sigFiles))+'_recoJet' ].Clone() )
        plotSimpleComparison( dataHistos[ivar+'_data_recoJet'].Clone(), 'data', allMCHisto, 'allBkgs', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles)), rebinX=variables[ivar][0], version=args.version  )

        ######## Cross check: plotting response matrix
        tdrStyle.SetPadRightMargin(0.12)
        print '|------> Cross check: plotting response matrix for signal'
        can2D = TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[ivar+"_"+next(iter(sigFiles))+'_respJet'].Draw("colz")
        can2D.SaveAs('Plots/'+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))+'_responseMatrix.png')

        #print '----> Using RooUnfold'
        #R = ROOT.RooUnfoldResponse( signalHistos['recoJet'].Clone(), signalHistos['genJet'].Clone(), signalHistos['respJet'].Clone()  )
        #unfold = ROOT.RooUnfoldInvert( R, dataHisto )
        #hUnf = unfold.Hreco()

        ######## Creating datacard for combine
        datacardName = 'datacard_'+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))
        print '|------> Creating datacard: ', datacardName

        datacard = open( datacardName+'.txt', 'w')
        datacard.write("* imax\n")
        datacard.write("* jmax\n")
        datacard.write("* kmax\n")
        datacard.write("----------------\n")
        datacard.write("bin ")
        for ibin in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue  ### remove low stat bin
            datacard.write("Reco_"+str(ibin)+" ")
        datacard.write("\n")

        datacard.write("observation ")
        for ibin in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue
            datacard.write( str( dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin) )+" " )
        datacard.write("\n")
        datacard.write("----------------\n")

        datacard.write("bin ")
        for ibinReco in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue
            for ibinGen in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
                if signalHistos[ivar+'_'+next(iter(sigFiles))+'_respJet'].GetBinContent(ibinGen,ibinReco)<1: continue
                datacard.write("Reco_"+str(ibinReco)+' ')    ### one for ibinGen and one for ibinReco
            if args.process.startswith('data'): datacard.write("Reco_"+str(ibinReco)+' ')    ### one for all Bkgs
        datacard.write("\n")

        datacard.write("process ")
        for ibinReco in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue
            for ibinGen in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
                if signalHistos[ivar+'_'+next(iter(sigFiles))+'_respJet'].GetBinContent(ibinGen,ibinReco)<1: continue
                datacard.write("Gen_"+str(ibinGen)+' ')
            if args.process.startswith('data'): datacard.write("Bkg ")  ## bkg
        datacard.write("\n")

        datacard.write("process ")
        for ibinReco in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue
            for ibinGen in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
                if signalHistos[ivar+'_'+next(iter(sigFiles))+'_respJet'].GetBinContent(ibinGen,ibinReco)<1: continue
                datacard.write(" -"+str(ibinGen))   ## 0 -1, -2 --> for signal
            if args.process.startswith('data'): datacard.write(" 1 ")                   ## bkg >0 for bkg
        datacard.write("\n")

        datacard.write("rate ")
        for ibinReco in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
            if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibin)<1: continue
            for ibinGen in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1 ):
                if signalHistos[ivar+'_'+next(iter(sigFiles))+'_respJet'].GetBinContent(ibinGen,ibinReco)<1: continue
                datacard.write( str( round(signalHistos[ivar+'_'+next(iter(sigFiles))+'_respJet'].GetBinContent(ibinGen,ibinReco),2))+" " )
            if args.process.startswith('data'): datacard.write( str(round(allBkgHisto.GetBinContent(ibinReco),2))+' ' )
        datacard.write("\n")
        datacard.write("----------------\n")

        #### Creating the combine command
        po=' '.join(["--PO map='.*Gen_%d:r_bin%d[1,-20,20]'"%(ibinGen,ibinGen) for ibinGen in range(1, dataHistos[ivar+'_data_recoJet'].GetNbinsX()+1) if dataHistos[ivar+'_data_recoJet'].GetBinContent(ibinGen)>0 ])
        cmdText2Workspace = "text2workspace.py "+('--X-allow-no-background ' if args.process.startswith('MC') else ' ' )+datacardName+".txt -o "+datacardName+".root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+po
        datacard.write("### RUN WITH COMMANDS: ####\n")
        datacard.write("# "+cmdText2Workspace+"\n")
        cmdCombine = 'combine -M MultiDimFit --algo singles -d '+datacardName+'.root -t 0 -n '+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles))
        datacard.write("# Run with command: "+cmdCombine+" \n")
        datacard.write("############################\n")
        print '|------> To run combine: \n', cmdText2Workspace + '\n' + cmdCombine


##########################################################################
def loadHistograms( samples, variables, isMC=True ):
    """docstring for loadHistograms"""

    allHistos = {}
    for var in variables:
        for isam in samples:
            for ih in [ 'recoJet', 'genJet', 'respJet' ]:
                allHistos[var+'_'+isam+'_'+ih] = samples[isam][1].Get( 'jetObservables/'+ih+var )
                if isMC:
                    MCScale = checkDict( samples[isam][0], dictSamples )['XS'] * args.lumi / checkDict( samples[isam][0], dictSamples )['2016']['nevents']
                    allHistos[var+'_'+isam+'_'+ih].Scale( MCScale )

                if not ih.startswith('resp'):
                    if len(variables[var])==1: allHistos[var+'_'+isam+'_'+ih].Rebin( variables[var][0] )
                    else: allHistos[var+'_'+isam+'_'+ih] = allHistos[var+'_'+isam+'_'+ih].Rebin( len(variables[var])-1, allHistos[var+'_'+isam+'_'+ih].GetName()+"_Rebin", array( 'd', variables[var] ) )
                else:
                    if len(variables[var])==1: allHistos[var+'_'+isam+'_'+ih].Rebin2D( variables[var][0], variables[var][0] )
                    else:
                        #### fancy way to create variable binning TH2D
                        tmpHisto = TH2F( allHistos[var+'_'+isam+'_'+ih].GetName()+isam+"_Rebin", allHistos[var+'_'+isam+'_'+ih].GetName()+isam+"_Rebin", len(variables[var])-1, array( 'd', variables[var]), len(variables[var])-1, array( 'd', variables[var]) )

                        tmpArrayContent = np.zeros((len(variables[var]), len(variables[var])))
                        tmpArrayError = np.zeros((len(variables[var]), len(variables[var])))

                        for biny in range( 1, allHistos[var+'_'+isam+'_'+ih].GetNbinsY()+1 ):
                            by = allHistos[var+'_'+isam+'_'+ih].GetYaxis().GetBinCenter( biny )
                            for binx in range( 1, allHistos[var+'_'+isam+'_'+ih].GetNbinsX()+1 ):
                                bx = allHistos[var+'_'+isam+'_'+ih].GetXaxis().GetBinCenter(binx)
                                for iY in range( len(variables[var])-1 ):
                                    for iX in range( len(variables[var])-1 ):
                                        if (by<variables[var][iY+1] and by>variables[var][iY]) and (bx<variables[var][iX+1] and bx>variables[var][iX]):
                                            jbin = allHistos[var+'_'+isam+'_'+ih].GetBin(binx,biny)
                                            tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + allHistos[var+'_'+isam+'_'+ih].GetBinContent( jbin )
                                            tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + TMath.Power( allHistos[var+'_'+isam+'_'+ih].GetBinError( jbin ), 2 )

                        for biny in range( 1, tmpHisto.GetNbinsY()+1 ):
                            for binx in range( 1, tmpHisto.GetNbinsX()+1 ):
                                tmpHisto.SetBinContent( tmpHisto.GetBin(binx,biny), tmpArrayContent[binx-1][biny-1] )

                        allHistos[var+'_'+isam+'_'+ih] = tmpHisto

    return allHistos



###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="data", help="Process to unfold: data or MC." )
    ##parser.add_argument("-r", "--runCombine", action='store_true', dest="runCombine", help="Run combine (true)" )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=35920., help='Luminosity, example: 1.' )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dataFile = {}
    if args.process.startswith('MC'): dataFile['data'] = [ 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root' ), 'Data', 'kBlack' ]
    else: dataFile['data'] = [ 'data', TFile( 'Rootfiles/jetObservables_histograms_SingleMuon2016ALL.root' ), 'Data', 'kBlack' ]

    sigFiles = {}
    sigFiles['ttbarMadgraph'] = [ 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]
    #sigFiles['ttbarPowheg'] = [ 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]

    bkgFiles = {}
    bkgFiles['ST_s'] = [ 'ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8', TFile( 'Rootfiles/jetObservables_histograms_ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t_antitop'] = [ 'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1', TFile( 'Rootfiles/jetObservables_histograms_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t_top'] = [ 'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1', TFile( 'Rootfiles/jetObservables_histograms_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_antitop'] = [ 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4', TFile( 'Rootfiles/jetObservables_histograms_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_top'] = [ 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4', TFile( 'Rootfiles/jetObservables_histograms_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['WJets'] = [ 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8', TFile( 'Rootfiles/jetObservables_histograms_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root' ), 'WJets', 'kCyan' ]
    ##bkgFiles[] = [ '', TFile( 'Rootfiles/jetObservables_histograms_'+ibkg+'.root' ), '', 'kMagenta' ]

    variables = {}
    variables[ 'Tau21' ] = [ 4 ]
    #variables[ 'Tau21' ] = [ 0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1. ]

    if args.process.endswith('plot'):
        for ivar in variables:
            plotUnfoldCombined(sigFiles['ttbarMadgraph'], ivar, args.lumi, variables[ivar], ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+next(iter(sigFiles)), 'Madgraph' )
    else: createDatacards( dataFile, sigFiles, bkgFiles, variables )


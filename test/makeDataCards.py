#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from ROOT import *
#import CMS_lumi as CMS_lumi
#18 import tdrstyle as tdrstyle
####gReset()
gROOT.SetBatch()
#gROOT.ForceStyle()
#tdrstyle.setTDRStyle()
#gStyle.SetOptStat(0)


def createDatacards( dataFile, MCFile, histo ):
    """docstring for convertRootToNumpy"""

    ### Getting input histos
    dataHisto = dataFile.Get( 'jetObservables/recoJet'+histo )
    dataHisto.Rebin(2)

    MCHisto = {}
    for ih in [ 'recoJet', 'genJet', 'respJet' ]:
        MCHisto[ih] = MCFile.Get( 'jetObservables/'+ih+histo )
        if not ih.startswith('resp'):
            #MCHisto[ih].Scale( 1./MCHisto[ih].Integral() )
            MCHisto[ih].Rebin(2)
        else: MCHisto[ih].Rebin2D(2,2)

    MCHisto['recoJetScaled'] = MCHisto['recoJet'].Clone()
    MCHisto['recoJetScaled'].Scale( dataHisto.Integral() )

    if args.process.startswith('datacard'):
        datacard = open('datacard.txt', 'w')
        datacard.write("* imax\n")
        datacard.write("* jmax\n")
        datacard.write("* kmax\n")
        datacard.write("----------------\n")
        datacard.write("bin ")
        for ireco in range(1, dataHisto.GetNbinsX()+1 ):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            datacard.write("Reco_"+str(ireco)+" ")
        datacard.write("\n")

        datacard.write("observation ")
        for ireco in range(1, dataHisto.GetNbinsX()+1 ):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            datacard.write( str( dataHisto.GetBinContent(ireco))+" " )
        datacard.write("\n")
        datacard.write("----------------\n")

        cleanup=True
        datacard.write("bin ")
        for ireco in range(1, MCHisto['recoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['genJet'].GetNbinsX()+1):
                # remove un-necessary processes
                if cleanup and MCHisto['respJet'].GetBinContent(ireco,igen)<0.01: continue
                datacard.write("Reco_%d "%ireco) ##sig igen, in reco ireco
            #datacard.write("Reco_%d "%ireco)## bkg
        datacard.write("\n")

        datacard.write("process ")
        for ireco in range(1, MCHisto['recoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['genJet'].GetNbinsX()+1):
                if cleanup and MCHisto['respJet'].GetBinContent(ireco,igen)<0.01: continue
                datacard.write("Gen_%d "%igen) ##sig igen, in reco ireco
            #datacard.write("Bkg ")## bkg
        datacard.write("\n")

        datacard.write("process ")
        for ireco in range(1, MCHisto['recoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['genJet'].GetNbinsX()+1):
                if cleanup and MCHisto['respJet'].GetBinContent(ireco,igen)<0.01: continue
                datacard.write("%d "%(-igen)) ## 0 -1, -2 --> for signal
            #datacard.write("1 ")## bkg >0 for bkg
        datacard.write("\n")

        datacard.write("rate ")
        for ireco in range(1, MCHisto['recoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<0.01: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['genJet'].GetNbinsX()+1):
                if cleanup and MCHisto['respJet'].GetBinContent(ireco,igen)<0.01: continue
                datacard.write("%.2f "% ( MCHisto['respJet'].GetBinContent(ireco,igen)))
            #datacard.write("1 ")
            #datacard.write("%.2f "%(MCHisto['recoJetScaled'].GetBinContent(ireco) if MCHisto['recoJetScaled'].GetBinContent(ireco)>0 else 1))##
        datacard.write("\n")
        datacard.write("----------------\n")

        po=' '.join(["--PO map='.*/Gen_%d:r_bin%d[1,-20,20]'"%(igen,igen) for igen in range(1, dataHisto.GetNbinsX()+1) if dataHisto.GetBinContent(igen)>0 ])
        datacard.write("### RUN WITH COMMANDS: ####\n")
        datacard.write("# text2workspace.py --X-allow-no-background datacard.txt -o datacard.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+po+"\n")
        datacard.write("# Run with command: combine -M MultiDimFit --algo singles -d datacard.root -t 0 \n")
        datacard.write("############################\n")

    ###test
    if args.process.startswith('plot'):

        combineFile = TFile('higgsCombineTest.MultiDimFit.mH120.root')
        combineTree = combineFile.Get('limit')
        MCHisto['genJet'].Scale( dataHisto.Integral() )
        test = MCHisto['genJet'].Clone()
        test.Reset()
        tmp = [ 1, +1.128, -0.185, +0.740, -0.015, +1.980, +1.035, -0.030, +3.746, -0.826, +2.570, +0.590, +2.297, +0.607, +1.991, -1.773, -1.354, -1.814, 1, 1 ]
        for ireco in range(1, test.GetNbinsX()):
            if ireco>1:
                tmpHisto = TH1F('tmp'+str(ireco), 'tmp'+str(ireco), 80, -20, 20)
                combineTree.Draw("r_bin"+str(ireco)+'>>tmp'+str(ireco))
                res = tmpHisto.GetMean()
            else: res = 1
            test.SetBinContent( ireco, MCHisto['genJet'].GetBinContent(ireco)*res )
            #test.SetBinError( MCHisto['genJet'].GetBinError(ireco)*tmp[ireco] )

        test.SetMarkerStyle(8)

        can = TCanvas('can', 'can', 750, 800 )
        MCHisto['genJet'].Draw("hist")
        test.Draw("ep same")
        MCHisto['genJet'].SetMaximum( 1.2*max(test.GetMaximum(), MCHisto['genJet'].GetMaximum()) )
        can.SaveAs('test.png')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="datacard", help="Process: datacard or plot." )
    parser.add_argument("-f", "--inputFolder", action='store', dest="inputFolder", default="", help="Path of the folder where the root files are located (wihtout the root://blah//)." )
    parser.add_argument("-o", "--outputFolder", action='store', dest="outputFolder", default="test", help="Name of output folder" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    MCFile = TFile( 'histograms_TTJets.root' )
    dataFile = TFile( 'histograms_TT.root' )

    histo = 'Tau21'

    createDatacards( dataFile, MCFile, histo )


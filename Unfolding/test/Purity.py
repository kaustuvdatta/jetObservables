#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import ROOT
from array import array
from collections import OrderedDict
from variables import nSubVariables, nSubVariables_WSel, nSubVariables_topSel
from runTUnfold import loadHistograms
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
sys.path.insert(0,'../../Skimmer/test/')
from datasets import checkDict, dictSamples

####gReset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)

colorPallete = [ 0, 2, 4, 8, 12, 28, 30 ]

def runPurity( name, bkgFiles, variables, sel ):
    """docstring for createDataCards"""

    numBins = 0
    numBinsList = [0]
    dictHistos = OrderedDict()
    for ivar in variables:
        if not ivar.startswith( name ): continue
        ### Getting input histos
        signalHistos = loadHistograms( bkgFiles, ivar, sel, sysUnc=[], lumi=1, year=args.year, process='MC', variables=variables )
        outputDir=args.outputFolder+'Plots/'+args.selection.split('_')[1]+'/Purity/'+args.year+'/'+ivar+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)

        legend=ROOT.TLegend(0.15,0.80,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)
        legend.SetNColumns( 4 )

        dictHistos[ 'acceptanceGraph_'+ivar ] = signalHistos[signalLabel+'_accepgen'+ivar+args.selection].Clone()
        dictHistos[ 'acceptanceGraph_'+ivar ].SetLineWidth(2)
        dictHistos[ 'acceptanceGraph_'+ivar ].SetLineColor(colorPallete[1])
        dictHistos[ 'acceptanceGraph_'+ivar ].Divide(  signalHistos[signalLabel+'_gen'+ivar+args.selection] )
        legend.AddEntry( dictHistos[ 'acceptanceGraph_'+ivar ], 'Acceptance', 'l' )

        dictHistos[ 'fakeGraph_'+ivar ] = signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+args.selection+'_genBin'].Clone()
        dictHistos[ 'fakeGraph_'+ivar ].SetLineWidth(2)
        dictHistos[ 'fakeGraph_'+ivar ].SetLineColor(colorPallete[2])
        dictHistos[ 'fakeGraph_'+ivar ].Divide(  signalHistos[signalLabel+'_reco'+ivar+'_nom'+args.selection+'_genBin'] )
        legend.AddEntry( dictHistos[ 'fakeGraph_'+ivar ], 'Fakerate', 'l' )

        #### create diagonal
        tmpHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].Clone()
        tmpHisto.Rebin2D( 1, 2 )
        dictHistos[ 'purityGraph_'+ivar ] = signalHistos[signalLabel+'_reco'+ivar+'_nom'+args.selection+'_genBin'].Clone()
        dictHistos[ 'purityGraph_'+ivar ].Reset()
        dictHistos[ 'stabilityGraph_'+ivar ] = dictHistos[ 'purityGraph_'+ivar ].Clone()
        for ibin in range( 1, tmpHisto.GetNbinsX()+1 ):
            diagReco = tmpHisto.GetBinContent( tmpHisto.GetBin(ibin, ibin) )
            dictHistos[ 'purityGraph_'+ivar ].SetBinContent( ibin, diagReco )
            diagGen = tmpHisto.GetBinContent( tmpHisto.GetBin(ibin, ibin) )
            dictHistos[ 'stabilityGraph_'+ivar ].SetBinContent( ibin, diagGen )

        dictHistos[ 'purityGraph_'+ivar ].SetLineWidth(2)
        dictHistos[ 'purityGraph_'+ivar ].SetLineColor(colorPallete[3])
        dictHistos[ 'purityGraph_'+ivar ].Divide(  signalHistos[signalLabel+'_truereco'+ivar+'_nom'+args.selection+'_genBin'] )
        legend.AddEntry( dictHistos[ 'purityGraph_'+ivar ], 'Purity', 'l' )

        dictHistos[ 'stabilityGraph_'+ivar ].SetLineWidth(2)
        dictHistos[ 'stabilityGraph_'+ivar ].SetLineColor(colorPallete[4])
        dictHistos[ 'stabilityGraph_'+ivar ].Divide(  signalHistos[signalLabel+'_accepgen'+ivar+args.selection] )
        legend.AddEntry( dictHistos[ 'stabilityGraph_'+ivar ], 'Stability', 'l' )

        ##### removing tau21 and tau32 from total plots
        if not ivar.endswith(('21', '32')):
            numBins = numBins + dictHistos[ 'purityGraph_'+ivar ].GetNbinsX()
            numBinsList.append(numBins)

        multigraph = ROOT.THStack()
        for i in dictHistos:
            if i.endswith(ivar): multigraph.Add( dictHistos[i]  )

        ###### Plot unfolding results
        ROOT.gStyle.SetPadRightMargin(0.05)
        canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)

        multigraph.Draw("hist nostack")
        multigraph.GetXaxis().SetTitle( variables[ivar]['label'] )
        multigraph.GetYaxis().SetTitle( 'Percentage' )
        multigraph.SetMaximum( 1.1 )
        multigraph.GetYaxis().SetTitleOffset(0.8)

        legend.Draw()
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2017+2018' if args.year.startswith('all') else args.year )
        CMS_lumi.relPosX = 0.11
        CMS_lumi.CMS_lumi(canvas, 4, 0)

        textBox=ROOT.TLatex()
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        textBox.SetTextFont(62) ### 62 is bold, 42 is normal
        #textBox.DrawLatex(0.65, 0.75, sel.split('Sel')[0].split('_')[1]+' Selection' )

        canvas.SaveAs(outputDir+ivar+'_'+signalLabel+sel+'_Purity_'+args.version+'.'+args.ext)
        if args.ext.startswith('pdf'):
            canvas.SaveAs(outputDir+ivar+'_'+signalLabel+sel+'_Purity_'+args.version+'.png')


        ########## Plot response matrix
        print ('|------> Cross check: plotting response matrix for signal')
        ROOT.gStyle.SetPadRightMargin(0.15)
        #ROOT.gStyle.SetPalette(ROOT.kGistEarth)
        #ROOT.TColor.InvertPalette()
        can2D = ROOT.TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetXaxis().SetTitle('Accepted Gen '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitle('True Reco '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitleOffset( 0.8 )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].Draw("colz")
        CMS_lumi.relPosX = 0.12
        CMS_lumi.CMS_lumi(can2D, 4, 0)
        can2D.SaveAs(outputDir+ivar+'_'+signalLabel+sel+'_responseMatrix'+args.version+'.'+args.ext)
        if args.ext.startswith('pdf'):
            can2D.SaveAs(outputDir+ivar+'_'+signalLabel+sel+'_responseMatrix'+args.version+'.png')

    if not args.only:

        outputFileName = 'reco'+name+'_purity_'+signalLabel+'_combinePlots_'+args.version+'.'+args.ext
        print('Processing.......', outputFileName)

        legend=ROOT.TLegend(0.20,0.80,0.60,0.90)
        legend.SetNColumns(2)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.06)

        dictHistos['purity'] = ROOT.TH1F('purity', 'purity', numBinsList[-1], 0, numBinsList[-1])
        legend.AddEntry( dictHistos[ 'purity' ], 'Purity', 'lep' )
        dictHistos['stability'] = ROOT.TH1F('stability', 'stability', numBinsList[-1], 0, numBinsList[-1])
        legend.AddEntry( dictHistos[ 'stability' ], 'Stability', 'lep' )

        tmpNbin = 0
        Xlabels = []
        for ivar in variables:
            if not ivar.startswith( name ): continue
            if ivar.startswith( ('Jet', 'sdJet' ) ) and not ivar.endswith(('21', '32')):
                Xlabels.append( '#'+nSubVariables[ivar]['label'].split('#')[1] )
                for ibin in range(1, dictHistos['purityGraph_'+ivar].GetNbinsX()+1):
                    tmpNbin = tmpNbin+1
                    dictHistos['purity'].SetBinContent( tmpNbin, dictHistos['purityGraph_'+ivar].GetBinContent(ibin) )
                    dictHistos['stability'].SetBinContent( tmpNbin, dictHistos['stabilityGraph_'+ivar].GetBinContent(ibin) )

        ROOT.gStyle.SetPadRightMargin(0.05)
        #ROOT.gStyle.SetPadLeftMargin(0.08)
        canvasoutputFileName = ROOT.TCanvas('c1purity', 'c1purity', 1400, 500 )
        ROOT.gStyle.SetPadTickX(0)
        #pad1.SetLogy()
        dictHistos['purity'].GetXaxis().SetNdivisions(100)
        dictHistos['purity'].GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#tau_{X}}' )
        dictHistos['purity'].GetYaxis().SetTitleSize( 0.06 )
        dictHistos['purity'].GetYaxis().SetTitleOffset( 0.8 )
        #dictHistos['purity'].SetMaximum( dictHistos['stability'].GetMaximum()*1.2 )
        dictHistos['purity'].SetMaximum( 0.99 )
        dictHistos['purity'].SetMinimum( 0.26 )
        dictHistos['purity'].SetLineColor( ROOT.kBlack )
        dictHistos['purity'].SetLineWidth( 2 )

        dictHistos['stability'].SetLineColor( ROOT.kMagenta )
        dictHistos['stability'].SetLineWidth( 2 )
        dictHistos['purity'].Draw('hist')
        dictHistos[ 'stability' ].Draw('hist same')

        ### division lines
        lines = {}
        for i in numBinsList[1:-1]:
            lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
            lines[i].SetLineColor(ROOT.kGray)
            lines[i].Draw('same')

        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.relPosX = 0.07
        CMS_lumi.CMS_lumi( canvasoutputFileName, 4, 0)
        legend.Draw()

        textBox=ROOT.TLatex()
        #textBox.SetNDC()
        textBox.SetTextSize(0.06)
        textBox.SetTextAlign(12)
        textBoxList = {}
        for i in range(1, len(numBinsList)):
            textBoxList[i] = textBox.Clone()
            textBoxList[i].DrawLatex(numBinsList[i-1]+(numBinsList[i]-numBinsList[i-1])/2., 0.22, Xlabels[i-1] )

        outputDir=args.outputFolder+'Plots/'+args.selection.split('_')[1]+'/Purity/'+args.year+'/'
        canvasoutputFileName.SaveAs( outputDir+'/'+outputFileName )
        if args.ext.startswith('pdf'):
            canvasoutputFileName.SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
        ROOT.gStyle.SetPadTickX(1)

###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--selection", action='store', dest="selection", default="_dijetSel", help="Selection to unfold: _dijetSel, _WSel, _topSel" )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )
    parser.add_argument("--outputFolder", action='store', dest="outputFolder", default="", help="Output folder" )
    parser.add_argument("--inputFolder", action='store', dest="inputFolder", default="", help="Input folder" )
    parser.add_argument("--main", action='store', dest="main", choices=['Ptbin', 'HTbin', 'herwig'], default='Ptbin', help="For dijet sel, main signal QCD" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
    #args.inputFolder = os.environ['CMSSW_BASE']+'/src/jetObservables/Unfolding/test/Rootfiles/'

    if args.selection.startswith('_dijet'):
        if args.main.startswith('Ptbin'):
            signalLabelBegin = 'QCD_Pt_'
            signalLabel = 'QCD_Pt_3200toInf'
        elif args.main.startswith('HTbin'):
            signalLabelBegin = 'QCD_HT'
            signalLabel = 'QCD_HT2000toInf'
        elif args.main.startswith('herwig'):
            signalLabelBegin = 'QCD_Pt-'
            signalLabel = 'QCD_Pt-150to3000'
    else:
        signalLabel = 'TTToSemiLeptonic'
        sigPlotLabel = 'Powheg+Pythia8'
        signalLabelBegin = 'TTToSemiLeptonic'
        
        altSignalLabelBegin = 'TTJets'
        altSigPlotLabel = 'aMC@NLO-FXFX+Pythia8'
        altSignalLabel = 'TTJets'

    bkgFiles = {}
    for iy in ( ['2017', '2018'] if args.year.startswith('all') else [ args.year ] ):
        for isam in dictSamples:
            if not checkDict( isam, dictSamples )[iy]['skimmerHisto'].endswith('root'): continue
            if args.selection.startswith('_dijet') and ( 'MuEnriched' in isam ): continue
            if isam.startswith(signalLabelBegin):
                bkgFiles[isam.split('_Tune')[0]] = [
                                ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[iy]['skimmerHisto'] ),
                                checkDict( isam, dictSamples )
                            ]
    if args.selection.startswith('_WSel'): nSubVariables = nSubVariables_WSel
    elif args.selection.startswith('_topSel'): nSubVariables = nSubVariables_topSel
    print (bkgFiles.keys())
    #### define variables
    if args.only:
        filterVariables = { k:v for (k,v) in nSubVariables.items() if k.endswith(args.only)  }
        if len(filterVariables)>0 : variables = filterVariables
        else:
            print('|------> Variable not found. Have a nice day')
            sys.exit(0)
    else: variables = nSubVariables
    #if args.selection.startswith('_dijet'): variables = { k:v for (k,v) in variables.items() if k.startswith(('Jet1', 'Jet2')) }
    #else: variables = { k:v for (k,v) in variables.items() if k.startswith('Jet_') }

    plotList = [ 'Jet1', 'Jet2', 'sdJet1', 'sdJet2' ] if args.selection.startswith('_dijet') else [ 'Jet' ]

    for i in plotList:
        runPurity( i, bkgFiles, variables, args.selection )


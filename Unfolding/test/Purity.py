#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import ROOT
from datasets import checkDict, dictSamples
from array import array
import numpy as np
from DrawHistogram import plotSimpleComparison
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle

####gReset()
ROOT.gROOT.SetBatch()
ROOT.gROOT.ForceStyle()
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)

colorPallete = [ 0, 2, 4, 8, 12, 28, 30 ]

def runPurity( bkgFiles, variables, sel ):
    """docstring for createDataCards"""

    ### Getting input histos
    signalHistos = loadHistograms( bkgFiles, variables, sel )

    for ivar in variables:

        outputDir='Plots/Unfold/'+args.year+'/'+ivar+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)

        dictGraph = {}
        legend=ROOT.TLegend(0.15,0.80,0.90,0.90)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.03)
        legend.SetNColumns( len(variables[ivar][0]) )
        multigraph = ROOT.TMultiGraph()

        print ('|------> Computing purity '+ivar )
        dummy=1
        for ihisto in signalHistos:
            if ivar in ihisto:
                matrix = np.zeros((signalHistos[ihisto].GetNbinsX(),signalHistos[ihisto].GetNbinsY()))
                binCenters = []
                binWidths = []
                for biny in range( 1, signalHistos[ihisto].GetNbinsY()+1 ):
                    binCenters.append( signalHistos[ihisto].GetYaxis().GetBinCenter(biny) )
                    binWidths.append( signalHistos[ihisto].GetYaxis().GetBinWidth(biny)/2. )
                    for binx in range( 1, signalHistos[ihisto].GetNbinsX()+1 ):
                        jbin = signalHistos[ihisto].GetBin(binx,biny)
                        binCont = signalHistos[ihisto].GetBinContent(jbin)
                        matrix[binx-1][biny-1] = binCont
                recoBins = np.sum( matrix, axis=0 )
                genBins = np.sum( matrix, axis=1 )

                purity = np.nan_to_num( genBins/(recoBins+genBins) )         #### it is flipped to make it percentage
                stability = np.nan_to_num( recoBins/(recoBins+genBins) )     #### it is flipped to make it percentage


                dictGraph[ 'purGraph_'+ihisto ] = ROOT.TGraphErrors( signalHistos[ihisto].GetNbinsX(), array( 'd', binCenters ), array( 'd', purity ), array( 'd', binWidths ), array( 'd', [0]*len(binWidths) ) )
                dictGraph[ 'purGraph_'+ihisto ].SetLineWidth(2)
                dictGraph[ 'purGraph_'+ihisto ].SetLineColor(colorPallete[dummy])
                legend.AddEntry( dictGraph[ 'purGraph_'+ihisto ], 'Purity Bin '+str(signalHistos[ihisto].GetXaxis().GetBinWidth(1)), 'l' )

                dictGraph[ 'staGraph_'+ihisto ] = ROOT.TGraphErrors( signalHistos[ihisto].GetNbinsX(), array( 'd', binCenters ), array( 'd', stability ), array( 'd', binWidths ), array( 'd', [0]*len(binWidths) ) )
                dictGraph[ 'staGraph_'+ihisto ].SetLineWidth(2)
                dictGraph[ 'staGraph_'+ihisto ].SetLineStyle(2)
                dictGraph[ 'staGraph_'+ihisto ].SetLineColor(colorPallete[dummy])
                legend.AddEntry( dictGraph[ 'staGraph_'+ihisto ], 'Stability Bin '+str(signalHistos[ihisto].GetXaxis().GetBinWidth(1)), 'l' )

                dummy=dummy+1

        for i in dictGraph: multigraph.Add( dictGraph[i]  )

        ###### Plot unfolding results
        canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)

        multigraph.Draw("AP")
        multigraph.GetXaxis().SetTitle( variables[ivar][1] )
        multigraph.GetYaxis().SetTitle( 'Percentage' )
        multigraph.GetYaxis().SetRangeUser( 0, 1.3 )
        multigraph.GetYaxis().SetTitleOffset(0.8)

        legend.Draw()
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+args.year
        CMS_lumi.relPosX = 0.11
        CMS_lumi.CMS_lumi(canvas, 4, 0)

        textBox=ROOT.TLatex()
        textBox.SetNDC()
        textBox.SetTextSize(0.04)
        textBox.SetTextFont(62) ### 62 is bold, 42 is normal
        textBox.DrawLatex(0.65, 0.75, sel.split('Sel')[0].split('_')[1]+' Selection' )

        canvas.SaveAs(outputDir+ivar+sel+'_Purity_'+args.version+'.'+args.ext)

##########################################################################
def loadHistograms( samples, variables, sel ):
    """docstring for loadHistograms"""

    allHistos = {}
    for var in variables:
        for isam in samples:
            ih = 'resp'+var+'_nom'+sel
            tmpHisto = samples[isam][0].Get( 'jetObservables/'+ih )
            for ibin in variables[var][0]:
                if isinstance(ibin, int):
                    allHistos[isam+'_'+ih+str(ibin)] = tmpHisto.Clone()
                    allHistos[isam+'_'+ih+str(ibin)].SetName('Rebin_'+str(ibin))
                    allHistos[isam+'_'+ih+str(ibin)].Rebin2D( ibin, ibin )
                else:
                    #### fancy way to create variable binning TH2D
                    ##### IT IS NOT USED IN THIS SCRIPT YET
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


    return allHistos



###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="data", help="Process to unfold: data or MC." )
    parser.add_argument("-s", "--selection", action='store', dest="selection", default="_dijetSel", help="Selection to unfold: _dijetSel, _WSel, _topSel" )
    ##parser.add_argument("-r", "--runCombine", action='store_true', dest="runCombine", help="Run combine (true)" )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    inputFolder='Rootfiles/'+args.version

    bkgFiles = {}
    bkgFiles['QCDHerwig'] = [
                            ROOT.TFile.Open( checkDict( 'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7', dictSamples )[args.year]['skimmerHisto'] ),
                            checkDict( 'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7', dictSamples )
            ]

    variables = {}
    for ijet in [ ('Jet1', 'Outer'), ('Jet2', 'Central') ]:
        variables[ ijet[0]+'_tau21' ] = [ [ 20, 10, 5, 2 ], ijet[1]+' AK8 jet #tau_{21}' ]  #### original bin width is 0.1, rebinning means 0.1 times number in list
        variables[ ijet[0]+'_tau32' ] = [ [ 20, 10, 5, 2 ], ijet[1]+' AK8 jet #tau_{32}' ]
        variables[ ijet[0]+'_tau_0p5_1' ] = [ [ 40, 20, 10, 4 ], ijet[1]+' AK8 jet #tau_{1}^{0.5}' ]
        variables[ ijet[0]+'_tau_0p5_2' ] = [ [ 100, 50, 25, 10 ], ijet[1]+' AK8 jet #tau_{2}^{0.5}' ]
        variables[ ijet[0]+'_tau_0p5_3' ] = [ [ 100, 50, 25, 10 ], ijet[1]+' AK8 jet #tau_{3}^{0.5}' ]
        variables[ ijet[0]+'_tau_0p5_4' ] = [ [ 100, 50, 25, 10 ], ijet[1]+' AK8 jet #tau_{4}^{0.5}' ]
        variables[ ijet[0]+'_tau_1_1' ] = [ [ 40, 20, 10, 4 ], ijet[1]+' AK8 jet #tau_{1}^{1}' ]
        variables[ ijet[0]+'_tau_1_2' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{2}^{1}' ]
        variables[ ijet[0]+'_tau_1_3' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{3}^{1}' ]
        variables[ ijet[0]+'_tau_1_4' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{4}^{1}' ]
        variables[ ijet[0]+'_tau_2_1' ] = [ [ 40, 20, 10, 4 ], ijet[1]+' AK8 jet #tau_{1}^{2}' ]
        variables[ ijet[0]+'_tau_2_2' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{2}^{2}' ]
        variables[ ijet[0]+'_tau_2_3' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{3}^{2}' ]
        variables[ ijet[0]+'_tau_2_4' ] = [ [ 200, 100, 50, 25 ], ijet[1]+' AK8 jet #tau_{4}^{2}' ]

    runPurity( bkgFiles, variables, args.selection )


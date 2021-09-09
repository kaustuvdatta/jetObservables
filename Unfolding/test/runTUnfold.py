#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import ROOT
from array import array
import numpy as np
import yoda
from collections import OrderedDict
from multiprocessing import Process
from DrawHistogram import plotSimpleComparison, plotSysComparison
from variables import nSubVariables, nSubVariables_WSel, nSubVariables_topSel
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

colors = [ 2, 4, 9, 8, 28, 30, 42, 13, 46 ]
##########################################################################
def runTUnfold( dataFile, sigFiles, bkgFiles, variables, sel, sysUncert ):
    """docstring for createDataCards"""

    for ivar in variables:

        outputDir=args.outputFolder+'Plots/'+sel.split('_')[1]+'/Unfold/'+args.year+'/'+ivar+'/'+args.process+'/'
        if not os.path.exists(outputDir): os.makedirs(outputDir)

        if args.year.startswith('all'):
            print ('resp'+ivar+'_nom'+sel+signalLabel+'_Rebin')
            signalHistos = {
                    signalLabel+'_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+'_nom'+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel),
                    signalLabel+'_accepgen'+ivar+sel : dataFile[ivar+'_2017'].Get(signalLabel+'_accepgen'+ivar+sel),
                    signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin': dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin')
                    }
            if not args.process.startswith('MCCrossClosure'):
                signalHistos[ signalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+'_nom'+sel) )
                signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel) )
                signalHistos[ signalLabel+'_accepgen'+ivar+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_accepgen'+ivar+sel) )
                signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin') )
            
            for sys in sysUncert:
                for upDown in [ 'Up', 'Down' ]:
                    if sys.startswith(('_model', '_hdamp', '_Tune')): continue
                    signalHistos[ signalLabel+'_reco'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(signalLabel+'_reco'+ivar+sys+upDown+sel)
                    if not args.process.startswith('MCCrossClosure'): signalHistos[ signalLabel+'_reco'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+sys+upDown+sel) )
                    signalHistos[ signalLabel+'_resp'+ivar+sys+upDown+sel ] = dataFile[ivar+'_2017'].Get(signalLabel+'_resp'+ivar+sys+upDown+sel)
                    if not args.process.startswith('MCCrossClosure'): signalHistos[ signalLabel+'_resp'+ivar+sys+upDown+sel ].Add( dataFile[ivar+'_2018'].Get(signalLabel+'_resp'+ivar+sys+upDown+sel) )

            altSignalHistos = {
                altSignalLabel+'_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel),
                altSignalLabel+'_recoJetfrom_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(altSignalLabel+'_recoJetfrom_resp'+ivar+'_nom'+sel),
                altSignalLabel+'_genJetfrom_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(altSignalLabel+'_genJetfrom_resp'+ivar+'_nom'+sel),
                }
           
            if not args.process.startswith('MCCrossClosure'):
                altSignalHistos[ altSignalLabel+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_resp'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_recoJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_recoJetfrom_resp'+ivar+'_nom'+sel) )
                altSignalHistos[ altSignalLabel+'_genJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(altSignalLabel+'_genJetfrom_resp'+ivar+'_nom'+sel) )
            '''
            varSignalHistos = {
                            varSignalLabels[0]+'_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(varSignalLabels[0]+'_resp'+ivar+'_nom'+sel),
                            varSignalLabels[0]+'_recoJetfrom_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(varSignalLabels[0]+'_recoJetfrom_resp'+ivar+'_nom'+sel),
                            varSignalLabels[0]+'_genJetfrom_resp'+ivar+'_nom'+sel : dataFile[ivar+'_2017'].Get(varSignalLabels[0]+'_genJetfrom_resp'+ivar+'_nom'+sel),
                            }
                       
            if not args.process.startswith('MCCrossClosure'):
                varSignalHistos[ varSignalLabels[0]+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(varSignalLabels[0]+'_resp'+ivar+'_nom'+sel) )
                varSignalHistos[ varSignalLabels[0]+'_recoJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(varSignalLabels[0]+'_recoJetfrom_resp'+ivar+'_nom'+sel) )
                varSignalHistos[ varSignalLabels[0]+'_genJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(varSignalLabels[0]+'_genJetfrom_resp'+ivar+'_nom'+sel) )

            for v in varSignalLabels[1:]:
                varSignalHistos[ v+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(v+'_resp'+ivar+'_nom'+sel) )
                varSignalHistos[ v+'_recoJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(v+'_recoJetfrom_resp'+ivar+'_nom'+sel) )
                varSignalHistos[ v+'_genJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2017'].Get(v+'_genJetfrom_resp'+ivar+'_nom'+sel) )
                if not args.process.startswith('MCCrossClosure'):
                    varSignalHistos[ v+'_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(v+'_resp'+ivar+'_nom'+sel) )
                    varSignalHistos[ v+'_recoJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(v+'_recoJetfrom_resp'+ivar+'_nom'+sel) )
                    varSignalHistos[ v+'_genJetfrom_resp'+ivar+'_nom'+sel ].Add( dataFile[ivar+'_2018'].Get(v+'_genJetfrom_resp'+ivar+'_nom'+sel) )
            '''
            allHistos = {
                    'dataMinusBkgs' : dataFile[ivar+'_2017'].Get( 'dataMinusBkgs' ),
                    'dataMinusBkgsGenBin' : dataFile[ivar+'_2017'].Get( 'dataMinusBkgsGenBin' )
                    }
            allHistos[ 'dataMinusBkgs' ].Add( dataFile[ivar+'_2018'].Get( 'dataMinusBkgs' ) )
            allHistos[ 'dataMinusBkgsGenBin' ].Add( dataFile[ivar+'_2018'].Get( 'dataMinusBkgsGenBin' ) )

            if args.process.startswith('MCCrossClosure'):
                allHistos[ 'dataMinusBkgs' ] = dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel)
                allHistos[ 'dataMinusBkgsGenBin' ] = dataFile[ivar+'_2018'].Get(signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin')


            dataHistos = { }
            bkgHistos = { }

        else:
            print('|-------> Running single year '+args.year)
            ### Getting input histos
            mainSigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(signalLabelBegin)  }
            signalHistos = loadHistograms( mainSigFiles, ivar, sel, sysUnc=[ isys for isys in  sysUncert if not isys.startswith(('_model', '_hdamp', '_Tune')) ], lumi=args.lumi, year=args.year, process=args.process, variables=variables )
            print ("Signal histograms:", signalHistos)
            tmp2SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(altSignalLabelBegin)  }
            altSignalHistos = loadHistograms( tmp2SigFiles, ivar, sel, sysUnc=[], respOnly=True, lumi=args.lumi, year=args.year, process=args.process, variables=variables)
            tmp3SigFiles = { k:v for (k,v) in sigFiles.items() if k.startswith(varSignalLabelBegin)  }
            varSignalHistos = loadHistograms( tmp3SigFiles, ivar, sel, sysUnc=[], respOnly=True, lumi=args.lumi, year=args.year, process=args.process, variables=variables)

            if args.process.startswith("MC"):
                if args.process.startswith('MCSelfClosure'):
                    dataHistos = { 'data_reco'+k.split(('_reco' if sel.startswith('_dijet') else 'Leptonic'))[1] : v for (k,v) in signalHistos.items() if '_reco' in k }
                else:
                    dataHistos = loadHistograms( tmp2SigFiles, ivar, sel, isMC=True, addGenInfo=False, lumi=args.lumi, year=args.year, process=args.process, variables=variables )
                    dataHistos = { 'data_reco'+k.split('_reco')[1] : v for (k,v) in dataHistos.items() if '_reco' in k }
                    
            else:
                dataHistos = loadHistograms( dataFile, ivar, sel, isMC= False, lumi=args.lumi, year=args.year, process=args.process, variables=variables )
            print("Data histos", dataHistos.keys())

            bkgHistos = loadHistograms( bkgFiles, ivar, sel, sysUnc=[], lumi=args.lumi, year=args.year, process=args.process, variables=variables) if args.process.startswith('data') else {}

            ####### Fix fake and true reco
            signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel]= signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel] = signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone()
            signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel].Add( signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel], -1 )

            signalHistos[signalLabel+'_accepgen'+ivar+sel] = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionX()
            signalHistos[signalLabel+'_missgen'+ivar+sel] = signalHistos[signalLabel+'_gen'+ivar+sel].Clone()
            signalHistos[signalLabel+'_missgen'+ivar+sel].Add( signalHistos[signalLabel+'_accepgen'+ivar+sel], -1 )

            print '|------> Unfolding '+ivar

            ######## Cross check: plotting data vs all MC Scaled
            print ('|------> Cross check: plotting data vs all MC')
            allHistos = {}
            print (dataHistos.keys())
            allHistos[ 'allBkgHisto' ] = dataHistos['data_reco'+ivar+'_nom'+sel].Clone()
            allHistos[ 'allBkgHisto' ].Reset()
            allHistos[ 'allBkgHistoGenBin' ] = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Clone()
            allHistos[ 'allBkgHistoGenBin' ].Reset()
            for ibkg in bkgHistos:
                if ibkg.endswith('_reco'+ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( bkgHistos[ibkg].Clone() )
                if ibkg.endswith('_reco'+ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( bkgHistos[ibkg].Clone() )
            allHistos[ 'allMCHisto' ] = allHistos[ 'allBkgHisto' ].Clone()
            allHistos[ 'allMCHisto' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel ].Clone() )
            allHistos[ 'allMCHistoGenBin' ] = allHistos[ 'allBkgHistoGenBin' ].Clone()
            allHistos[ 'allMCHistoGenBin' ].Add( signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone() )
            ### For dijet, scale QCD to data
            scaleFactor =1
            if sel.startswith('_dijet'):
                scaleFactor = dataHistos['data_reco'+ivar+'_nom'+sel].Integral() / allHistos[ 'allMCHisto' ].Integral()
                scaleFactorGenBin = dataHistos['data_reco'+ivar+'_nom'+sel+'_genBin'].Integral() / allHistos[ 'allMCHistoGenBin' ].Integral()
                allHistos[ 'allMCHisto' ].Scale( scaleFactor )
                allHistos[ 'allBkgHisto' ].Scale( scaleFactor )
                allHistos[ 'allBkgHistoGenBin' ].Scale( scaleFactorGenBin )
                for ihsig in signalHistos:
                    if ihsig.endswith(sel):
                        print (ihsig)
                        signalHistos[ihsig].Scale( scaleFactor )
                    if ihsig.endswith('genBin'): signalHistos[ihsig].Scale( scaleFactorGenBin )
                    #if 'resp' in ihsig: signalHistos[ihsig].Scale( scaleFactor )   ### dont notice the difference

            plotSimpleComparison( dataHistos['data_reco'+ivar+'_nom'+sel].Clone(), 'data', allHistos[ 'allMCHisto' ], 'allBkgs', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+"_nom", rebinX=variables[ivar]['bins'], version=sel+'_'+args.version, outputDir=outputDir )

#            ######## Adding missgen to overflow and underflow
#            for ibin in range( 1, signalHistos[signalLabel+'_missgen'+ivar+sel].GetNbinsX()+1 ):
#                signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].SetBinContent(  signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetBin( 0, ibin ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )
#                signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].SetBinError(  signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetBin( 0, ibin ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinError(ibin) )
#                signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].SetBinContent(  signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetBin( -1, ibin ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )
#                signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].SetBinError(  signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetBin( -1, ibin ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinError(ibin) )
#                altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].SetBinContent(  altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].GetBin( 0, ibin ), altSignalHistos[altSignalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )
#                altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].SetBinError(  altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].GetBin( 0, ibin ), altSignalHistos[altSignalLabel+'_missgen'+ivar+sel].GetBinError(ibin) )
#                altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].SetBinContent(  altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].GetBin( -1, ibin ), altSignalHistos[altSignalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )
#                altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].SetBinError(  altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].GetBin( -1, ibin ), altSignalHistos[altSignalLabel+'_missgen'+ivar+sel].GetBinError(ibin) )
##                #for sys in sysUncert:
##                #    for upDown in [ 'Up', 'Down' ]:
##                #        signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel].SetBinContent(  signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel].GetBin( ibin, 0 ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )
##                #        signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel].SetBinContent(  signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel].GetBin( ibin, -1 ), signalHistos[signalLabel+'_missgen'+ivar+sel].GetBinContent(ibin) )

            ####### Cross check response matrix
            tmpGenHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionX()
            plotSimpleComparison( tmpGenHisto, 'projection', signalHistos[signalLabel+'_accepgen'+ivar+sel].Clone(), 'Regular AccepGen', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionGen", rebinX=1, version=sel+'_'+args.version, outputDir=outputDir )
            tmpRecoHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            plotSimpleComparison( tmpRecoHisto, 'projection', signalHistos[signalLabel+'_truereco'+ivar+'_nom'+sel].Clone(), 'Regular TrueReco', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+"_TestProjectionReco", rebinX=1, version=sel+'_'+args.version, outputDir=outputDir )
            tmpRecoHisto.Scale( scaleFactor )

            ####### Removing bkgs from data
            fakeHistos = { k:v for (k,v) in signalHistos.items() if 'fakereco' in k }
            for ih in fakeHistos:
                if ih.endswith(ivar+'_nom'+sel): allHistos[ 'allBkgHisto' ].Add( fakeHistos[ih] )
                if ih.endswith(ivar+'_nom'+sel+'_genBin'): allHistos[ 'allBkgHistoGenBin' ].Add( fakeHistos[ih] )
            plotSimpleComparison( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone(), 'data', allHistos[ 'allBkgHisto' ], 'Bkg+fakes', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataBkgFakes", rebinX=1, version=sel+'_'+args.version, outputDir=outputDir )

            allHistos[ 'dataMinusBkgs' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone()
            allHistos[ 'dataMinusBkgs' ].Add( allHistos[ 'allBkgHisto' ].Clone(), -1 )
            allHistos[ 'dataMinusBkgs' ].Scale( 1/allHistos[ 'dataMinusBkgs' ].Integral() )
            allHistos[ 'dataMinusBkgsGenBin' ] = dataHistos[ 'data_reco'+ivar+'_nom'+sel+'_genBin' ].Clone()
            allHistos[ 'dataMinusBkgsGenBin' ].Add( allHistos[ 'allBkgHistoGenBin' ].Clone(), -1 )
            allHistos[ 'dataMinusBkgsGenBin' ].Scale( 1/allHistos[ 'dataMinusBkgsGenBin' ].Integral() )

            tmpHisto = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].ProjectionY()
            tmpHisto.Scale( 1/tmpHisto.Integral() )
            plotSimpleComparison( allHistos[ 'dataMinusBkgs' ], 'data-Bkgs', tmpHisto, 'signal true reco', ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+"_TestDataMinusBkgs", rebinX=1, version=sel+'_'+args.version, outputDir=outputDir )


        ######## Cross check: plotting response matrix
        print ('|------> Cross check: plotting response matrix for signal')
        ROOT.gStyle.SetPadRightMargin(0.15)
        #ROOT.gStyle.SetPalette(ROOT.kGistEarth)
        #ROOT.TColor.InvertPalette()
        can2D = ROOT.TCanvas(ivar+'can2D', ivar+'can2D', 750, 500 )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetXaxis().SetTitle('Accepted Gen '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitle('True Reco '+variables[ivar]['label'])
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetTitleOffset( 0.8 )
        #signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].GetYaxis().SetRange( variables[ivar]['bins'][0], variables[ivar]['bins'][-1] )
        signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel].Draw("colz")
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2017+2018' if args.year.startswith('all') else args.year )
        CMS_lumi.relPosX = 0.12
        CMS_lumi.CMS_lumi(can2D, 4, 0)
        can2D.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+sel+'_responseMatrix'+args.version+'.'+args.ext)

        if args.runMLU: runMaxLikelihoodUnfold( dataHistos, signalHistos, bkgHistos, ivar, sel, signalLabel, sysUncert, outputDir )
        else:
            ######## TUnfold part
            print '|------> TUnfolding starts:'

            ##### Defining options for TUnfold
            tunfolder = ROOT.TUnfoldDensity(
                                                signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel], ### response matrix. According to TUnfold, this distribution does NOT have to be normalized
                                                ROOT.TUnfold.kHistMapOutputHoriz,  #### kHistMapOutputVert if x->reco and y->gen, kHistMapOutputHoriz if x->gen and y->reco
                                                ROOT.TUnfold.kRegModeCurvature,   ##### Regularization Mode : ROOT.TUnfold.kRegModeCurvature regularizes based on the 2nd derivative of the output. More information wrt the other options can be gained from reading the source code
                                                ROOT.TUnfold.kEConstraintNone,    ##### Constraint : TUnfold.kEConstraintNone meaning we do not constrain further, the other option is to force constraint of area. (Need to look into this!!)
                                                ROOT.TUnfoldDensity.kDensityModeBinWidth  ##### Density Mode: ROOT.TUnfoldDensity.kDensityModeBinWidth uses the bin width to normalize the event rate in a given bin, accounting for non-uniformity in bin widths as discussed in section 7.2.1 of the TUnfold paper
                                                )

            ##### Defining input (data recoJet )
            print '|------> TUnfolding adding input:'
            #tunfolder.SetInput( dataHistos[ 'data_reco'+ivar+'_nom'+sel ].Clone() )
            tunfolder.SetInput( allHistos[ 'dataMinusBkgs' ] )

            ###### Removing bkgs from data using TUnfold. Better to subtract bkgs beforehand
            #for ibkg in bkgHistos:
            #    if ibkg.endswith('_recoJet'+ivar+'_nom'+sel+'_Normalized'):
            #        print '|--------> Removing this bkg: ', ibkg
            #        tunfolder.SubtractBackground( bkgHistos[ibkg], ibkg )
            #####tunfolder.SubtractBackground( allHistos[ 'allBkgHisto' ], 'bkg', 1 )

            ###### Adding SYS unc
            if len(sysUncert)>0 :
                print '|------> TUnfolding adding uncert:'
                dictUncHistos = {}
                for sys in sysUncert:
                    if sys.startswith(('_puWeight','_jesTotal','_jer')):
                        dictUncHistos[sys+' Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                        dictUncHistos[sys+' Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                        for upDown in [ 'Up', 'Down' ]:
                            print sys+upDown
                            tunfolder.AddSysError(
                                                signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel],
                                                sys+upDown,
                                                ROOT.TUnfold.kHistMapOutputHoriz,
                                                ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                                )
                            can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+sys+upDown, ivar+'can2DNorm'+sys+upDown, 750, 500 )
                            signalHistos[signalLabel+'_resp'+ivar+sys+upDown+sel].Draw("colz")
                            can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+signalLabel+sel+sys+upDown+'Normalized_responseMatrix'+args.version+'.'+args.ext)

                    #### adding model uncertainty
                    elif sys.startswith(('_model')):
                        #if not args.process.startswith('MCSelfClosure'):
                        print('modelUnc')
                        tunfolder.AddSysError(
                                            altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel],
                                            'modelUncTotal',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNormAltSignal', ivar+'can2DNormAltSignal', 750, 500 )
                        altSignalHistos[altSignalLabel+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+altSignalLabel+sel+'Normalized_alt_responseMatrix'+args.version+'.'+args.ext)

                    #### adding model variation uncertainties
                    #varSignalLabels = ['varTTToSemileptonic_TuneCP5_erdON', 'varTTToSemileptonic_TuneCP5CR1', 'varTTToSemileptonic_TuneCP5CR2', 
                    #       'varTTToSemileptonic_TuneCP5Up','varTTToSemileptonic_TuneCP5Down','varTTToSemileptonic_hdampUp', 'varTTToSemileptonic_hdampDown', 
                    #       'varTTToSemileptonic_mtop166p5', 'varTTToSemileptonic_mtop169p5', 'varTTToSemileptonic_mtop171p5', 'varTTToSemileptonic_mtop173p5', 'varTTToSemileptonic_mtop175p5']
                    
                    elif sys.startswith('_hdamp'):
                        #if not args.process.startswith('MCSelfClosure'):
                        print('hdampUnc')
                        #dictUncHistos[sys+' Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                        #dictUncHistos[sys+' Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_hdampUp'+'_resp'+ivar+'_nom'+sel],
                                            'hdampUp',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'hdampUp', ivar+'can2DNorm'+'hdampUp', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_hdampUp'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampUp'+'_resp'+sel+'Normalized_responseMatrix'+args.version+'.'+args.ext)
                        
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_hdampDown'+'_resp'+ivar+'_nom'+sel],
                                            'hdampDown',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'hdampDown', ivar+'can2DNorm'+'hdampDown', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_hdampDown'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_hdampDown'+'_resp'+sel+'Normalized_responseMatrix'+args.version+'.'+args.ext)

                    elif sys.startswith('_Tune'):
                        #if not args.process.startswith('MCSelfClosure'):
                        print('TuneUnc')
                        #dictUncHistos[sys+' Up'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Up'+sel].Clone()
                        #dictUncHistos[sys+' Down'] = signalHistos[signalLabel+'_reco'+ivar+sys+'Down'+sel].Clone()
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_resp'+ivar+'_nom'+sel],
                                            'TuneCP5Up',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5Up', ivar+'can2DNorm'+'TuneCP5Up', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_TuneCP5Up'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Up'+'_resp'+sel+'Normalized_responseMatrix'+args.version+'.'+args.ext)
                        
                        tunfolder.AddSysError(
                                            varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_resp'+ivar+'_nom'+sel],
                                            'TuneCP5Down',
                                            ROOT.TUnfold.kHistMapOutputHoriz,
                                            ROOT.TUnfoldSys.kSysErrModeMatrix, #### kSysErrModeMatrix the histogram sysError corresponds to an alternative response matrix. kSysErrModeShift the content of the histogram sysError are the absolute shifts of the response matrix. kSysErrModeRelative the content of the histogram sysError specifies the relative uncertainties
                                            )
                        can2DNorm = ROOT.TCanvas(ivar+'can2DNorm'+'TuneCP5Down', ivar+'can2DNorm'+'TuneCP5Down', 750, 500 )
                        varSignalHistos['varTTToSemileptonic_TuneCP5Down'+'_resp'+ivar+'_nom'+sel].Draw("colz")
                        can2DNorm.SaveAs(outputDir+ivar+'_from'+('Data' if args.process.startswith('data') else 'MC')+'_'+'TTToSemileptonic_TuneCP5Down'+'_resp'+sel+'Normalized_responseMatrix'+args.version+'.'+args.ext)
                        

                #### Making unc plot
                plotSysComparison( signalHistos[signalLabel+'_reco'+ivar+'_nom'+sel].Clone(),
                                    dictUncHistos,
                                    ivar+'_'+signalLabel+'_allSys',
                                    labelX=variables[ivar]['label'],
                                    version=sel+'_'+args.version,
                                    year= ( '2017+2018' if args.year.startswith('all') else args.year ),
                                    outputDir=outputDir
                                    )

            ###### Running the unfolding
            print '|------> TUnfolding doUnfold:'
            tunfolder.DoUnfold(0)

            ###### Regularization
    #        nScan=50
    #        tauMin=0.0
    #        tauMax=0.0
    #        iBest=0
    #
    #        logTauX = ROOT.MakeNullPointer(ROOT.TSpline)
    #        logTauY = ROOT.MakeNullPointer(ROOT.TSpline)
    #        lCurve = ROOT.MakeNullPointer(ROOT.TGraph)
    #        ## this method scans the parameter tau and finds the kink in the L curve finally, the unfolding is done for the best choice of tau
    #        tunfolder.ScanLcurve(nScan,tauMin,tauMax,lCurve,logTauX,logTauY)
            #########################

            ##### Get output of unfolding
            allHistos [ 'unfoldHisto'+ivar ] = tunfolder.GetOutput("unfoldHisto"+ivar).Clone()
            allHistos [ 'foldHisto'+ivar ] = tunfolder.GetFoldedOutput("folded"+ivar).Clone()
            
            #### Get Probability matrix
            allHistos[ 'probaMatrix'+ivar ] = tunfolder.GetProbabilityMatrix('probaMatrix'+ivar).Clone()

            #### Get various covariances
            print '|------> TUnfolding covariances'
            allHistos[ 'cov'+ivar ] = tunfolder.GetEmatrixTotal("cov"+ivar, "Covariance Matrix")
            allHistos[ 'cov_uncorr_'+ivar ] = tunfolder.GetEmatrixSysUncorr("cov_uncorr"+ivar, "Covariance Matrix from Uncorrelated Uncertainties")
            allHistos[ 'cov_uncorr_data_'+ivar ] = tunfolder.GetEmatrixInput("cov_uncorr_data"+ivar, "Covariance Matrix from Stat Uncertainties of Input Data")
            #### cov = cov_uncorr + cov_uncorr_data + other uncertaitnies
            #### stat = cov_uncorr + cov_uncorr_data
            
            uncerUnfoldHisto = OrderedDict()
            allHistos[ 'unfoldHistowoUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone()        # Unfolding and stat unc
            allHistos[ 'unfoldHistoStatUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone("unfoldHistoStatUnc")        # Unfolding and stat unc
            allHistos[ 'unfoldHistoTotUnc'+ivar ] = allHistos[ 'unfoldHisto'+ivar ].Clone("unfoldHistoTotUnc")          # Total uncertainty
            uncerUnfoldHisto[ivar+'_StatTotal'] = allHistos[ 'unfoldHistoTotUnc'+ivar ].Clone(ivar+'_StatTotal')
            uncerUnfoldHisto[ivar+'_StatTotal'].Reset()
            uncerUnfoldHisto[ivar+'_Data statTotal'] = uncerUnfoldHisto[ivar+'_StatTotal'].Clone(ivar+'_InputStatTotal')
            uncerUnfoldHisto[ivar+'_Resp. Matrix statTotal'] = uncerUnfoldHisto[ivar+'_StatTotal'].Clone(ivar+'_UncorUncerTotal')
            allHistos[ 'unfoldHistoOnlyStatUnc'+ivar ] = uncerUnfoldHisto[ivar+'_StatTotal'].Clone('unfoldHistoOnlyStatUnc'+ivar)  # only err for ratio
            for ibin in range( 0, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX()+1 ):
                unc_tot = ROOT.TMath.Sqrt( allHistos[ 'cov'+ivar ].GetBinContent(ibin,ibin) )
                allHistos[ 'unfoldHistoTotUnc'+ivar ].SetBinContent(ibin, unc_tot )
                allHistos[ 'unfoldHisto'+ivar ].SetBinError(ibin, unc_tot )
                stat_tot = ROOT.TMath.Sqrt(allHistos[ 'cov_uncorr_data_'+ivar ].GetBinContent(ibin,ibin) + allHistos[ 'cov_uncorr_'+ivar ].GetBinContent(ibin,ibin))
                uncerUnfoldHisto[ivar+'_StatTotal'].SetBinContent(ibin, stat_tot)
                allHistos[ 'unfoldHistoOnlyStatUnc'+ivar ].SetBinContent( ibin, 1 )
                allHistos[ 'unfoldHistoOnlyStatUnc'+ivar ].SetBinError( ibin, stat_tot/allHistos[ 'unfoldHistoOnlyStatUnc'+ivar ].GetBinWidth(ibin) )
                uncerUnfoldHisto[ivar+'_Data statTotal'].SetBinContent(ibin, ROOT.TMath.Sqrt(allHistos[ 'cov_uncorr_data_'+ivar ].GetBinContent(ibin,ibin)))
                uncerUnfoldHisto[ivar+'_Resp. Matrix statTotal'].SetBinContent(ibin, ROOT.TMath.Sqrt(allHistos[ 'cov_uncorr_'+ivar ].GetBinContent(ibin,ibin)))
                
            ##### Get systematic shifts of output
            if len(sysUncert)>0 :
                print '|------> TUnfolding uncertainties:'

                for sys in sysUncert:
                    if not sys.startswith('_model'):
                        for upDown in [ 'Up', 'Down' ]:
                            print sys+upDown
                            uncerUnfoldHisto[ivar+sys+upDown] = tunfolder.GetDeltaSysSource(sys+upDown, "unfoldHisto_"+ivar+sys+upDown+"shift", "-1#sigma")
                            try: uncerUnfoldHisto[ivar+sys+upDown].SetLineStyle(1)
                            except ReferenceError: uncerUnfoldHisto.pop( ivar+sys+upDown, None )

                        # Create total uncertainty and sys uncertainty plots.
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'] = allHistos[ 'unfoldHisto'+ivar ].Clone()       # Syst uncertainty
                        uncerUnfoldHisto[ivar+sys.upper()+'Total'].Reset()
                        for i in xrange( 0, allHistos[ 'unfoldHisto'+ivar ].GetNbinsX() + 1):
                            try: yup = abs( uncerUnfoldHisto[ivar+sys+'Up'].GetBinContent(i) )
                            except KeyError: yup = 0
                            try: ydn = abs( uncerUnfoldHisto[ivar+sys+'Down'].GetBinContent(i) )
                            except KeyError: ydn = 0
                            dy = ROOT.TMath.Sqrt( (yup**2 + ydn**2) )
                            uncerUnfoldHisto[ivar+sys.upper()+'Total'].SetBinContent( i, dy )

                    elif sys.startswith('_model'):
                        #if not args.process.startswith('MCSelfClosure'):
                        uncerUnfoldHisto[ivar+'_Physics ModelTotal'] = allHistos[ 'unfoldHisto'+ivar ].Clone(ivar+'_modelUncTotal')
                        uncerUnfoldHisto[ivar+'_Physics ModelTotal'].Reset()
                        tmpModelHisto = tunfolder.GetDeltaSysSource('modelUncTotal', "unfoldHisto_"+ivar+"modelUncshift", "-1#sigma")
                        for i in xrange( 0, tmpModelHisto.GetNbinsX() + 1):
                            uncerUnfoldHisto[ivar+'_Physics ModelTotal'].SetBinContent( i, abs(tmpModelHisto.GetBinContent(i)) )
                        print('_modelUncTotal')



            ###### Plot unfolding results
            print '|------> Drawing unfold plot:'
            drawUnfold( ivar, allHistos[ 'dataMinusBkgsGenBin' ].Clone(),
                            signalHistos[ signalLabel+'_accepgen'+ivar+sel ].Clone(),
                            allHistos[ 'unfoldHisto'+ivar ].Clone(),
                            allHistos[ 'unfoldHistowoUnc'+ivar ].Clone(),
                            tunfolder.GetFoldedOutput("folded"+ivar).Clone(),
                            signalHistos[ signalLabel+'_reco'+ivar+'_nom'+sel+'_genBin' ].Clone(),
                            allHistos[ 'unfoldHistoOnlyStatUnc'+ivar ].Clone(),
                            variables[ivar]['label'],
                            variables[ivar]['bins'][-1],
                            variables[ivar]['alignLeg'],
                            outputDir+ivar+sel+'_from'+('Data' if args.process.startswith('data') else 'MC')+signalLabel+(''.join(sysUncert))+'_Tunfold_'+args.version+'.'+args.ext
                            )

            ######### Plotting Uncertainties
            print '|------> Drawing unfold uncertainty plot:'
            drawUncertainties(ivar, allHistos[ 'unfoldHistoTotUnc'+ivar ],
                            uncerUnfoldHisto,
                            variables[ivar]['label'],
                            variables[ivar]['alignLeg'],
                            outputDir+ivar+sel+'_from'+('Data' if args.process.startswith('data') else 'MC')+(''.join(sysUncert))+'_Tunfold_UNC_'+args.version+'.'+args.ext
                            )

            ######### Saving Histograms
            def renamingHistos( dictHistos ):
                for isam, ihis in dictHistos.items():
                    ihis.SetName(isam)
                    ihis.SetTitle(isam)
                    ihis.Write()

            outputRootName = outputDir+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root'
            print '|------> Saving histograms in rootfile: ', outputRootName
            outputRoot = ROOT.TFile.Open( outputRootName, 'recreate' )
            renamingHistos( signalHistos )
            renamingHistos( altSignalHistos )
            renamingHistos( dataHistos )
            renamingHistos( bkgHistos )
            renamingHistos( allHistos )
            renamingHistos( uncerUnfoldHisto )
            tunfolder.Write()
            outputRoot.Close()

            print '|------> Saving histograms in yodafile: ', outputRootName.replace('.root', '.yoda')
            histToYoda = [  yoda.root.to_yoda( allHistos [ 'unfoldHisto'+ivar ] ) ]
            yoda.writeYODA( histToYoda, outputRootName.replace('.root', '.yoda') )


##########################################################################
def loadHistograms( samples, var, sel, sysUnc=[], isMC=True, addGenInfo=True, respOnly=False, lumi=1., variables={}, year='2017', process='data' ):
    """docstring for loadHistograms"""

    SYSUNC = [ '_nom' ] + [ s+u for u in ['Up', 'Down'] for s in sysUnc ]

    allHistos = {}
    for isam in samples:
        tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC ]
        if isMC and addGenInfo: tmpList = tmpList + [ 'gen'+var+sel ] + [ 'resp'+var+syst+sel for syst in SYSUNC ]
        #if isMC and addGenInfo: tmpList = tmpList + [ 'gen'+var+sel, 'missgen'+var+sel, 'accepgen'+var+sel, 'truereco'+var+'_nom'+sel, 'fakereco'+var+'_nom'+sel ] + [ 'resp'+var+syst+sel for syst in SYSUNC ]
        if respOnly: tmpList = [ 'resp'+var+syst+sel for syst in SYSUNC ] #+ [ 'missgen'+var+sel ]
        for ih in tmpList:
            print 'Processing '+isam+' '+ih
            if isMC:
                allHistos[isam+'_'+ih] = samples[isam][0].Get( 'jetObservables/'+ih )
                tmpIsam = 'TT' if isam.startswith('data') else isam
                MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1][year][('nGenWeights' if process.startswith('data') else 'nGenWeights') ]
                allHistos[isam+'_'+ih].Scale( MCScale )
            else:
                if sel.startswith('_dijet') and process.startswith('data'):
                    tmpdataHistos = {}
                    for it in checkDict( 'JetHT', dictSamples )[year]['triggerList']:
                        tmpdataHistos[ it ] = samples[isam][0].Get( 'jetObservables/'+ih.replace( sel, '_'+it+sel ) )
                        tmpdataHistos[ it ].Scale( checkDict( 'JetHT', dictSamples )[year]['triggerList'][it] )
                    allHistos[ isam+'_'+ih ] = tmpdataHistos[next(iter(tmpdataHistos))].Clone()
                    allHistos[ isam+'_'+ih ].Reset()
                    for i in tmpdataHistos: allHistos[isam+'_'+ih].Add( tmpdataHistos[i] )
                else:
                    allHistos[isam+'_'+ih] = samples[isam][0].Get( 'jetObservables/'+ih )

            if len(variables[var]['bins'])==1:
                genBin = variables[var]['bins'][0]
                recoBin = variables[var]['bins'][0]/2
            else:
                genBin = variables[var]['bins']
                recoBin = np.sort(np.append( variables[var]['bins'], np.array([ (variables[var]['bins'][i]+variables[var]['bins'][i+1])/2 for i in range(len(variables[var]['bins'])-1) ]  ) ))

            if not ih.startswith('resp'):
                if len(variables[var]['bins'])==1:
                    if ih.startswith(('reco','fake', 'true')):
                        allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih].Clone()
                        allHistos[isam+'_'+ih+'_genBin'].Rebin( genBin )
                        allHistos[isam+'_'+ih].Rebin( recoBin )
                    else: allHistos[isam+'_'+ih].Rebin( genBin )
                else:
                    if ih.startswith(('reco','fake','true')):
                        allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih].Clone()
                        allHistos[isam+'_'+ih+'_genBin'] = allHistos[isam+'_'+ih+'_genBin'].Rebin( len(genBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin_genBin", array( 'd', genBin ) )
                        allHistos[isam+'_'+ih] = allHistos[isam+'_'+ih].Rebin( len(recoBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin", array( 'd', recoBin ) )
                    else:
                        allHistos[isam+'_'+ih] = allHistos[isam+'_'+ih].Rebin( len(genBin)-1, allHistos[isam+'_'+ih].GetName()+"_Rebin", array( 'd', genBin ) )
            else:
                if len(variables[var]['bins'])==1: allHistos[isam+'_'+ih].Rebin2D( genBin, recoBin )
                else:

                    #### fancy way to create variable binning TH2D
                    tmpHisto = ROOT.TH2F( allHistos[isam+'_'+ih].GetName()+isam+"_Rebin", allHistos[isam+'_'+ih].GetName()+isam+"_Rebin", len(genBin)-1, array( 'd', genBin), len(recoBin)-1, array( 'd', recoBin) )

                    tmpArrayContent = np.zeros((len(genBin), len(recoBin)))
                    tmpArrayError = np.zeros((len(genBin), len(recoBin)))

#                    for biny in range( 1, allHistos[isam+'_'+ih].GetNbinsY()+1 ):
#                        for binx in range( 1, allHistos[isam+'_'+ih].GetNbinsX()+1 ):
#                            tmpHisto.Fill( allHistos[isam+'_'+ih].GetXaxis().GetBinCenter(binx), allHistos[isam+'_'+ih].GetYaxis().GetBinCenter(biny), allHistos[isam+'_'+ih].GetBinContent( binx, biny ) )

                    for biny in range( 1, allHistos[isam+'_'+ih].GetNbinsY()+1 ):
                        by = allHistos[isam+'_'+ih].GetYaxis().GetBinCenter( biny )
                        for binx in range( 1, allHistos[isam+'_'+ih].GetNbinsX()+1 ):
                            bx = allHistos[isam+'_'+ih].GetXaxis().GetBinCenter(binx)
                            for iX in range( len(genBin)-1 ):
                                for iY in range( len(recoBin)-1 ):
                                    if (bx<genBin[iX+1] and bx>genBin[iX]) and (by<recoBin[iY+1] and by>recoBin[iY]):
                                        jbin = allHistos[isam+'_'+ih].GetBin(binx,biny)
                                        tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + allHistos[isam+'_'+ih].GetBinContent( jbin )
                                        tmpArrayError[iX][iY] = tmpArrayError[iX][iY] + ROOT.TMath.Power( allHistos[isam+'_'+ih].GetBinError( jbin ), 2 )

                    for biny in range( 1, tmpHisto.GetNbinsY()+1 ):
                        for binx in range( 1, tmpHisto.GetNbinsX()+1 ):
                            #if (binx <= ((biny+1)/2.)+4) and (binx >= ((biny+1)/2.)-4):
                            tmpHisto.SetBinContent( tmpHisto.GetBin(binx,biny), tmpArrayContent[binx-1][biny-1] )
                            tmpHisto.SetBinError( tmpHisto.GetBin(binx,biny), ROOT.TMath.Sqrt(tmpArrayError[binx-1][biny-1] ) )

                    tmpHisto.Sumw2()
                    allHistos[isam+'_'+ih] = tmpHisto

                #if isMC: allHistos[isam+'_'+ih].Scale( MCScale )
                ##### For tests, projections directly from 2D
                allHistos[isam+'_genJetfrom_'+ih] = allHistos[isam+'_'+ih].ProjectionY()
                allHistos[isam+'_recoJetfrom_'+ih] = allHistos[isam+'_'+ih].ProjectionX()


    tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k }
    for ih in tmpHistos:
        for jh in allHistos:
            if (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ):
                tmpHistos[ih].Add( allHistos[jh].Clone() )
    if len(tmpHistos)>0: allHistos = tmpHistos

    return allHistos

##########################################################################
def drawUnfold( ivar, dataJetHisto, genJetHisto, unfoldHisto, unfoldHistowoUnc, foldHisto, recoJetHisto, ratioUncHisto, labelX, maxX, tlegendAlignment, outputName ):
    """docstring for drawUnfold"""

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 750, 750 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.207,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()

    pad1.cd()

    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.65,0.90,0.88)
    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    dataIntegral = dataJetHisto.Integral()
    dataJetHisto.Scale(1, 'width')  ### divide by bin width
    #######dataJetHisto.Scale(1/dataJetHisto.Integral(), 'width')  ### divide by bin width
    dataJetHisto.SetLineWidth(2)
    dataJetHisto.SetLineStyle(2)
    dataJetHisto.SetLineColor(ROOT.kMagenta)
    legend.AddEntry( dataJetHisto, ('Data-Bkgs' if args.process.startswith('data') else 'MC '+('Self-' if args.process.startswith('MCSelfClosure') else '')+'Closure' ), 'l' )

    #genJetHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    genJetHisto.Scale(dataIntegral/genJetHisto.Integral(), 'width')  ### divide by bin width
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(1)
    legend.AddEntry( genJetHisto, 'Accepted Gen', 'l' )

    unfoldHisto.Scale(1, 'width')  ### divide by bin width
    #######unfoldHisto.Scale(dataIntegral/unfoldHisto.Integral(), 'width')  ### divide by bin width
    unfoldHisto.SetMarkerStyle(4)
    unfoldHisto.SetMarkerSize(2)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    legend.AddEntry( unfoldHisto, 'Unfolded, total unc', 'pl' )

    #unfoldHistowoUnc.Scale(1, 'width')  ### divide by bin width
    #unfoldHistowoUnc.Scale(dataIntegral/unfoldHisto.Integral(), 'width')  ### divide by bin width
    unfoldHistowoUnc.SetMarkerStyle(0)
    unfoldHistowoUnc.SetMarkerColor(ROOT.kRed)
    unfoldHistowoUnc.SetLineColor(ROOT.kRed-3)
    unfoldHistowoUnc.SetLineWidth(2)
    #legend.AddEntry( unfoldHistowoUnc, 'Unfolded, stat+unf unc', 'l' )

    #if args.selfClosure:
    foldHisto.Rebin( 2 )
    #foldHisto.Scale(1, 'width')  ### divide by bin width
    foldHisto.Scale(dataIntegral/foldHisto.Integral(), 'width')  ### divide by bin width
    foldHisto.SetLineWidth(2)
    foldHisto.SetLineStyle(2)
    foldHisto.SetLineColor(8)
    legend.AddEntry( foldHisto, 'Folded', 'l' )

    #recoJetHisto.Scale(1, 'width')  ### divide by bin width
    recoJetHisto.Scale(dataIntegral/recoJetHisto.Integral(), 'width')  ### divide by bin width
    recoJetHisto.SetLineWidth(2)
    recoJetHisto.SetLineStyle(2)
    recoJetHisto.SetLineColor(ROOT.kBlue)
    legend.AddEntry( recoJetHisto, 'True Reco level MC', 'l' )

    unfoldHisto.GetYaxis().SetTitle( '#frac{1}{d#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( 1.2*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum(), foldHisto.GetMaximum(), dataJetHisto.GetMaximum() ] )  )

    unfoldHisto.Draw( "histe")
    genJetHisto.Draw( "esame")
    #unfoldHistowoUnc.Draw( "e1 same")
    foldHisto.Draw( "histe same")
    #if args.selfClosure: foldHisto.Draw( "histe same")
    dataJetHisto.Draw( "histe same")
    recoJetHisto.Draw( "histe same")
    #if args.process.startswith('data'): recoJetHisto.Draw( "hist same")

    legend.Draw()
    if args.process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if args.selection.startswith('dijet') else '')+str( round( (args.lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if args.year.startswith('all') else ", "+( '2017+2018' if args.year.startswith('all') else args.year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2017+2018' if args.year.startswith('all') else args.year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)

    pad2.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    tmpPad2.GetXaxis().SetTitle( labelX )
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
    hRatioUp.Divide( genJetHisto, unfoldHisto, 'pois' )
    #hRatioUp.SetLineColor(kRed-4)
    #hRatioUp.SetLineWidth(2)
    hRatioUp.SetMarkerStyle(8)
    hRatioUp.Draw('P0')
    #hRatioDown.Draw('P same')
    ratioUncHisto.SetFillColor(ROOT.kBlack)
    ratioUncHisto.SetFillStyle(3004)
    ratioUncHisto.Draw('E2 same')

    can.SaveAs(outputName)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)

###########################################################################
def drawUncertainties( ivar, unfoldHistoTotUnc, uncerUnfoldHisto, labelX, tlegendAlignment, outputName ):
    """docstring for drawUncUncertainties"""

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 750, 500 )
    ##canUnc.SetLogy()

    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.65,0.65,0.90,0.88)
    else: legend=ROOT.TLegend(0.20,0.65,0.40,0.88)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.AddEntry( unfoldHistoTotUnc, 'Total Unc', 'l' )

    #unfoldHistoTotUnc.SetMinimum(0.00000000001)
    unfoldHistoTotUnc.SetLineWidth(2)
    uncScaleFactor = 1/unfoldHistoTotUnc.Integral()
    unfoldHistoTotUnc.Scale( uncScaleFactor )
    #unfoldHistoTotUnc.Scale( 1, 'width' )
    unfoldHistoTotUnc.GetXaxis().SetTitle(labelX)
    unfoldHistoTotUnc.GetYaxis().SetTitle('Fractional Uncertainty')
    #unfoldHistoTotUnc.GetYaxis().SetTitleOffset( 0.8 )
    unfoldHistoTotUnc.SetMarkerStyle(4)
    unfoldHistoTotUnc.Draw('hist')

    dummy=0
    for k in uncerUnfoldHisto:
        if k.endswith('Total') and not k.endswith(('StatTotal')):
            print k
            legend.AddEntry( uncerUnfoldHisto[k], k.split('_')[-1].split('Total')[0].replace('TOTAL', '').replace('WEIGHT', ''), 'l' )
            uncerUnfoldHisto[k].SetLineColor(colors[dummy])
            uncerUnfoldHisto[k].SetLineWidth(2)
            uncerUnfoldHisto[k].Scale( uncScaleFactor )
            #uncerUnfoldHisto[k].Scale( 1, 'width' )
            uncerUnfoldHisto[k].Draw("hist same")
            dummy=dummy+1

    ##### this is a test
    tmp = OrderedDict()
    for i in xrange( 0, unfoldHistoTotUnc.GetNbinsX() + 1):
        tmp[i] = 0
        for k in uncerUnfoldHisto:
            if k.endswith('Total') and not k.endswith(('StatTotal')):
                tmp[i] = tmp[i] + ( uncerUnfoldHisto[k].GetBinContent( i )**2 )
                #print(i, k, tmp[i], uncerUnfoldHisto[k].GetBinContent( i ), ( uncerUnfoldHisto[k].GetBinContent( i )**2 ))
    tmpHisto = unfoldHistoTotUnc.Clone()
    tmpHisto.Reset()
    for i,j in tmp.items():
        tmpHisto.SetBinContent( i, ROOT.TMath.Sqrt( j ) )
    tmpHisto.SetLineColor(ROOT.kBlack)
    #tmpHisto.Scale( 1, 'width' )
    tmpHisto.Draw("hist same")   ### this should be the same as TotalUnc

    legend.Draw()
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2017+2018' if args.year.startswith('all') else args.year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)

    canUnc.SaveAs(outputName)

###########################################################################
def runMaxLikelihoodUnfold( dataHistos, signalHistos, bkgHistos, ivar, sel, signalLabel, sysUncert, outputDir ):
    """docstring for runMaxLikelihoodUnfold"""

    print('|-----> RUNNING MAX LIKELIHOOD UNFOLD (combine)')

    textDict = OrderedDict()
    textDict['bin'] = [ 'bin' ]
    textDict['obs'] = [ 'observation' ]
    textDict['bin2'] = [ 'bin' ]
    textDict['proc'] = [ 'process' ]
    textDict['proc2'] = [ 'process' ]
    textDict['rate'] = [ 'rate' ]
    po = []
    setParameters=[]

    ### renaming histos to avoid confusion
    dataHisto = dataHistos['data_reco'+ivar+'_nom'+sel]
    dataHisto.Rebin(2)
    respMatrix = signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel]
    respMatrix.Rebin2D( 1, 2 )
    fakeHisto = signalHistos[signalLabel+'_fakereco'+ivar+'_nom'+sel]
    fakeHisto.Rebin(2)

    dataHisto.SetLineColor(ROOT.kRed)
    tmpMC = respMatrix.ProjectionY()
    tmpMC.Add( fakeHisto )
    tmpMC.SetLineColor(ROOT.kBlue)
    plotSimpleComparison( dataHisto, 'data', tmpMC, 'respPlusFake', ivar+'_'+signalLabel+"_TestRespFakeVsRecoCombine", rebinX=1, version=sel, outputDir=outputDir )

    for ibinReco in range(1, respMatrix.GetNbinsY()+1 ):
        if dataHisto.GetBinContent(ibinReco)<1: continue
        dataBin = dataHisto.GetBinContent(ibinReco)
        textDict['bin'].append('Reco_'+str(ibinReco))
        textDict['obs'].append( str( dataBin ) )
        po.append("--PO map='.*/Gen_"+str(ibinReco)+':r_bin'+str(ibinReco)+"[1,-10,10]'")
        setParameters.append("r_bin"+str(ibinReco)+"=1")

        binRecoSum=0
        for ibinGen in range(1, respMatrix.GetNbinsX()+1 ):
            #if respMatrix.GetBinContent(ibinGen,ibinReco) < (dataBin*0.01): continue
            ###print(ibinReco, ibinGen, respMatrix.GetBinContent(ibinGen,ibinReco))
            binRecoSum=binRecoSum+respMatrix.GetBinContent(ibinGen,ibinReco)
            textDict['bin2'].append("Reco_"+str(ibinReco))    ### one for ibinGen and one for ibinReco
            textDict['proc'].append("Gen_"+str(ibinGen))
            textDict['proc2'].append("-"+str(ibinGen))   ## 0 -1, -2 --> for signal
            textDict['rate'].append( str( round(respMatrix.GetBinContent(ibinGen,ibinReco),2)) )
        if(binRecoSum==0):
            print 'removing rbin'+str(ibinReco)
            if not textDict['bin'][-1].startswith('bin'): textDict['bin'].pop()
            if not textDict['obs'][-1].startswith('observation'): textDict['obs'].pop()
            if not textDict['bin2'][-1].startswith('bin'): textDict['bin2'].pop()
            if not textDict['proc'][-1].startswith('process'): textDict['proc'].pop()
            if not textDict['proc2'][-1].startswith('process'): textDict['proc2'].pop()
            if not textDict['rate'][-1].startswith('rate'): textDict['rate'].pop()
            po.pop()

        if fakeHisto.GetBinContent(ibinReco)<1: continue
        textDict['bin2'].append("Reco_"+str(ibinReco))    ### one for ibinGen and one for ibinReco
        textDict['proc'].append("Bkg")
        textDict['proc2'].append("1" ) #str(ibinReco))   ## 0 -1, -2 --> for signal
        textDict['rate'].append( str( round(fakeHisto.GetBinContent(ibinReco),2)) )

    ######## Creating datacard for combine
    datacardName = 'datacard_'+ivar+'_from'+signalLabel
    print '|------> Creating datacard: ', outputDir+'/'+datacardName

    datacard = open( outputDir+'/'+datacardName+'.txt', 'w')
    datacard.write("imax *\n")
    datacard.write("jmax *\n")
    datacard.write("kmax *\n")
    datacard.write("----------------\n")
    datacard.write("shapes * * FAKE\n")
    datacard.write("----------------\n")
    for x in textDict:
        datacard.write(' '.join(textDict[x])+'\n')
    datacard.write("----------------\n")
    datacard.close()

    #### Creating the combine command
    runScript = open( outputDir+'/runCombine'+signalLabel+'.sh', 'w')
    po = ' '.join(po)
    cmdText2Workspace = "text2workspace.py "+datacardName+".txt -o "+datacardName+".root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+po
    runScript.write("### RUN WITH COMMANDS: ####\n")
    runScript.write(cmdText2Workspace+"\n")
    setParameters= ','.join(setParameters)
    cmdCombine = 'combine -M MultiDimFit -d '+datacardName+'.root -n '+ivar+'_from'+signalLabel+' --algo singles --cminDefaultMinimizerStrategy=0'#' --setParameters='+setParameters
    runScript.write(cmdCombine+" \n")
    runScript.close()
    print '|------> To run combine: \n', cmdText2Workspace + '\n' + cmdCombine #+ '\n' + cmdCombine2


###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="data", help="Process to unfold: data or MC." )
    parser.add_argument("-s", "--selection", action='store', dest="selection", default="_dijetSel", help="Selection to unfold: _dijetSel, _WSel, _topSel" )
    parser.add_argument("--only", action='store', dest="only", default="", help="Submit only one variable" )
    parser.add_argument("--main", action='store', dest="main", choices=['Ptbin', 'HTbin', 'herwig'], default='Ptbin', help="For dijet sel, main signal QCD" )
    parser.add_argument("--alt", action='store', dest="alt", choices=['Ptbin', 'HTbin', 'herwig'], default='herwig', help="For dijet sel, alternative signal QCD" )
    parser.add_argument("--noUnc", action='store_true', dest="noUnc", default=False, help="Run without Unc " )
    parser.add_argument("-r", "--runMLU", action='store_true', default=False, help="Run max likelihood unfold (combine)" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016, 2017, 2018.' )
    parser.add_argument('-e', '--ext', action='store', default='png', help='Extension of plots.' )
    parser.add_argument('--plotOnly', action='store_true', default=False, dest='plotOnly',  help='Plot only.' )
    parser.add_argument("--inputFolder", action='store', dest="inputFolder", default="os.environ['CMSSW_BASE']+'/src/jetObservables/Unfolding/test/Samples/", help="input folder" )
    parser.add_argument("--outputFolder", action='store', dest="outputFolder", default="os.environ['CMSSW_BASE']+'/src/jetObservables/Unfolding/test/Results/", help="Output folder" )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=0., help='Luminosity, example: 1.' )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    #### define variables
    if args.selection=='_WSel': nSubVariables = nSubVariables_WSel
    elif args.selection=='_topSel': nSubVariables = nSubVariables_topSel

    if args.only:
        filterVariables = { k:v for (k,v) in nSubVariables.items() if k.endswith(args.only)  }
        if len(filterVariables)>0 : variables = filterVariables
        else:
            print('|------> Variable not found. Have a nice day')
            sys.exit(0)
    else: variables = nSubVariables
    if args.selection.startswith('_dijet'): variables = { k:v for (k,v) in variables.items() if k.startswith(('Jet1', 'Jet2', 'sdJet1', 'sdJet2')) }
    else: variables = { k:v for (k,v) in variables.items() if k.startswith('Jet_') }

    #### syst
    if args.selection.startswith('_dijet'): sysUncert = [] if (args.noUnc or args.process.startswith('MCSelfClosure')) else [ '_jesTotal', '_jer', '_puWeight', '_isrWeight', '_fsrWeight', '_model' ]
    elif args.selection.startswith('_W'): sysUncert = [] if (args.noUnc or args.process.startswith('MCSelfClosure')) else [ '_model', '_hdamp', '_Tune' ] #'_jesTotal', '_jer', '_puWeight', '_isrWeight', '_fsrWeight',
    elif args.selection.startswith('_top'): sysUncert = [] if (args.noUnc or args.process.startswith('MCSelfClosure')) else [ '_model']#, '_hdamp', '_Tune', '_CR', '_mtop' ] #'_jesTotal', '_jer', '_puWeight', '_isrWeight', '_fsrWeight',
    
    #### Define main and alt signals
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
        if args.alt.startswith('Ptbin'):
            altSignalLabelBegin = 'QCD_Pt_'
            altSignalLabel = 'QCD_Pt_3200toInf'
        elif args.alt.startswith('HTbin'):
            altSignalLabelBegin = 'QCD_HT'
            altSignalLabel = 'QCD_HT2000toInf' 
        elif args.alt.startswith('herwig'):
            altSignalLabelBegin = 'QCD_Pt-'
            altSignalLabel = 'QCD_Pt-150to3000'
    else:
        signalLabel = 'TTToSemiLeptonic'
        signalLabelBegin = 'TTToSemiLeptonic'
        varSignalLabelBegin = 'varTTToSemileptonic_'
        varSignalLabels = ['varTTToSemileptonic_TuneCP5Up','varTTToSemileptonic_TuneCP5Down','varTTToSemileptonic_hdampUp', 'varTTToSemileptonic_hdampDown', 
                            #'varTTToSemileptonic_TuneCP5_erdON', 'varTTToSemileptonic_TuneCP5CR1', 'varTTToSemileptonic_TuneCP5CR2', 
                            #'varTTToSemileptonic_mtop166p5', 'varTTToSemileptonic_mtop169p5', 'varTTToSemileptonic_mtop171p5', 'varTTToSemileptonic_mtop173p5', 'varTTToSemileptonic_mtop175p5'
                          ]
        altSignalLabelBegin = 'TTJets'
        altSignalLabel = 'TTJets'

    #### Files
    dataFile = {}
    sigFiles = {}
    bkgFiles = {}
    #### All data runs from single year files
    if args.year.startswith('all'):
        print('|------> Running combination of years')
        if not args.inputFolder:
            print('|------> For running all the data you must specify the input folder.\n Have a nice day :)')
            sys.exit(0)

        for ivar in variables.keys():
            dataFile[ivar+'_2017'] = ROOT.TFile.Open( args.inputFolder+'/Plots/'+args.selection.split('_')[1]+'/Unfold/2017/'+ivar+'/'+('MCClosure' if 'Cross' in args.process else args.process )+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root' )
            dataFile[ivar+'_2018'] = ROOT.TFile.Open( args.inputFolder+'/Plots/'+args.selection.split('_')[1]+'/Unfold/2018/'+ivar+'/'+('MCClosure' if 'Cross' in args.process else args.process )+'/outputHistograms_main_'+signalLabel+'_alt_'+altSignalLabel+'.root' )

    #### Single year unfolding runs on skimmers
    else:
        print('|------> Running single year')
        args.inputFolder = os.environ['CMSSW_BASE']+'/src/jetObservables/Unfolding/test/Samples/'
        if args.process.startswith('data'): dataFile['data'] = [ ROOT.TFile.Open(args.inputFolder+checkDict( ('JetHT' if args.selection.startswith("_dijet") else 'SingleMuon'), dictSamples )[args.year]['skimmerHisto']) ]
        for iy in ( ['2017', '2018'] if args.year.startswith('all') else [ args.year ] ):
            args.lumi = args.lumi + checkDict( ( 'JetHT' if args.selection.startswith('dijet') else 'SingleMuon' ), dictSamples )[iy]['lumi']

        if args.selection.startswith(('_W', '_top')):
            for isam in dictSamples:
                if isam.startswith(('ST', 'W', 'Z', 'TTTo2L2Nu', 'QCD_MuEnriched')):
                    bkgFiles[isam.split('_Tune')[0]] = [
                                    ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[args.year]['skimmerHisto'] ),
                                    checkDict( isam, dictSamples )]
            #print ("Background Files:", bkgFiles)                        

        if args.selection.startswith('_dijet'):
            for isam in dictSamples:
                if 'MuEnriched' in isam: continue
                if not checkDict( isam, dictSamples )[args.year]['skimmerHisto'].endswith('root'): continue
                if isam.startswith( signalLabelBegin ):
                    sigFiles[isam.split('_Tune')[0]] = [
                                    ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[args.year]['skimmerHisto'] ),
                                    checkDict( isam, dictSamples )
                                ]
                if isam.startswith( altSignalLabelBegin ) and not ( altSignalLabel.startswith( signalLabel ) ):
                    sigFiles[isam.split('_Tune')[0]] = [
                                    ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[args.year]['skimmerHisto'] ),
                                    checkDict( isam, dictSamples )
                                ]
        else:
            sigFiles['TTToSemiLeptonic'] = [
                            ROOT.TFile.Open( args.inputFolder+checkDict( 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8', dictSamples )[args.year]['skimmerHisto'] ),
                            checkDict( 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8', dictSamples )
                        ]
            sigFiles['TTJets'] = [
                            ROOT.TFile.Open( args.inputFolder+checkDict( 'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8', dictSamples )[args.year]['skimmerHisto'] ),
                            checkDict( 'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8', dictSamples )
                        ]

            for isam in dictSamples:
                if isam.startswith( varSignalLabelBegin ) :
                    sigFiles[isam.split('_13TeV')[0]] = [
                                        ROOT.TFile.Open( args.inputFolder+checkDict( isam, dictSamples )[args.year]['skimmerHisto'] ),
                                        checkDict( isam, dictSamples )
                                    ]

            #print ("Signal files:", sigFiles.keys())
            #print ("Background files:", bkgFiles.keys())

    #print ( dataFile, sigFiles, bkgFiles, variables, args.selection, sysUncert )

    p = Process( target=runTUnfold, args=( dataFile, sigFiles, bkgFiles, variables, args.selection, sysUncert ) )
    p.start()
    p.join()


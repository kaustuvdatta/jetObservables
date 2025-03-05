from collections import OrderedDict
import copy
import pprint 
import ROOT
import numpy as np
import array
from array import array
import bisect
import scipy
from scipy import stats
import gc

#from legend import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
####gReset()
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

ROOT.TH1.StatOverflows(ROOT.kTRUE)
ROOT.TH2.StatOverflows(ROOT.kTRUE)

#from root_numpy import array2hist, hist2array
#import histoHelpers
from histoHelpers import *#th2_to_ndarray
import os
import glob
import sys
import math

sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
import tdrstyle as tdrstyle
ROOT.gROOT.ForceStyle()
#tdrstyle.setTDRStyle()

lumi=0.
canvas = {}
textBox=ROOT.TLatex()
textBox.SetTextSize(0.10)
textBox.SetTextAlign(12)

##########################################################################
######################## Helpers for unfoldings ##########################
##########################################################################
def bottomLineTest( ivar, dataHisto, dataHistoLabel, MCHisto, covMatrix, varInfo, outputLabel, outputDir, rebin=1,
                    ext='png', version='Vn', selection='_dijet', process='data', no_null_bins=True):#, ignore_UF=True):
    #based on https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/my_unfolder.py#L2501
    
    print(f"Data nbins: {dataHisto.GetNbinsX()}",f"MC nbins: {MCHisto.GetNbinsX()}",f"cov nbinsX: {covMatrix.GetNbinsX()}",f"cov nbinsY: {covMatrix.GetNbinsY()}")
    if (dataHisto.GetNbinsX()!=MCHisto.GetNbinsX()):
        print("!!! ERROR:(dataHisto.GetNbinsX()!=MCHisto.GetNbinsX()) !!!" )
        dataBins=[]
        MCbins=[]
        for i in range(1,dataHisto.GetNbinsX()+2):#!=MCHisto.GetNbinsX()):
            dataBins.append(dataHisto.GetBinLowEdge(i))
        for i in range(1,MCHisto.GetNbinsX()+2):
            MCbins.append(MCHisto.GetBinLowEdge(i))
        
        print("Data bins", dataBins )
        print("MC bins", MCbins )
            
    #print(dataHisto.GetNbinsX(),MCHisto.GetNbinsX(),covMatrix.GetNbinsX(),covMatrix.GetNbinsY())
    
    
    if not(rebin==1):
        MCHisto.Rebin( rebin )  
        dataHisto.Rebin( rebin )
        covMatrix.Rebin2D( rebin, rebin )
    
    #if isinstance(covMatrix,ROOT.TH2): 
    cov_arr,_ = th2_to_ndarray(covMatrix.Clone())
    #if isinstance(dataHisto,ROOT.TH1): 
    data_arr,_ = th1_to_ndarray(dataHisto.Clone())
    #if isinstance(MCHisto,ROOT.TH1): 
    mc_arr,_ = th1_to_ndarray(MCHisto.Clone())
    print(cov_arr.shape,data_arr.shape,mc_arr.shape)
    
    if no_null_bins:
        
        null_bins = get_null_bins(cov_arr)
        print(f'null bins {dataHistoLabel} space:', null_bins)
        
        cov_arr = remove_null_bins(cov_arr, null_bins)
        mc_arr = remove_null_bins(mc_arr, null_bins)
        data_arr = remove_null_bins(data_arr, null_bins)    
    
    try:
        delta = data_arr - mc_arr
        #print(f'Data b.c.: {data_arr}, MC b.c.: {mc_arr}, Residuals:{delta}')        
    except ValueError:
        print(f'Data b.c.: {data_arr}, MC b.c.: {mc_arr}')
        print("chi2 cannot be calculated since something is wrong with the data mc event th2->ndarrays, please check what's going on")
        return 1,1,1
    ##### computing chi2 and inverted matrix
    try: 
        v_inv = np.linalg.inv(cov_arr)
    except np.linalg.LinAlgError:
        print("Using pseudo-inverse instead since true inv. operation via np.linalg.inv() failed")
        v_inv = np.linalg.pinv(cov_arr, rcond=1E-30)
        
    inter = v_inv.dot(delta.T)
    chi2 = delta.dot(inter)[0][0]
    
    data_nonzero = len([i for i in range(data_arr[0].shape[0]) if data_arr[0][i]!=0])
    mc_nonzero = len([i for i in range(mc_arr[0].shape[0]) if mc_arr[0][i]!=0])
    
    
    ndof = max(data_nonzero,mc_nonzero)#delta.shape[1]# # only consider n bins where at least one has data - if both 0, don't count it
    print(1.-scipy.stats.chi2.cdf(chi2, int(ndof)))
    
    p = 1.-scipy.special.gammainc(chi2/2.,ndof/2.)#1.-scipy.stats.chi2.cdf(chi2, int(ndof))
        
    print(f'chi2 for {dataHistoLabel}, ndf, p, chi2/ndf = ', np.round(chi2,5), ndof, p, np.round(chi2/ndof,5))
    
    return chi2, ndof, p

def calc_chi2_stats(one_hist, other_hist, cov_matrix):
    one_vec = one_hist #, one_err = th1_to_ndarray(one_hist, False)
    # print(one_err)
    other_vec = other_hist#, _ = th1_to_ndarray(other_hist, False)
    delta = one_vec - other_vec
    if isinstance(cov_matrix, ROOT.TH2):
        v, _ = th2_to_ndarray(cov_matrix)
    else:
        v = cov_matrix
    # print("delta:", delta)
    # v = np.diag(np.diag(v))  
    # print("v:", v)
    try:
        v_inv = np.linalg.inv(v)
    except np.linalg.LinAlgError:
        print("Trying pseudo-inverse instead")
        v_inv = np.linalg.pinv(v, rcond=1E-30)
    inter = v_inv.dot(delta.T)
    # print("parts:", delta * inter.T)
    chi2 = delta.dot(inter)[0][0]
    ndof = delta.shape[1]
    p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
    return chi2, ndof, p


def fold_generator_level(hist_truth, probability_matrix, bins_reco, oflow=False):
    
    # Convert ROOT TH1 to vector
    gen_vec, gen_vec_err = th1_to_ndarray(hist_truth, oflow_x=oflow)
    
    #convert TH2 to ndarray for PM
    probaM, _ = th2_to_ndarray(probability_matrix, oflow)

    # Multiply array (vector) with PM (transpose from row vec to column vec as necessary)
    folded_vec = probaM.dot(gen_vec.T)

    # Convert array to TH1
    folded_mc_truth = ndarray_to_th1(folded_vec.T, has_oflow_x=oflow, offset=0., bins=bins_reco)

    # Err. prop.: if y = Ax, with covariance matrices Vyy and Vxx,respectively, then Vyy = (A*Vxx)*A^T
    vxx, _ = th2_to_ndarray(make_diag_cov_hist_from_errors(hist_truth, inverse=False), oflow)
    result = probaM.dot(vxx)
    
    folded_covariance = result.dot(probaM.T)
    print(folded_covariance.shape)
    folded_errors = make_hist_from_diagonal_errors(folded_covariance,bins=bins_reco)
    update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_mc_truth)
    return folded_mc_truth


def DoUnfolding(Response,Reco):
    tunfolder = ROOT.TUnfoldDensity(Response,
                                    ROOT.TUnfold.kHistMapOutputHoriz,
                                    ROOT.TUnfold.kRegModeCurvature, 
                                    ROOT.TUnfold.kEConstraintNone, 
                                    ROOT.TUnfoldDensity.kDensityModeBinWidth)
    tunfolder.SetInput(Reco)
    tunfolder.DoUnfold(0.)
    return tunfolder.GetOutput("MC_unfolded")

def get_folded_unfolded(folded, unfolded, cov_tot, probaM, oflow=True):
    # don't use getfoldedoutput, because it doesn't have the updated errors from the total error matrix
    # so we'll have to do it ourselves
    # 1. Make unfolded hist into TVector/TMatrix

    # 2. Make response 2d hist into matrix

    # 3. Multiply the two, convert to TH1

    # Get the TUnfold one for reference, although its errors will be wrong
    
    #print(folded.GetNbinsX(),unfolded.GetNbinsX())
    
    
    #probability matrix is simply the normalised RM, 
    #so need to consider that there is an UF bin on the y-axis for miss(gen)-corrections
    probaM, _ = th2_to_ndarray(probaM, oflow_x=False, oflow_y=False, uflow_x=False, uflow_y=False)
    #print(probaM.shape, "with UF on y axis")
    
    folded_unfolded_tunfold = folded.Clone()

    # Get unfolded results as array
    unfolded_vector, _ = th1_to_ndarray(unfolded, oflow_x=False)
    #print(unfolded_vector.shape)
    
    # Multiply
    # Note that we need to transpose from row vec to column vec
    folded_vec = probaM.dot(unfolded_vector.T)
    #print(folded_vec.shape)
    
    
    bins = array('d',[folded.GetBinLowEdge(i) for i in range(1,folded.GetNbinsX()+2)])

    
    # Convert vector to TH1
    folded_unfolded = ndarray_to_th1(folded_vec.T, has_oflow_x=False, offset=0., bins=bins)

    # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
    # respectively, then Vyy = (A*Vxx)*A^T
    unfolded_covariance_matrix, _ = th2_to_ndarray(cov_tot , oflow_x=False, oflow_y=False, uflow_x=False, uflow_y=False )
    result = probaM.dot(unfolded_covariance_matrix)
    #print(result.shape)
    
    folded_covariance = result.dot(probaM.T)
    #print(folded_covariance.shape)
    
    
    folded_errors = make_hist_from_diagonal_errors(folded_covariance,bins=bins)
    
    #print(folded_errors.GetBinLowEdge(4))
    #print(folded_unfolded.GetBinLowEdge(4))
    
    update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_unfolded)

    return folded_unfolded.Clone(folded.GetName()+'errorFixed')

def CrossClosure(response1,reco1,response2,reco2):
    unf11=DoUnfolding(response1,reco1)
    unf12=DoUnfolding(response2,reco1)
    unf21=DoUnfolding(response1,reco2)
    unf22=DoUnfolding(response2,reco2)
    return unf11,unf12,unf21,unf22

def SelfClosure(response1,reco2,response2,reco1):
    unf21=DoUnfolding(response1,reco2)
    unf12=DoUnfolding(response2,reco1)
    return unf21,unf12




#############################################################################
###################### Helpers for plotting unfoldings ######################
#############################################################################

def plotSimpleComparison( inFile1, sample, inFile2, sample2, name, rebinX=1, xmin='', xmax='', labX=0.92, labY=0.50, axisX='', axisY='', log=False, ext='png', Norm=False, version='', outputDir='Plots/' ):
    """"Take two root files, make simple comparison plot"""

    outputFileName = name+'_'+sample+sample2+'_simpleComparisonPlot'+version+'.'+ext
    #print('Processing.......', outputFileName)

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
    
    if Norm:
        if 'data' in histo.GetName().lower():
            norm = histo.Integral()#Clone()
        elif 'data' in histo2.GetName().lower():
            norm = histo2.Integral()#Clone()
        else:
            norm=1.
    else:
        norm=1.
        

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
    histo.SetLineColorAlpha(ROOT.kRed,0.7)
    histo.SetMarkerColor(ROOT.kRed)
    histo.SetLineWidth(2)
    histo2.SetLineColor(ROOT.kBlue)
    #if Norm:
    histo.Scale(1./norm,'width')
    histo.Draw('histE')
    histo2.Scale(1./norm,'width')
    histo2.Draw('histE same')
    #else:
    ##    histo.Scale(1./norm,'width')
    #    histo2.Scale(1./norm,'width')

    #    histo.Draw('histe')
    #    histo2.Draw('histe same')
    if not axisY: histo.GetYaxis().SetTitle( 'Events / '+str(binWidth) )
    if axisX: histo.GetXaxis().SetTitle( axisX )

    legend.Draw()

    canvas[name].SaveAs( outputDir+outName )
    if ext.startswith('pdf'):
        canvas[name].SaveAs( outputDir+outName.replace('pdf', 'png') )

def plotSysComparison( nomHisto, dictUncHistos, outputName, labelX='', log=False, version='', ext='png', year='2017', outputDir='Plots/' ): #from Alejandro

    colors = [ 2, 4,  9, 8, 28, 30, 42, 13, 12, 40, 46, 3, 24, 26, 219, 92, 48, 49, 37, 38, 33, 17, 50, 205, 225, 94, 221, 16,  225, 128]
    
    
    outputFileName = outputName+'_'+version+'.'+ext
    #print ('Processing.......', outputFileName)

    binWidth = nomHisto.GetBinWidth(1)

    legend=ROOT.TLegend(0.35,0.6,0.80,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.02)
    legend.SetBorderSize(0)

    multiGraph = ROOT.TMultiGraph()
    gnom = ROOT.TGraphAsymmErrors()
    gnom.Divide( nomHisto.Clone(), nomHisto.Clone(), 'pois' )
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetMarkerStyle(1)
    gnom.SetLineWidth(2)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    multiGraph.Add( gnom )

    dictgraph = {}
    col_counter=0
    
    #print (dictUncHistos)
    for ih in dictUncHistos:
        #print(ih, col_counter, len(colors))
        dictgraph[ih] = ROOT.TGraphAsymmErrors()
        dictgraph[ih].Divide( dictUncHistos[ih], nomHisto, 'pois' )
        if not col_counter==len(colors)-1:
            dictgraph[ih].SetLineColor( colors[col_counter] )
            dictgraph[ih].SetLineStyle( 2 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        else:
            col_counter = col_counter-len(colors)+2
            dictgraph[ih].SetLineColor( colors[col_counter] )
            dictgraph[ih].SetLineStyle( 3 )
            if 'jes' in ih: dictgraph[ih].SetLineStyle( 1 )
            dictgraph[ih].SetMarkerStyle(1)
            dictgraph[ih].SetLineWidth( 1 )
        if 'jes' in ih and ('2017' in ih or '2018' in ih): legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        else: legend.AddEntry( dictgraph[ih], ih.split('_')[1] , 'l' )
        multiGraph.Add( dictgraph[ih] )
        col_counter=col_counter+1

    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 750, 500 )
    if log: canUnc.SetLogy()
    multiGraph.GetYaxis().SetTitle( 'Ratio Unc/Nominal' )
    multiGraph.GetXaxis().SetTitle( labelX )
    multiGraph.SetMaximum( 3. )
    multiGraph.SetMinimum( -1.)
    multiGraph.Draw('ALP')

    CMS_lumi.cmsTextOffset = 0.0
    CMS_lumi.relPosX = 0.13
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+year
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    if ext.startswith('pdf'):
        canUnc.SaveAs( outputDir + outputFileName.replace('pdf', 'png') )
    del canUnc

    
def drawDataMCReco( ivar, selection, year, lumi, process,
                    dataJetHisto, nominal_recoJetHisto, 
                    alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto,
                    labelX, jetType, maxX, tlegendAlignment, outputName,log=False ):
    """docstring for drawDataMCReco (dijets)"""
    print ("Drawing Data/MC")
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    #can.cd()
    if log==True:
        #print(ivar,outputName)
        outputName=outputName.replace("DataMC_","DataMC_logPlot_")
        #dataJetHisto.Rebin(2)
        #nominal_recoJetHisto.Rebin(2)
        #alt0_recoJetHisto.Rebin(2)
        #alt1_recoJetHisto.Rebin(2)
        #alt2_recoJetHisto.Rebin(2)
        #can.SetLogy()
        #print(ivar,outputName)
        
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    pad1.Draw()
    pad2.Draw()
    can.cd()
    pad1.cd()
    
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    
    if tlegendAlignment.startswith('right'): 
        legend=ROOT.TLegend(0.64,0.6,0.9,0.79)
        #legend=ROOT.TLegend(0.65,0.6,0.90,0.8)

    else: 
        #legend=ROOT.TLegend(0.20,0.6,0.40,0.8)
        legend=ROOT.TLegend(0.19,0.6,0.40,0.79)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03 )
    legend.SetBorderSize(0)
    print(dataJetHisto,nominal_recoJetHisto,alt0_recoJetHisto,alt1_recoJetHisto,alt2_recoJetHisto)
    
    dataHisto = dataJetHisto.Clone()#normalise_hist
    recoHisto = nominal_recoJetHisto.Clone()#normalise_hist
    alt0recoHisto = alt0_recoJetHisto.Clone()#normalise_hist
    alt1recoHisto = alt1_recoJetHisto.Clone()#normalise_hist
    alt2recoHisto = alt2_recoJetHisto.Clone()#normalise_hist
    
    if 'tau' in ivar: dataHisto.Scale(1, 'width')  ### divide by bin width
    #dataHisto.Scale(1/dataHisto.Integral(), 'width')  ### divide by bin width
    dataHisto.SetMarkerStyle(8)
    dataHisto.SetMarkerSize(2)
    dataHisto.SetMarkerColor(ROOT.kBlack)
    dataHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( dataHisto, 'Data', 'pe' )
    
    if 'tau' in ivar: recoHisto.Scale(1, 'width')
    #genJetHisto.Scale(scaleFactor)
    #genJetHisto.Scale(1/genJetHisto.Integral(), 'width')  ### divide by bin width
    recoHisto.SetLineWidth(1)
    recoHisto.SetLineColor(colors[0])
    recoHisto.SetMarkerColor(colors[0])
    recoHisto.SetMarkerStyle(25)
    recoHisto.SetMarkerSize(2)
    legend.AddEntry( recoHisto, 'MG5-MLM+Pythia8', 'lp' )

    print(labelX)

    if 'tau' in ivar or '#' in labelX:
        dataHisto.GetYaxis().SetTitle( 'Events')#'#frac{dN}{d#'+labelX.split('#')[1]+'}' )##frac{1}{dN} , +'    [A.U.]'
        #print( '#frac{dN}{d#'+labelX.split('#')[1]+'}')
    else:
        bw = np.round(recoHisto.GetXaxis().GetBinLowEdge(2)-recoHisto.GetXaxis().GetBinLowEdge(1),3)
        label=None
        if 'pt' in ivar:
            label = 'p_{T}'
        elif 'mass'in ivar:
            label = 'm'
        elif 'softdrop' in ivar:
            label = 'm_{SD}'
        elif 'y' in ivar:
            label = 'y'
        elif 'npv' in ivar.lower():
            label= 'N_{PV}'
        else:
            pass
        #if label:
        #    dataHisto.GetYaxis().SetTitle( '#frac{dN}{d'+label+f'}} [Events/{bw} GeV]')# if 'pt' in ivar or 'mass' in ivar or 'softdrop' in ivar else '') ) 
        #else:
        dataHisto.GetYaxis().SetTitle( f' Events/{bw} '+ ( '[GeV]' if 'pt' in ivar or 'mass' in ivar or 'softdrop' in ivar else '') ) 
    
    dataHisto.GetYaxis().SetTitleOffset(1.15)
    dataHisto.GetYaxis().CenterTitle()
    dataHisto.GetYaxis().SetTitleSize(0.05)
    dataHisto.GetYaxis().SetLabelSize(0.05)
    dataHisto.GetXaxis().SetTitleSize(0.0)
    dataHisto.GetXaxis().SetLabelSize(0.0)
    dataHisto.GetXaxis().SetTickLength(0.)
    
    dataHisto.SetMaximum( 1.8*max([ recoHisto.GetMaximum(), dataHisto.GetMaximum()] ) if not('pt') in ivar else 40.*max([ recoHisto.GetMaximum(), dataHisto.GetMaximum()] )  )
    dataHisto.SetMinimum(0. if not log else 0.01)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), dataHisto.GetMaximum()] ) )
    
    dataHisto.SetTitle('')
    can.SetTitle('')
    ROOT.TGaxis.SetMaxDigits(4)#,'y')
    #ROOT.TGaxis.SetExponentOffset(-5,0,'y')
        
    dataHisto.Draw( "E1")

    
    #alt1recoHisto.Scale(1, 'width')  ### divide by bin width
    #alt1recoHisto.SetLineWidth(2)
    #alt1recoHisto.SetLineColor(ROOT.kCyan+3)
    #alt1recoHisto.SetMarkerColor(ROOT.kCyan+3)
    #alt1recoHisto.SetMarkerStyle(25)
    #alt1recoHisto.SetMarkerSize(2)
    #legend.AddEntry( alt1recoHisto, 'MG5-MLM+Pythia8', 'lp' )
    #alt1recoHisto.Draw("histe1 same")
    recoHisto.Draw( "histe1 same")

    if 'tau' in ivar: alt0recoHisto.Scale(1, 'width')  ### divide by bin width
    alt0recoHisto.SetLineWidth(1)
    alt0recoHisto.SetLineColor(colors[1])#ROOT.kBlue)
    alt0recoHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    alt0recoHisto.SetMarkerStyle(25)
    alt0recoHisto.SetMarkerSize(2)
    legend.AddEntry( alt0recoHisto, 'MG5-MLM+Herwig7', 'lp' )
    alt0recoHisto.Draw("histe1 same")
    
    
    if 'tau' in ivar: alt2recoHisto.Scale(1, 'width')  ### divide by bin width
    alt2recoHisto.SetLineWidth(1)
    alt2recoHisto.SetLineColor(colors[2])#ROOT.kGray+4)
    alt2recoHisto.SetMarkerColor(colors[2])#ROOT.kGray+4)
    alt2recoHisto.SetMarkerStyle(25)
    alt2recoHisto.SetMarkerSize(2)
    legend.AddEntry( alt2recoHisto, 'Pythia8', 'lp' )
    alt2recoHisto.Draw("histe1 same")
    if log: ROOT.gPad.SetLogy()
    else: ROOT.gPad.SetLogy(0)
    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = ( 'Central Dijet' if 'Central' in jetType  else 'Forward Dijet' )#+' dijet region'
    elif selection.startswith("_W"): seltext = 'Boosted W region'
    elif selection.startswith("_top"): seltext = 'Boosted top region'
    selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.88, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.04)

    selText.SetNDC()
    
    if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_dijet") and 'Forward' in jetType : seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )

    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ("#leq" if selection.startswith("_dijet") else " ")+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+", " +( '2016+2017+2018' if year.startswith('all') else year ) 
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.08
    CMS_lumi.CMS_lumi(pad1, 4, 10 if tlegendAlignment.startswith('right') else 3) 
    
    can.cd()
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    
    tmpPad2= pad2.DrawFrame( recoHisto.GetXaxis().GetBinLowEdge(1), 0., maxX, 1.9 )
    #print (labelX)
    tmpPad2.GetYaxis().SetTitle( "Data/Sim." )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(0.1, 2.1 )
    tmpPad2.GetYaxis().CenterTitle()
    
       
    if 'tau' in ivar: tmpPad2.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    else: tmpPad2.GetXaxis().SetTitle(labelX)
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.1, 'y')
    tmpPad2.SetTitleSize(0.1, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_nominal = ROOT.TGraphAsymmErrors()#len(l_bins)-1,x_bins,y_vals)
    ratio_nominal.Divide( dataHisto,recoHisto, 'pois' )
    ratio_nominal.SetLineColor(colors[0])
    ratio_nominal.SetMarkerColor(colors[0])
    ratio_nominal.SetMarkerSize(2)
    ratio_nominal.GetXaxis().SetNdivisions(505)
    #ratio_nominal.GetYaxis().SetNdivisions(505)
    ratio_nominal.SetMarkerStyle(25)
    ratio_nominal.Draw('PE1 ')
    
    ratio_alt0MC = ROOT.TGraphAsymmErrors()
    ratio_alt0MC.Divide(  dataHisto, alt0recoHisto,'pois' )
    ratio_alt0MC.SetLineColor(colors[1])
    ratio_alt0MC.SetMarkerColor(colors[1])
    ratio_alt0MC.SetMarkerStyle(25)
    ratio_alt0MC.SetMarkerSize(2)
    ratio_alt0MC.Draw('PE1 same')
    
    #ratio_alt1MC = ROOT.TGraphAsymmErrors()
    #ratio_alt1MC.Divide(  dataHisto,alt1recoHisto, 'pois' )
    #ratio_alt1MC.SetLineColor(ROOT.kCyan+3)
    #ratio_alt1MC.SetMarkerColor(ROOT.kCyan+3)
    #ratio_alt1MC.SetMarkerStyle(25)
    #ratio_alt1MC.SetMarkerSize(2)
    #ratio_alt1MC.Draw('PE1 same')
    
    ratio_alt2MC = ROOT.TGraphAsymmErrors()
    ratio_alt2MC.Divide(  dataHisto, alt2recoHisto, 'pois' )
    ratio_alt2MC.SetLineColor(colors[2])#ROOT.kGray+4)
    ratio_alt2MC.SetMarkerColor(colors[2])#ROOT.kGray+4)
    ratio_alt2MC.SetMarkerStyle(25)
    ratio_alt2MC.SetMarkerSize(2)
    ratio_alt2MC.Draw('PE1 same')
    
    ratioLegend=ROOT.TLegend(0.20,0.85,0.8,0.95)
    ratioLegend.SetTextSize(0.06)
    ratioLegend.SetNColumns(4)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_nominal, 'MG5-MLM+P8', 'lp' )
    ratioLegend.AddEntry( ratio_alt0MC, 'MG5-MLM+H7', 'lp' )
    #ratioLegend.AddEntry( ratio_alt1MC, 'MG5-MLM+P8', 'lp' )
    ratioLegend.AddEntry( ratio_alt2MC, 'P8+P8', 'lp' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    pad2.Update()
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    print(outputName,png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)
                    


    
def drawUnfold(ivar, selection, process, year, lumi,
               dataJetHisto, genJetHisto, 
               unfoldHistowoUnc, altMCHisto, foldHisto, recoJetHisto,
               cov_tot, cov_datastat_tot, labelX, maxX, tlegendAlignment, outputName,
               altMC1Histo = None, altMC2Histo = None, altMC1Histo_label = None, altMC2Histo_label = None, 
               nomMCHisto_label = None, altMCHisto_label = None,
               extraMC=False, includeFSR = False, fsrUpHisto = None, fsrDownHisto=None, noNorm=False
              ):
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.60,0.61,0.89,0.89)

    else: legend=ROOT.TLegend(0.16,0.61,0.45,0.89)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)

    
    unfoldHisto = unfoldHistowoUnc.Clone('unfoldHisto'+ivar)
    unfoldHistoStatUnc = unfoldHistowoUnc.Clone('unfoldHistoStatUnc'+ivar)
    
    dataJetHisto.SetTitle("")
    print("data(minus bkgs).Integral()",dataJetHisto.Integral())
    
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    unfoldHistoStatUnc.SetTitle("")
    
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    
    recoJetHisto.SetTitle("")
    print("(RM proj.Y )recoJetHisto.Integral()",recoJetHisto.Integral())
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())
    
    
    dataScaling = unfoldHisto.Integral()
    print (dataScaling)
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the transformed covariance matrix, from the unnormalised -> the normalised space
    
    cov_normTot_np, normed_covTot = get_normalised_cov(unfoldHisto.Clone(), 
                                                       cov_tot.Clone())
    cov_norm_dataStat_np, normed_cov_dataStat = get_normalised_cov(unfoldHistoStatUnc.Clone(), 
                                                                   cov_datastat_tot.Clone())
    
    
    
    unfoldHistoDataStatErr=unfoldHistoStatUnc.Clone()

    unfoldHistoDataStatErr.Sumw2()
    unfoldHisto.Sumw2()
    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()

    unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')

    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.))
    if not(noNorm): get_th1_normedCovErrors(unfoldHistoDataStatErr, cov_norm_dataStat_np)
    unfoldHistoDataStatErr.Scale(1.,'width')

    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.))
    if not(noNorm): get_th1_normedCovErrors(unfoldHisto, cov_normTot_np)
    unfoldHisto.Scale(1.,'width')

    
    dataJetHisto.Scale(1./(dataJetHisto.Integral() if not(noNorm) else 1.),'width')
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.),'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.),'width')
    foldHisto.Scale(1./(foldHisto.Integral() if not(noNorm) else 1.),'width')
    recoJetHisto.Scale(1./(recoJetHisto.Integral() if not(noNorm) else 1.),'width')
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        fsrUpHisto.Scale(1./(fsrUpHisto.Integral() if not(noNorm) else 1.),'width')
        fsrDownHisto.Sumw2()
        fsrDownHisto.Scale(1./(fsrDownHisto.Integral() if not(noNorm) else 1.),'width')
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        altMC1Histo.Scale(1./(altMC1Histo.Integral() if not(noNorm) else 1.),'width')
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo.Scale(1./(altMC2Histo.Integral() if not(noNorm) else 1.),'width')
            
            altMC2Histo.SetTitle("")

    
    
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(2)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(2)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(2)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(22)


        fsrDownHisto.SetMarkerSize(2)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(23)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )

   
    if 'tau' in labelX: 
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( (1.6 if '21' in ivar or '32' in ivar else 1.56)*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )

    unfoldHisto.Draw( "E1")

    #altMCHisto.Scale(1, 'width')  ### divide by bin width
    altMCHisto.SetLineWidth(2)
    altMCHisto.SetMarkerSize(2)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR up", 'pe')

        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR down", 'pe')
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'PWHG+H7','lpe')#
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            #altMC2Histo.Scale(1, 'width')  ### divide by bin width
            altMC2Histo.SetLineWidth(2)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(2)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            #altMC1Histo.Scale(1, 'width')  ### divide by bin width
            altMC1Histo.SetLineWidth(2)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(2)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            
            legend.AddEntry( altMC1Histo, altMC1Histo_label.replace('-FXFX',''),'lpe')#'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1 same")

        
    genJetHisto.Draw( "histE1 same")
    altMCHisto.Draw("histE1 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1 same")
        fsrDownHisto.Draw("PE1 same")

    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.044)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.80, seltext )
    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        if year=='all': 
            if 'dijet' in selection:
                CMS_lumi.lumi_13TeV = ('#leq 135' if 'dijet' in selection else '138')+" fb^{-1} (13 TeV)"+('' if year.startswith('all') else ", "+( '' if year.startswith('all') else year ) )
        else:
            CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistoDataStatErr.Clone()
    ratio_datastatUnc.Divide(unfoldHistowoUnc)
    ratio_totalUnc = unfoldHisto.Clone()
    ratio_totalUnc.Divide(unfoldHistowoUnc)
    
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    print (labelX)
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.50 )
    #tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.13, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )

    ratio_totalUnc.GetYaxis().CenterTitle()
    ratio_totalUnc.GetXaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetXaxis().SetTitleSize(0.13)

    ratio_totalUnc.GetYaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetYaxis().SetTitleSize(0.12)
    ratio_totalUnc.GetXaxis().SetNdivisions(505)
    ratio_totalUnc.GetYaxis().SetNdivisions(505)
    
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)

    ratio_totalUnc.SetFillColorAlpha(14,0.8)
    ratio_totalUnc.SetLineColor(14)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3354)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.5,0.5)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerStyle(22)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerStyle(23)
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
    
    hRatio.SetMarkerSize(2)
    hRatio.Draw('P0 same')
    
    hRatio2.SetMarkerSize(2)
    hRatio2.Draw('P0 same')
    
    hRatio5.SetMarkerSize(2)
    hRatio5.Draw('P0 same')
    
    if includeFSR:
        hRatio3.SetMarkerSize(2)
        hRatio3.Draw('P0 same')

        hRatio4.SetMarkerSize(2)
        hRatio4.Draw('P0 same')
    
    
    ratioLegend=ROOT.TLegend(0.15,0.85,0.7,0.95)
    ratioLegend.SetTextSize(0.088)
    ratioLegend.SetNColumns(3)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_totalUnc, 'Data total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)    

def drawUnfold_normedCovErr(   ivar, selection, process, year, lumi,
                               dataJetHisto, genJetHisto, unfoldHisto, unfoldHistoStatUnc, unfoldHistowoUnc, altMCHisto, foldHisto, recoJetHisto,
                               cov_tot, cov_datastat_tot, labelX, maxX, tlegendAlignment, outputName,
                               altMC1Histo = None, altMC2Histo = None, altMC1Histo_label = None, altMC2Histo_label = None, 
                               nomMCHisto_label = None, altMCHisto_label = None,
                               extraMC=False, includeFSR = False, fsrUpHisto = None, fsrDownHisto=None, noNorm=False
                              ):
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    
    dataJetHisto.SetTitle("")
    print("data(minus bkgs).Integral()",dataJetHisto.Integral())
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    unfoldHistoStatUnc.SetTitle("")
    #print("unfoldHistoStatUnc.Integral()",unfoldHistoStatUnc.Integral())
    #unfoldHistowoUnc.SetTitle("")
    #print("unfoldHistowoUnc.Integral()",unfoldHistowoUnc.Integral())
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    #foldHisto.SetTitle("")
    #print("foldHisto.Integral()",foldHisto.Integral())
    recoJetHisto.SetTitle("")
    print("(RM proj.Y )recoJetHisto.Integral()",recoJetHisto.Integral())
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())

            
    
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.68,0.61,0.90,0.89)

    else: legend=ROOT.TLegend(0.16,0.61,0.38,0.89)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    
    #bins = variables[ivar]['bins']

    unfoldHistoTot = unfoldHisto.Clone()
    dataScaling = unfoldHisto.Integral()
    
    print (dataScaling)
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the covariance matrix, from the normalise -> the unnormalised space
    #normed_cov_tot_matrix, normed_cov_tot = GetNormalizedTMatrixandTH2(cov_tot.Clone(),"normed_cov_tot", unfoldHisto.Clone())
    
    #normed_cov_datastat_tot_matrix, normed_cov_datastat_tot = GetNormalizedTMatrixandTH2(cov_datastat_tot.Clone(),"normed_cov_dastat_tot", unfoldHisto.Clone())
    unfoldHistoDataStatErr=unfoldHistoStatUnc.Clone()
    unfoldHistoDataStatErr.Sumw2()
    unfoldHisto.Sumw2()
    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()
    
    unfoldHistoDataStatErr.Scale(1./(unfoldHistoDataStatErr.Integral() if not(noNorm) else 1.),'width')
    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not(noNorm) else 1.),'width')
    dataJetHisto.Scale(1./(dataJetHisto.Integral() if not(noNorm) else 1.),'width')
    genJetHisto.Scale(1./(genJetHisto.Integral() if not(noNorm) else 1.),'width')
    unfoldHistowoUnc.Scale(1./(unfoldHistowoUnc.Integral() if not(noNorm) else 1.),'width')
    altMCHisto.Scale(1./(altMCHisto.Integral() if not(noNorm) else 1.),'width')
    foldHisto.Scale(1./(foldHisto.Integral() if not(noNorm) else 1.),'width')
    recoJetHisto.Scale(1./(recoJetHisto.Integral() if not(noNorm) else 1.),'width')
    
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        fsrUpHisto.Scale(1./(fsrUpHisto.Integral() if not(noNorm) else 1.),'width')
        fsrDownHisto.Sumw2()
        fsrDownHisto.Scale(1./(fsrDownHisto.Integral() if not(noNorm) else 1.),'width')
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        altMC1Histo.Scale(1./(altMC1Histo.Integral() if not(noNorm) else 1.),'width')
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection and altMC2Histo:
            altMC2Histo.Sumw2()
            altMC2Histo.Scale(1./(altMC2Histo.Integral() if not(noNorm) else 1.),'width')
            
            altMC2Histo.SetTitle("")

    
    
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(2)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(2)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(2)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(22)


        fsrDownHisto.SetMarkerSize(2)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(23)
    
    legend.AddEntry( genJetHisto, nomMCHisto_label, 'lpe' )

   
    if 'tau' in labelX: 
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( (1.6 if '21' in ivar or '32' in ivar else 1.56)*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )

    unfoldHisto.Draw( "E1")

    #altMCHisto.Scale(1, 'width')  ### divide by bin width
    altMCHisto.SetLineWidth(2)
    altMCHisto.SetMarkerSize(2)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR up", 'pe')

        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR down", 'pe')
        
    legend.AddEntry( altMCHisto, altMCHisto_label, 'lp' )#'PWHG+H7','lpe')#
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            #altMC2Histo.Scale(1, 'width')  ### divide by bin width
            altMC2Histo.SetLineWidth(2)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(2)
            
            legend.AddEntry( altMC2Histo, altMC2Histo_label, 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            #altMC1Histo.Scale(1, 'width')  ### divide by bin width
            altMC1Histo.SetLineWidth(2)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(2)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            legend.AddEntry( altMC1Histo, altMC1Histo_label,'lpe')#'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1 same")

        
    genJetHisto.Draw( "histE1 same")
    altMCHisto.Draw("histE1 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1 same")
        fsrDownHisto.Draw("PE1 same")

    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.040)

    selText.SetNDC()
    
    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.80, seltext )
    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistoDataStatErr.Clone()
    ratio_datastatUnc.Divide(unfoldHistowoUnc)
    ratio_totalUnc = unfoldHisto.Clone()
    ratio_totalUnc.Divide(unfoldHistowoUnc)
    
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    print (labelX)
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.50 )
    #tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.13, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    ratio_totalUnc.GetXaxis().SetTitle( '#'+labelX.split('#')[1] )
    ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )

    ratio_totalUnc.GetYaxis().CenterTitle()
    ratio_totalUnc.GetXaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetXaxis().SetTitleSize(0.13)

    ratio_totalUnc.GetYaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetYaxis().SetTitleSize(0.12)
    ratio_totalUnc.GetXaxis().SetNdivisions(505)
    ratio_totalUnc.GetYaxis().SetNdivisions(505)
    
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)

    ratio_totalUnc.SetFillColorAlpha(14,0.8)
    ratio_totalUnc.SetLineColor(14)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3354)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.5,0.5)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(2)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(2)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio3.SetMarkerStyle(22)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(2)
        hRatio4.SetMarkerStyle(23)
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(2)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
    
    hRatio.SetMarkerSize(2)
    hRatio.Draw('P0 same')
    
    hRatio2.SetMarkerSize(2)
    hRatio2.Draw('P0 same')
    
    hRatio5.SetMarkerSize(2)
    hRatio5.Draw('P0 same')
    
    if includeFSR:
        hRatio3.SetMarkerSize(2)
        hRatio3.Draw('P0 same')

        hRatio4.SetMarkerSize(2)
        hRatio4.Draw('P0 same')
    
    
    ratioLegend=ROOT.TLegend(0.15,0.85,0.7,0.95)
    ratioLegend.SetTextSize(0.088)
    ratioLegend.SetNColumns(3)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_totalUnc, 'Data total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)         
    
def drawClosures(ivar, selection, process, year, lumi, genJetHisto, genJetHistoCross, unfoldHisto, unfoldHistoCross,
                 ratioUncHisto, ratiototUncHisto, ratiosystUncHisto, labelX, maxX, tlegendAlignment, 
                 outputName, nomMCHisto_label = None, altMCHisto_label = None, noNorm = False ):
    
    if process.startswith('MCCrossClosure'):
    
        genJetHistoCross.SetTitle("") 
        unfoldHistoCross.SetTitle("")
        
    else:

        genJetHisto.SetTitle("") 
        unfoldHisto.SetTitle("")

    
    """docstring for drawClosures"""
    print ("Drawing unfolding closure")
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 1500, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad1.Draw()
    
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.54,0.6,0.86,0.9)
    else: legend=ROOT.TLegend(0.20,0.65,0.45,0.9)
        
    legend.SetFillStyle(0)
    legend.SetTextSize(0.032)
    legend.SetBorderSize(0)
    
    
    if process.startswith('MCCrossClosure'): 
        print(genJetHisto.Integral(), unfoldHisto.Integral(),unfoldHistoCross.Integral(),genJetHistoCross.Integral())
    
    unfoldHisto.Scale(1./(unfoldHisto.Integral() if not noNorm else 1.),'width')
    unfoldHisto.SetMarkerStyle(4)
    unfoldHisto.SetMarkerColor(ROOT.kRed)
    unfoldHisto.SetLineColor(ROOT.kRed)
    unfoldHisto.SetLineWidth(2)
    
    legend.AddEntry( unfoldHisto, (f'{nomMCHisto_label} (self-closure)' if process.startswith('MCSelfClosure') else f'{nomMCHisto_label} unf. w/ {nomMCHisto_label}'), 'pe' )
    
    
    genJetHisto.Scale(1./(genJetHisto.Integral() if not noNorm else 1.),'width')
    genJetHisto.SetLineWidth(2)
    genJetHisto.SetLineColor(ROOT.kBlue)
    genJetHisto.SetMarkerStyle(0)
    genJetHisto.SetLineStyle(2)
    legend.AddEntry( genJetHisto, f'{nomMCHisto_label} (gen)', 'lp' )
    
    
    if 'tau' in labelX:
        unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d#'+labelX.split('#')[1]+'}' )
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
            
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.Draw()
    unfoldHisto.SetMaximum( 1.6*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    can.Update()
    can.Modified()
    unfoldHisto.Draw( "E")
    genJetHisto.Draw( "histe same")
    
    if not process.startswith('MCSelfClosure'):
        
        unfoldHistoCross.Scale(1./(unfoldHistoCross.Integral() if not noNorm else 1.),'width')
        
        unfoldHistoCross.SetMarkerStyle(26)
        #unfoldHistoCross.SetMarkerSize(2)
        unfoldHistoCross.SetMarkerColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineColor(ROOT.kRed+4)
        unfoldHistoCross.SetLineWidth(2)
        legend.AddEntry( unfoldHistoCross, f'{nomMCHisto_label} unf. w/ {altMCHisto_label}', 'pe')

        genJetHistoCross.Scale(1./(genJetHistoCross.Integral() if not noNorm else 1.),'width')
        
        genJetHistoCross.SetLineWidth(2)
        genJetHistoCross.SetLineColor(ROOT.kMagenta)
        genJetHistoCross.SetMarkerStyle(0)
        genJetHistoCross.SetLineStyle(2)
        legend.AddEntry( genJetHistoCross, f'{nomMCHisto_label} (gen)', 'lp')
        
        unfoldHistoCross.Draw( "E same")
        genJetHistoCross.Draw( "histe same")

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.045)

    selText.SetNDC()

    if selection.startswith("_dijet"): seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
    elif selection.startswith("_W"): seltext = ' Boosted W region'
    elif selection.startswith("_top"): seltext = ' Boosted top region'
    selText.DrawLatex( ( 0.2 if tlegendAlignment.startswith('right') else 0.68 ), 0.87, seltext )

    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if selection.startswith('dijet') else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    can.cd()
    
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1)#;
    
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()   
    pad2.cd()
    
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    tmpPad2.GetXaxis().SetTitle( labelX )
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Unf.}" if 'Self' in process else "#frac{Unf. with alt. MC}{ Unf. with nom. MC}"  )
    tmpPad2.GetYaxis().SetTitleOffset( 0.5 )
    tmpPad2.GetYaxis().SetRangeUser(0.7,1.4 )
    tmpPad2.GetXaxis().SetRangeUser(unfoldHisto.GetBinLowEdge(1),unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2) )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.12, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    if 'Self' in process:

        hRatioUp = ROOT.TGraphAsymmErrors()
        hRatioUp.Divide( genJetHisto, unfoldHisto, 'pois' )
        hRatioUp.SetLineColor(ROOT.kBlack)
        hRatioUp.SetMarkerColor(ROOT.kBlack)
        hRatioUp.SetLineWidth(2)
        hRatioUp.SetMarkerStyle(25)
        #hRatioUp.GetXaxis().SetLimits(0.,unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        hRatioUp.Draw('P0')

    else:
        hRatioUp2 = ROOT.TGraphAsymmErrors()
        hRatioUp2.Divide( unfoldHistoCross, unfoldHisto, 'pois' )
        hRatioUp2.SetLineColor(ROOT.kBlack)
        hRatioUp2.SetMarkerColor(ROOT.kBlack)
        hRatioUp2.SetLineWidth(2)
        hRatioUp2.SetMarkerStyle(25)
        #hRatioUp2.GetXaxis().SetLimits(0.,unfoldHisto.GetBinLowEdge(unfoldHisto.GetNbinsX()+2))
        hRatioUp2.Draw('P0')
    
    
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)
    
    
    
def get_condition_number(histo):
    
    ## based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L41
    Nx = histo.GetNbinsX()
    Ny = histo.GetNbinsY()
    RMx = histo.ProjectionX( 'RMx', 0, -1 )

    m = ROOT.TMatrixD( Ny, Nx )   ### need to swap the axes
    for ibin in range(1, Nx+1):
        normalization = RMx.GetBinContent(ibin)
        if (normalization>0):
            for jbin in range( 1, Ny+1 ):
                m[jbin-1][ibin-1] = histo.GetBinContent(ibin,jbin) / normalization
    svd = ROOT.TDecompSVD(m)
    v = ROOT.TVectorD( svd.GetSig() )
    Min = v[0]
    Max = v[0]
    for ibin in range( 0, Nx ):
        if (abs(v[ibin]) < 1e-5 ): break
        Min = v[ibin]
    conditionNumber = round( Max/Min, 2 ) if Min > 0 else 1000000
    
    return conditionNumber
    
def draw2D( ivar, histo, varInfo, outputDir, outputLabel='data', addCorrelation=False, addCondition=False, addInvertedMatrix=False,ext='pdf',selection='_dijetSel',version='vNew',year='2017',pngToo=False, outputName=None,  ):

    if not os.path.exists(outputDir): os.makedirs(outputDir)
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_'+version+'.'+ext


    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetPadTopMargin(0.08)

    can2D = ROOT.TCanvas(ivar+'can2D'+histo.GetName(), ivar+'can2D'+histo.GetName(), 750, 500 ) if not('body' in ivar) else ROOT.TCanvas(ivar+'can2D'+histo.GetName(), ivar+'can2D'+histo.GetName(), 800, 600 )
    histo.GetXaxis().SetTitle('Accepted Gen '+varInfo['label'])
    histo.GetYaxis().SetTitle('True Reco '+varInfo['label'])
    histo.GetYaxis().SetTitleOffset( 0.8 )
    #signalHistos[signalLabel+'_resp'+ivar+'_nom'+sel]covHisto.GetYaxis().SetRange( variables[ivar]['bins'][0], variables[ivar]['bins'][-1] )
    histo.Draw("colz")
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(can2D, 4, 0)

    textBox.SetTextSize(0.04)
    if addCorrelation:
        print('|-----> Correlation: ', histo.GetCorrelationFactor())
        textBoxCorr = textBox.Clone()
        textBoxCorr.DrawLatex( 0.06 if not('body' in ivar) else 0.1, varInfo['bins'][-1]-( .05*(varInfo['bins'][-1]-varInfo['bins'][0]) ), '#color[8]{Corr. Factor = '+str(round(histo.GetCorrelationFactor(),3))+'}' )

    if addCondition:   
                
        conditionNumber = get_condition_number(histo.Clone())
        if conditionNumber<=10.:
            print('|-----> Condition Number: ', conditionNumber)
        else:
            print('##############################################')
            print (f' WARNING: Condition Number>10: {conditionNumber} ')
            print('##############################################')
            #print('|-----> Condition Number: ', conditionNumber)
        
        textBoxCond = textBox.Clone()
        textBoxCond.DrawLatex( 0.05, varInfo['bins'][-1]-( .1*(varInfo['bins'][-1]-varInfo['bins'][0]) ), '#color[8]{Cond. Number = '+str(conditionNumber)+'}' )

    can2D.SaveAs(outputName)
    if ext.startswith('pdf') and pngToo:
        can2D.SaveAs( outputName.replace('pdf', 'png') )
    del(can2D)
    gc.collect()
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12)


def combinePlots( name, dictHistos, numBins, mainHistoLabel, otherHisto, otherHistoLabel, outputLabel, variables, ext, process, log, year, runMLU, version,  axisX='', outputDir='Plots/', ratioOnly=False):

    """docstring for combinePlots"""

    outputFileName = name+'_'+outputLabel+'_combinePlots_'+version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.10,0.80,0.60,0.90)
    legend.SetNColumns(3 if runMLU else 2)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.06)

    dictHistos['combData'] = ROOT.TH1F('combData', 'combData', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combData' ], mainHistoLabel, 'lep' )
    dictHistos['combUnfold'] = ROOT.TH1F('combUnfold', 'combUnfold', numBins[-1], 0, numBins[-1])
    legend.AddEntry( dictHistos[ 'combUnfold' ], otherHistoLabel, 'lep' )
    if runMLU:
        dictHistos['combMLU'] = dictHistos['combUnfold'].Clone()
        legend.AddEntry( dictHistos[ 'combMLU' ], 'MLU', 'lep' )

    tmpNbin = 0
    Xlabels = []
    for ivar,ih in dictHistos.items():
        if ivar.startswith('Jet') and not ivar.endswith(('21', '32')):
            Xlabels.append( '#'+variables[ivar]['label'].split('#')[1] )
            for ibin in range(1, ih['data'].GetNbinsX()+1):
                tmpNbin = tmpNbin+1
                dictHistos['combData'].SetBinContent( tmpNbin, ih[otherHisto].GetBinContent(ibin) )
                dictHistos['combData'].SetBinError( tmpNbin, ih[otherHisto].GetBinError(ibin) )
                dictHistos['combUnfold'].SetBinContent( tmpNbin, ih['unfold'].GetBinContent(ibin) )
                dictHistos['combUnfold'].SetBinError( tmpNbin, ih['unfold'].GetBinError(ibin) )
                if runMLU:
                    dictHistos['combMLU'].SetBinContent( tmpNbin, ih['MLU'].GetBinContent(ibin) )
                    dictHistos['combMLU'].SetBinError( tmpNbin, ih['MLU'].GetBinError(ibin) )


    if not ratioOnly:
        canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 750 )
        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetPadLeftMargin(0.08)
        ROOT.gStyle.SetPadTickX(0)
        pad1 = ROOT.TPad(ivar+'1', "Fit",0.,0.330,1.00,1.00,-1)
        pad2 = ROOT.TPad(ivar+'2', "Pull",0,0.00,1.00,0.40,-1);
        pad1.Draw()
        pad2.Draw()

        pad1.cd()
        #pad1.SetLogy()
        dictHistos['combData'].GetXaxis().SetNdivisions(100)
        dictHistos['combData'].GetYaxis().SetTitle( '#frac{d#sigma}{d#tau_{X}}' )
        dictHistos['combData'].GetYaxis().SetTitleSize( 0.06 )
        dictHistos['combData'].GetYaxis().SetTitleOffset( 0.6 )
        dictHistos['combData'].SetMaximum( dictHistos['combUnfold'].GetMaximum()*1.2 )
        dictHistos['combData'].SetMinimum( 0.001 )
        dictHistos['combData'].SetLineColor( ROOT.kBlack )
        dictHistos['combData'].SetLineWidth( 2 )

        dictHistos['combUnfold'].SetLineColor( ROOT.kMagenta )
        dictHistos['combUnfold'].SetLineWidth( 2 )
        dictHistos['combData'].Draw('E')
        dictHistos[ 'combUnfold' ].Draw('E same')
        if runMLU:
            dictHistos['combMLU'].SetLineColor( 8 )
            dictHistos['combMLU'].SetLineWidth( 2 )
            dictHistos[ 'combMLU' ].Draw('E same')

        ### division lines
        lines = {}
        for i in numBins[1:-1]:
            lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
            lines[i].SetLineColor(ROOT.kGray)
            lines[i].Draw('same')

        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.relPosX = 0.07
        CMS_lumi.CMS_lumi( pad1, 4, 0)
        legend.Draw()

        pad2.cd()
        pad2.SetGridy()
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.2)

        tmppad= pad2.DrawFrame(0,0.2,numBins[-1],1.8)
        tmppad.GetYaxis().SetTitle( "Data/Unfold" )
        tmppad.GetXaxis().SetTitle(  axisX )
        tmppad.GetYaxis().SetTitleOffset( 0.4 )
        tmppad.GetYaxis().CenterTitle()
        tmppad.SetLabelSize(0., 'x')
        tmppad.SetTitleSize(0., 'x')
        tmppad.SetLabelSize(0.10, 'y')
        tmppad.SetTitleSize(0.10, 'y')
        tmppad.SetNdivisions(100, 'x')
        tmppad.SetNdivisions(505, 'y')
        pad2.Modified()
        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dictHistos['combUnfold'], dictHistos[ 'combData' ], 'pois' )
        hRatio.SetMarkerStyle(8)
        hRatio.Draw('P0')
        if runMLU:
            hRatioMLU = ROOT.TGraphAsymmErrors()
            hRatioMLU.Divide( dictHistos['combMLU'], dictHistos[ 'combData' ], 'pois' )
            hRatioMLU.SetMarkerStyle(8)
            hRatioMLU.SetMarkerColor( 8 )
            hRatioMLU.Draw('P0 same')
            hRatio.Draw('P0 same')

        for i in lines: lines[i].Draw('same')

        textBox.SetTextSize(0.10)
        textBox.SetTextAlign(12)
        textBoxList = {}
        
        for i in range(1, len(numBins)):
            textBoxList[i] = textBox.Clone()
            textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., 0., Xlabels[i-1] )
        textBox.SetTextSize(0.04)
    else:

        legend=ROOT.TLegend(0.60,0.80,0.90,0.90)
        legend.SetNColumns(3 if runMLU else 2)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.06)

        ROOT.gStyle.SetPadRightMargin(0.05)
        ROOT.gStyle.SetPadLeftMargin(0.08)
        ROOT.gStyle.SetPadTickX(0)
        canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 500 )
        canvas[outputFileName].SetGridy()

        hRatio = ROOT.TGraphAsymmErrors()
        hRatio.Divide( dictHistos[ 'combData' ], dictHistos['combUnfold'], 'pois' )
        if process.startswith('MCCrossClosure'): labelLegend = 'MCSelfClosure Ind.Sample'
        elif process.startswith('MCClosure'): labelLegend = 'MCCrossClosure'
        else: labelLegend = process
        legend.AddEntry( hRatio, labelLegend, 'lep' )
        hRatio.SetMarkerStyle(8)

        hRatio.GetXaxis().SetNdivisions(100)
        hRatio.GetYaxis().SetTitle( 'Sim. / Data'+(' (MC)' if process.startswith('MC') else '') )
        hRatio.GetYaxis().SetTitleSize( 0.06 )
        hRatio.GetYaxis().SetTitleOffset( 0.6 )
        hRatio.GetXaxis().SetLimits( 0., numBins[-1] )
        hRatio.SetMaximum( 2. )
        hRatio.SetMinimum( 0. )

        hRatio.Draw('AP')

        ### division lines
        lines = {}
        for i in numBins[1:-1]:
            lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
            lines[i].SetLineColor(ROOT.kGray)
            lines[i].Draw('same')

        CMS_lumi.lumiTextSize = 0.6
        CMS_lumi.relPosX = 0.07
        CMS_lumi.CMS_lumi( canvas[outputFileName], 4, 0)
        legend.Draw()

        aBox = textBox.Clone()
        aBox.SetTextSize(0.06)
        #aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Central Jet' if name.startswith('recoJet2') else 'Outer Jet')+"}{Dijet Selection}}"  )
        if 'WSel' in selection: aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Leading jet')+"}{Boosted W selection}}"  )
        elif 'topSel'in selection: aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Leading jet')+"}{Boosted top selection}}"  )
        textBox.SetTextAlign(12)
        textBoxList = {}
        for i in range(1, len(numBins)):
            textBoxList[i] = textBox.Clone()
            textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., -0.1, Xlabels[i-1] )

    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]
    ROOT.gStyle.SetPadTickX(1)

def combineRatioPlots( name, ratioDicts, numBins, outputLabel,process, ext, log, version, selection, axisX='', outputDir='Plots/'):
    """docstring for combineRatioPlots"""

    outputFileName = name+'_'+outputLabel+'_combineBLTPlots_'+version+'.'+ext
    if log: outputFileName = outputFileName.replace('Plots','Plots_Log')
    print('Processing.......', outputFileName)

    legend=ROOT.TLegend(0.5,0.8,0.8,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.055)
    legend.SetBorderSize(0)
    
    tmpNbin = 0
    Xlabels = []
    yvaluesGen = np.array( [0]*numBins[-1], 'd' )
    yErrLowValuesGen = np.array( [0]*numBins[-1], 'd' )
    yErrHighValuesGen = np.array( [0]*numBins[-1], 'd' )
    yvaluesReco = np.array( [0]*numBins[-1], 'd' )
    yErrLowValuesReco = np.array( [0]*numBins[-1], 'd' )
    yErrHighValuesReco = np.array( [0]*numBins[-1], 'd' )
    for ivar,ih in ratioDicts.items():
        if ivar.startswith('Jet') and not ivar.endswith(('21', '32')):
            Xlabels.append( '#'+variables[ivar]['label'].split('#')[1] )
            for ibin in range( ih[0].GetN() ):
                a = ctypes.c_double(0.)
                b = ctypes.c_double(0.)
                ih[0].GetPoint(ibin, a, b)
                yvaluesGen[tmpNbin] = b.value
                yErrLowValuesGen[tmpNbin] = ih[0].GetErrorYlow(ibin)
                yErrHighValuesGen[tmpNbin] = ih[0].GetErrorYhigh(ibin)
                c = ctypes.c_double(0.)
                d = ctypes.c_double(0.)
                ih[1].GetPoint(ibin, c, d )
                yvaluesReco[tmpNbin] = d.value
                yErrLowValuesReco[tmpNbin] = ih[1].GetErrorYlow(ibin)
                yErrHighValuesReco[tmpNbin] = ih[1].GetErrorYhigh(ibin)
                tmpNbin = tmpNbin+1

    xvalues = np.array( range(1,numBins[-1]+1), 'd' )
    xvalues = xvalues - 0.5

    dictHistos = {}
    dictHistos['combRecoRatio'] = ROOT.TGraphAsymmErrors(len(xvalues), np.array(xvalues, 'd'), yvaluesReco, np.array([0.5]*len(xvalues), 'd' ), np.array([0.5]*len(xvalues), 'd' ), yErrLowValuesReco, yErrHighValuesReco )
    legend.AddEntry( dictHistos[ 'combRecoRatio' ], 'detector level (stat. unc. only)', 'lep' )
    dictHistos['combGenRatio'] = ROOT.TGraphAsymmErrors(len(xvalues), np.array(xvalues, 'd'), yvaluesGen, np.array([0]*len(xvalues), 'd' ), np.array([0]*len(xvalues), 'd' ), yErrLowValuesGen, yErrHighValuesGen )
    legend.AddEntry( dictHistos[ 'combGenRatio' ], 'hadron level (all uncorr. unc.)', 'lep' )


    canvas[outputFileName] = ROOT.TCanvas('c1'+name, 'c1'+name, 1400, 500 )
    canvas[outputFileName].SetGridy()
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.08)
    ROOT.gStyle.SetPadTickX(0)
    dictHistos['combRecoRatio'].SetLineColor( ROOT.kBlack )
    dictHistos['combRecoRatio'].SetMarkerColor( ROOT.kBlack )
    dictHistos['combRecoRatio'].SetLineWidth( 1 )

    dictHistos['combGenRatio'].SetLineColor( ROOT.kMagenta )
    dictHistos['combGenRatio'].SetMarkerColor( ROOT.kMagenta )
    dictHistos['combGenRatio'].SetLineWidth( 1 )

    multiGraph = ROOT.TMultiGraph()
    multiGraph.Add( dictHistos['combRecoRatio'] )
    multiGraph.Add( dictHistos['combGenRatio'] )

    multiGraph.GetXaxis().SetNdivisions(100)
    multiGraph.GetYaxis().SetTitle( 'Sim. / Data'+(' (MC)' if process.startswith('MC') else '') )
    multiGraph.GetYaxis().SetTitleSize( 0.06 )
    multiGraph.GetYaxis().SetTitleOffset( 0.6 )
    multiGraph.GetXaxis().SetLimits( 0., xvalues[-1]+1 )
    multiGraph.SetMaximum( 4.5 )
    multiGraph.SetMinimum( 0.5)
    multiGraph.Draw('AP')

    aBox = textBox.Clone()
    aBox.SetTextSize(0.06)
    #aBox.DrawLatex( 5, 1.80, '#bf{#splitline{'+('Central Jet' if name.startswith('recoJet2') else 'Outer Jet')+"}{Dijet Selection}}"  )
    if 'WSel' in selection: aBox.DrawLatex( 5, 4.10, '#bf{Boosted W selection}'  )
    elif 'topSel'in selection: aBox.DrawLatex( 5, 4.10, '#bf{Boosted top selection}'  )

    ### division lines
    lines = {}
    for i in numBins[1:-1]:
        lines[i] = ROOT.TGraph(2, array('d', [i,i]), array('d', [0, 200]) )
        lines[i].SetLineColor(ROOT.kGray)
        lines[i].Draw('same')

    CMS_lumi.lumiTextSize = 0.6
    CMS_lumi.relPosX = 0.07
    CMS_lumi.CMS_lumi( canvas[outputFileName], 4, 0)
    legend.Draw()

    textBox.SetTextAlign(12)
    textBoxList = {}
    
    print (len(numBins),len(Xlabels), numBins, Xlabels)

    for i in range(1, len(numBins)):
        textBoxList[i] = textBox.Clone()
        textBoxList[i].DrawLatex(numBins[i-1]+(numBins[i]-numBins[i-1])/2., -0.1, Xlabels[i-1] )
    
    canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName )
    if ext.startswith('pdf'):
        canvas[outputFileName].SaveAs( outputDir+'/'+outputFileName.replace('pdf', 'png') )
    del canvas[outputFileName]
    ROOT.gStyle.SetPadTickX(1)

def set_dynamic_y_range(histList, extra_margin=1.3):
    """
    Given a list of histograms (or a single histogram), determine the min and max values
    to adjust the Y-axis range and require adding some buffer space above the max for legends.

    Parameters:
        histList (list or ROOT.TH1): List of histograms or a single histogram.
        extra_margin (float): Factor by which to multiply the max for spacing.

    Returns:
        None (the first histogram in histList will have its Y range updated)
    """
    if not isinstance(histList, list):
        histList = [histList]

    minVal = 1e9
    maxVal = -1e9

    for h in histList:
        for ibin in range(1, h.GetNbinsX()+1):
            val = h.GetBinContent(ibin)
            if val < minVal: minVal = val
            if val > maxVal: maxVal = val

    hMain = histList[0]
    hMain.GetYaxis().SetRangeUser(minVal * 0.95, maxVal * extra_margin)    
    
    
def set_dynamic_y_range_errRatioHist(histList, extra_margin=1.4, bottom_margin=0.9):
    """
    Given a list of error bar histograms (ie, hist w. errors divide by same hist with no y-errors), 
    determine the min and max values thereof to adjust the Y-axis range and require
    adding some buffer space above the max for legends.

    Parameters:
        histList (list or ROOT.TH1): List of histograms or a single histogram.
        extra_margin (float): Factor by which to multiply the max for spacing.

    Returns:
        None (the first histogram in histList will have its Y range updated)
    """
    if not isinstance(histList, list):
        histList = [histList]

    minVal = 1e9
    maxVal = -1e9

    for h in histList:
        for ibin in range(1, h.GetNbinsX()+1):
            val = 1.-h.GetBinError(ibin)
            if val < minVal: 
                minVal = val
            val = 1.+h.GetBinError(ibin)
            if val > maxVal: 
                maxVal = val
            #print(val)

    hMain = histList[0]
    hMain.GetYaxis().SetRangeUser(minVal * bottom_margin, maxVal * extra_margin)        
    
def get_colour_palette_as_list(palette_requested=None):
    if palette_requested=='vf_10' or palette_requested==None:
        hex_list = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
    elif palette_requested=='vf_8':
        hex_list = ["#1845fb", "#ff5e02", "#c91f16", "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"]
    elif palette_requested=='vf_6':
        hex_list = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
    else:
        hex_list = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
    colour_list = []
    for i in hex_list:
        colour_list.append(ROOT.TColor.GetColor(i))
    return colour_list


def set_palette_from_list(color_list):
    arr = array('i', color_list)
    ROOT.gStyle.SetPalette(len(color_list), arr)
    
def drawUncertainties_from_err_shifts_theoryVariations(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, 
                                                       cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, norming=True ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    colors_cr = [ROOT.kMagenta+2,  ROOT.kBlue-4, 433]  
    colors_syst = get_colour_palette_as_list('vf_8')[1:]#list(reversed())#[ROOT.kBlue+1, ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kGreen+2]
    
    #colors = get_colour_palette_as_list('vf_10')
    #[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    #ROOT.gStyle.SetPalette(len(colors),array('i', colors))
        
    upstyles =   [20,21,34,29,22,23,29,47,33,43,39,41,39,45,117,114,48]
    downstyles = [24,25,28,30,26,32,30,46,27,42,37,40,37,44,38 ,60 , 5 ]
    
    modelVariations = True if 'VariationUNC' in outputName else False
    
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    #otherUncs = OrderedDict()
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)
    
    #if tlegendAlignment.startswith('right'): 
    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoNoNorm.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    

    up_counter=0#1
    down_counter=0#1
    col_counter=0#1
    col_counter_jes=0
    

    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    '''
    if not('dijet' in selection): 
        totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    else:
        if 'all' in year:
            totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        
        if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
    '''
    set_dynamic_y_range_errRatioHist(totalErrHist,1.3)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')

    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    #rmStatErrHist.Draw('L same ')
    #h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    #bkgSubErrHist.Draw('L same ')
    #h2.Draw("L same")
    

    cr_histos = OrderedDict()

    CR1_key=None
    CR2_key=None
    erdOn_key=None
    #set_palette_from_list(colors_cr)

    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            #print (k, col_counter)
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            #cr_histos[k] = normalise_hist(cr_histos[k].Clone())
            cr_histos[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
            if 'cr1' in k.lower():
                CR1_key=k
            elif 'cr2' in k.lower():
                CR2_key=k
            elif 'erd' in k.lower():
                erdOn_key=k

            cr_histos[k].SetLineStyle(1)
            cr_histos[k].SetLineWidth(2)
            cr_histos[k].SetMarkerSize(0)
            #cr_histos[k].SetLineColor(colors[col_counter])
            #col_counter+=1
            cr_histos[k].SetFillColor(0)

    cr_keys = [CR1_key, CR2_key, erdOn_key]
    for i, key in enumerate(cr_keys):
        if key:
            color = colors_cr[i % len(colors_cr)]
            cr_histos[key].SetLineColor(color)
            cr_histos[key].SetMarkerColor(color)
            cr_histos[key].Draw('L same')

    if CR1_key: legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')
    if CR2_key: legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')
    if erdOn_key: legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')
    
    #cr_histos[CR1_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')

    #cr_histos[CR2_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')

    #cr_histos[erdOn_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')

   
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                up_counter=up_counter+1
    
    col_counter=0
    col_counter_jes=0
    
    syst_sources_up = list(normeduncerUnfoldHistoshiftsUp.keys())
    for i, k in enumerate(syst_sources_up):
        
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsUp[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(color)
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
            elif 'CP5' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'UE tune (CP5)', 'p' )
            elif 'mtop' in k:
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'm_{top}', 'p' )
            

            else: 
                print('else in th. syst comp maker', k)
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
                
    
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                down_counter=down_counter+1
    
    #set_palette_from_list(colors_syst)
    syst_sources_down = list(normeduncerUnfoldHistoshiftsDown.keys())
    col_counter=0
    for i, k in enumerate(syst_sources_down):
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsDown[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(color)

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k or 'CP5' in k or 'mtop' in k:
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")

            else: 
                print('else in th. syst comp maker', k)
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    #legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    #legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
        

def drawUncertainties_from_err_shifts_theoryVariations_unitNorm(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, 
                                                       cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, norming=True ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    colors_cr = [ROOT.kMagenta+2,  ROOT.kBlue-4, 433]  
    colors_syst = get_colour_palette_as_list('vf_8')[1:]#list(reversed())#[ROOT.kBlue+1, ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kGreen+2]
    
    #colors = get_colour_palette_as_list('vf_10')
    #[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    #ROOT.gStyle.SetPalette(len(colors),array('i', colors))
        
    upstyles =   [20,21,34,29,22,23,29,47,33,43,39,41,39,45,117,114,48]
    downstyles = [24,25,28,30,26,32,30,46,27,42,37,40,37,44,38 ,60 , 5 ]
    
    modelVariations = True if 'VariationUNC' in outputName else False
    
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    #otherUncs = OrderedDict()
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)
    
    #if tlegendAlignment.startswith('right'): 
    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoNoNorm.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    

    up_counter=0#1
    down_counter=0#1
    col_counter=0#1
    col_counter_jes=0
    

    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    '''
    if not('dijet' in selection): 
        totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    else:
        if 'all' in year:
            totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        
        if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
            totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
        else:
            totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
    '''
    set_dynamic_y_range_errRatioHist(totalErrHist,1.3)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')

    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    #rmStatErrHist.Draw('L same ')
    #h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    #bkgSubErrHist.Draw('L same ')
    #h2.Draw("L same")
    

    cr_histos = OrderedDict()

    CR1_key=None
    CR2_key=None
    erdOn_key=None
    #set_palette_from_list(colors_cr)

    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            #print (k, col_counter)
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            #cr_histos[k] = normalise_hist(cr_histos[k].Clone())
            cr_histos[k].Scale(1./(cr_histos[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
            if 'cr1' in k.lower():
                CR1_key=k
            elif 'cr2' in k.lower():
                CR2_key=k
            elif 'erd' in k.lower():
                erdOn_key=k

            cr_histos[k].SetLineStyle(1)
            cr_histos[k].SetLineWidth(2)
            cr_histos[k].SetMarkerSize(0)
            #cr_histos[k].SetLineColor(colors[col_counter])
            #col_counter+=1
            cr_histos[k].SetFillColor(0)

    cr_keys = [CR1_key, CR2_key, erdOn_key]
    for i, key in enumerate(cr_keys):
        if key:
            color = colors_cr[i % len(colors_cr)]
            cr_histos[key].SetLineColor(color)
            cr_histos[key].SetMarkerColor(color)
            cr_histos[key].Draw('L same')

    if CR1_key: legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')
    if CR2_key: legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')
    if erdOn_key: legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')
    
    #cr_histos[CR1_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')

    #cr_histos[CR2_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')

    #cr_histos[erdOn_key].Draw('L same PLC PMC')
    #legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')

   
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(normeduncerUnfoldHistoshiftsUp[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                up_counter=up_counter+1
    
    col_counter=0
    col_counter_jes=0
    
    syst_sources_up = list(normeduncerUnfoldHistoshiftsUp.keys())
    for i, k in enumerate(syst_sources_up):
        
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsUp[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(color)
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
            elif 'CP5' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'UE tune (CP5)', 'p' )
            elif 'mtop' in k:
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'm_{top}', 'p' )
            

            else: 
                print('else in th. syst comp maker', k)
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
                
    
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue


            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(normeduncerUnfoldHistoshiftsDown[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                #normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(3)
                #normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                #normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                #if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                #col_counter=col_counter+1    
                down_counter=down_counter+1
    
    #set_palette_from_list(colors_syst)
    syst_sources_down = list(normeduncerUnfoldHistoshiftsDown.keys())
    col_counter=0
    for i, k in enumerate(syst_sources_down):
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsDown[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(color)

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k or 'CP5' in k or 'mtop' in k:
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")

            else: 
                print('else in th. syst comp maker', k)
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    #legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    #legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
        
def drawUncertainties_from_err_shifts(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, with_modelUnc=True, norming=False ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
    colors = colors_list#[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
        
    upstyles =   [20,21,22,29,23,34,47,33,43, 117,114,48]  #39,41, 45,
    downstyles = [24,25,26,30,32,28,46,27,42, 38 ,60 , 5 ]  #37,40, 44,
    
        
    modelkey=None
    JES_key=None
    JER_key=None
    btag_key=None
    btagUncIncluded=False
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)

    #if tlegendAlignment.startswith('right'): 
    #    legend=ROOT.TLegend(0.2,0.65,0.8,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone('jesHistoUpMax')
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone('jesHistoDownMax')
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    if 'all' in year:
        jerHistoUpMax = unfoldHistoTotUnc.Clone('jerHistoUpMax')
        jerHistoUpMax.Reset()
        jerHistoDownMax = unfoldHistoTotUnc.Clone('jerHistoDownMax')
        jerHistoDownMax.Reset()

        jerHistoUpMax.Sumw2()
        jerHistoDownMax.Sumw2()
    
    
    #print(uncerUnfoldHisto.keys())
    for k in uncerUnfoldHisto:
        if 'modeltotal'in k.lower() and 'shifthist' in k.lower() and (modelkey==None) and with_modelUnc: 
            modelkey=k
        elif 'jes' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower() and not('const' in k.lower()) and (JES_key==None):
            JES_key=k
            print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            jesHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            jesHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            jerHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            jerHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            if not('dijet' in selection):
                btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
                btagHistoUpMax.Reset()
                btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
                btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            btagHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            btagHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #    print(k)
            #    btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                if 'tau_2_2' in k: print (k,text, up_counter, col_counter,upstyles[up_counter],colors[col_counter])
                col_counter=col_counter+1    
                up_counter=up_counter+1
            
    #up_counter=1
    #down_counter=1
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
                
            if 'btag' in k.lower(): continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()  if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        #modelUnc = normalise_hist(modelUnc.Clone())
        modelUnc.Scale(1./(unftot if norming else 1.),'width')#
        modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
        modelUnc.SetLineStyle(1)
        modelUnc.SetLineWidth(2)
        modelUnc.SetMarkerSize(0)
        modelUnc.SetLineColor(28)
        #col_counter+=1
        modelUnc.SetFillColor(0)
    
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    #if not('dijet' in selection): 
    #    totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    #else:
    #    if 'all' in year:
    #        totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    
    #    if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
   
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.4,0.95)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')
    
    
    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    rmStatErrHist.Draw('L same ')
    h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    bkgSubErrHist.Draw('L same ')
    h2.Draw("L same")
    if with_modelUnc: modelUnc.Draw('L same')

    
    
    for k in otherUncs:
        if ('cr' in k.lower() or 'erd' in k.lower()): continue
        #print(k)
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print ("OtherUncs loop", text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        
        upmax_ibin = 1.
        downmax_ibin = 0.
        diff = 0.
        
        upmax_ibin = jesHistoUpMax.GetBinContent(ibin)
        diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
        
        if (diff>=1. or diff<0.):
            print(f'WARNING: JES total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
            upmax_ibin = 0.
            downmax_ibin = 0.
            diff = 0.
        else:
            downmax_ibin = 1.-diff
        print ("JES total",ibin, upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if 'all' in year:

        for ibin in range(1,jerHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = jerHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: JER total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("JER total",ibin, upmax_ibin,downmax_ibin)
            jerHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            jerHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if not('dijet' in selection) and btagUncIncluded:

        for ibin in range(1,btagHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = btagHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: b-tagging total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("b-tagging total",ibin, upmax_ibin,downmax_ibin)
            btagHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            btagHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)
       
    col_counter=9
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(colors[col_counter])
    jesHistoDownMax.SetLineColor(colors[col_counter])
    jesHistoUpMax.SetMarkerColor(colors[col_counter])
    jesHistoDownMax.SetMarkerColor(colors[col_counter])
    jesHistoUpMax.SetMarkerStyle(39)#upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(37)#downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(2.0)
    jesHistoDownMax.SetMarkerSize(2.0)
    up_counter+=1
    down_counter+=1
    col_counter+=1

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(2.0)
        jerHistoDownMax.SetMarkerSize(2.0)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:# and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(2.0)
            btagHistoDownMax.SetMarkerSize(2.0)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #    legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if ('jes' in k.lower() and not('const' in k.lower())) or 'model'in k.lower() or 'tag' in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower() or ('all' in year and 'jer' in k.lower()): 
            continue
        else:
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'l1' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'L1 prefiring', 'p' )
            elif 'unclus' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'MET uncl. en.', 'p' )

            elif 'const' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'neut' in k.lower():
                    text="Neutral ES"
                elif 'charg' in k.lower():
                    text="Charged ES"
                elif 'photon' in k.lower():
                    text="Photon ES"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )


            else: 
                print('else in non-th. syst comp maker', k)
                
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'Shower & hadronization', 'l' )    
    legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
    
def drawUncertainties_from_err_shifts_unitNorm(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, year, unftot, selection, with_modelUnc=True, norming=False ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
    colors = colors_list#[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
        
    upstyles =   [20,21,22,29,23,34,47,33,43, 117,114,48]  #39,41, 45,
    downstyles = [24,25,26,30,32,28,46,27,42, 38 ,60 , 5 ]  #37,40, 44,
    
        
    modelkey=None
    JES_key=None
    JER_key=None
    btag_key=None
    btagUncIncluded=False
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)

    #if tlegendAlignment.startswith('right'): 
    #    legend=ROOT.TLegend(0.2,0.65,0.8,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone('jesHistoUpMax')
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone('jesHistoDownMax')
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    if 'all' in year:
        jerHistoUpMax = unfoldHistoTotUnc.Clone('jerHistoUpMax')
        jerHistoUpMax.Reset()
        jerHistoDownMax = unfoldHistoTotUnc.Clone('jerHistoDownMax')
        jerHistoDownMax.Reset()

        jerHistoUpMax.Sumw2()
        jerHistoDownMax.Sumw2()
    
    
    #print(uncerUnfoldHisto.keys())
    for k in uncerUnfoldHisto:
        if 'modeltotal'in k.lower() and 'shifthist' in k.lower() and (modelkey==None) and with_modelUnc: 
            modelkey=k
        elif 'jes' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower() and not('const' in k.lower()) and (JES_key==None):
            JES_key=k
            print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            jesHistoUpMax.Scale(1./(jesHistoUpMax.Integral() if norming else 1.),'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            jesHistoDownMax.Scale(1./(jesHistoDownMax.Integral() if norming else 1.),'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            jerHistoUpMax.Scale(1./(jerHistoUpMax.Integral() if norming else 1.),'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            jerHistoDownMax.Scale(1./(jerHistoDownMax.Integral() if norming else 1.),'width')
            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            if not('dijet' in selection):
                btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
                btagHistoUpMax.Reset()
                btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
                btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            btagHistoUpMax.Scale(1./(btagHistoUpMax.Integral() if norming else 1.),'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            btagHistoDownMax.Scale(1./(btagHistoDownMax.Integral() if norming else 1.),'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #    print(k)
            #    btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            normeduncerUnfoldHistoshiftsUp[k].Scale(1./(normeduncerUnfoldHistoshiftsUp[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                if 'tau_2_2' in k: print (k,text, up_counter, col_counter,upstyles[up_counter],colors[col_counter])
                col_counter=col_counter+1    
                up_counter=up_counter+1
            
    #up_counter=1
    #down_counter=1
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
                
            if 'btag' in k.lower(): continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()  if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            normeduncerUnfoldHistoshiftsDown[k].Scale(1./(normeduncerUnfoldHistoshiftsDown[k].Integral() if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        #modelUnc = normalise_hist(modelUnc.Clone())
        modelUnc.Scale(1./(modelUnc.Integral() if norming else 1.),'width')#
        modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
        modelUnc.SetLineStyle(1)
        modelUnc.SetLineWidth(2)
        modelUnc.SetMarkerSize(0)
        modelUnc.SetLineColor(28)
        #col_counter+=1
        modelUnc.SetFillColor(0)
    
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    #if not('dijet' in selection): 
    #    totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    #else:
    #    if 'all' in year:
    #        totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    
    #    if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
   
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.4,0.95)
    
    totalErrHist.GetXaxis().SetTitle('#'+labelX.split('#')[1])
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')
    
    
    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    rmStatErrHist.Draw('L same ')
    h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    bkgSubErrHist.Draw('L same ')
    h2.Draw("L same")
    if with_modelUnc: modelUnc.Draw('L same')

    
    
    for k in otherUncs:
        if ('cr' in k.lower() or 'erd' in k.lower()): continue
        #print(k)
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print ("OtherUncs loop", text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        
        upmax_ibin = 1.
        downmax_ibin = 0.
        diff = 0.
        
        upmax_ibin = jesHistoUpMax.GetBinContent(ibin)
        diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
        
        if (diff>=1. or diff<0.):
            print(f'WARNING: JES total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
            upmax_ibin = 0.
            downmax_ibin = 0.
            diff = 0.
        else:
            downmax_ibin = 1.-diff
        print ("JES total",ibin, upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if 'all' in year:

        for ibin in range(1,jerHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = jerHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: JER total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("JER total",ibin, upmax_ibin,downmax_ibin)
            jerHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            jerHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if not('dijet' in selection) and btagUncIncluded:

        for ibin in range(1,btagHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = btagHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: b-tagging total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            print ("b-tagging total",ibin, upmax_ibin,downmax_ibin)
            btagHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            btagHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)
       
    col_counter=9
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(colors[col_counter])
    jesHistoDownMax.SetLineColor(colors[col_counter])
    jesHistoUpMax.SetMarkerColor(colors[col_counter])
    jesHistoDownMax.SetMarkerColor(colors[col_counter])
    jesHistoUpMax.SetMarkerStyle(39)#upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(37)#downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(2.0)
    jesHistoDownMax.SetMarkerSize(2.0)
    up_counter+=1
    down_counter+=1
    col_counter+=1

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(2.0)
        jerHistoDownMax.SetMarkerSize(2.0)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:# and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(2.0)
            btagHistoDownMax.SetMarkerSize(2.0)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #    legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if ('jes' in k.lower() and not('const' in k.lower())) or 'model'in k.lower() or 'tag' in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower() or ('all' in year and 'jer' in k.lower()): 
            continue
        else:
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'l1' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'L1 prefiring', 'p' )
            elif 'unclus' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'MET uncl. en.', 'p' )

            elif 'const' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'neut' in k.lower():
                    text="Neutral ES"
                elif 'charg' in k.lower():
                    text="Charged ES"
                elif 'photon' in k.lower():
                    text="Photon ES"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )


            else: 
                print('else in non-th. syst comp maker', k)
                
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'Shower & hadronization', 'l' )    
    legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    
def plotSysComparison2( nomHisto, dictUncHistos, outputName, labelX='', 
                       log=False, version='', ext='pdf', year='2017', 
                       outputDir='Plots/',
                       #sys_pref = 'jes',
                       mode = 'onlyJES', sysList=[],pngToo=False,
                      ):
    """Draw comparison of input systematic uncertainties/variations vs. nominal"""
    colors = [ 95, 38, 6, 7, 8, 42, 50, 218, 225, 30, 16,51, 61, 67, 89, 133, 142, 208, 36, 2, 144, 225, 198, 190, 83, 167, 207, 209, 212, 216,  227, 150, 93, 40]

    if 'jes' in mode.lower():
        jes_uncorr_list = [
                                '_jesAbsolute_2016', '_jesBBEC1_2016', '_jesEC2_2016', '_jesHF_2016', '_jesRelativeSample_2016',
                                '_jesAbsolute_2017', '_jesBBEC1_2017', '_jesEC2_2017', '_jesHF_2017', '_jesRelativeSample_2017',
                                '_jesAbsolute_2018', '_jesBBEC1_2018', '_jesEC2_2018', '_jesHF_2018', '_jesRelativeSample_2018'
                              ] 
    elif 'jer' in mode.lower() and year=='all':
        jes_uncorr_list = [
                                '_jer_2016', '_jer_2017', '_jer_2018'
                            ] 
        print(sysList)
    jes_corr_list = ['_jesAbsolute', '_jesBBEC1', '_jesEC2', '_jesFlavorQCD', '_jesHF', '_jesRelativeBal']
    
    if '18' in year and '_jesHEMIssue' in sysList: 
        jes_corr_list.append('_jesHEMIssue')
        
    elif 'all' in year and '_jesHEMIssue_2018' in sysList and 'jes' in mode.lower(): 
        jes_uncorr_list.append('_jesHEMIssue_2018')
        #print(sysList)
        #print(jes_uncorr_list)
        
    
    
    outputFileName = outputName+'_'+version+'.'+ext
    print ('Processing plots for sys comparisons......', outputFileName)
    
    legend=ROOT.TLegend(0.2,0.7,0.9,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.022)
    legend.SetBorderSize(0)

    multiGraph = ROOT.TMultiGraph()
    nomHisto.Sumw2()
    
    
    nomHisto=nomHisto.Clone()
    #print("Before scaling, integral:", nomHisto.Integral())

    nomHisto.Scale(1.,'width')
    
    #print("After normalising and scaling to b.w., nom integral:", nomHisto.Integral())
    
    gnom = ROOT.TGraphAsymmErrors() #nomHisto.Clone()
    gnom.Divide(nomHisto,nomHisto,'pois')
    
    gnom.SetLineColor(ROOT.kBlack)
    gnom.SetLineWidth(3)
    legend.AddEntry( gnom, 'Nominal' , 'l' )
    multiGraph.Add( gnom )
    
    
    dictShifts = {}
    col_counter=0
    colUp_counter=0
    colDown_counter=0
    dictGraphs = {}
    print(dictUncHistos.keys())
    for ih in dictUncHistos.keys():
        
        dictUncHistos[ih].Sumw2()
        if (('jes' in mode.lower()) and not('jes' in ih.lower())) or (('jer' in mode.lower()) and not('jer' in ih.lower())):
            continue
        elif not('jes' in mode.lower()) and ('jes' in ih.lower()) and not('const' in ih.lower()):
            continue
        #print(ih)
        #print("Before scaling, integral:", ih,dictUncHistos[ih].Integral())
        dictShifts[ih] = dictUncHistos[ih].Clone()

        dictShifts[ih].Scale(1.,'width')
        dictGraphs[ih] = ROOT.TGraphAsymmErrors()
        dictGraphs[ih].Divide(dictShifts[ih],nomHisto,'pois')
        
        if 'model' in 'ih' or 'erd' in ih or 'mtop' in ih or 'CR' in ih:
            col_counter = colUp_counter
            colUp_counter+=1
            colDown_counter+=1
            dictGraphs[ih].SetLineStyle( 2 )#dictShifts[ih]
            
        else:
            if 'up' in ih.lower():
                col_counter = colUp_counter
                dictGraphs[ih].SetLineStyle( 3 )#dictShifts[ih]
                #colors= colors_up
                colUp_counter+=1
                
            elif 'down' in ih.lower():
                col_counter = colDown_counter
                dictGraphs[ih].SetLineStyle( 2 )#dictShifts[ih]
                #colors = colors_down
                colDown_counter+=1

        if not col_counter==len(colors)-1:
            dictGraphs[ih].SetLineColor( colors[col_counter] )#dictShifts[ih]
        else:
            col_counter = col_counter-len(colors)+2
            if 'up' in ih.lower(): 
                colUp_counter = col_counter-len(colors)+2
            if 'down' in ih.lower(): 
                colDown_counter = col_counter-len(colors)+2
            dictGraphs[ih].SetLineColor( colors[col_counter] )#dictShifts[ih]


        dictGraphs[ih].SetMarkerStyle(0)
        dictGraphs[ih].SetLineWidth( 2 )
        
        y=year if not('all' in year or '+' in year) else ''
        #print(ih)
        
        stringtocheck=ih.split('_' )[1] if not('const' in ih) else (ih.split('_' )[1]+ih.split('_' )[2]).replace('JES', ' ES ')#if (year=='all') else ih#+'_'+(year if not('+' in year) else '_fullRunII')
            
        if ('jes' in ih and not('const' in ih)) or ( 'jer' in ih and year=='all'):
            flagUncorr=False
            for j in jes_uncorr_list:
                #print(j)
                if j in ih:
                    stringtocheck=ih[1:]
                    #print("Sys. comp plot for uncorr sources, working on:", j,ih,stringtocheck)
                    flagUncorr=True
                    break
            
           
        legend.AddEntry( dictGraphs[ih], stringtocheck, 'l' )#+'_'+y
        multiGraph.Add(dictGraphs[ih])
        
        
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    canUnc = ROOT.TCanvas('canUnc', 'canUnc',  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    
    multiGraph.GetYaxis().SetTitle( 'Variation/Nominal' )
    multiGraph.GetXaxis().SetTitle( labelX )
    if not('jes' in mode.lower() or 'jer' in mode.lower() or 'const' in outputFileName.lower()):
        multiGraph.SetMaximum( 1.7 )
        multiGraph.SetMinimum( 0.5 )
    
    else:
        if not('const' in outputFileName.lower() or 'jes' in mode.lower()):
            multiGraph.SetMaximum( 1.10 )
            multiGraph.SetMinimum( 0.96 )
        elif 'const' in outputFileName.lower():
            multiGraph.SetMaximum( 1.1 )
            multiGraph.SetMinimum( 0.95 )
        else:
            
            multiGraph.SetMaximum( 1.06 )
            multiGraph.SetMinimum( 0.97 )
        
        
    
    multiGraph.Draw('ALP')

    
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    canUnc.Update()
    #canUnc.cd()
    legend.Draw()

    canUnc.SaveAs( outputDir + outputFileName )
    if 'pdf' in outputFileName and pngToo: canUnc.SaveAs( outputDir + outputFileName.replace('png', 'pdf') )
    del canUnc

    
    
def getFilesInDictSamples_fromROOT(labelBegin='', year_list=['2017','2018','all']):
    Files={}
    for iy in year_list:
        for isam in dictSamples:
            if not checkDict( isam, dictSamples )[iy]['skimmerHisto'].endswith('root'): continue
            if isam.startswith(labelBegin):
                Files[isam] = [
                                ROOT.TFile.Open(inputFolder+checkDict( isam, dictSamples )[iy]['skimmerHisto'] ),
                                checkDict( isam, dictSamples )
                            ]
    return Files

def createRecoBins(genBins):
    recoBins=[genBins[0]]
    for i in range(1,len(genBins)):
        recoBins.append(np.round(genBins[i-1]+(genBins[i]-genBins[i-1])/2.,3))
        recoBins.append(genBins[i])
        #recoBins.append(b)
    return recoBins
"""
def makePSplot_simple(purity,stability,variables,var,outputDir,ext='pdf',dictHistos=OrderedDict(),bins=[0.,1.],sel='_dijetSel',year='2017',signalLabelBegin='QCD_HT_MG5+P8'):
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    colorPallete = [ 0, 2, 4, 8, 12, 28, 30 ]

    ROOT.gStyle.SetPadRightMargin(0.05)
    canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)
    #canvas.cd()
    p=ROOT.TH1D("Purity",";;",len(bins)-1,array('d',bins))
    s=ROOT.TH1D("Stability",";;",len(bins)-1,array('d',bins))
    
    for i in range(len(purity)):
        p.SetBinContent(i+1,purity[i])
        s.SetBinContent(i+1,stability[i])
    
    p.SetLineWidth(2)
    p.SetLineColor(colorPallete[3])
    #purity.SetLineStyle(1)
    s.SetLineWidth(2)
    s.SetLineColor(colorPallete[4])
    
    legend=ROOT.TLegend(0.15,0.15,0.90,0.35)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.04)
    #legend.SetNColumns( 4 )
    
    
    legend.AddEntry( p, 'Purity', 'l' )
    legend.AddEntry( s,  'Stability' , 'l' )
    
    p.SetMaximum( 1. )
    p.SetMinimum( 0.2 )
    p.SetLineColor( ROOT.kBlack )
    p.SetMarkerColor( ROOT.kBlack )
    p.SetMarkerSize( 0.5 )
    p.SetLineWidth( 2 )

    s.SetLineColor( ROOT.kMagenta )
    s.SetMarkerColor( ROOT.kMagenta )
    s.SetMarkerSize( 0.5 )
    s.SetLineWidth( 2 )
    p.GetXaxis().SetTitle( variables[var]['label'] )
    p.Draw('hist')
    s.Draw('hist same')
    
    dictHistos[ 'purityGraph_'+var ] = p.Clone()
    dictHistos[ 'stabilityGraph_'+var ] = s.Clone()

    '''
    numBins = numBins + dictHistos[ 'purityGraph_'+var ].GetNbinsX()
    numBinsList.append(numBins)
    multigraph = ROOT.THStack()
    for i in dictHistos:
        if i.endswith(var): multigraph.Add( dictHistos[i]  )

    ROOT.gStyle.SetPadRightMargin(0.05)
    canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)

    multigraph.Draw("hist nostack")
    multigraph.GetXaxis().SetTitle( variables[var]['label'] )
    multigraph.GetYaxis().SetTitle( 'Percentage' )
    multigraph.SetMaximum( 1.1 )
    multigraph.GetYaxis().SetTitleOffset(0.8)
    '''
    #canvas.PlaceLegend()
    legend.Draw()
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canvas, 4, 0)
    '''
    textBox=ROOT.TLatex()
    textBox.SetNDC()
    textBox.SetTextSize(0.04)
    textBox.SetTextFont(62) ### 62 is bold, 42 is normal
    '''
    #textBox.DrawLatex(0.65, 0.75, sel.split('Sel')[0].split('_')[1]+' Selection' )
    canvas.Update()
    canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.'+ext)
    if ext.startswith('pdf'):
        canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.png')
    return 1#numBins
"""
def makePSplot_simple(purity,stability,
                      accepGen,#=None,
                      fakeReco,#=None,
                      variables,
                      var,outputDir, ext,
                      dictHistos=OrderedDict(),
                      
                      bins=[0.,1.],
                      year='2017',
                      #ext,#='pdf',
                      sel='_WSel',
                      signalLabelBegin='TTToSemiLeptonic'):

    if not os.path.exists(outputDir): os.makedirs(outputDir)

    ROOT.gStyle.SetPadRightMargin(0.05)
    canvas = ROOT.TCanvas('canvas', 'canvas', 750, 500)
    p=ROOT.TH1D("Purity",";;",len(bins)-1,array('d',bins))
    s=ROOT.TH1D("Stability",";;",len(bins)-1,array('d',bins))
    a=accepGen.Clone()
    f=fakeReco.Clone()
    for i in range(len(purity)):
        p.SetBinContent(i+1,purity[i])
        s.SetBinContent(i+1,stability[i])

    legend=ROOT.TLegend(0.15,0.15,0.90,0.35)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.04)
    #legend.SetNColumns( 4 )



    p.SetMaximum( 1. )
    p.SetMinimum( 0.2 )
    p.SetLineColor( ROOT.kBlack )
    p.SetMarkerColor( ROOT.kBlack )
    p.SetMarkerSize( 0.5 )
    p.SetLineWidth( 2 )

    s.SetLineColor( ROOT.kMagenta )
    s.SetMarkerColor( ROOT.kMagenta )
    s.SetMarkerSize( 0.5 )
    s.SetLineWidth( 2 )
    p.GetXaxis().SetTitle( variables[var]['label'] )
    legend.AddEntry( p, 'Purity', 'l' )
    legend.AddEntry( s,  'Stability' , 'l' )

    #a.SetLineStyle(2)
    #f.SetLineStyle(2)
    #a.SetLineWidth(2)
    #f.SetLineWidth(2)
    #a.SetLineColor(ROOT.kBlue)
    #f.SetLineColor(ROOT.kRed)

    #legend.AddEntry( a,  'Acceptance' , 'l' )
    #legend.AddEntry( f, 'Fake rate', 'l' )

    p.Draw('hist')
    s.Draw('hist same')
    #a.Draw('hist same')
    #f.Draw('hist same')

    dictHistos[ 'purityGraph_'+var ] = p.Clone()
    dictHistos[ 'stabilityGraph_'+var ] = s.Clone()
    #dictHistos[ 'acceptanceRateGraph_'+var ] = a.Clone()
    #dictHistos[ 'fakeRateGraph_'+var ] = f.Clone()


    legend.Draw()
    CMS_lumi.extraText = "Simulation"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canvas, 4, 0)

    canvas.Update()
    canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+year+'.'+ext)
    if ext.startswith('pdf'):
        canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+'_'+year+'.png')
    return 1#numBins
def createCanvasPads():
    c = ROOT.TCanvas("c", "canvas", 800, 600)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()
 
    return c, pad1, pad2

#below from https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/print_bottom_line_test.py#L112C1-L135C15
def get_null_bins(h):
    null_bins = []
    if isinstance(h, (ROOT.TH1, ROOT.TH2)):
        if isinstance(h, ROOT.TH2):
            h_proj = h.ProjectionX()
        else:
            h_proj = h
        for ix in range(1, h_proj.GetNbinsX()+1):
            if h_proj.GetBinContent(ix) == 0:
                null_bins.append(ix)
        return null_bins
    else:
        proj = h.sum(axis=0)
        return np.where(proj == 0)[0]


def remove_null_bins(arr, null_bins):
    if len(arr.shape) > 1:
        for ax in range(len(arr.shape)):
            if arr.shape[ax] > 1:
                arr = np.delete(arr, null_bins, axis=ax)
    else:
        arr = np.delete(arr, null_bins)
    return arr


"""
def bottomLineTest( ivar, dataHisto, mainHistoLabel, datarecoHisto, recoHisto, genHisto, datacovMatrix, hadcovMatrix, varInfo, outputLabel, outputDir, runMLU=False, ext='pdf',selection='_dijetSel',version='vNew' ):
    #based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L74

    outputDir=outputDir+'/'+ivar+'/'+process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    
    
    ##### computing chi2 and inverted matrix
    vector = []
    ndf = 0
    for ibin in range(1, datacovMatrix.GetNbinsX()+1):
        if (datacovMatrix.GetBinContent(ibin, ibin) > 0.):
            ndf = ndf + 1
            vector.append( datarecoHisto.GetBinContent( ibin ) - recoHisto.GetBinContent( ibin ) )
            
    print (ndf,datacovMatrix.GetNbinsX(),datacovMatrix.GetNbinsY())
    #assert ndf==covMatrix.GetNbinsX()
    matrix = np.eye( ndf, ndf )
    
    for ibin in range(1, datacovMatrix.GetNbinsX()+1):
        if (datacovMatrix.GetBinContent(ibin, ibin) > 0.):
            for jbin in range(1, datacovMatrix.GetNbinsY()+1):
                matrix[ibin-1][jbin-1] = datacovMatrix.GetBinContent( ibin, jbin )
    
    vector = np.array( vector )
    invMatrix = np.linalg.inv( matrix )
    chi2 = np.dot( vector, np.dot( invMatrix, vector ) )
    print('Detector level: chi2, ndf, chi2/ndf = ', chi2, ndf, chi2/ndf)
    
    
    vector = []
    ndf = 0
    for ibin in range(1, hadcovMatrix.GetNbinsX()+1):
        if (hadcovMatrix.GetBinContent(ibin, ibin) > 0.):
            ndf = ndf + 1
            vector.append( dataHisto.GetBinContent( ibin ) - genHisto.GetBinContent( ibin ) )
            
    print (ndf,hadcovMatrix.GetNbinsX(),hadcovMatrix.GetNbinsY())
    
    matrix = np.eye( ndf, ndf)
    
    for ibin in range(1, hadcovMatrix.GetNbinsX()+1):
        if (hadcovMatrix.GetBinContent(ibin, ibin) > 0.):
            for jbin in range(1, hadcovMatrix.GetNbinsY()+1):
                matrix[ibin-1][jbin-1] = hadcovMatrix.GetBinContent( ibin, jbin )

    #assert ndf==covMatrix.GetNbinsX()
    vector = np.array( vector )
    invMatrix = np.linalg.inv( matrix )
    chi2 = np.dot( vector, np.dot( invMatrix, vector ) )
    print('Hadron level: chi2, ndf, chi2/ndf = ', chi2, ndf, chi2/ndf)
    
    


    #### plotting inverted matrix
    invertedMatrix = hadcovMatrix.Clone()
    invertedMatrix.Reset()
    for ibin in range(1, hadcovMatrix.GetNbinsX()):
        for jbin in range(1, hadcovMatrix.GetNbinsY()):
            invertedMatrix.SetBinContent( ibin, jbin, invMatrix[ibin-1][jbin-1]  )
    draw2D( ivar, invertedMatrix, varInfo, outputLabel=outputLabel+'_invertedMatrix', outputDir=outputDir.split(ivar)[0],selection=selection,version=version )
    
    ##### plotting ratios together
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_bottomLineTest_'+version+'.'+ext

    canRatio = ROOT.TCanvas('canRatio'+ivar, 'canRatio'+ivar,  10, 10, 750, 500 )

    genRatio = ROOT.TGraphAsymmErrors()
    genRatio.Divide( dataHisto, genHisto, 'pois' )
    recoRatio = ROOT.TGraphAsymmErrors()
    recoRatio.Divide( datarecoHisto, recoHisto, 'pois' )

    legend=ROOT.TLegend(0.15,0.70,0.40,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.AddEntry( genRatio, 'Hadron level (w/ stat + unf. unc.)', 'pl' )
    legend.AddEntry( recoRatio, 'Detector level (stat unc.)', 'pl' )
    legend.SetBorderSize(0)
    
    genRatio.SetLineWidth(2)
    genRatio.GetXaxis().SetTitle(varInfo['label'])
    genRatio.GetXaxis().SetLimits( varInfo['bins'][0], varInfo['bins'][-1] )
    genRatio.GetYaxis().SetTitle('Data / Simulation')
    genRatio.GetYaxis().SetTitleOffset( 0.8 )
    genRatio.SetMarkerStyle(8)
    genRatio.GetYaxis().SetRangeUser(-5.,10.)
    genRatio.Draw('AP0')
    recoRatio.SetLineColor(ROOT.kRed)
    recoRatio.SetLineWidth(2)
    recoRatio.SetMarkerStyle(4)
    recoRatio.Draw('P0 same')

    lineOne = ROOT.TGraph(2, array('d', [0, 1]), array('d', [1, 1]))
    lineOne.Draw('same')

    legend.Draw()
    CMS_lumi.extraText = " Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canRatio, 4, 0)
    canRatio.SaveAs(outputName)

    return [genRatio, recoRatio]

def bottomLineTest2(ivar, dataHisto, mainHistoLabel, 
                    datarecoHisto, recoHisto, genHisto, 
                    datacovMatrix, hadcovMatrix, 
                    varInfo, outputLabel, outputDir, 
                    selection = '_dijetSel',  ext='pdf' ):
    #based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L74

    outputDir=outputDir+'/'+ivar+'/'+process+'/'
    if not os.path.exists(outputDir): os.makedirs(outputDir)

    #recoHisto.Rebin( 2 )  ### because data and covMatrix have less number of bins
    #datarecoHisto.Rebin(2)
    #recoHisto.Scale( dataHisto.Integral()/recoHisto.Integral() )
    
    #genHisto.Scale( dataHisto.Integral()/genHisto.Integral() )

    #dataHisto.Scale(datarecoHisto.Integral()/dataHisto.Integral())

    #datarecoHisto.Scale( 1 )
    
    
    ##### computing chi2 and inverted matrix
    vector = []
    ndf = 0
    chi2_1 = 0.
    for ibin in range(1, datacovMatrix.GetNbinsX()):
        if (datacovMatrix.GetBinContent(ibin, ibin) > 0):
            chi2_1+=pow((recoHisto.GetBinContent( ibin )-datarecoHisto.GetBinContent( ibin )), 2. ) / np.sqrt(datacovMatrix.GetBinContent(ibin,ibin))
            #vector.append(  )
            
    chi2_1 /= recoHisto.GetNbinsX()
    print('Detector level: chi2, ndf, chi2/ndf = ', chi2_1, recoHisto.GetNbinsX(), chi2_1/recoHisto.GetNbinsX())
    
    
    chi2_2 = 0.
    for ibin in range(1, hadcovMatrix.GetNbinsX()):
        if (hadcovMatrix.GetBinContent(ibin, ibin) > 0):
            chi2_2+=pow(( dataHisto.GetBinContent( ibin ) - genHisto.GetBinContent( ibin ) ), 2) / np.sqrt(hadcovMatrix.GetBinContent(ibin,ibin))
            
    print (ndf,hadcovMatrix.GetNbinsX(),hadcovMatrix.GetNbinsY())
    
    chi2_2/=dataHisto.GetNbinsX()
    print('Hadron level: chi2, ndf, chi2/ndf = ', chi2_2, dataHisto.GetNbinsX(), chi2_2/dataHisto.GetNbinsX())
    
    ##### plotting ratios together
    outputName = outputDir+ivar+'_'+selection+'_'+outputLabel+'_bottomLineTest_'+version+'.'+ext

    canRatio = ROOT.TCanvas('canRatio'+ivar, 'canRatio'+ivar,  10, 10, 750, 500 )

    genRatio = ROOT.TGraphAsymmErrors()
    genRatio.Divide( dataHisto, genHisto, 'pois' )
    recoRatio = ROOT.TGraphAsymmErrors()
    recoRatio.Divide( datarecoHisto, recoHisto, 'pois' )

    legend=ROOT.TLegend(0.15,0.70,0.40,0.90)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.AddEntry( genRatio, 'Hadron level (w/ stat + unf. unc.)', 'pl' )
    legend.AddEntry( recoRatio, 'Detector level (stat unc.)', 'pl' )
    legend.SetBorderSize(0)
    
    genRatio.SetLineWidth(2)
    genRatio.GetXaxis().SetTitle(varInfo['label'])
    genRatio.GetXaxis().SetLimits( varInfo['bins'][0], varInfo['bins'][-1] )
    genRatio.GetYaxis().SetTitle('Data / Simulation')
    genRatio.GetYaxis().SetTitleOffset( 0.8 )
    genRatio.SetMarkerStyle(8)
    genRatio.GetYaxis().SetRangeUser(-5.,10.)
    genRatio.Draw('AP0')
    recoRatio.SetLineColor(ROOT.kRed)
    recoRatio.SetLineWidth(2)
    recoRatio.SetMarkerStyle(4)
    recoRatio.Draw('P0 same')

    lineOne = ROOT.TGraph(2, array('d', [0, 1]), array('d', [1, 1]))
    lineOne.Draw('same')

    legend.Draw()
    CMS_lumi.extraText = " Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.11
    CMS_lumi.CMS_lumi(canRatio, 4, 0)
    canRatio.SaveAs(outputName)

    return [genRatio, recoRatio]

def bottomLineTest( ivar, dataHisto, dataHistoLabel, MCHisto, covMatrix, varInfo, outputLabel, outputDir, rebin=1,
                    ext='png', version='V2023', selection='_dijet', process='data'):
    '''based on https://gitlab.cern.ch/DasAnalysisSystem/InclusiveJet/-/blob/master/UnfoldingSampleND/bin/unfold.cc#L74'''
    
    print("Data nbins","MC nbins","cov nbinsX","cov nbinsY")
        
    if not(rebin==1):
        MCHisto.Rebin( rebin )  ### because data and covMatrix have less number of bins
        dataHisto.Rebin( rebin )
        covMatrix.Rebin2D( rebin, rebin )
    #print(dataHisto.GetNbinsX(),MCHisto.GetNbinsX(),covMatrix.GetNbinsX(),covMatrix.GetNbinsY())
     
    ##### computing chi2 and inverted matrix
    vector = []
    ndf = 0
    
    for ibin in range(1, covMatrix.GetNbinsX()+1):
        if (covMatrix.GetBinContent(ibin, ibin) > 0):# and dataHisto.GetBinContent( ibin ) >0 and MCHisto.GetBinContent( ibin ) > 0 :
            ndf = ndf + 1
            vector.append( dataHisto.GetBinContent( ibin ) - MCHisto.GetBinContent( ibin ) )
            #for jbin in range(1, covMatrix.GetNbinsY()+1):
            #    matrix[ibin-1][jbin-1] = covMatrix.GetBinContent( ibin, jbin )
        else: continue
    
    matrix = np.eye( ndf, ndf) #covMatrix.GetNbinsX(), covMatrix.GetNbinsX() )
    
    for ibin in range(1, ndf+1):
        if (covMatrix.GetBinContent(ibin, ibin) > 0):# and dataHisto.GetBinContent( ibin ) >0 and MCHisto.GetBinContent( ibin ) > 0 :
            #ndf = ndf + 1
            #vector.append( dataHisto.GetBinContent( ibin ) - MCHisto.GetBinContent( ibin ) )
            for jbin in range(1, ndf+1):#covMatrix.GetNbinsY()+1):
                matrix[ibin-1][jbin-1] = covMatrix.GetBinContent( ibin, jbin )
                
        else:
            #ndf = ndf + 1
            continue
            #vector.append( 0.)#dataHisto.GetBinContent( ibin ) - recoHisto.GetBinContent( ibin ) )
            #print(f"WARNING: empty bin in cov or data or MC histo; bin contents of these in bin {ibin} listed as follows: {covMatrix.GetBinContent( ibin ), dataHisto.GetBinContent( ibin ), MCHisto.GetBinContent( ibin )}")
            #for jbin in range(1, covMatrix.GetNbinsY()+1):
            #    matrix[ibin-1][jbin-1] = machineEps()#covMatrix.GetBinContent( ibin, jbin ) #hack for now, forces crash  if condition not met since matrix will be singular and non-invertible

    #assert ndf==covMatrix.GetNbinsX()
    vector = np.array( vector )
    invMatrix = np.linalg.inv( matrix )
    chi2 = np.dot( vector, np.dot( invMatrix, vector ) )
    print(f'chi2 for {dataHistoLabel}, ndf, chi2/ndf = ', chi2, ndf, chi2/ndf)
    
    return chi2, ndf

"""

############# n dim helpers ################

def plot_combined_MCCrossClosure(combined_alt0_truth_hist,
                                 unfolded_combined,
                                 unfolded_combined_alt,
                                 combined_truth_hist,
                                 selection='_dijetSel',
                                 labelX='',
                                 outputDir = '../Plots_January25_dijetSel_NDim/dijetSel/',
                                 outputFilename = 'combined_MCCrossClosure',
                                 tlegendAlignment='right',                                 
                                 year='all',
                                 ext='.pdf',
                                 maxYFactor=1.04
                                ):

    #combined_alt0_truth_hist = altSignalHistos['combined_H7MLMQCD_HT2000toInf_gen6bodyOC_nom_dijetSel'].Clone()
    

    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    can = ROOT.TCanvas('can'+'CrossClosure', 'can'+'CrossClosure',  10, 10, 2000, 1500 )
    pad1 = ROOT.TPad("pad1"+'CrossClosure', "Main",0,0.3,1.00,1.00,-1)

    pad1.Draw()
    
    can.cd()
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)

    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.55,0.61,0.77,0.89)

    else: legend=ROOT.TLegend(0.16,0.61,0.38,0.89)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.033)
    legend.SetBorderSize(0)

    unfold_integral=unfolded_combined.Integral()
    #combined_truth_hist.Scale(1./unfold_integral)

    #unfolded_combined.Scale(1./unfold_integral)

    combined_truth_hist.SetLineColor(ROOT.kBlue)
    #combined_alt0_truth_hist.Scale(1./unfold_integral)
    combined_truth_hist.SetLineStyle(1)
    unfolded_combined.SetLineStyle(2)
    #unfolded_combined_alt.Scale(1./unfold_integral)
    legend.AddEntry( unfolded_combined, 'MG5-MLM+P8 unf. w/ MG5-MLM+P8' if 'dijet' in selection.lower() else 'PWHG+P8 unf. w/ PWHG+P8', 'pe' )
    legend.AddEntry( combined_truth_hist, 'MG5-MLM+P8 (gen)' if 'dijet' in selection.lower() else 'PWHG+P8 (gen)' , 'lp' )
    unfolded_combined.SetMarkerColor(ROOT.kRed)
    unfolded_combined.SetLineColor(ROOT.kRed)
    unfolded_combined.GetYaxis().SetTitleSize(0.05)
    #unfolded_combined.Draw()
    #unfolded_combined.GetYaxis().SetRangeUser(0., ( 2.2*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] )  ))
    #combined_truth_hist.SetMaximum(  1.56*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] )  )
    #combined_truth_hist.GetYaxis().SetRangeUser(0., (  maxYFactor*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] ) ))
    combined_truth_hist.SetMaximum( ( maxYFactor*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum(),
                                                       combined_alt0_truth_hist.GetMaximum(), unfolded_combined_alt.GetMaximum(),
                                                     ] ) ))
    combined_truth_hist.SetMinimum( 0.)
    #pad1.Modified()
    #pad1.Update()
    #can.Modified()
    #can.Update()
    
    
    #can.Update()
    #can.Modified()
    combined_truth_hist.GetYaxis().SetTitle('N_{events}')
    combined_truth_hist.GetYaxis().SetTitleSize(0.05)
    
    combined_truth_hist.Draw('histE')
    ROOT.TGaxis.SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0.005, "y")
    
    unfolded_combined.Draw('histE same')

    unfolded_combined_alt.SetMarkerStyle(26)
    unfolded_combined_alt.SetMarkerColor(ROOT.kRed+4)
    unfolded_combined_alt.SetLineColor(ROOT.kRed+4)
    unfolded_combined_alt.SetLineWidth(1)

    combined_alt0_truth_hist.SetLineWidth(1)
    combined_alt0_truth_hist.SetLineColor(ROOT.kMagenta)
    combined_alt0_truth_hist.SetMarkerStyle(0)
    combined_alt0_truth_hist.SetLineStyle(2)

    legend.AddEntry( unfolded_combined_alt, 'MG5-MLM+P8 unf. w/ MG5-MLM+H7' if 'dijet' in selection.lower() else 'PWHG+P8 unf. w/ PWHG+P7', 'pe' )
    legend.AddEntry( combined_alt0_truth_hist, 'MG5-MLM+H7 (gen)' if 'dijet' in selection.lower() else 'PWHG+H7 (gen)' , 'lp' )

    unfolded_combined_alt.Draw('E same')
    combined_alt0_truth_hist.Draw('histE same')

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.040)

    selText.SetNDC()


    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.80, seltext )
    legend.Draw()
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    can.cd()
    pad2 = ROOT.TPad("pad2"+'CrossClosure', "Ratio",0,0.00,1.00,0.30,-1)#;
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio = ROOT.TGraphAsymmErrors()
    #ratio.Divide(combined_truth_hist,unfolded_combined,'pois')
    ratio.Divide(unfolded_combined_alt,unfolded_combined,'pois')
    
    
    tmpPad2= pad2.DrawFrame( combined_truth_hist.GetXaxis().GetBinLowEdge(1), 0., combined_truth_hist.GetXaxis().GetBinLowEdge(combined_truth_hist.GetNbinsX()), 1.9 )
    #print (labelX)
    tmpPad2.GetYaxis().SetRangeUser(0.6, 1.4 )
    tmpPad2.GetXaxis().SetTitleOffset( 0.9 )
    
    tmpPad2.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')# '#'+labelX.split('#')[1] )#.SetTitle(f'{}N-subjettiness basis')
    tmpPad2.GetYaxis().SetTitleOffset( 0.50 )
    #tmpPad2.GetYaxis().SetTitle( "#frac{Unf.}{Sim.}" )
    tmpPad2.GetYaxis().SetTitle( "#frac{Alt RM unf.}{Nom. RM. unf.}" )   
    tmpPad2.GetYaxis().CenterTitle()
    
    tmpPad2.SetLabelSize(0.13, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    #pad2.cd()

    #ratio = ROOT.TRatioPlot(combined_truth_hist,unfolded_combined,)
    #ratio.GetLowerPad().etYaxis().SetRangeUser(0.8, 1.2)
    
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetLineWidth(1)
    ratio.SetMarkerStyle(25)
    ratio.SetMarkerSize(1.5)
    #set_dynamic_y_range_errRatioHist(ratio,1.5,0.5)
    
    ratio.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')# '#'+labelX.split('#')[1] )
    ratio.GetXaxis().SetTitleOffset( 0.9 )
    ratio.GetYaxis().SetTitle( "#frac{Alt RM unf.}{Nom. RM. unf.}" )
    ratio.GetYaxis().SetTitleOffset( 0.50 )
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleOffset( 0.50 ) 
    ratio.GetXaxis().SetLabelSize(0.12)
    ratio.GetXaxis().SetTitleSize(0.13)

    ratio.GetYaxis().SetLabelSize(0.12)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.Draw('PE1')
    #can.Draw()
    can.SaveAs(outputDir+outputFilename+ext)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12) 
    #return can, pad1, pad2

def drawUnfold_Ndim(ivar, process, 
                    lumi, 
                    dataJetHisto, 
                    genJetHisto, 
                    unfoldHisto, 
                    unfoldHistoStatUnc, 
                    unfoldHistowoUnc,
                    altMCHisto, 
                    foldHisto, 
                    recoJetHisto, 
                    cov_tot, 
                    cov_datastat_tot, 
                    labelX, 
                    maxX, 
                    tlegendAlignment, 
                    outputName,
                    year='all',
                    selection='_dijetSel',
                    altMC1Histo = None, 
                    altMC2Histo = None, 
                    altMC1Histo_label = None, 
                    altMC2Histo_label = None, 
                    extraMC=False,
                    includeFSR = False, fsrUpHisto = None, fsrDownHisto=False):
    
    """docstring for drawUnfold"""
    print ("Drawing unfolding for:",ivar)
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    
    colors = [ROOT.TColor.GetColor("#e42536"),ROOT.TColor.GetColor("#5790fc"),ROOT.TColor.GetColor("#f89c20")]
    
    dataJetHisto.SetTitle("")
    print("data(minus bkgs).Integral()",dataJetHisto.Integral())
    genJetHisto.SetTitle("")
    print("genJetHisto.Integral()",genJetHisto.Integral())
    unfoldHisto.SetTitle("")
    print("unfoldHisto.Integral()",unfoldHisto.Integral())
    unfoldHistoStatUnc.SetTitle("")
    #print("unfoldHistoStatUnc.Integral()",unfoldHistoStatUnc.Integral())
    #unfoldHistowoUnc.SetTitle("")
    #print("unfoldHistowoUnc.Integral()",unfoldHistowoUnc.Integral())
    altMCHisto.SetTitle("")
    print("altMCHisto.Integral()",altMCHisto.Integral())
    #foldHisto.SetTitle("")
    #print("foldHisto.Integral()",foldHisto.Integral())
    recoJetHisto.SetTitle("")
    print("(RM proj.Y )recoJetHisto.Integral()",recoJetHisto.Integral())
    if includeFSR: 
        fsrUpHisto.SetTitle("")
        print("fsrUpHisto.Integral()",fsrUpHisto.Integral())
        fsrDownHisto.SetTitle("")
        print("fsrDownHisto.Integral()",fsrDownHisto.Integral())

            
    
    can = ROOT.TCanvas('can'+ivar, 'can'+ivar,  10, 10, 2000, 1500 )
    pad1 = ROOT.TPad("pad1"+ivar, "Main",0,0.3,1.00,1.00,-1)
    pad1.Draw()
    
    can.cd()
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.68,0.61,0.90,0.89)

    else: legend=ROOT.TLegend(0.16,0.61,0.38,0.89)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    
    #bins = variables[ivar]['bins']

    unfoldHistoTot = unfoldHisto.Clone()
    dataScaling = unfoldHisto.Integral()
    
    print (dataScaling)
    #use unnormed unfold histo to build the jacobian for the correct propagation of errors
    #via the covariance matrix, from the normalise -> the unnormalised space
    #normed_cov_tot_matrix, normed_cov_tot = GetNormalizedTMatrixandTH2(cov_tot.Clone(),"normed_cov_tot", unfoldHisto.Clone())
    
    #normed_cov_datastat_tot_matrix, normed_cov_datastat_tot = GetNormalizedTMatrixandTH2(cov_datastat_tot.Clone(),"normed_cov_dastat_tot", unfoldHisto.Clone())
    unfoldHistoDataStatErr=unfoldHistoStatUnc.Clone()
    unfoldHistoDataStatErr.Sumw2()
    unfoldHisto.Sumw2()
    dataJetHisto.Sumw2()
    genJetHisto.Sumw2()
    unfoldHistowoUnc.Sumw2()
    altMCHisto.Sumw2()
    foldHisto.Sumw2()
    recoJetHisto.Sumw2()
    
    #unfoldHistoDataStatErr.Scale(1./dataScaling, 'width')#normalise_hist(unfoldHistoStatUnc.Clone())    
    #unfoldHisto.Scale(1./dataScaling, 'width')#normalise_hist(unfoldHisto.Clone())    
    #dataJetHisto.Scale(1./dataScaling, 'width')#normalise_hist(dataJetHisto.Clone())
    #genJetHisto.Scale(1./dataScaling, 'width')#normalise_hist(genJetHisto.Clone())
    #unfoldHistowoUnc.Scale(1./dataScaling, 'width')#normalise_hist(unfoldHistowoUnc.Clone())#_divide_bin_width
    #altMCHisto.Scale(1./dataScaling, 'width')#normalise_hist(altMCHisto.Clone())
    #foldHisto.Scale(1./dataScaling, 'width')#normalise_hist(foldHisto.Clone())
    #recoJetHisto.Scale(1./dataScaling, 'width')#normalise_hist(recoJetHisto.Clone())
    
    
    
    
    if includeFSR: 
        fsrUpHisto.Sumw2()
        #fsrUpHisto.Scale(1./dataScaling, 'width')#normalise_hist(fsrUpHisto.Clone())
        fsrDownHisto.Sumw2()
        #fsrDownHisto.Scale(1./dataScaling, 'width')#normalise_hist(fsrDownHisto.Clone())    
        
        
    
    
    if extraMC:

        altMC1Histo.Sumw2()
        #altMC1Histo.Scale(1./dataScaling, 'width')#normalise_hist(altMC1Histo.Clone())
        
        altMC1Histo.SetTitle("")
        if 'dijet' in selection:
            altMC2Histo.Sumw2()
            #altMC2Histo.Scale(1./dataScaling, 'width')#normalise_hist(altMC2Histo.Clone())
            
            altMC2Histo.SetTitle("")

    
    
    
    unfoldHisto.SetMarkerStyle(8)
    unfoldHisto.SetMarkerSize(1.5)
    unfoldHisto.SetMarkerColor(ROOT.kBlack)
    unfoldHisto.SetLineColor(ROOT.kBlack)
    legend.AddEntry( unfoldHisto, 'Data', 'pe' )
    
    
    genJetHisto.SetLineWidth(1)
    genJetHisto.SetLineColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerColor(colors[0])#ROOT.kRed)
    genJetHisto.SetMarkerSize(1.5)
    genJetHisto.SetMarkerStyle(25)
    if includeFSR: 
        fsrUpHisto.SetMarkerSize(1.5)
        fsrUpHisto.SetLineColor(46)
        fsrUpHisto.SetMarkerColor(46)
        fsrUpHisto.SetMarkerStyle(22)


        fsrDownHisto.SetMarkerSize(1.5)
        fsrDownHisto.SetLineColor(46)
        fsrDownHisto.SetMarkerColor(46)
        fsrDownHisto.SetMarkerStyle(23)
    
    legend.AddEntry( genJetHisto, 'MG5-MLM+P8' if 'dijet' in selection else 'PWHG+P8', 'lpe' )

   
    if 'body' in labelX: 
        unfoldHisto.GetYaxis().SetTitle( 'N_{events}')#frac{1}{#sigma} #frac{d#sigma}{d#(6-body_OC)'+'}' )#labelX.split('#')[1]+
    else:
        label=None
        if 'pt' in labelX:
            label = 'p_T'
        elif 'mass'in labelX:
            label = 'm'
        elif 'softdrop' in labelX:
            label = 'm_SD'
        else:
            pass
        if label: unfoldHisto.GetYaxis().SetTitle( '#frac{1}{#sigma} #frac{d#sigma}{d'+label+'}' )
    #unfoldHisto.GetYaxis().SetTitleOffset(0.95)
    unfoldHisto.GetYaxis().SetTitleSize(0.05)
    unfoldHisto.SetMaximum( (1.6 if '21' in ivar or '32' in ivar else 1.56)*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] )  )
    unfoldHisto.SetMinimum(0.)
    #pad1.GetYaxis().SetRangeUser(0,1.5*max([ genJetHisto.GetMaximum(), unfoldHisto.GetMaximum()] ) )

    unfoldHisto.Draw( "E1")
    ROOT.TGaxis.SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0.005, "y")
    
    
    #altMCHisto.Scale(1, 'width')  ### divide by bin width
    altMCHisto.SetLineWidth(1)
    altMCHisto.SetMarkerSize(1.5)
    altMCHisto.SetLineColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerColor(colors[1])#ROOT.kBlue)
    altMCHisto.SetMarkerStyle(25)
    
    if includeFSR: 

        legend.AddEntry(fsrUpHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR up", 'pe')

        legend.AddEntry(fsrDownHisto, #('MG5-MLM+P8, ' if 'dijet' in selection else 'PWHG+P8, ') + 
                        "FSR down", 'pe')
        
    legend.AddEntry( altMCHisto, 'MG5-MLM+H7' if 'dijet' in selection else'PWHG+H7','lpe')# 'aMC@NLO+Pythia8', 'lp' )
    
    
    
    if extraMC:
        
        
        if 'dijet' in selection: 
        
            #altMC2Histo.Scale(1, 'width')  ### divide by bin width
            altMC2Histo.SetLineWidth(1)
            altMC2Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC2Histo.SetMarkerStyle(25)
            altMC2Histo.SetMarkerSize(1.5)
            
            legend.AddEntry( altMC2Histo, 'P8+P8' if 'Pt' in altMC2Histo_label else 'MG5+P8', 'lpe' )
        
            altMC2Histo.Draw("histE1 same")
        else:
            #altMC1Histo.Scale(1, 'width')  ### divide by bin width
            altMC1Histo.SetLineWidth(1)
            altMC1Histo.SetLineColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerColor(colors[2])#ROOT.kGray+4)
            altMC1Histo.SetMarkerStyle(25)
            altMC1Histo.SetMarkerSize(1.5)
            #print("altMC1Histo.Integral()",altMC1Histo.Integral())
            legend.AddEntry( altMC1Histo, 'aMC@NLO-FxFx+P8', 'lpe' )
            altMC1Histo.Draw("histE1 same")

        
    genJetHisto.Draw( "histE1 same")
    altMCHisto.Draw("histE1 same")
    if includeFSR: 
        fsrUpHisto.Draw( "PE1 same")
        fsrDownHisto.Draw("PE1 same")

    
    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.040)

    selText.SetNDC()


    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.80, seltext )
    
    legend.Draw()
    if process.startswith('data'):
        CMS_lumi.extraText = "Preliminary"
        CMS_lumi.lumi_13TeV = ('#leq' if 'dijet' in selection else '')+str( round( (lumi/1000.), 2 ) )+" fb^{-1}, 13 TeV"+('' if year.startswith('all') else ", "+( '2016+2017+2018' if year.startswith('all') else year ) )
    else:
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    
    
    can.cd()
    pad2 = ROOT.TPad("pad2"+ivar, "Ratio",0,0.00,1.00,0.30,-1);
    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio_datastatUnc = unfoldHistoDataStatErr.Clone()
    ratio_datastatUnc.Divide(unfoldHistowoUnc)
    ratio_totalUnc = unfoldHisto.Clone()
    ratio_totalUnc.Divide(unfoldHistowoUnc)
    
    tmpPad2= pad2.DrawFrame( 0, 0., maxX, 1.9 )
    #print (labelX)
    tmpPad2.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    tmpPad2.GetYaxis().SetTitleOffset( 0.50 )
    #tmpPad2.GetYaxis().SetRangeUser(0.3,1.9 )
    
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.13, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    
    
    ratio_datastatUnc.SetFillColorAlpha(ROOT.kAzure+7,0.7)
    ratio_datastatUnc.SetLineColor(ROOT.kAzure+7)#,0.5)
    ratio_datastatUnc.SetLineColor(0)
    ratio_datastatUnc.SetLineWidth(0)
    ratio_datastatUnc.SetFillStyle(3245)
    ratio_totalUnc.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')# '#'+labelX.split('#')[1] )
    ratio_totalUnc.GetXaxis().SetTitleOffset( 0.9 )
    ratio_totalUnc.GetYaxis().SetTitle( "#frac{Sim.}{Data}" )
    ratio_totalUnc.GetYaxis().SetTitleOffset( 0.50 )

    ratio_totalUnc.GetYaxis().SetRangeUser(0.3,1.9 )

    ratio_totalUnc.GetYaxis().CenterTitle()
    ratio_totalUnc.GetXaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetXaxis().SetTitleSize(0.13)

    ratio_totalUnc.GetYaxis().SetLabelSize(0.12)
    ratio_totalUnc.GetYaxis().SetTitleSize(0.12)
    ratio_totalUnc.GetXaxis().SetNdivisions(505)
    ratio_totalUnc.GetYaxis().SetNdivisions(505)
    
    ratio_datastatUnc.SetMarkerStyle(0)
    ratio_datastatUnc.SetMarkerSize(0)

    ratio_totalUnc.SetFillColorAlpha(14,0.8)
    ratio_totalUnc.SetLineColor(14)
    ratio_totalUnc.SetLineColor(0)
    ratio_totalUnc.SetLineWidth(0)
    ratio_totalUnc.SetFillStyle(3354)
    ratio_totalUnc.SetMarkerStyle(0)
    ratio_totalUnc.SetMarkerSize(0)
    set_dynamic_y_range_errRatioHist(ratio_totalUnc,1.5,0.5)
    ratio_totalUnc.Draw('E2')
    ratio_datastatUnc.Draw('E2 SAME')
    
   

    hRatio = ROOT.TGraphAsymmErrors()
    hRatio.Divide( genJetHisto, unfoldHisto, 'pois' )
    hRatio.SetLineColor(colors[0])#ROOT.kRed)
    hRatio.SetMarkerColor(colors[0])#ROOT.kRed)
    #hRatio.SetLineWidth(1)
    hRatio.SetMarkerStyle(25)
    
    
    hRatio2 = ROOT.TGraphAsymmErrors()
    hRatio2.Divide( altMCHisto, unfoldHisto, 'pois' )
    hRatio2.SetLineColor(colors[1])#ROOT.kBlue)
    hRatio2.SetMarkerColor(colors[1])#ROOT.kBlue)
    #hRatio.SetLineWidth(1)
    hRatio2.SetMarkerStyle(25)
    if includeFSR: 
        hRatio3 = ROOT.TGraphAsymmErrors()
        hRatio3.Divide( fsrUpHisto, unfoldHisto, 'pois' )
        hRatio3.SetLineColor(46)
        hRatio3.SetMarkerColor(46)
        #hRatio.SetLineWidth(1)
        hRatio3.SetMarkerStyle(22)


        hRatio4 = ROOT.TGraphAsymmErrors()
        hRatio4.Divide( fsrDownHisto, unfoldHisto, 'pois' )
        hRatio4.SetLineColor(46)
        hRatio4.SetMarkerColor(46)
        #hRatio.SetLineWidth(1)
        hRatio4.SetMarkerStyle(23)
    
    if extraMC:
        

        hRatio5 = ROOT.TGraphAsymmErrors()
        hRatio5.Divide( altMC2Histo if 'dijet' in selection else altMC1Histo, unfoldHisto, 'pois' )
        hRatio5.SetLineColor(colors[2])#ROOT.kGray+4)
        hRatio5.SetMarkerColor(colors[2])#ROOT.kGray+4)
        #hRatio4.SetLineWidth(1)
        hRatio5.SetMarkerStyle(25)
        #hRatio5.Draw('P0 same')
    
    hRatio.SetMarkerSize(1.5)
    hRatio.Draw('P0 same')
    
    hRatio2.SetMarkerSize(1.5)
    hRatio2.Draw('P0 same')
    
    hRatio5.SetMarkerSize(1.5)
    hRatio5.Draw('P0 same')
    
    if includeFSR:
        hRatio3.SetMarkerSize(1.5)
        hRatio3.Draw('P0 same')

        hRatio4.SetMarkerSize(1.5)
        hRatio4.Draw('P0 same')
    
    
    ratioLegend=ROOT.TLegend(0.15,0.85,0.7,0.95)
    ratioLegend.SetTextSize(0.088)
    ratioLegend.SetNColumns(3)
    ratioLegend.SetFillColorAlpha(10,0.6)
    ratioLegend.SetBorderSize(0)
    #ratioLegend.SetTextSize(0.1)
    ratioLegend.AddEntry( ratio_totalUnc, 'Data total unc.', 'f' )
    ratioLegend.AddEntry( ratio_datastatUnc, 'Data stat. unc.', 'f' )
    #ratioLegend.AddEntry( ratiosystUncHisto, 'Syst.', 'f' )
    ratioLegend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    can.SaveAs(outputName)
    can.SaveAs(png)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12) 

    

def plot_combined_MCSelfClosure( unfolded_combined,
                                 combined_truth_hist,
                                 outputDir = '../Plots_January25_dijetSel_NDim/dijetSel/',
                                 outputFilename = 'combined_MCSelfClosure',
                                 selection='_dijetSel',
                                 labelX='',
                                 ext='.pdf',
                                 year='all',
                                 tlegendAlignment='right',
                                 process='MCSelfClosure',
                                 maxYFactor=1.04
                                ):
    
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    
    can = ROOT.TCanvas('can'+'SelfClosure', 'can'+'SelfClosure',  10, 10, 2000, 1500 )
    pad1 = ROOT.TPad("pad1"+'SelfClosure', "Main",0,0.3,1.00,1.00,-1)

    pad1.Draw()
    
    can.cd()
    pad1.cd()
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.02)
    if tlegendAlignment.startswith('right'): legend=ROOT.TLegend(0.66,0.61,0.88,0.89)

    else: legend=ROOT.TLegend(0.16,0.61,0.38,0.89)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.033)
    legend.SetBorderSize(0)


    unfold_integral=unfolded_combined.Integral()
    combined_truth_hist.SetLineColor(ROOT.kBlue)
    #combined_truth_hist.Scale(1./unfold_integral)
    combined_truth_hist.SetLineStyle(1)
    unfolded_combined.SetLineStyle(2)
    #unfolded_combined.Scale(1./unfold_integral)
    #legend.AddEntry( unfolded_combined, ('MG5-MLM+P8 (self-closure)' if process.startswith('MCSelfClosure') else 'MG5-MLM+P8 unf. w/ MG5-MLM+P8'), 'pe' )
    
    legend.AddEntry( unfolded_combined, 'MG5-MLM+P8 (self-closure)' if 'dijet' in selection.lower() else 'PWHG+P8 (self-closure)', 'pe' )
    legend.AddEntry( combined_truth_hist, 'MG5-MLM+P8 (gen)' if 'dijet' in selection.lower() else 'PWHG+P8 (gen)' , 'lp' )
    unfolded_combined.SetMarkerColor(ROOT.kRed)
    unfolded_combined.SetLineColor(ROOT.kRed)
    unfolded_combined.GetYaxis().SetTitleSize(0.05)
    #unfolded_combined.Draw()
    #unfolded_combined.GetYaxis().SetRangeUser(0., ( 2.2*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] )  ))
    #combined_truth_hist.SetMaximum(  1.56*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] )  )
    #combined_truth_hist.GetYaxis().SetRangeUser(0., (  maxYFactor*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] ) ))
    combined_truth_hist.SetMaximum( ( maxYFactor*max([ combined_truth_hist.GetMaximum(), unfolded_combined.GetMaximum()] ) ))
    combined_truth_hist.SetMinimum( 0.)
    #pad1.Modified()
    #pad1.Update()
    #can.Modified()
    #can.Update()
    
    
    #can.Update()
    #can.Modified()
    combined_truth_hist.GetYaxis().SetTitle('N_{events}')
    combined_truth_hist.GetYaxis().SetTitleSize(0.05)
    
    combined_truth_hist.Draw('histE')
    ROOT.TGaxis.SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0.005, "y")
    
    unfolded_combined.Draw('histE same')


    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.042)

    selText.SetNDC()
    
    dijetOffset = 0
    
    if selection.startswith("_dijet"): 
        seltext = 'Central Dijet'#( 'Central' if 'Central' in labelX  else 'Outer' )+' dijet region'
        dijetOffset = 0.15
    elif selection.startswith("_W"): seltext = 'Boosted W-enriched'
    elif selection.startswith("_top"): seltext = 'Boosted top-enriched'
    
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.87, seltext )

    selText = textBox.Clone()
    selText.SetTextFont(42)
    selText.SetTextSize(0.040)

    selText.SetNDC()


    #if selection.startswith("_dijet") and 'Central' in jetType : seltext = 'p_{T}>200 GeV' 
    if selection.startswith("_dijet"): seltext = 'p_{T}>200 GeV' 
    elif selection.startswith("_W"): seltext = 'p_{T}>200 GeV, 65<m_{jet}<125 GeV' 
    elif selection.startswith("_top"): seltext = 'p_{T}>400 GeV, 140<m_{jet}<300 GeV'
    #selText.DrawLatex( ( 0.65 if tlegendAlignment.startswith('right') else 0.2 ), 0.83, seltext )
    selText.DrawLatex( ( 0.19 if tlegendAlignment.startswith('right') else 0.55+dijetOffset ), 0.80, seltext )
    legend.Draw()
    CMS_lumi.extraText = "Simulation Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(pad1, 4, 0)
    can.cd()
    pad2 = ROOT.TPad("pad2"+'SelfClosure', "Ratio",0,0.00,1.00,0.30,-1);

    ROOT.gStyle.SetOptFit(1)
    pad2.SetGrid()
    pad2.SetTopMargin(0.)
    pad2.SetBottomMargin(0.3)
    pad2.Draw()
    pad2.cd()
    
    ratio = ROOT.TGraphAsymmErrors()
    ratio.Divide(combined_truth_hist,unfolded_combined,'pois')
    #set_dynamic_y_range_errRatioHist(ratio,1.5,0.5)
    
    
    tmpPad2= pad2.DrawFrame( 0, 0., combined_truth_hist.GetXaxis().GetBinLowEdge(combined_truth_hist.GetNbinsX()), 1.9 )
    #print (labelX)
    if 'dijet' in selection:
        tmpPad2.GetYaxis().SetRangeUser(0.9, 1.1 )
    else:
        
        tmpPad2.GetYaxis().SetRangeUser(0.7, 1.3 )
    

    tmpPad2.GetXaxis().SetTitleOffset( 0.9 )
    
    tmpPad2.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')# '#'+labelX.split('#')[1] )#.SetTitle(f'{}N-subjettiness basis')
    tmpPad2.GetYaxis().SetTitleOffset( 0.50 )
    #tmpPad2.GetYaxis().SetTitle( "#frac{Unf.}{Sim.}" )
    
    tmpPad2.GetYaxis().SetTitle( "#frac{Unf.}{Sim.}" )
    tmpPad2.GetYaxis().CenterTitle()
    tmpPad2.SetLabelSize(0.13, 'x')
    tmpPad2.SetTitleSize(0.12, 'x')
    tmpPad2.SetLabelSize(0.12, 'y')
    tmpPad2.SetTitleSize(0.12, 'y')
    tmpPad2.SetNdivisions(505, 'x')
    tmpPad2.SetNdivisions(505, 'y')
    pad2.Modified()
    pad2.Update()
    pad2.Draw()
    can.Update()
    #pad2.cd()

    #ratio = ROOT.TRatioPlot(combined_truth_hist,unfolded_combined,)
    #ratio.GetLowerPad().etYaxis().SetRangeUser(0.8, 1.2)
    
    ratio.SetMarkerSize(1.5)
    
    ratio.GetXaxis().SetTitle(f'{labelX} N-subjettiness basis')# '#'+labelX.split('#')[1] )
    ratio.GetXaxis().SetTitleOffset( 0.9 )
    ratio.GetYaxis().SetTitle( "#frac{Unf.}{Sim.}" )
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleOffset( 0.50 ) 
    ratio.GetXaxis().SetLabelSize(0.12)
    ratio.GetXaxis().SetTitleSize(0.13)

    ratio.GetYaxis().SetLabelSize(0.12)
    ratio.GetYaxis().SetTitleSize(0.12)
    
    ratio.Draw('PE1')
    
    can.SaveAs(outputDir+outputFilename+ext)
    ROOT.gStyle.SetPadRightMargin(0.09)     ## reseating
    ROOT.gStyle.SetPadLeftMargin(0.12) 
    
def drawUncertainties_from_err_shifts_theoryVariations_Ndim(ivar, 
                                                            unfoldHistoTotUnc, unfoldHistowoUnc, 
                                                            unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, 
                                                            unfoldHistoBkgSubUnc, uncerUnfoldHisto, 
                                                            cov_tot, cov_datastat_tot, 
                                                            cov_rmstat_tot, cov_bkg_tot, labelX, 
                                                            tlegendAlignment, 
                                                            outputName, unftot, selection, 
                                                            norming=True,
                                                            year='all'
                                                           ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing theory/model variation uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    colors_cr = [ROOT.kMagenta+2,  ROOT.kBlue-4, 433]  
    colors_syst = get_colour_palette_as_list('vf_8')[1:]#list(reversed())#[ROOT.kBlue+1, ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kGreen+2]
    
    #colors = get_colour_palette_as_list('vf_10')
    #[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    #ROOT.gStyle.SetPalette(len(colors),array('i', colors))
        
    upstyles =   [20,21,34,29,22,23,29,47,33,43,39,41,39,45,117,114,48]
    downstyles = [24,25,28,30,26,32,30,46,27,42,37,40,37,44,38 ,60 , 5 ]
    
    modelVariations = True if 'VariationUNC' in outputName else False
    
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    #otherUncs = OrderedDict()
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)
    
    #if tlegendAlignment.startswith('right'): 
    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    #unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoNoNorm.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    

    up_counter=0#1
    down_counter=0#1
    col_counter=0#1
    col_counter_jes=0
    

    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    
    set_dynamic_y_range_errRatioHist(totalErrHist,1.3)
    
    totalErrHist.GetXaxis().SetTitle(f'{labelX}')#6-body_OC')
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')

    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
   
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    

    cr_histos = OrderedDict()

    CR1_key=None
    CR2_key=None
    erdOn_key=None

    for k in uncerUnfoldHisto:

        if ('cr1' in k.lower() or 'cr2' in k.lower() or 'erd' in k.lower()) and '_shifthist' in k.lower():
            cr_histos[k] = uncerUnfoldHisto[k].Clone()
            cr_histos[k].Sumw2()
            cr_histos[k] = convert_syst_shift_to_error_ratio_hist(cr_histos[k].Clone(),
                                                                  unfoldHistoTotUnc.Clone())
            if 'cr1' in k.lower():
                CR1_key=k
            elif 'cr2' in k.lower():
                CR2_key=k
            elif 'erd' in k.lower():
                erdOn_key=k

            cr_histos[k].SetLineStyle(1)
            cr_histos[k].SetLineWidth(2)
            cr_histos[k].SetMarkerSize(0)
            cr_histos[k].SetFillColor(0)

    cr_keys = [CR1_key, CR2_key, erdOn_key]
    for i, key in enumerate(cr_keys):
        if key:
            color = colors_cr[i % len(colors_cr)]
            cr_histos[key].SetLineColor(color)
            cr_histos[key].SetMarkerColor(color)
            cr_histos[key].Draw('L same')

    if CR1_key: legend.AddEntry(cr_histos[CR1_key],'CR1', 'l')
    if CR2_key: legend.AddEntry(cr_histos[CR2_key],'CR2', 'l')
    if erdOn_key: legend.AddEntry(cr_histos[erdOn_key],'ERD on', 'l')
    
   
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),
                                                                                       unfoldHistoTotUnc.Clone())                            
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                up_counter=up_counter+1
    
    col_counter=0
    col_counter_jes=0
    
    syst_sources_up = list(normeduncerUnfoldHistoshiftsUp.keys())
    for i, k in enumerate(syst_sources_up):
        
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsUp[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(color)
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'h_{damp}', 'p' )
            elif 'CP5' in k: 
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'UE tune (CP5)', 'p' )
            elif 'mtop' in k:
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'm_{top}', 'p' )

            else: 
                print('else in th. syst comp maker', k)
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
                
    
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            
            if 'cr' in text.lower() or 'erd' in text.lower(): continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
                                                                                       
            if 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                down_counter=down_counter+1
    
    #set_palette_from_list(colors_syst)
    syst_sources_down = list(normeduncerUnfoldHistoshiftsDown.keys())
    col_counter=0
    for i, k in enumerate(syst_sources_down):
        if 'model'in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower(): 
            continue
        else:
            color = colors_syst[i % len(colors_syst) + (1 if len(colors_syst)>col_counter>0 else 0)]
            col_counter+=1
            normeduncerUnfoldHistoshiftsDown[k].SetLineColor(color)
            normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(color)

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'damp' in k or 'CP5' in k or 'mtop' in k:
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")

            else: 
                print('else in th. syst comp maker', k)
                normeduncerUnfoldHistoshiftsDown[k].Draw("P same")
        
    #legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    #legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
        

    
        
def drawUncertainties_from_err_shifts_Ndim(ivar, unfoldHistoTotUnc, unfoldHistowoUnc, unfoldHistoDataStatUnc, unfoldHistoRMStatUnc, unfoldHistoBkgSubUnc, uncerUnfoldHisto, cov_tot, cov_datastat_tot, cov_rmstat_tot, cov_bkg_tot, labelX, tlegendAlignment, outputName, unftot, selection, with_modelUnc=True, norming=False,year='all' ):
    
    #print('All uncertainty keys from uncerUnfoldHisto', uncerUnfoldHisto.keys())
    
    print (f'|------> Procesing uncertainty plot for {ivar} {"with" if norming else "without"} norming of err_shift_hists by unfolding total={unftot} ')
    
    colors_list = list(reversed(get_colour_palette_as_list('vf_10')[1:]))+[ROOT.TColor.GetColor('#c849a9'),61,30]
    
    colors = colors_list#[ 95, 7, 6, 38, 8, 42, 50, 218, 225, 30, 16, 51, 83, 61, 167, 207, 209, 212, 216, 198, 190, 67, 89, 133, 142, 208, 36, 2, 144, 225, 227, 150, 93, 40]
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.15)
        
    upstyles =   [20,21,22,29,23,34,47,33,43, 117,114,48]  #39,41, 45,
    downstyles = [24,25,26,30,32,28,46,27,42, 38 ,60 , 5 ]  #37,40, 44,
    
        
    modelkey=None
    JES_key=None
    JER_key=None
    btag_key=None
    btagUncIncluded=False
    
    normeduncerUnfoldHistoshiftsUp = OrderedDict()
    normeduncerUnfoldHistoshiftsDown = OrderedDict()
    otherUncs = OrderedDict()
    
    
    
    canUnc = ROOT.TCanvas('canUnc'+ivar, 'canUnc'+ivar,  10, 10, 1500, 1000 )
    canUnc.SetTopMargin(0.08)

    legend=ROOT.TLegend(0.2,0.65,0.9,0.9)

    #if tlegendAlignment.startswith('right'): 
    #    legend=ROOT.TLegend(0.2,0.65,0.8,0.9)
    #else: 
    #    legend=ROOT.TLegend(0.35,0.65,0.95,0.9)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    
    unfoldHistoNoNorm = unfoldHistoTotUnc.Clone()
    
    #unfoldHistowoUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoTotUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoDataStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoRMStatUnc.Scale(1./(unftot if norming else 1.),'width')#
    #unfoldHistoBkgSubUnc.Scale(1./(unftot if norming else 1.),'width')#
    
    unfoldHistoNoNorm.SetTitle("")
    unfoldHistowoUnc.SetTitle("")
    unfoldHistoTotUnc.SetTitle("")
    unfoldHistoDataStatUnc.SetTitle("")
    unfoldHistoRMStatUnc.SetTitle("")
    unfoldHistoBkgSubUnc.SetTitle("")
    
    unfoldHistoNoNorm.Sumw2()
    unfoldHistoTotUnc.Sumw2()
    unfoldHistoDataStatUnc.Sumw2()
    unfoldHistoRMStatUnc.Sumw2()
    unfoldHistoBkgSubUnc.Sumw2()
    unfoldHistowoUnc.Sumw2()
    
    jesHistoUpMax = unfoldHistoTotUnc.Clone('jesHistoUpMax')
    jesHistoUpMax.Reset()
    jesHistoDownMax = unfoldHistoTotUnc.Clone('jesHistoDownMax')
    jesHistoDownMax.Reset()
    
    jesHistoUpMax.Sumw2()
    jesHistoDownMax.Sumw2()
    
    if 'all' in year:
        jerHistoUpMax = unfoldHistoTotUnc.Clone('jerHistoUpMax')
        jerHistoUpMax.Reset()
        jerHistoDownMax = unfoldHistoTotUnc.Clone('jerHistoDownMax')
        jerHistoDownMax.Reset()

        jerHistoUpMax.Sumw2()
        jerHistoDownMax.Sumw2()
    
    
    #print(uncerUnfoldHisto.keys())
    for k in uncerUnfoldHisto:
        if 'modeltotal'in k.lower() and 'shifthist' in k.lower() and (modelkey==None) and with_modelUnc: 
            modelkey=k
        elif 'jes' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower() and not('const' in k.lower()) and (JES_key==None):
            JES_key=k
            print(JES_key)
            jesHistoUpMax = uncerUnfoldHisto[k].Clone()
            jesHistoUpMax.Sumw2()
            #jesHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoUpMax = convert_syst_shift_to_error_ratio_hist(jesHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jesHistoDownMax = uncerUnfoldHisto[k].Clone()
            jesHistoDownMax.Sumw2()
            #jesHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jesHistoDownMax = convert_syst_shift_to_error_ratio_hist(jesHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
        elif ('jer' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and ('all' in year) and (JER_key==None):
            JER_key=k
            print(JER_key)
            jerHistoUpMax = uncerUnfoldHisto[k].Clone()
            jerHistoUpMax.Sumw2()
            #jerHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoUpMax = convert_syst_shift_to_error_ratio_hist(jerHistoUpMax.Clone(), 
                                                                   unfoldHistoTotUnc.Clone())
            jerHistoDownMax = uncerUnfoldHisto[k].Clone()
            jerHistoDownMax.Sumw2()
            #jerHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            jerHistoDownMax = convert_syst_shift_to_error_ratio_hist(jerHistoDownMax.Clone(), 
                                                                     unfoldHistoTotUnc.Clone())
            
        elif ('btag' in k.lower() and 'shifthist' in k.lower() and 'total' in k.lower()) and (btag_key==None):
            btag_key = k
            btagUncIncluded=True

            if not('dijet' in selection):
                btagHistoUpMax = unfoldHistoTotUnc.Clone('btagHistoUpMax')
                btagHistoUpMax.Reset()
                btagHistoDownMax = unfoldHistoTotUnc.Clone('btagHistoDownMax')
                btagHistoDownMax.Reset()
                
            btagHistoUpMax = uncerUnfoldHisto[k].Clone()
            btagHistoUpMax.Sumw2()
            #btagHistoUpMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoUpMax = convert_syst_shift_to_error_ratio_hist(btagHistoUpMax.Clone(), 
                                                                    unfoldHistoTotUnc.Clone())
            btagHistoDownMax = uncerUnfoldHisto[k].Clone()
            btagHistoDownMax.Sumw2()
            #btagHistoDownMax.Scale(1./(unftot if norming else 1.),'width')
            btagHistoDownMax = convert_syst_shift_to_error_ratio_hist(btagHistoDownMax.Clone(), 
                                                                      unfoldHistoTotUnc.Clone())
    
    up_counter=0
    down_counter=0
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
        
        if ('shifthist' in k.lower() and 'up' in k.lower()):# and not k.endswith(('TotalUnc', 'SystTotal', 'StatTotal')) and not 'CM' in k:
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
            if 'btag' in k.lower(): 
            #    print(k)
            #    btagUncIncluded = True 
                continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsUp[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsUp[k].Sumw2()
            #normeduncerUnfoldHistoshiftsUp[k] = normalise_hist(normeduncerUnfoldHistoshiftsUp[k].Clone())
            #normeduncerUnfoldHistoshiftsUp[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsUp[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsUp[k].Clone(),                            
                                                                                       unfoldHistoTotUnc.Clone())
            
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#'BTAG' in text or 'LEPTON' in text 
                normeduncerUnfoldHistoshiftsUp[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsUp[k].SetMarkerStyle(upstyles[up_counter])
                if 'tau_2_2' in k: print (k,text, up_counter, col_counter,upstyles[up_counter],colors[col_counter])
                col_counter=col_counter+1    
                up_counter=up_counter+1
            
    #up_counter=1
    #down_counter=1
    col_counter=0
    col_counter_jes=0
    
    for k in uncerUnfoldHisto:
           
        if ('shifthist' in k.lower() and 'down' in k.lower()):
            
            if '_jes' in k.lower() or (('all' in year) and 'jer' in k.lower()):
                continue
                
            if 'btag' in k.lower(): continue
            #print(k)
            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper()  if not('ALL' in text.upper()) else text.upper().replace('ALL','')
            if 'cr' in text.lower() or 'erd' in text.lower() or 'model' in text.lower() or 'DAMP' in text or 'MTOP' in text or 'TUNE' in text:
                continue

            normeduncerUnfoldHistoshiftsDown[k] = uncerUnfoldHisto[k].Clone()
            normeduncerUnfoldHistoshiftsDown[k].Sumw2()
            #normeduncerUnfoldHistoshiftsDown[k] = normalise_hist(normeduncerUnfoldHistoshiftsDown[k].Clone())
            #normeduncerUnfoldHistoshiftsDown[k].Scale(1./(unftot if norming else 1.),'width')#./(unftot if norming else 1.)
            normeduncerUnfoldHistoshiftsDown[k] = convert_syst_shift_to_error_ratio_hist(normeduncerUnfoldHistoshiftsDown[k].Clone(),
                                                                                         unfoldHistoTotUnc.Clone())
              
            if 'ISR' in text or 'L1' in text or 'FSR' in text or ('JER' in text and not('all' in year)) or ('PU' in text and not('DAMP' in text)) or 'PDF' in text or 'const' in text.lower() or 'unclus' in text.lower():#r 'BTAG' in text or 'LEPTON' in text
                normeduncerUnfoldHistoshiftsDown[k].SetLineStyle(2 if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetLineColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerColor(colors[col_counter])
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerSize(2.0)# if not('L1' in text) else 1)
                normeduncerUnfoldHistoshiftsDown[k].SetMarkerStyle(downstyles[down_counter])
                if 'tau_2_2' in k: print (k,text, down_counter, col_counter,downstyles[down_counter],colors[col_counter])
                down_counter=down_counter+1
                col_counter=col_counter+1 
            
    
          
    #print ("Other uncs' keys", modelkey,btag_key)#,lepton_key)
    if with_modelUnc:
        modelUnc = uncerUnfoldHisto[modelkey].Clone()
        modelUnc.Sumw2()
        #modelUnc = normalise_hist(modelUnc.Clone())
        #modelUnc.Scale(1./(unftot if norming else 1.),'width')#
        modelUnc = convert_syst_shift_to_error_ratio_hist(modelUnc.Clone(), unfoldHistoTotUnc.Clone())
        modelUnc.SetLineStyle(1)
        modelUnc.SetLineWidth(2)
        modelUnc.SetMarkerSize(0)
        modelUnc.SetLineColor(28)
        #col_counter+=1
        modelUnc.SetFillColor(0)
    
    dataStatErrHist = unfoldHistoDataStatUnc.Clone()
    dataStatErrHist.Sumw2()
    rmStatErrHist = unfoldHistoRMStatUnc.Clone()
    rmStatErrHist.Sumw2()
    bkgSubErrHist = unfoldHistoBkgSubUnc.Clone()
    bkgSubErrHist.Sumw2()
    totalErrHist = unfoldHistoTotUnc.Clone()
    totalErrHist.Sumw2()

    dataStatErrHist.Divide(unfoldHistowoUnc)
    totalErrHist.Divide(unfoldHistowoUnc)
    
    totalErrHist.SetLineWidth(0)
    totalErrHist.GetYaxis().SetTitle('Variation/nominal')
    totalErrHist.GetYaxis().SetTitleSize(0.05)
    #if not('dijet' in selection): 
    #    totalErrHist.GetYaxis().SetRangeUser(0.3,1.8)
    #else:
    #    if 'all' in year:
    #        totalErrHist.GetYaxis().SetRangeUser(0.5,1.5)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    
    #    if '_2_3' in ivar or '_2_4' in ivar or '_2_5' in ivar or '_1p5_3' in ivar or '_1p5_4' in ivar or '_1p5_5' in ivar:
    #        totalErrHist.GetYaxis().SetRangeUser(0.4,1.6)
    #    else:
    #        totalErrHist.GetYaxis().SetRangeUser(0.7,1.45)
   
    set_dynamic_y_range_errRatioHist(totalErrHist,1.25 if ('dijet' in selection) else 1.4,0.95)
    
    totalErrHist.GetXaxis().SetTitle(f'{labelX}')#'6-body_OC')
    totalErrHist.SetLineWidth(0)
    totalErrHist.SetLineStyle(2)
    totalErrHist.SetFillColorAlpha(14,0.7)#ROOT.kGray+3
    totalErrHist.SetMarkerSize(0)
    totalErrHist.SetFillStyle(3354)
    totalErrHist.SetLineColor(14)#ROOT.kGray+3)
    totalErrHist.Draw(' E2')
    
    dataStatErrHist.SetLineWidth(0)
    dataStatErrHist.SetLineStyle(2)
    dataStatErrHist.SetMarkerSize(0)
    dataStatErrHist.SetFillStyle(3245)
    dataStatErrHist.SetFillColorAlpha(ROOT.kAzure+7,0.6)#ROOT.kAzure+7)
    dataStatErrHist.SetLineColor(ROOT.kAzure+7)#ROOT.kAzure+7)
    dataStatErrHist.Draw('E2 same')
    
    
    h1 = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),-1)
    rmStatErrHist = convert_error_bars_to_error_ratio_hist(rmStatErrHist.Clone(),1)

    rmStatErrHist.SetLineWidth(2)
    h1.SetLineWidth(2)
    rmStatErrHist.SetLineStyle(9)
    h1.SetLineStyle(9)
    h1.SetLineColor(1)
    rmStatErrHist.SetLineColor(1)
    h1.SetMarkerSize(0)
    rmStatErrHist.SetMarkerSize(0)
    rmStatErrHist.Draw('L same ')
    h1.Draw("L same")
    #h.Delete()
    
    h2 = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),-1)
    bkgSubErrHist = convert_error_bars_to_error_ratio_hist(bkgSubErrHist.Clone(),1)
    
    bkgSubErrHist.SetLineWidth(2)
    h2.SetLineWidth(2)
    bkgSubErrHist.SetLineStyle(7)
    h2.SetLineStyle(7)
    h2.SetLineColor(50)
    bkgSubErrHist.SetLineColor(50)
    h2.SetMarkerSize(0)
    bkgSubErrHist.SetMarkerSize(0)
    bkgSubErrHist.Draw('L same ')
    h2.Draw("L same")
    if with_modelUnc: modelUnc.Draw('L same')

    
    
    for k in otherUncs:
        if ('cr' in k.lower() or 'erd' in k.lower()): continue
        #print(k)
        text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
        print ("OtherUncs loop", text, k)
        #h0 = 0
        h0 = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),-1)
        otherUncs[k] = convert_error_bars_to_error_ratio_hist(otherUncs[k].Clone(),1)
        otherUncs[k].Draw('L same')
        h0.Draw('L same')
    
    
    
    for ibin in range(1,jesHistoUpMax.GetNbinsX()+1):
        
        upmax_ibin = 1.
        downmax_ibin = 0.
        diff = 0.
        
        upmax_ibin = jesHistoUpMax.GetBinContent(ibin)
        diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
        
        if (diff>=1. or diff<0.):
            print(f'WARNING: JES total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
            upmax_ibin = 0.
            downmax_ibin = 0.
            diff = 0.
        else:
            downmax_ibin = 1.-diff
        #print ("JES total",ibin, upmax_ibin,downmax_ibin)
        jesHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
        jesHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if 'all' in year:

        for ibin in range(1,jerHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = jerHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: JER total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            #print ("JER total",ibin, upmax_ibin,downmax_ibin)
            jerHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            jerHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)

    if not('dijet' in selection) and btagUncIncluded:

        for ibin in range(1,btagHistoUpMax.GetNbinsX()+1):
            upmax_ibin = 1.
            downmax_ibin = 0.
            diff = 0.
            upmax_ibin = btagHistoUpMax.GetBinContent(ibin)
            diff = upmax_ibin - 1. if upmax_ibin>1 else 1. - upmax_ibin
            
            if (diff>=1. or diff<0.):
                print(f'WARNING: b-tagging total contrib, diff.: {upmax_ibin,diff} is >=1 or <0 in bin {ibin}, setting to 0')
                upmax_ibin = 0.
                downmax_ibin = 0.
                diff = 0.
            else:
                downmax_ibin = 1.-diff
            #print ("b-tagging total",ibin, upmax_ibin,downmax_ibin)
            btagHistoUpMax.SetBinContent(ibin,1.+diff)#ibin,upmax_ibin)
            btagHistoDownMax.SetBinContent(ibin,1.-diff)#ibin,downmax_ibin)
       
    col_counter=9
    legend.AddEntry(jesHistoUpMax,'JES', 'p')
    jesHistoUpMax.SetLineColor(colors[col_counter])
    jesHistoDownMax.SetLineColor(colors[col_counter])
    jesHistoUpMax.SetMarkerColor(colors[col_counter])
    jesHistoDownMax.SetMarkerColor(colors[col_counter])
    jesHistoUpMax.SetMarkerStyle(39)#upstyles[up_counter])
    jesHistoDownMax.SetMarkerStyle(37)#downstyles[down_counter])
    jesHistoUpMax.SetMarkerSize(2.0)
    jesHistoDownMax.SetMarkerSize(2.0)
    up_counter+=1
    down_counter+=1
    col_counter+=1

    
    if 'all' in year:

        legend.AddEntry(jerHistoUpMax,'JER', 'p')
        jerHistoUpMax.SetLineColor(ROOT.kCyan+3)
        jerHistoDownMax.SetLineColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoDownMax.SetMarkerColor(ROOT.kCyan+3)
        jerHistoUpMax.SetMarkerStyle(41)#upstyles[up_counter])
        jerHistoDownMax.SetMarkerStyle(40)#downstyles[down_counter])
        jerHistoUpMax.SetMarkerSize(2.0)
        jerHistoDownMax.SetMarkerSize(2.0)
        
        up_counter+=1
        down_counter+=1
    
    
    if not('dijet' in selection):
        
        
        if btagUncIncluded:# and not(btag_key!=None):
            btagHistoUpMax.SetLineColor(colors[col_counter])
            btagHistoDownMax.SetLineColor(colors[col_counter])
            btagHistoUpMax.SetMarkerColor(colors[col_counter])
            btagHistoDownMax.SetMarkerColor(colors[col_counter])
            btagHistoUpMax.SetMarkerStyle(45)#upstyles[up_counter])
            btagHistoDownMax.SetMarkerStyle(44)#downstyles[down_counter])

            btagHistoUpMax.SetMarkerSize(2.0)
            btagHistoDownMax.SetMarkerSize(2.0)
            btagHistoUpMax.Draw('P same')
            btagHistoDownMax.Draw('P same')
            up_counter=up_counter+1
            down_counter=down_counter+1
            col_counter+=1
               
        
        #if not(lepton_key==None): 
        #    legend.AddEntry(leptonUp,'Lepton wt.', 'p')
    
    for k in normeduncerUnfoldHistoshiftsUp:
        if ('jes' in k.lower() and not('const' in k.lower())) or 'model'in k.lower() or 'tag' in k.lower() or 'bkg' in k.lower() or 'cr' in k.lower() or 'erd' in k.lower() or ('all' in year and 'jer' in k.lower()): 
            continue
        else:
            normeduncerUnfoldHistoshiftsUp[k].Draw("P same")

            text = (k.split('_shiftHist')[0].replace('Up','').replace('Down','').replace('Weight', '')).split(ivar+'_')[1]
            text=text.upper() if not('ALL' in text.upper()) else text.upper().replace('ALL','')

            if 'l1' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'L1 prefiring', 'p' )
            elif 'unclus' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], 'MET uncl. en.', 'p' )
                

            elif 'const' in k.lower():
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'neut' in k.lower():
                    text="Neutral ES"
                elif 'charg' in k.lower():
                    text="Charged ES"
                elif 'photon' in k.lower():
                    text="Photon ES"
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )


            else: 
                print('else in non-th. syst comp maker', k)
                
                normeduncerUnfoldHistoshiftsDown[k.replace('Up', 'Down')].Draw("P same")
                if 'asandpdf' in k.lower():
                    text = f"Scale and PDF"
                
                legend.AddEntry( normeduncerUnfoldHistoshiftsUp[k], text, 'p' )
        
        #print (text)
    jesHistoUpMax.Draw('P same')
    jesHistoDownMax.Draw('P same')
    if 'all' in year:
        jerHistoUpMax.Draw('P same')
        jerHistoDownMax.Draw('P same')

    if not('dijet' in selection):
        if btagUncIncluded: legend.AddEntry(btagHistoUpMax,'b-tagging', 'p')
    
    if with_modelUnc: legend.AddEntry( modelUnc, 'Shower & hadronization', 'l' )    
    legend.AddEntry( bkgSubErrHist, 'Background stat.', 'l' )    
    legend.AddEntry( rmStatErrHist, 'Response matrix stat.', 'l' )    
    legend.AddEntry( dataStatErrHist, 'Data stat.', 'f' )    
    legend.AddEntry( totalErrHist, 'Total uncertainty', 'f' )   
    
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_13TeV = "13 TeV, "+ ( '2016+2017+2018' if year.startswith('all') else year )
    CMS_lumi.relPosX = 0.10
    CMS_lumi.CMS_lumi(canUnc, 4, 0)
    #ROOT.gROOT.ForceStyle()
    #tdrstyle.setTDRStyle()
    #canUnc.SetLogy()
    canUnc.Update()
    
    legend.Draw()
    png = outputName.split('.pdf')[0]+'.png'
    canUnc.SaveAs(outputName)
    canUnc.SaveAs(png)
    


def combine_hist_2D_withUF(
                            aBlankTH2 = None,
                            RM_lists_dict=[],
                            gen_lists_dict=[],
                            truereco_lists_dict=[],
                
                            ):
    combined_response_matrix = aBlankTH2.Clone()
    gen_bin_offset = 0
    reco_bin_offset = 0
    
    if len(gen_lists_dict)==0:
        for i in range(len(RM_lists_dict)):
            gen_lists_dict.append(RM_lists_dict[i].ProjectionX(f'{i}_projX')) #just need for bin counting

    if len(truereco_lists_dict)==0:
        for i in range(len(RM_lists_dict)):
            truereco_lists_dict.append(RM_lists_dict[i].ProjectionY(f'{i}_projY')) #just need for bin counting
    
    for i in range(len(RM_lists_dict)):
        rm = RM_lists_dict[i]
        gen_bins = gen_lists_dict[i].GetNbinsX()+1
        reco_bins = truereco_lists_dict[i].GetNbinsX()+2

        for g in range(1,gen_bins):
            #print(g)#,gb,rb)
            for r in range(1,reco_bins):
                combined_response_matrix.SetBinContent(gen_bin_offset + g, reco_bin_offset + r, rm.GetBinContent(g, r))
                combined_response_matrix.SetBinError(gen_bin_offset + g, reco_bin_offset + r, rm.GetBinError(g, r))
            #print(r)    
            #include reco underflows for each observable (corrections for misreconstruction rate)
            combined_response_matrix.SetBinContent(gen_bin_offset + g, 0, rm.GetBinContent(g, 0))
            combined_response_matrix.SetBinError(gen_bin_offset + g, 0, rm.GetBinError(g, 0))

            #rb+=1
            #gb+=1
        #print(g)#,gb,rb)
        gen_bin_offset += gen_bins
        reco_bin_offset += reco_bins
    gc.collect()    
    return combined_response_matrix

def combine_all_histogram_types(
    allVarsDict, 
    varList,
    histTypes = [
        "reco", "truereco", "fakereco",
        "reco_genBin", "truereco_genBin", "fakereco_genBin",
        "gen", "accepgen", "missgen",
        "respWithMiss"
    ],
    samplePrefLabel="MLMQCD_HT2000toInf",  
    sel="_dijetSel",
    sysName = '_nom',
    extraRecoGap=2,     
    extraGenGap=1,      
    verbose=False,
    combined_var_name = "6bodyOC"
                                ):
    """
    Combines all sub-variable histograms (one dict entry per obs) into a set of 'combined' histograms plus a binMap for each category.

    Parameters
    ----------
    allVarsDict : dict[str -> dict[str->TH1 or TH2]]
        Outer dict: keys are obs names. 
        Value: a dict mapping histogramName->histogramObject for that obs.
        Example: allVars_signalHistos or allVars_dataHistos, etc.
    varList : list of str
        Which obs to combine (keys in allVarsDict).
    histTypes : list of str
        The histogram types we look foruse (provided in allVars dicts): reco, gen, respWithMiss, etc.
    samplePrefLabel : str
        The prefix used for sample naming, e.g. "MLMQCD_HT2000toInf" or "data" or "sysMLMQCD_..." 
        so that for each obs, we have keys like "MLMQCD_HT2000toInf_recoJet_tau_0p25_1_nom_dijetSel", etc.
    sel : str ---  "_dijetSel, _WSel, _topSel"
    extraRecoGap : int
        Number of empty bins to insert between sub-vars in the combined 'reco' dimension/axis.
    extraGenGap : int
        Number of empty bins to insert between sub-vars in the combined 'gen' dimension/axis.
    verbose : bool
        Print debug if True.

    Returns
    -------
    A dictionary with:
      {
        "combined_{histType}": TH1F (if found),
        "binMap_{histType}": { ivar: (startBin,endBin) },
        also, e.g.,
        "binMap_respWithMiss": { ivar: (genStart, genEnd, recoStart, recoEnd) }
      }
    """

    outDict = {}
    #Store one combined histogram (if found) per type of histogram in input dict.
    #Keep a bin map for each type, a la binMap_reco[ivar] = (startBin, endBin) for  1D histos
    #or for the 2D (resp. matrix), binMap_respWithMiss[ivar] = (genStart, genEnd, recoStart, recoEnd).

    #Keep track of how many bins each type needs in the reco dimension vs gen dimension for generalisation's sake
    #(currently hard-set to 2x finer reco bins in same global bin range as gen bins which have unit width, ie, reco bin widths uniformly 0.5)
    #For 1D reco-scheme-binned histograms, increment totalRecoBins_{type}.
    #For 1D gen-scheme-binned histograms, increment totalGenBins_{type}.
    #For 2D resp. increment on both axes.
    totalRecoBins = {}
    totalGenBins  = {}

    for t in histTypes:
        totalRecoBins[t] = 0
        totalGenBins[t]  = 0

    #Build record of which histName for each (ivar, type) to combine them afterwards.
    foundHists = {}
    hasResp = 0 
    for ivar in varList:
        subDict = allVarsDict[ivar]  
        for t in histTypes:
            # build the expected name pattern for histos with a small helper

            def build_expected_name(t, ivar, sel):
                """Return the suffix for the histo naming convention."""
                
                if "_genBin" in t:
                    # e.g. "MLMQCD_HT2000toInf_recoJet_tau_0p25_1_nom_dijetSel_genBin"
                    mainT = t.replace("_genBin", "")  
                    return f"_{mainT}{ivar}{sysName}{sel}_genBin"
                else:
                    # e.g. "MLMQCD_HT2000toInf_recoJet_tau_0p25_1_nom_dijetSel", "MLMQCD_HT2000toInf_respWithMissJet_tau_0p25_1_nom_dijetSel"
                    return f"_{t}{ivar}{sysName}{sel}"

            suffix = build_expected_name(t, ivar, sel)
            #print(suffix,samplePrefLabel)
            
            for histKey in subDict.keys():
                if histKey.startswith(samplePrefLabel) and suffix in histKey:
                    #print(histKey)
                    hObj = subDict[histKey]
                    foundHists[(ivar, t)] = hObj
                    if t == "respWithMiss":
                        hasResp+=1
                        # 2D
                        nX = hObj.GetNbinsX()
                        nY = hObj.GetNbinsY()
                        
                        totalGenBins[t]  += (nX + extraGenGap)
                        totalRecoBins[t] += (nY + extraRecoGap)
                    else:
                        # 1D
                        nBins = hObj.GetNbinsX()
                        if ("gen" in t.lower()) or (t in ["gen", "accepgen", "missgen"]):
                            totalGenBins[t] += (nBins + extraGenGap)
                        else:
                            totalRecoBins[t] += (nBins + extraRecoGap)

                    break  

    if verbose:
        print("[combine_all_histogram_types] Summed bin counts:")
        for t in histTypes:
            if totalRecoBins[t] or totalGenBins[t]:
                print(f"  type={t}, totalRecoBins={totalRecoBins[t]}, totalGenBins={totalGenBins[t]}")

    
    binMaps = {}
    for t in histTypes:
        if not any((ivar, t) in foundHists for ivar in varList):
            continue
        if t == "respWithMiss":
            nx = totalGenBins[t]
            ny = totalRecoBins[t]
            if nx < 1 or ny < 1:
                continue
            h2 = ROOT.TH2F(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}",
                           nx, 0, nx,
                           ny, 0, ny // 2 if ny>2 else ny)  # or just ny
            h2.Sumw2()
            outDict[f"combined_{t}"] = h2
            binMaps[t] = {}
        else:
            
            if ("gen" in t) or (t in ["gen", "accepgen", "missgen"]):
                nb = totalGenBins[t]
                h1 = ROOT.TH1F(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}", nb, 0, nb)
                h1.Sumw2()
                outDict[f"combined_{t}"] = h1
                binMaps[t] = {}
            else:
                nb = totalRecoBins[t]
                h1 = ROOT.TH1F(f"combined_{t}+{samplePrefLabel+sysName}", f"combined_{t}+{samplePrefLabel+sysName}", nb, 0, nb // 2 if nb>2 else nb)
                h1.Sumw2()
                outDict[f"combined_{t}"] = h1
                binMaps[t] = {}

    # Fill all combined histos, for each obs find the hist of a certain type, t, and if it exists, offset bins, copy contents/errors from original histos
    offsets_reco = {t:0 for t in histTypes}
    offsets_gen  = {t:0 for t in histTypes}
    #rm_gen_bin_offset = 0 
    #rm_reco_bin_offset = 0
    
    
    for ivar in varList:
        for t in histTypes:
            if (ivar, t) not in foundHists:
                continue
            hObj = foundHists[(ivar, t)]
            # copy bins to outDict
            if t == "respWithMiss":
                #h2_comb = outDict[f"combined_{t}"]
                #nX = hObj.GetNbinsX()+1
                #nY = hObj.GetNbinsY()+2
                
                #gxOff = rm_gen_bin_offset#offsets_gen[t]
                #ryOff = rm_reco_bin_offset#offsets_reco[t]
                #binMaps[t][ivar] = (gxOff+1, gxOff+nX, ryOff+1, ryOff+nY)
                
                
                #for gx in range(1, nX):
                #    for ry in range(1, nY):
                #        c = hObj.GetBinContent(gx, ry)
                #        e = hObj.GetBinError(gx, ry)
                #        h2_comb.SetBinContent(gxOff + gx, ryOff + ry, c)
                #        h2_comb.SetBinError(gxOff + gx, ryOff + ry, e)
                # handle misreconstructed gen in reco UF
                #h2_comb.SetBinContent(gxOff + gx, 0, hObj.GetBinContent(g, 0))
                #h2_comb.SetBinError(gxOff + gx, 0, hObj.GetBinError(g, 0))

                #offsets_gen[t]  += (nX + extraGenGap)
                #offsets_reco[t] += (nY + extraRecoGap)
                #rm_gen_bin_offset += nX
                #rm_reco_bin_offset += nY
                continue
            else:
                # 1D
                h1_comb = outDict[f"combined_{t}"]
                nBins   = hObj.GetNbinsX()
                if ("gen" in t) or (t in ["gen", "accepgen", "missgen"]):
                    baseOff = offsets_gen[t]
                    binMaps[t][ivar] = (baseOff+1, baseOff + nBins)
                    for iBin in range(1, nBins+1):
                        c = hObj.GetBinContent(iBin)
                        e = hObj.GetBinError(iBin)
                        h1_comb.SetBinContent(baseOff + iBin, c)
                        h1_comb.SetBinError(baseOff + iBin, e)
                    # Insert gap bins
                    for g in range(nBins+1, nBins+1 + extraGenGap):
                        h1_comb.SetBinContent(baseOff + g, 0)
                        h1_comb.SetBinError(baseOff + g, 0)
                    offsets_gen[t] += (nBins + extraGenGap)
                else:
                    baseOff = offsets_reco[t]
                    binMaps[t][ivar] = (baseOff+1, baseOff + nBins)
                    for iBin in range(1, nBins+1):
                        c = hObj.GetBinContent(iBin)
                        e = hObj.GetBinError(iBin)
                        h1_comb.SetBinContent(baseOff + iBin, c)
                        h1_comb.SetBinError(baseOff + iBin, e)
                    for g in range(nBins+1, nBins+1 + extraRecoGap):
                        h1_comb.SetBinContent(baseOff + g, 0)
                        h1_comb.SetBinError(baseOff + g, 0)
                    offsets_reco[t] += (nBins + extraRecoGap)
    #2D
    if hasResp>1:
        RM_lists_dict = []
        gen_lists_dict = []
        truereco_lists_dict = [] 

        for ivar in varList:
            for t in histTypes:
                if (ivar, t) not in foundHists:
                    continue
                if '_gen' in t and not ('genBin' in t):
                    gen_lists_dict.append(foundHists[(ivar, t)].Clone())
                elif '_truereco' in t and not ('genBin' in t):
                    truereco_lists_dict.append(foundHists[(ivar, t)].Clone())
                elif 'respWithMiss' in t:
                    RM_lists_dict.append(foundHists[(ivar, t)].Clone())


        outDict[f"combined_{t}"]  = combine_hist_2D_withUF( outDict[f"combined_{t}"].Clone(outDict[f"combined_{t}"].GetName()+'_blankClone'),
                                                            RM_lists_dict,
                                                            gen_lists_dict,
                                                            truereco_lists_dict,
                                                          
                                                          )
    
    # Return combined 1- and/or 2-D histos
    
    finalDict = {}
    for t in histTypes:
        if f"combined_{t}" in outDict:
            
            if not('genBin' in t): 
                finalDict[f"combined_{samplePrefLabel}_{t}{combined_var_name}{sysName}{sel}"] = outDict[f"combined_{t}"]
            else:
                finalDict[f"combined_{samplePrefLabel}_{t.split('_genBin')[0]}{combined_var_name}{sysName}{sel}_genBin"] = outDict[f"combined_{t}"]

    for t in binMaps:
        
        if binMaps[t]: 
            if not('genBin' in t):
                finalDict[f"binMap_{samplePrefLabel}_{t}{combined_var_name}{sysName}{sel}"] = binMaps[t]
            else:
                finalDict[f"binMap_{samplePrefLabel}_{t.split('_genBin')[0]}{combined_var_name}{sysName}{sel}_genBin"] = binMaps[t]
    return finalDict

def build_combined_covariance_matrix(
    hist_list,
    correlation_matrix,
    use_off_diag_corr=False
):
    """
    Takes:
      hist_list           : list of 1D ROOT histograms (TH1F, etc.),
      correlation_matrix  : 2D numpy array describing correlation among obs corresponding to input hists in list,
      use_off_diag_corr   : bool; if False, off-diagonal blocks are set to zero instead
                            of correlation_matrix[i, j].

    Returns:
      A TH2D (combined_cov_hist) representing the combined covariance matrix
      for all histograms in 'hist_list' after combining them in one big 1D histo.
      Diagonal (blocks) contains each 1D histogram's bin variances.
      Off-diagonal blocks incorporate correlations between observables or are zeroed on use_off_diag_corr input value (default=False),
      If zero, just one big diagonal input covariance is returned.
    """

    # For each histogram, build its diagonal (co)variance array
    # (ie, bin error^2 per bin along diagonal of new combined cov). 
    cov_matrices = []
    n_bins_total = 0

    for reco_hist in hist_list:
        n_bins = reco_hist.GetNbinsX() + 2
        n_bins_total += n_bins
        
        #for blocks on a given diagonal
        cov_matrix = np.zeros((n_bins, n_bins))
        for i in range(1, n_bins):
            error = reco_hist.GetBinError(i)
            cov_matrix[i-1, i-1] = error**2
        cov_matrices.append(cov_matrix)
        
    print(n_bins_total)
    
    # Create square, combined covariance as a TH2 with dimension n_bins_total 
    # over range of global bins [0, n_bins_total//2] used for reco axes in other
    # 1-/2-D combined hists
    combined_cov_hist = ROOT.TH2D(
        "combined_cov_matrix",
        "Combined Covariance Matrix",
        n_bins_total, 0, n_bins_total/2,
        n_bins_total, 0, n_bins_total/2
    )
    combined_cov_hist.Sumw2()

    # Fill the diagonal blocks from each individual observables' 1-D, bin-wise variances
    bin_offset = 0
    for cov_matrix in cov_matrices:
        n_bins = cov_matrix.shape[0]
        for i in range(n_bins):
            for j in range(n_bins):
                combined_cov_hist.SetBinContent(
                    bin_offset + i,
                    bin_offset + j,
                    cov_matrix[i, j]
                )
        bin_offset += n_bins

    # Fill the off-diagonal blocks, and entries in blocks, using the correlation_matrix for the observables
    # unless 'use_off_diag_corr' is False, then corr is set to 0.
    bin_offset_i = 0
    for i in range(len(cov_matrices)):
        cov_matrix_i = cov_matrices[i]
        n_bins_i = cov_matrix_i.shape[0]
        bin_offset_j = 0
        for j in range(len(cov_matrices)):
            cov_matrix_j = cov_matrices[j]
            n_bins_j = cov_matrix_j.shape[0]

            if i != j:
                for k in range(n_bins_i):
                    for l in range(n_bins_j):
                        if use_off_diag_corr:
                            corr = correlation_matrix[i, j]
                        else:
                            corr = 0.0
                        combined_cov_value = corr * np.sqrt(cov_matrix_i[k, k] * cov_matrix_j[l, l])
                        # Apply +1 offset in bin indexing for ROOT
                        combined_cov_hist.SetBinContent(
                            bin_offset_i + k + 1,
                            bin_offset_j + l + 1,
                            combined_cov_value
                        )
            bin_offset_j += n_bins_j
        bin_offset_i += n_bins_i

    return combined_cov_hist


import os
from PIL import Image
from PyPDF2 import PdfMerger, PdfReader, PdfWriter

def png_to_pdf(png_files, output_file):
    images = []
    for file in png_files:
        img = Image.open(file)
        if img.mode == 'RGBA':
            img = img.convert('RGB')
        images.append(img)

    images[0].save(output_file, save_all=True, append_images=images[1:])

def merge_pdfs(pdf_files, output_file):
    merger = PdfMerger()
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(output_file)
    merger.close()

def combine_images_and_pdfs(input_files, output_file):
    png_files = [file for file in input_files if file.lower().endswith('.png')]
    pdf_files = [file for file in input_files if file.lower().endswith('.pdf')]

    temp_pdf_files = []
    if png_files:
        temp_pdf = "temp_images.pdf"
        png_to_pdf(png_files, temp_pdf)
        temp_pdf_files.append(temp_pdf)
    
    all_pdf_files = pdf_files + temp_pdf_files
    if all_pdf_files:
        merge_pdfs(all_pdf_files, output_file)

    for temp_file in temp_pdf_files:
        os.remove(temp_file)

def generate_latex_table(observables, data, column_titles, filename=None):
    """
    Generates a LaTeX table based on the input data.

    Args:
    observables (list of str): List of row identifiers (first column, like "Observable").
    data (list of lists): List of lists containing the data for each row (each list represents a row).
    column_titles (list of str): List of column titles (like 'Observable', '$\chi^2$', 'ndf', '$\frac{\chi^2}{ndf}$').
    filename (str): If provided, saves the LaTeX table to a file with the given name. Otherwise, it prints it.

    Returns:
    str: LaTeX table code.
    """
    
    latex_table = "\\begin{table}[ht!]\n\\centering\n\\begin{tabular}{|" + "c|" * len(column_titles) + "}\n\\hline\n"
    
    latex_table += " & ".join(column_titles) + " \\\\\n\\hline\n"
    
    for i, observable in enumerate(observables):
        row_data = " & ".join(map(str, data[i]))
        latex_table += f"{observable} & {row_data} \\\\\n"
    
    latex_table += "\\hline\n\\end{tabular}\n\\caption{Your table caption here.}\n\\label{tab:your_label}\n\\end{table}"

    if filename:
        with open(filename, 'w') as file:
            file.write(latex_table)
        print(f"LaTeX table saved to {filename}.")
    else:
        print(latex_table)
    
    return latex_table


def extendTH1(h, extendUF=True, extendOF=True):
    """
    Given a TH1 (e.g. TH1F), return a new TH1 with extra bins
    that include the underflow (UF) and/or overflow (OF) entries.
    
    The extra bin(s) will have the same width as the nominal first
    (for underflow) and last (for overflow) bins.
    
    Parameters:
      h         : The original TH1 histogram.
      extendUF  : If True, add an extra (leftmost) bin for underflow.
      extendOF  : If True, add an extra (rightmost) bin for overflow.
    
    Returns:
      A new TH1 histogram with the “extended” x‐axis.
    """
    # Number of nominal bins
    n = h.GetNbinsX()
    axis = h.GetXaxis()
    x_min = axis.GetXmin()  # lower edge of first nominal bin
    x_max = axis.GetXmax()  # upper edge of last nominal bin
    first_bin_width = axis.GetBinWidth(1)
    last_bin_width = axis.GetBinWidth(n)
    
    # Get the original bin edges.
    # If the histogram was created with nonuniform binning, GetXbins() returns
    # a TArrayD of bin edges (of size n+1). Otherwise it is empty.
    bins_array = axis.GetXbins()
    if bins_array.GetSize() > 0:
        # Non-uniform binning: get edges for bins 1..n+1.
        orig_edges = [axis.GetBinLowEdge(i) for i in range(1, n+2)]
    else:
        # Uniform binning.
        orig_edges = [x_min + i*(x_max - x_min)/n for i in range(0, n+1)]
    
    # Build the new edge array.
    new_edges = []
    if extendUF:
        new_edges.append(x_min - first_bin_width)
    new_edges.extend(orig_edges)
    if extendOF:
        new_edges.append(x_max + last_bin_width)
    new_edges_arr = array('d', new_edges)
    
    # Create new histogram (here TH1F is assumed; for TH1D use TH1D).
    h_ext = ROOT.TH1F(h.GetName() + "_ext", h.GetTitle() + " (extended)", len(new_edges_arr)-1, new_edges_arr)
    h_ext.Sumw2()  # preserve Sumw2
    
    # Offset in bin numbering: if we extended UF then new bin 1 is the UF bin.
    offset = 1 if extendUF else 0
    
    # Copy underflow if extended.
    if extendUF:
        h_ext.SetBinContent(1, h.GetBinContent(0))
        h_ext.SetBinError(1, h.GetBinError(0))
    else:
        h_ext.SetBinContent(0, h.GetBinContent(0))
        h_ext.SetBinError(0, h.GetBinError(0))
    
    # Copy the nominal bins.
    for i in range(1, n+1):
        new_bin = i + offset
        h_ext.SetBinContent(new_bin, h.GetBinContent(i))
        h_ext.SetBinError(new_bin, h.GetBinError(i))
    
    # Copy overflow if extended.
    if extendOF:
        new_bin = n + 1 + offset
        h_ext.SetBinContent(new_bin, h.GetBinContent(n+1))
        h_ext.SetBinError(new_bin, h.GetBinError(n+1))
    else:
        h_ext.SetBinContent(h_ext.GetNbinsX()+1, h.GetBinContent(n+1))
        h_ext.SetBinError(h_ext.GetNbinsX()+1, h.GetBinError(n+1))
    
    return h_ext
'''
def extendTH2(h, extendUF_x=True, extendOF_x=True, extendUF_y=True, extendOF_y=True):
    """
    Given a TH2 (e.g. TH2F), return a new TH2 histogram whose binning in x and y is either
    “extended” (i.e. the original underflow and/or overflow bins are inserted as extra visible bins)
    or, if not requested, left in the conventional underflow/overflow locations.
    
    Parameters:
      h          : The original TH2 histogram.
      extendUF_x : If True, include the original x-axis underflow as the first visible bin.
                   If False, leave it as an underflow (non‐visible) bin.
      extendOF_x : If True, include the original x-axis overflow as the last visible bin.
                   If False, leave it as an overflow (non‐visible) bin.
      extendUF_y : If True, include the original y-axis underflow as the first visible bin.
                   If False, leave it as an underflow (non‐visible) bin.
      extendOF_y : If True, include the original y-axis overflow as the last visible bin.
                   If False, leave it as an overflow (non‐visible) bin.
    
    Returns:
      A new TH2 histogram with the modified binning and with bin contents (and errors) copied appropriately.
    """
    from array import array

    n_x = h.GetNbinsX()
    axisX = h.GetXaxis()
    x_min = axisX.GetXmin()
    x_max = axisX.GetXmax()
    first_bin_width_x = axisX.GetBinWidth(1)
    last_bin_width_x = axisX.GetBinWidth(n_x)
    binsX = axisX.GetXbins()
    if binsX.GetSize() > 0:
        orig_edges_x = [axisX.GetBinLowEdge(i) for i in range(1, n_x+2)]
    else:
        orig_edges_x = [x_min + i*(x_max - x_min)/n_x for i in range(0, n_x+1)]
    new_edges_x = []
    # If we want to “extend” the underflow, add an extra bin on the left.
    if extendUF_x:
        new_edges_x.append(x_min - first_bin_width_x)
    new_edges_x.extend(orig_edges_x)
    # Similarly for overflow.
    if extendOF_x:
        new_edges_x.append(x_max + last_bin_width_x)
    new_edges_x_arr = array('d', new_edges_x)
    new_n_x = len(new_edges_x_arr) - 1  

    n_y = h.GetNbinsY()
    axisY = h.GetYaxis()
    
    y_min = axisY.GetXmin()
    y_max = axisY.GetXmax()
    first_bin_width_y = axisY.GetBinWidth(1)
    last_bin_width_y = axisY.GetBinWidth(n_y)
    binsY = axisY.GetXbins()
    if binsY.GetSize() > 0:
        orig_edges_y = [axisY.GetBinLowEdge(i) for i in range(1, n_y+2)]
    else:
        orig_edges_y = [y_min + i*(y_max - y_min)/n_y for i in range(0, n_y+1)]
    new_edges_y = []
    if extendUF_y:
        new_edges_y.append(y_min - first_bin_width_y)
    new_edges_y.extend(orig_edges_y)
    if extendOF_y:
        new_edges_y.append(y_max + last_bin_width_y)
    new_edges_y_arr = array('d', new_edges_y)
    new_n_y = len(new_edges_y_arr) - 1  # number of visible y bins

    h_ext = ROOT.TH2F(h.GetName() + "_ext", h.GetTitle() + " (extended)",
                      new_n_x, new_edges_x_arr,
                      new_n_y, new_edges_y_arr)
    h_ext.Sumw2()

    
    for i_new in range(1, new_n_x+1):
        # Determine the original x bin index corresponding to new bin i_new.
        if extendUF_x:
            if i_new == 1:
                orig_i = 0           # first visible bin is original underflow
            elif extendOF_x and i_new == new_n_x:
                orig_i = n_x + 1     # last visible bin is original overflow
            else:
                orig_i = i_new - 1   # in between: shift by 1 because UF was added
        else:
            # If not extending underflow, the visible bins are just the nominal ones,
            # except that if overflow is extended, the last visible bin comes from the original OF.
            if extendOF_x and i_new == new_n_x:
                orig_i = n_x + 1
            else:
                orig_i = i_new
        for j_new in range(1, new_n_y+1):
            # Determine the original y bin index.
            if extendUF_y:
                if j_new == 1:
                    orig_j = 0
                elif extendOF_y and j_new == new_n_y:
                    orig_j = n_y + 1
                else:
                    orig_j = j_new - 1
            else:
                if extendOF_y and j_new == new_n_y:
                    orig_j = n_y + 1
                else:
                    orig_j = j_new
            new_bin = h_ext.GetBin(i_new, j_new)
            h_ext.SetBinContent(new_bin, h.GetBinContent(orig_i, orig_j))
            h_ext.SetBinError(new_bin, h.GetBinError(orig_i, orig_j))

    # --- Copy non-extended underflow/overflow bins ---
    # For any axis that was NOT extended, copy the original underflow/overflow bins
    # into the corresponding non-visible bins of h_ext.

    # For x-axis:
    if not extendUF_x:
        for j in range(0, h_ext.GetNbinsY()+2):
            h_ext.SetBinContent(0, j, h.GetBinContent(0, j))
            h_ext.SetBinError(0, j, h.GetBinError(0, j))
    if not extendOF_x:
        for j in range(0, h_ext.GetNbinsY()+2):
            h_ext.SetBinContent(h_ext.GetNbinsX()+1, j, h.GetBinContent(n_x+1, j))
            h_ext.SetBinError(h_ext.GetNbinsX()+1, j, h.GetBinError(n_x+1, j))
    # For y-axis:
    if not extendUF_y:
        for i in range(0, h_ext.GetNbinsX()+2):
            h_ext.SetBinContent(i, 0, h.GetBinContent(i, 0))
            h_ext.SetBinError(i, 0, h.GetBinError(i, 0))
    if not extendOF_y:
        for i in range(0, h_ext.GetNbinsX()+2):
            h_ext.SetBinContent(i, h_ext.GetNbinsY()+1, h.GetBinContent(i, n_y+1))
            h_ext.SetBinError(i, h_ext.GetNbinsY()+1, h.GetBinError(i, n_y+1))

    return h_ext
'''

def extendTH2(h, extendUF_x=True, extendOF_x=True, extendUF_y=True, extendOF_y=True):
    """
    Create a new TH2F whose binning in x and y optionally extends underflow/overflow
    into the visible range. Corner bins (x-flow, y-flow) follow the chosen setting on
    each axis independently:
      - If an axis is extended for underflow/overflow, that flow bin becomes visible
        along that axis.
      - If an axis is NOT extended, that flow bin remains in underflow/overflow for
        that axis.
    Thus, a partially extended corner (e.g. x=overflow is extended, y=underflow is not)
    ends up in the new histogram at (x=lastVisibleBin, y=underflowBin=0), without merging
    into the first visible y bin.

    Parameters:
      h          : The original TH2 histogram (TH2F assumed).
      extendUF_x : If True, x underflow becomes the first visible bin in x.
      extendOF_x : If True, x overflow becomes the last visible bin in x.
      extendUF_y : If True, y underflow becomes the first visible bin in y.
      extendOF_y : If True, y overflow becomes the last visible bin in y.

    Returns:
      A new TH2F with the extended axes and contents/errors correctly placed, including
      partial-flow corners.
    """
    import math
    from array import array

    n_x = h.GetNbinsX()
    n_y = h.GetNbinsY()

    axisX = h.GetXaxis()
    axisY = h.GetYaxis()

    # Original bin edges (x)
    x_min = axisX.GetXmin()
    x_max = axisX.GetXmax()
    binsX = axisX.GetXbins()
    if binsX.GetSize() > 0:
        orig_edges_x = [axisX.GetBinLowEdge(i) for i in range(1, n_x+2)]
    else:
        # uniform binning
        orig_edges_x = [x_min + i*(x_max - x_min)/n_x for i in range(0, n_x+1)]
    first_bin_width_x = axisX.GetBinWidth(1)
    last_bin_width_x  = axisX.GetBinWidth(n_x)

    # Original bin edges (y)
    y_min = axisY.GetXmin()
    y_max = axisY.GetXmax()
    binsY = axisY.GetXbins()
    if binsY.GetSize() > 0:
        orig_edges_y = [axisY.GetBinLowEdge(i) for i in range(1, n_y+2)]
    else:
        # uniform binning
        orig_edges_y = [y_min + i*(y_max - y_min)/n_y for i in range(0, n_y+1)]
    first_bin_width_y = axisY.GetBinWidth(1)
    last_bin_width_y  = axisY.GetBinWidth(n_y)

    # Build the new bin-edge arrays
    new_edges_x = []
    if extendUF_x:
        new_edges_x.append(orig_edges_x[0] - first_bin_width_x)
    new_edges_x.extend(orig_edges_x)
    if extendOF_x:
        new_edges_x.append(orig_edges_x[-1] + last_bin_width_x)

    new_edges_y = []
    if extendUF_y:
        new_edges_y.append(orig_edges_y[0] - first_bin_width_y)
    new_edges_y.extend(orig_edges_y)
    if extendOF_y:
        new_edges_y.append(orig_edges_y[-1] + last_bin_width_y)

    new_edges_x_arr = array('d', new_edges_x)
    new_edges_y_arr = array('d', new_edges_y)

    # Number of visible bins in new histogram
    n_new_x = len(new_edges_x_arr) - 1
    n_new_y = len(new_edges_y_arr) - 1

    # Create new histogram
    h_ext = ROOT.TH2F(h.GetName()+"_ext", h.GetTitle()+" (extended)",
                      n_new_x, new_edges_x_arr,
                      n_new_y, new_edges_y_arr)
    h_ext.Sumw2()

    #------------------------------------------------------------
    # 1) Helper to map old bin index (0..n+1) -> new bin index
    #    (0..n_new+1). If that axis is extended for underflow,
    #    old underflow(0) -> new bin 1, else -> 0. Similarly for
    #    overflow (n+1). Nominal bins map to [1..n] or [2..n+1].
    #------------------------------------------------------------
    def map_axis_bin(old_bin, n, extendUF, extendOF):
        # old_bin can be 0..(n+1)
        # new_n = # of visible bins = n + (1 if UF extended) + (1 if OF extended)
        n_visible = n + (1 if extendUF else 0) + (1 if extendOF else 0)

        if old_bin == 0:   # underflow
            return 1 if extendUF else 0
        elif old_bin == n+1:  # overflow
            return n_visible if extendOF else (n_visible + 1)
        else:
            # nominal bin => shift by +1 if we extended the underflow
            offset = 1 if extendUF else 0
            newb = old_bin + offset
            return newb

    #------------------------------------------------------------
    # 2) Fill all bins (including corners) in ONE pass.
    #    We loop over old bins [0..n_x+1, 0..n_y+1] and add them
    #    into the new histogram bin that corresponds.
    #    This automatically handles partial corners without merging.
    #------------------------------------------------------------
    for old_i in range(0, n_x+2):  # x in [0..n_x+1]
        for old_j in range(0, n_y+2):  # y in [0..n_y+1]
            c  = h.GetBinContent(old_i, old_j)
            ce = h.GetBinError(old_i, old_j)
            if c == 0 and ce == 0:
                continue  # skip empty to save a bit of time

            # map to new bin indices
            i_new = map_axis_bin(old_i, n_x, extendUF_x, extendOF_x)
            j_new = map_axis_bin(old_j, n_y, extendUF_y, extendOF_y)

            # accumulate (in case multiple old bins map to the same new bin)
            old_c  = h_ext.GetBinContent(i_new, j_new)
            old_ce = h_ext.GetBinError(i_new, j_new)
            new_c  = old_c + c
            new_ce = math.sqrt(old_ce*old_ce + ce*ce)
            h_ext.SetBinContent(i_new, j_new, new_c)
            h_ext.SetBinError(i_new, j_new, new_ce)

    return h_ext





def extendTH2_old(h, extendUF_x=True, extendOF_x=True, extendUF_y=True, extendOF_y=True):
    """
    Given a TH2 (e.g. TH2F), return a new TH2 with extra bins along the x- and/or y-axes
    that include the underflow and/or overflow entries.
    
    For each axis the extra bin has the same width as the first (for underflow)
    or last (for overflow) nominal bin.
    
    Parameters:
      h          : The original TH2 histogram.
      extendUF_x : If True, add an extra x-bin at the left for underflow.
      extendOF_x : If True, add an extra x-bin at the right for overflow.
      extendUF_y : If True, add an extra y-bin at the bottom for underflow.
      extendOF_y : If True, add an extra y-bin at the top for overflow.
    
    Returns:
      A new TH2 histogram with extended x- and y-axes.
    """
    # --- Process the X axis ---
    n_x = h.GetNbinsX()
    axisX = h.GetXaxis()
    x_min = axisX.GetXmin()
    x_max = axisX.GetXmax()
    first_bin_width_x = axisX.GetBinWidth(1)
    last_bin_width_x = axisX.GetBinWidth(n_x)
    binsX = axisX.GetXbins()
    if binsX.GetSize() > 0:
        orig_edges_x = [axisX.GetBinLowEdge(i) for i in range(1, n_x+2)]
    else:
        orig_edges_x = [x_min + i*(x_max-x_min)/n_x for i in range(0, n_x+1)]
    new_edges_x = []
    if extendUF_x:
        new_edges_x.append(x_min - first_bin_width_x)
    new_edges_x.extend(orig_edges_x)
    if extendOF_x:
        new_edges_x.append(x_max + last_bin_width_x)
    new_edges_x_arr = array('d', new_edges_x)
    new_n_x = len(new_edges_x_arr) - 1

    # --- Process the Y axis ---
    n_y = h.GetNbinsY()
    axisY = h.GetYaxis()
    y_min = axisY.GetXmin()  # for Y axis, use GetXmin()/GetXmax() as well
    y_max = axisY.GetXmax()
    first_bin_width_y = axisY.GetBinWidth(1)
    last_bin_width_y = axisY.GetBinWidth(n_y)
    binsY = axisY.GetXbins()
    if binsY.GetSize() > 0:
        orig_edges_y = [axisY.GetBinLowEdge(i) for i in range(1, n_y+2)]
    else:
        orig_edges_y = [y_min + i*(y_max-y_min)/n_y for i in range(0, n_y+1)]
    new_edges_y = []
    if extendUF_y:
        new_edges_y.append(y_min - first_bin_width_y)
    new_edges_y.extend(orig_edges_y)
    if extendOF_y:
        new_edges_y.append(y_max + last_bin_width_y)
    new_edges_y_arr = array('d', new_edges_y)
    new_n_y = len(new_edges_y_arr) - 1

    # Create new TH2 histogram (here TH2F is assumed)
    h_ext = ROOT.TH2F(h.GetName() + "_ext", h.GetTitle() + " (extended)",
                      new_n_x, new_edges_x_arr,
                      new_n_y, new_edges_y_arr)
    h_ext.Sumw2()

    # Offsets in bin numbering: if underflow is extended then nominal bins start at new bin index 2.
    offset_x = 1 if extendUF_x else 0
    offset_y = 1 if extendUF_y else 0

    # Loop over the bins of the new histogram (which now includes the flow bins as actual bins)
    for i_new in range(1, new_n_x+1):
        # Map new bin i_new to original x bin index:
        if extendUF_x and i_new == 1:
            orig_i = 0      # original underflow bin in x
        elif extendOF_x and i_new == new_n_x:
            orig_i = n_x + 1  # original overflow bin in x
        else:
            orig_i = i_new - offset_x  # nominal bin (1..n_x)
        for j_new in range(1, new_n_y+1):
            if extendUF_y and j_new == 1:
                orig_j = 0      # original underflow bin in y
            elif extendOF_y and j_new == new_n_y:
                orig_j = n_y + 1  # original overflow bin in y
            else:
                orig_j = j_new - offset_y
            # Get the content and error from the original histogram.
            content = h.GetBinContent(orig_i, orig_j)
            err = h.GetBinError(orig_i, orig_j)
            new_bin = h_ext.GetBin(i_new, j_new)
            h_ext.SetBinContent(new_bin, content)
            h_ext.SetBinError(new_bin, err)
    return h_ext


def computeOOA(resp2D):
    """
    Computes per-bin out of acceptance correction factors of the form:
       (SB-only gen counts) / (nom+SB gen counts)
    by summing over all reco bins for each gen bin.
    
   
    """
    effCorr = []
    for i in range(0,resp2D.GetNbinsX()+1):
        # Numerator: integral over Y= nominal region only
        # Denominator: integral over Y= nominal + SB region
        
        num = resp2D.Integral(i, i+1, resp2D.GetNbinsY(), resp2D.GetNbinsY()+1)     # "nom" portion in Y 
        den = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY() + 1) # "nom + SB" portion

        if den != 0:
            effCorr.append(num / den)
        else:
            effCorr.append(1.0)  # or 0.0, depending on your preference
    return effCorr




def computeAcceptance(resp2D):
    """
    Computes per-bin acceptance correction factors of the form:
       (nom-only gen counts) / (nom+SB gen counts)
    by summing over all reco bins for each gen bin.
    
   
    """
    effCorr = []
    for i in range(0,resp2D.GetNbinsX()+1):
        # Numerator: integral over Y= nominal region only
        # Denominator: integral over Y= nominal + SB region
        
        num = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY())     # "nom" portion in Y 
        den = resp2D.Integral(i, i+1, 0, resp2D.GetNbinsY() + 1) # "nom + SB" portion

        if den != 0:
            effCorr.append(num / den)
        else:
            effCorr.append(1.0)  # or 0.0, depending on your preference
    return effCorr


def applyAcceptanceCorrection(hist1D, effCorr):
    """
    Scales each bin of 'hist1D' by the corresponding acceptance factor
    stored in 'effCorr'.  
    
    
    """
    
    for i in range(0, hist1D.GetNbinsX() + 2):
        scaleFactor = effCorr[i - 1]  # match bin i to effCorr[i-1]
        val  = hist1D.GetBinContent(i)
        err  = hist1D.GetBinError(i)
        hist1D.SetBinContent(i, val * scaleFactor)
        hist1D.SetBinError(i,  err * scaleFactor)
        
        
def applyInvAcceptanceCorrection(hist1D, effCorr):
    """
    Scales each bin of 'hist1D' by the corresponding acceptance factor
    stored in 'effCorr'.  
    
    
    """
    
    for i in range(0, hist1D.GetNbinsX() + 2):
        scaleFactor = effCorr[i - 1]  # match bin i to effCorr[i-1]
        val  = hist1D.GetBinContent(i)
        err  = hist1D.GetBinError(i)
        hist1D.SetBinContent(i, val * (1.-scaleFactor))
        hist1D.SetBinError(i,  err * (1.-scaleFactor))
        
        
def make_rel_uncertainty_plot(
                              dummy_unf_histo,
                              rel_unc_dict,
                              total_unc=None,
                              total_unc_label="Total",
                              outfilename="relative_uncertainties.pdf",
                              canvas_title="Relative Uncertainties",
                              x_axis_title="Jet Observable",#nSubVariables_dijetSel['Jet_tau_0p5_2']['label'],
                              y_axis_title="Relative uncertainty [%]",
                              y_max=100.0):
    """
    Produce a plot of relative uncertainties vs. bin center, 
    optionally include a band for the total unc. and lines 
    for individual sources.

    Parameters
    ----------
    
    rel_unc_dict : dict
        Dictionary of { "label": rel_unc_array }, where rel_unc_array has length Nbins
        and each element is the relative uncertainty for that bin. 
        These will be drawn as lines on the plot.
    total_unc : 1D array_like, optional
        If provided, length = Nbins, this will be drawn as a gray band 
        representing the total uncertainty. 
    outfilename : str, optional
        PDF file to save the plot into.
    canvas_title : str, optional
        Title for the top of the canvas.
    x_axis_title : str, optional
        X-axis label.
    y_axis_title : str, optional
        Y-axis label.
    y_max : float, optional
        Maximum on the y-axis for plotting. Adjust as needed.
    """
    #bin_edges = np.array(bin_edges, dtype=float)
    n_bins = dummy_unf_histo.GetNbinsX()
    
    n_bins = dummy_unf_histo.GetNbinsX()
    if total_unc is not None and len(total_unc) != n_bins:
        raise ValueError("Length of total_unc array does not match number of bins in hist_binning!")
    
    for label, arr in rel_unc_dict.items():
        if len(arr) != n_bins:
            raise ValueError(f"Length of array for '{label}' does not match number of bins in hist_binning!")

    
    c = ROOT.TCanvas("c","c",1500,1500)
    c.SetMargin(0.13,0.03,0.12,0.07)  # left, right, bottom, top

    
    frame_histo = dummy_unf_histo.Clone("frameHisto")#ROOT.TH1F("frame_histo", "", n_bins, bin_edges)
    frame_histo.Reset("ICE")
    frame_histo.SetTitle("")
    frame_histo.GetXaxis().SetTitle(x_axis_title)
    frame_histo.GetYaxis().SetTitle(y_axis_title)
    frame_histo.GetYaxis().SetRangeUser(0.0, y_max)
    frame_histo.Draw("AXIS")


    legend = ROOT.TLegend(0.65, 0.60, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    total_band = None
    if total_unc is not None:
        x_vals = []
        y_vals = []
        ex_vals = []
        ey_vals = []
        
        total_band = frame_histo.Clone("total_unc_band")
        total_band.Reset("ICE")
        
        for i in range(1, n_bins+1):
            val_percent = total_unc[i-1]*100.0
            total_band.SetBinContent(i, val_percent)
            

        '''
                    ROOT.TGraphErrors(n_bins, 
                                       np.array(x_vals, dtype=float), 
                                       np.array(y_vals, dtype=float),
                                       np.array(ex_vals, dtype=float),
                                       np.array(ey_vals, dtype=float))
        '''
        total_band.SetFillColor(14)
        
        total_band.SetFillColorAlpha(14, 0.5)
        total_band.SetLineColor(14)
        total_band.SetMarkerColor(14)
        
        total_band.SetLineWidth(1)
        #total_band.SetFillStyle(3254)  
        total_band.Draw("hist same")      

        legend.AddEntry(total_band, total_unc_label, "f")

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta+1,
              ROOT.kOrange+1, ROOT.kAzure+2, ROOT.kTeal+1, ROOT.kViolet+1, ROOT.kGray+2,
              ROOT.kMagenta
              
             ]
    styles = [1, 2, 3, 4, 5, 6, 7, 9, 10, 2, 3]  

    graphs = []
    color_index = 0
    style_index = 0

    for label, unc_array in rel_unc_dict.items():
        h = frame_histo.Clone(label.replace(' ',''))
        
        for i in range(1, n_bins+1):
            h.SetBinContent(i,unc_array[i-1]* 100.0)
            h.SetBinError(i,0.000001)
            

        h.SetMarkerStyle(1)
        h.SetMarkerColor(colors[color_index % len(colors)])
        h.SetLineColor(colors[color_index % len(colors)])
        h.SetLineStyle(styles[style_index % len(styles)])
        h.SetLineWidth(2)

        color_index += 1
        style_index += 1

        h.Draw("E same")
        legend.AddEntry(h, label, "l")
        graphs.append(h)

    legend.Draw()

    c.SaveAs(outfilename)

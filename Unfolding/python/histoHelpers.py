from collections import OrderedDict
import ROOT
import math
import numpy as np
import array
from array import array
import bisect
#from legend import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
####gReset()
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()
import os, glob, sys

from root_numpy import array2hist,hist2array
#from unfoldingPlottersAndHelpers import makeJacobian


########################################################################
###################### Histo modification helpers ######################
########################################################################
def makeJacobian(aTH1, aJac):
    N = aTH1.Integral(0,aTH1.GetNbinsX())
    for i in range(0,aTH1.GetNbinsX()):
        for j in range(0,aTH1.GetNbinsX()):
            if i==j: 
                aJac[i][j]=(N-aTH1.GetBinContent(i))/N/N 
                
            else: 
                aJac[i][j]=(-1.*aTH1.GetBinContent(i))/N/N 

                
def correlation_from_covariance(covariance,hist): #pythonic, need to rootify
    covariance = hist2array(covariance)
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return array2hist(correlation,hist)


def getRebinnedRescaled_TH1(hist_name, file_name, xrange, scale, isMC=True):
        
    FILE = ROOT.TFile.Open(file_name,'read')
    #tree = FILE.Get('jetObservables')
    hist = FILE.Get(hist_name)
    hist.SetDirectory(0)
    
    ROOT.TH1.AddDirectory(ROOT.kFALSE);
    
    bins = array('d', np.array(Bin(xrange)))
    
    newHist = ROOT.TH1F(hist.Rebin(len(bins)-1, "%s_rebin"%hist_name, bins))
    #newHist.SetStats(ROOT.kFALSE)
    if isMC==True: newHist.Scale(scale)
    newHist.SetDirectory(0)

    return newHist

def getRebinned_TH1(hist, xrange, scale, isMC=True):
    
    hist.SetDirectory(0)
    
    ROOT.TH1.AddDirectory(ROOT.kFALSE);
    
    bins = array('d', np.array(Bin(xrange)))
    
    newHist = ROOT.TH1F(hist.Rebin(len(bins)-1, "%s_rebin"%hist_name, bins))
    newHist.SetDirectory(0)
    
    return newHist


def rebin_RM_withUF(h_resp,genBin,recoBin):
    #### fancy way to create variable binning TH2D
    tmpHisto = ROOT.TH2F( h_resp.GetName()+"_Rebin", h_resp.GetName()+"_Rebin", len(genBin)-1, array( 'd', genBin), len(recoBin)-1, array( 'd', recoBin) )
    tmpHisto.Sumw2()                        

    for biny in range( 0, h_resp.GetNbinsY()+2 ):
        by = h_resp.GetYaxis().GetBinCenter( biny )
        for binx in range( 0, h_resp.GetNbinsX()+2 ):
            bx = h_resp.GetXaxis().GetBinCenter(binx)
            if not(binx==0 or biny==0 or binx==h_resp.GetNbinsX()+1 or biny==h_resp.GetNbinsY()+1):

                for iX in range( len(genBin)-1 ):
                    for iY in range( len(recoBin)-1 ):
                        if (bx<genBin[iX+1] and bx>genBin[iX]) and (by<recoBin[iY+1] and by>recoBin[iY]):
                            jbin = h_resp.GetBin(binx,biny)
                            
                            tmpHisto.SetBinContent( tmpHisto.GetBin(iX+1,iY+1), tmpHisto.GetBinContent(iX+1,iY+1)+h_resp.GetBinContent( jbin )) #tmpArrayContent[binx-1][biny-1] )
                            tmpHisto.SetBinError( tmpHisto.GetBin(iX+1,iY+1), np.sqrt((tmpHisto.GetBinError(iX+1,iY+1))**2.+(h_resp.GetBinError( jbin ))**2.))#tmpArrayError[binx-1][biny-1] ) )
            else:
                if not(binx==0 or biny==0):# or binx==h_resp.GetNbinsX()+1 or biny==h_resp.GetNbinsY()+1): 
                    continue

                jbin = h_resp.GetBin(binx,biny)

                if h_resp.IsBinUnderflow(jbin):

                    j=0 
                    for i in range(tmpHisto.GetNbinsX()+2):
                        if (bx>tmpHisto.GetXaxis().GetBinLowEdge(i) and bx<tmpHisto.GetXaxis().GetBinLowEdge(i+1)):
                            tmpHisto.SetBinContent(i,j, tmpHisto.GetBinContent(i,j)+h_resp.GetBinContent(jbin))
                            tmpHisto.SetBinError(i,j, np.sqrt((tmpHisto.GetBinError(i,j))**2.+(h_resp.GetBinError(jbin))**2.))
    #tmpHisto.Sumw2()
    tmpHisto.SetDirectory(0)
    #h_resp = copy.deepcopy(tmpHisto.Clone())

    return tmpHisto


#Helper functions used from Robin's old repo: https://github.com/raggleton/QGAnalysisPlotting

def renorm(arr2d, axis):
    # create version where each axis summed to 1
    # use where and out args to ensure nans are made into 0s
    summed = arr2d.sum(axis=axis, keepdims=True)
    
    return np.divide(arr2d, summed, where=summed!=0, out=np.zeros_like(arr2d))

def concat_row(arr2d, row_ind):
    # concat row row_ind + row_ind+1
    nrows, ncols = arr2d.shape
    if row_ind > nrows - 2:
        raise IndexError("Cannot concat row [%d] as only %d rows in matrix" % (row_ind, nrows))
    arr2d_new = np.zeros(shape=(nrows-1, ncols), dtype=float)
    new_row = arr2d[row_ind] + arr2d[row_ind+1]
    arr2d_new[row_ind] = new_row
    # fill in new matrix
    if row_ind > 0:
        # do that bit before the new row
        arr2d_new[:row_ind, ] = arr2d[:row_ind, ]
    if row_ind < nrows - 2:
        arr2d_new[row_ind+1:, :] = arr2d[row_ind+2:, :]
    return arr2d_new


def rebin_2d_hist(h2d, new_binning_x, new_binning_y):
    """Rebin a 2D histogram according to specific bin edges for x & y axes
    new_binning_x, new_binning_y are lists of tuple pairs of bin edges
    e.g. [(0, 1), (1, 4), (4, 10)]
    """
    #print("rebinning...")
    # convert pairs of bins to list of edges, including upper edge of last bin
    bin_edges_x = [b[0] for b in new_binning_x]
    bin_edges_x.append(new_binning_x[-1][1])

    bin_edges_y = [b[0] for b in new_binning_y]
    bin_edges_y.append(new_binning_y[-1][1])

    #print("rebin_2d_hist, new axes:", bin_edges_x, bin_edges_y)

    new_h2d = ROOT.TH2D(
        h2d.GetName()+"Rebin",
        ';'.join([h2d.GetTitle(), h2d.GetXaxis().GetTitle(), h2d.GetYaxis().GetTitle()]),
        len(new_binning_x),
        array('d', bin_edges_x),
        len(new_binning_y),
        array('d', bin_edges_y)
    )

    # Get original bin edges
    bins_x_orig = get_bin_edges(h2d, 'X')
    bins_y_orig = get_bin_edges(h2d, 'Y')

    # Get original bin contents
    #import pdb; pdb.set_trace()

    #print("getting OG BC")
    arr, err = th2_to_np_arr(h2d, errCalc=True)

    # Get map of old bin edges -> new bin edges
    # -1 to start at 0
    bin_groups_x = np.digitize(bins_x_orig, bin_edges_x) - 1
    bin_groups_y = np.digitize(bins_y_orig, bin_edges_y) - 1

    # Count cumulative number of entries in each group (remove last bin as upper edge)
    group_counts_x = np.bincount(bin_groups_x)[:-1].cumsum()
    group_counts_y = np.bincount(bin_groups_y)[:-1].cumsum()

    group_counts_x = np.insert(group_counts_x, 0, 0)
    group_counts_y = np.insert(group_counts_y, 0, 0)

    # Iterate over each group, sum, set new TH2 contents
    #print("setting TH2 contents/errs")
    for xind, (xl, xh) in enumerate(zip(group_counts_x[:-1], group_counts_x[1:]), 1):
        for yind, (yl, yh) in enumerate(zip(group_counts_y[:-1], group_counts_y[1:]), 1):
            new_bin_content = arr[yl:yh,xl:xh].sum()
            new_bin_err = np.sqrt(np.power(err[yl:yh,xl:xh], 2).sum())
            #if xind==yind: print(xind,yind,new_bin_content,new_bin_err)
            new_h2d.SetBinContent(xind, yind, new_bin_content)
            new_h2d.SetBinError(xind, yind, new_bin_err)

    #print("...done rebinning")
    return new_h2d

def get_bin_edge_pairs(new_bins,old_reco_bin_edges):
    bin_edge_pairs = [list(x) for x in zip(new_bins[:-1], new_bins[1:])]

    bin_edge_pairs[-1][1] = old_reco_bin_edges[-1]
    
    return bin_edge_pairs

def make_rebinned_2d_hist(h2d, new_binning, use_half_width_y=False):
    """Rebin 2D histogram using new binning.
    new_binning is list of tuple pairs of bin edges
    e.g. [(0, 1), (1, 4), (4, 10)]
    If use_half_width_y=False, uses new_binning for both x (gen) & y (reco) axes
    If True, creates bins that are half the width of new_binning for
    the y axes (reco)
    """
    if use_half_width_y:
        # create half width bins from new_binning
        reco_binning = []
        reco_bin_edges = get_bin_edges(h2d, 'Y')
        print("Creating half-bin width...")
        for s, e in new_binning:
            ideal_mid = (s+e)/2.
            # find the bin that closest matches this ideal_mid
            mid = reco_bin_edges[bisect.bisect_left(reco_bin_edges, ideal_mid)]
            reco_binning.append((s, mid))
            reco_binning.append((mid, e))
        print("...done")
        return rebin_2d_hist(h2d, new_binning, reco_binning)
    else:
        return rebin_2d_hist(h2d, new_binning, new_binning)
    
            
def get_bin_edges(hist, axis):
    """Get array of bin edges from hist. Must specify which axis to use."""
    axis = axis.lower()
    if axis not in ['x', 'y']:
        raise RuntimeError("get_bin_edges axis must be x or y")
    ax = hist.GetXaxis() if axis == "x" else hist.GetYaxis()
    bins = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins()+2)]
    return bins


def th2_to_np_arr(h,errCalc=True):
    #print("converting TH2 to np array for rebinning or calculations")
    array = np.zeros((h.GetNbinsY(), h.GetNbinsX()), dtype=float)
    if errCalc: errors = np.zeros((h.GetNbinsY(), h.GetNbinsX()), dtype=float)
    #print("converting TH2 to np array for rebinning or calculations")
        
    for ix in range(1, h.GetNbinsX() + 1):
        for iy in range(1, h.GetNbinsY() + 1):
            #if ix==iy:
            #    #print(ix,iy)
            #    #print(ix,iy,h.GetBinContent(ix, iy),h.GetBinError(ix, iy))
            array[iy-1][ix-1] = h.GetBinContent(ix, iy)
            if errCalc:errors[iy-1][ix-1] = h.GetBinError(ix, iy)
    #print("CONVERTED!")
    if errCalc: return array, errors
    else: return array,[]

def get_genORreco_bins_from_resp(resp, axis='gen'):
    return resp.sum(axis=0 if axis.lower() == 'gen' else 1)
           
            
#####################old, but some still used: taken from others as named below################################

#Taken from Christine/Ashley: 
#https://gitlab.cern.ch/asparker/QJetMass/-/blob/master/unfold/cppImplementation/do2DTUnfolding_Dec15lep.py#L1790
def tmatrixdsparse_to_ndarray(matrix):
    ndarr = np.zeros(shape=(matrix.GetNrows(), matrix.GetNcols()))

    rows_A = matrix.GetRowIndexArray()
    cols_A = matrix.GetColIndexArray()
    data_A = matrix.GetMatrixArray()
    for iy in range(matrix.GetNrows()):
        for indexA in range(rows_A[iy], rows_A[iy+1]):
            ix = cols_A[indexA]
            # print([x for x in self.GetXToHist()])
            # TODO: care about orientation?
            ndarr[iy, ix] = data_A[indexA]
    return ndarr
                
def GetNormalizedTMatrixandTH2( ath2, matrixname , aunfoldedth1 ):
    # should use matrixname  to give uique names to histograms created in loop

    J = ROOT.TMatrixD( aunfoldedth1.GetNbinsX(), aunfoldedth1.GetNbinsX())
    makeJacobian(aunfoldedth1 , J)

    cov_m,_ = th2_to_ndarray(ath2.Clone())
    covnormth2,_ = th2_to_ndarray(ath2.Clone())
    J = ROOT.TH2D(J.Clone())
    J, _ = th2_to_ndarray(J.Clone())

    covnorm = J@cov_m@J.T
    covnormth2 = ndarray_to_th2(covnorm)

    covnorm = ndarray_to_th2(covnorm)
    covnorm = th2_to_tmatrixd(covnorm)
    #covnorm_temp = ROOT.TMatrixD( cov_m, ROOT.TMatrixD.kMultTranspose, J)
    #print ("Transform total error matrix to normalized space using Jacobian ")
    #covnorm = ROOT.TMatrixD(J, ROOT.TMatrix.kMult, covnorm_temp )
    #print ("fill normalized th2 for later use...")

    #for xbin in range(0, ath2.GetNbinsX() +1):
    #    for ybin in range(0,ath2.GetNbinsY() +1  ):
    #        covnormth2.SetBinContent(xbin, ybin , covnorm[xbin][ybin])  
    #ErrorMatrixToHist(covnorm)
    #print ("th2 -> (TMatrixD ->normalized ->th2)-> [ tmatrix, th2 ]")
    #print (covnormth2)
    #print (matrixname)

    return [ covnorm.Clone() , covnormth2.Clone()  ]
    
#https://gitlab.cern.ch/asparker/QJetMass/-/blob/master/unfold/cppImplementation/do2DTUnfolding_Mar15_BottomLine.py#L119
def ConvertTH2toTMatrix( someth2 ):
    acov_m = ROOT.TMatrixD( someth2.GetNbinsX()+1, someth2.GetNbinsY()+1 )

    for xbin in range(0, someth2.GetNbinsX()+1 ):
        for ybin in range(0,someth2.GetNbinsY()+1  ):
            if xbin > someth2.GetNbinsX() :
                print (xbin)
                print (ybin)
            acov_m[xbin][ybin] = someth2.GetBinContent(  xbin, ybin)
    print (" finished filling TMatrixD from TH2 "  )
    return acov_m
                
def ConvertTH1toTMatrix(  someth2 ):
    acov_m1 = ROOT.TVectorD( someth2.GetNbinsX()+1 )

    for xbin in range(0, someth2.GetNbinsX()+1 ):
            #for ybin in range(0,someth2.GetNbinsY()+1  ):
            if xbin > someth2.GetNbinsX() :
                print (xbin)
                #print ybin
            acov_m1[xbin] = someth2.GetBinContent(  xbin)
    print (" finished filling TVectorD from TH1 "  )
    return acov_m1


def NormYourHisto(hist = None ) :
    ### Make a new empty TH1 so the OFL and UFL will be zero
    name = hist.GetTitle()+'_normed'
    histn = hist.Clone(name)
    integral=0
    histn.Reset()
    
    for h in range(0,hist.GetNbinsX()+1) :
        bc = hist.GetBinContent(h) 
        be = hist.GetBinError(h)
        histn.SetBinContent(h ,bc)    
        histn.SetBinError(h, be)  
        integral+=bc
    
    for h in range(0,histn.GetNbinsX()+1) :
        bc = histn.GetBinContent(h) 
        be = histn.GetBinError(h)
        #print bw
        if integral > 0. :
            bc = bc/integral
            be = be/integral            
        histn.SetBinContent(h , bc)    
        histn.SetBinError(h, be)  
    
    return histn


#Taken from Robin: https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/my_unfolder.py#L2358

def normalise_hist(h):
    if h.Integral() > 0:
        h.Scale(1./h.Integral())
    return h
def hist_divide_bin_width(h):
    """Create copy of hist, but each bin's contents is divide by the bin width"""
    h_new = h.Clone(h.GetName()+"DivideBinWidth")
    h_new.Scale(1., 'width')
    # if any bin has 0 entries before, it will now have nan, so we need to manually fix that
    for i in range(1, h_new.GetNbinsX()+1):
        if np.isnan(h_new.GetBinContent(i)):
            h_new.SetBinContent(i, 0)
            h_new.SetBinError(i, 0)
    h_new.SetEntries(h.GetEntries())  # needed as the default is to replace with integral
    return h_new

def normalise_hist_divide_bin_width(h):
    h_new = h.Clone(h.GetName()+"DivideBinWidth")
    normalise_hist(h_new)
    h_new = hist_divide_bin_width(h_new)
    return h_new

def get_syst_shifted_hist(syst_shift, unfolded=None):
    """Get histogram with systematic shift applied to bin contents
    Can specify starting hist, otherwise assumes unfolded w/no error
    """
    # TODO: syst_label -> ExpSystematic obj
    hist_shift = syst_shift.Clone(syst_shift.GetTitle()+'_shiftedbyNom')
    hist_shift.Add(unfolded)  # TODO what about errors?
    return hist_shift

def get_syst_error_hist(syst_shift, unfolded=None):
    """Get histogram with systematic shift as error bars
    Can specify starting hist, otherwise assumes unfolded w/no error
    """
    this_syst = self.get_exp_syst(syst_label)
    if this_syst.syst_error_bar is None:
        ref = unfolded or self.get_unfolded_with_no_errors()
        new_hist = self.convert_error_shift_to_error_bars(ref, self.get_syst_shift(syst_label))
        this_syst.syst_error_bar = new_hist
    return this_syst.syst_error_bar


def convert_error_bars_to_error_shift(h):
        """Create histogram with bin contents equal to error bar on h,
        and 0 error bars"""
        h_new = h.Clone(get_unique_str())
        for i in range(1, h.GetNbinsX()+1):
            h_new.SetBinContent(i, h.GetBinError(i))
            h_new.SetBinError(i, 0)
        return h_new


def convert_error_shift_to_error_bars(h_nominal, h_shift):
    """Create histogram with bin contents from h_nominal,
    and error bars from bin values of h_shift"""
    h = h_nominal.Clone(get_unique_str())
    for i in range(1, h_nominal.GetNbinsX()+1):
        h.SetBinError(i, h_shift.GetBinContent(i))
    return h

def convert_syst_shift_to_error_ratio_hist(h_syst, h_nominal):
    """Create h_syst / h_nominal without error bars"""
    h_new = h_syst.Clone(h_syst.GetName() + get_unique_str())
    h_new.Divide(h_nominal)
    for ix in range(1, h_new.GetNbinsX()+1):
        if h_nominal.GetBinContent(ix) == 0:
            h_new.SetBinContent(ix, 1)
        if h_syst.GetBinContent(ix) < 0:
            h_new.SetBinContent(ix, 1)
            print("Warning: _convert_syst_shift_to_error_ratio_hist bin", ix, "of h_syst, %s, < 0"%h_syst.GetName(),h_syst.GetBinContent(ix))
        h_new.SetBinError(ix, 0)
    return h_new
        
def convert_error_bars_to_error_ratio_hist(h, direction=1):
    """Create hist with bin content = (bin value ± bin error) / bin value, 0 error"""
    h_new = h.Clone(h.GetName() + get_unique_str())
    for ix in range(1, h_new.GetNbinsX()+1):
        if h.GetBinContent(ix) > 0:
            h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
        else:
            if h.GetBinContent(ix) < 0:
                h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                print("_convert_error_bars_to_error_ratio_hist() warning: bin %d content < 0!" % (ix))
            else:
                h_new.SetBinContent(ix, 1)
        h_new.SetBinError(ix, 0)
    return h_new


def get_th1_bin_centers(h):
    # TODO maintain same shape(1, n) as in th1_to_array?
    centers = np.array([h.GetBinLowEdge(i) + 0.5*h.GetBinWidth(i) for i in range(1, h.GetNbinsX()+1)])
    return centers

def scale_th2_bin_widths(h2d, bins):
    """Scale bins of a square TH2 by bin widths
    bins is a list of bin edges, must have 1 more value than the number of bins in h2d
    """
    if len(bins) != h2d.GetNbinsX()+1:
        print(bins)
        print(h2d.GetNbinsX())
        raise ValueError("Wrong number of bins to scale x axis")
    if len(bins) != h2d.GetNbinsY()+1:
        raise ValueError("Wrong number of bins to scale y axis")
    for ix, (binx_low, binx_high) in enumerate(zip(bins[:-1], bins[1:]), 1):
        for iy, (biny_low, biny_high) in enumerate(zip(bins[:-1], bins[1:]), 1):
            width_x = binx_high - binx_low
            width_y = biny_high - biny_low
            scale = width_x * width_y
            value = h2d.GetBinContent(ix, iy)
            err = h2d.GetBinError(ix, iy)
            h2d.SetBinContent(ix, iy, value / scale)
            h2d.SetBinError(ix, iy, err / scale)
                
def get_th1_bin_widths(h):
    # TODO maintain same shape(1, n) as in th1_to_array?
    widths = np.array([h.GetBinWidth(i) for i in range(1, h.GetNbinsX()+1)])
    return widths
def get_unique_str():
    return ROOT.TUUID().AsString()
def ndarray_to_th1(nd_array, has_oflow_x=False, offset=0.5, bins=None):
    """Convert numpy ndarray row vector to TH1, with shape (1, nbins)
    Use has_oflow_x to include the under/overflow bins
    """
    nbinsx = nd_array.shape[1]
    print(len(bins),nbinsx)
    nbins_hist = nbinsx
    if has_oflow_x:
        nbins_hist -= 2

    # need the 0.5 offset to match TUnfold
    if bins == None:
        h = ROOT.TH1F(get_unique_str(), "", nbins_hist, offset, nbins_hist+offset)
    else:
        if len(bins) != nbinsx+1:
            raise IndexError("len(bins) != nbinsx + 1")
        h = ROOT.TH1F(get_unique_str(), "", nbins_hist, bins)

    x_start = 1
    x_end = nbins_hist

    if has_oflow_x:
        x_start = 0
        x_end = nbins_hist+1

    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        h.SetBinContent(ix, nd_array[0][x_ind])
        h.SetBinError(ix, np.sqrt(abs(nd_array[0][x_ind])))
        #FIXME how to do errors
    return h

def th1_to_ndarray(hist_A, oflow_x=False, uflow_x=False):
    """Convert TH1 to numpy ndarray"""
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 1
    if uflow_x: 
        ncol += 1

    # makes column vectors
    contents = np.zeros(shape=(1, ncol), dtype=np.float64)
    errors = np.zeros(shape=(1, ncol), dtype=np.float64)

    # Get ROOT indices to loop over
    x_start = 0 if uflow_x else 1
    x_end = hist_A.GetNbinsX() +(1 if oflow_x else 0)
    #if oflow_x:
    #    x_end += 1

    # x_ind for numpy as always starts at 0
    # ix for ROOT
    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        #print("Convert TH1 to numpy ndarray", x_ind,ix)
        contents[0][x_ind] = hist_A.GetBinContent(ix)
        errors[0][x_ind] = hist_A.GetBinError(ix)

    return contents, errors

def th2_to_tmatrixd(hist, include_uflow=False, include_oflow=False):
    n_rows = hist.GetNbinsY()
    n_cols = hist.GetNbinsX()

    # ignore for now as too complicated
    # if include_uflow:
    #     n_rows += 1
    #     n_cols += 1
    # if include_oflow:
    #     n_rows += 1
    #     n_cols += 1

    # taken from https://root.cern.ch/doc/master/TH2_8cxx_source.html#l03739
    m = ROOT.TMatrixD(n_rows, n_cols)
    ilow = m.GetRowLwb()
    iup  = m.GetRowUpb()
    jlow = m.GetColLwb()
    jup  = m.GetColUpb()
    for i in range(ilow, iup+1):
        for j in range(jlow, jup+1):
            m[i,j] = hist.GetBinContent(j-jlow+1,i-ilow+1)
    return m

def th2_to_ndarray(hist_A, oflow_x=False, oflow_y=False, uflow_x=False, uflow_y=False):
    """Convert TH2 to numpy ndarray"""
    
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 1
    if uflow_x: 
        ncol += 1
    
    nrow = hist_A.GetNbinsY()
    if oflow_y:
        nrow += 1
    if uflow_y: 
        nrow += 1
        

    contents = np.zeros(shape=(nrow, ncol), dtype=np.float64)
    errors = np.zeros(shape=(nrow, ncol), dtype=np.float64)
    # access via contents[irow][icol]

    # Get ROOT indices to loop over
    y_start = 0 if uflow_y else 1
    y_end = hist_A.GetNbinsY() +(1 if oflow_y else 0)
    #if oflow_y:
    #    y_end += 1

    x_start = 0 if uflow_x else 1
    x_end = hist_A.GetNbinsX() +(1 if oflow_x else 0)
    #if oflow_x:
    #    x_end += 1

    # y_ind, x_ind for numpy as always starts at 0
    # iy, ix for ROOT
    for y_ind, iy in enumerate(range(y_start, y_end+1)):
        for x_ind, ix in enumerate(range(x_start, x_end+1)):
            contents[y_ind][x_ind] = hist_A.GetBinContent(ix, iy)
            errors[y_ind][x_ind] = hist_A.GetBinError(ix, iy)

    return contents, errors


def ndarray_to_th2(data, offset=0, binsx=None, binsy=None):
    nbinsy, nbinsx = data.shape
    bins_x = array('d', [x+offset for x in range(1, nbinsx+2)]) # e.g. offset = -0.5 for TUnfold
    bins_y = array('d', [x+offset for x in range(1, nbinsy+2)])
    if binsx is not None:
        if len(binsx) != nbinsx+1:
            raise IndexError("binsx wrong size")
        bins_x = binsx
    if binsy is not None:
        if len(binsy) != nbinsy+1:
            raise IndexError("binsy wrong size")
        bins_y = binsy
    h = ROOT.TH2D(get_unique_str(), "", nbinsx, bins_x, nbinsy, bins_y)
    for ix in range(nbinsx):
        for iy in range(nbinsy):
            h.SetBinContent(ix+1, iy+1, data[iy,ix])
            h.SetBinError(ix+1, iy+1, 0)
    return h

def make_hist_from_diagonal_errors(h2d, bins=None, do_sqrt=True, set_errors=True, offset=0.):
    """Make 1D hist, with errors or contents set to diagonal elements from h2d
    Can be TH2 or numpy.ndarray, cos we have to use both
    Yes that is majorly wack
    set_errors: True to set TH1 bin errors, otherwise sets bin contents
    offset is on bin edge from 0 (TUnfold is 0.5)
    """
    if isinstance(h2d, ROOT.TH2):
        nbins = h2d.GetNbinsX()
        hnew = ROOT.TH1D("h_diag" + get_unique_str(), "", len(bins)-1, bins)#nbins, offset, nbins+offset)
        for i in range(1, nbins+1):
            err = h2d.GetBinContent(i, i)
            if do_sqrt and err > 0:
                err = np.sqrt(err)
            if set_errors:
                hnew.SetBinContent(i, 0)
                hnew.SetBinError(i, err)
            else:
                hnew.SetBinContent(i, err)
                hnew.SetBinError(i, 0)
        return hnew
    elif isinstance(h2d, np.ndarray):
        nbins = h2d.shape[0]
        hnew = ROOT.TH1D("h_diag" + get_unique_str(), "", len(bins)-1, bins)#nbins, offset, nbins+offset)
        for i in range(1, nbins+1):
            err = h2d[i-1, i-1]
            if do_sqrt and err > 0:
                err = np.sqrt(err)
            if set_errors:
                hnew.SetBinContent(i, 0)
                hnew.SetBinError(i, err)
            else:
                hnew.SetBinContent(i, err)
                hnew.SetBinError(i, 0)
        return hnew

def update_hist_bin_error(h_orig, h_to_be_updated):
    """Change the errors in h_to_be_updated to those from h_orig"""
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


################# courtesy of Robin's code, not used currently #################
def get_unsmooth_bins(arr2d, axis='gen'):
    """Find bins that are not smooth, ie they are spikey"""
    # Get 1D hist, by summing over bins of other axis in 2D hist
    bins = get_1D_bins_from_arr2d(arr2d, axis)
    diffs = np.diff(bins)
    diff_signs = np.sign(np.diff(bins))
    bad_diff_bins = []
    for i in range(len(diff_signs)-2):
        if i == 0:
            # edge case: consider i, i+1
            # but avoid cases where binning is coarse - could be genuine
            # 1st bin and not a spike i.e. +ve gradient up to peak
            if diff_signs[0] != diff_signs[1] and diffs[0] < 0:
                print('>>>>! Bin', i, 'diff wrong sign:', diffs[0], diffs[1])
                bad_diff_bins.append(i)
        else:
            # consider i-1, i, i+1
            if (diff_signs[i] != diff_signs[i-1] and
                diff_signs[i] != diff_signs[i+1]):
                print('>>>>! Bin', i, 'diff wrong sign:', diffs[i-1], diffs[i], diffs[i+1])
                bad_diff_bins.append(i)
    return bad_diff_bins

def get_1D_bins_from_arr2d(arr2d, axis='gen'):
    return arr2d.sum(axis=0 if axis.lower() == 'gen' else 1)




from collections import OrderedDict
import ROOT
from ROOT import * 
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
def scale_2d_hist(data, xbins, ybins):
    """
    Scales a 2D histogram provided as either a NumPy array or a ROOT.TMatrixD by dividing
    each element by its corresponding bin area (binWidth_x * binWidth_y).

    Parameters:
      data (numpy.ndarray or ROOT.TMatrixD): 2D array or TMatrix containing histogram contents.
      xbins (list or array): Bin edges for the x-axis. Length should be (n_xbins+1).
      ybins (list or array): Bin edges for the y-axis. Length should be (n_ybins+1).

    Returns:
      Scaled data of the same type as the input. If the input is a NumPy array, a new scaled
      NumPy array is returned; if the input is a ROOT.TMatrixD, a new scaled TMatrixD is returned.
    """
    # Compute bin widths for each axis
    x_widths = [xbins[i+1] - xbins[i] for i in range(len(xbins)-1)]
    y_widths = [ybins[j+1] - ybins[j] for j in range(len(ybins)-1)]
    
    # Check if the input is a NumPy array.
    if isinstance(data, np.ndarray):
        scaled = data.copy()
        nrows, ncols = scaled.shape
        for i in range(nrows):
            for j in range(ncols):
                area = x_widths[i] * y_widths[j]
                # Protect against division by zero.
                if area != 0:
                    scaled[i, j] /= area
        return scaled

    # Check if the input is a ROOT.TMatrixD.
    elif isinstance(data, ROOT.TMatrixD):
        nrows = data.GetNrows()
        ncols = data.GetNcols()
        # Copy TMatrixD data into a NumPy array.
        arr = np.empty((nrows, ncols), dtype=float)
        for i in range(nrows):
            for j in range(ncols):
                arr[i, j] = data(i, j)
        # Scale the array by the bin areas.
        for i in range(nrows):
            for j in range(ncols):
                area = x_widths[i] * y_widths[j]
                if area != 0:
                    arr[i, j] /= area
        # Create a new TMatrixD from the flattened, scaled array.
        flat = arr.flatten()
        scaled_matrix = ROOT.TMatrixD(nrows, ncols, flat)
        return scaled_matrix

    else:
        raise TypeError("Input data type not supported. Expected numpy.ndarray or ROOT.TMatrixD.")


def scale_th2_by_bin_width(histIn):
    """
    Scales contents and errors of a TH2 histogram by dividing each bin by the bin's area,
    defined as (binWidth_x * binWidth_y). 
    """
    hist = histIn.Clone(histIn.GetName()+'_divByBW')
    nBinsX = hist.GetNbinsX()
    nBinsY = hist.GetNbinsY()
    
    # Loop over all bins (1-indexed)
    for ix in range(1, nBinsX + 1):
        binWidthX = hist.GetXaxis().GetBinWidth(ix)
        for iy in range(1, nBinsY + 1):
            binWidthY = hist.GetYaxis().GetBinWidth(iy)
            area = binWidthX * binWidthY
            # Avoid division by zero (shouldn't happen, but it's good practice)
            if area == 0:
                print(f"WARNING: error in scale_th2_by_bin_width(), bin area for bin i,j={ix,iy} is ZERO; skipping normalising it by BW ")
                continue
            # Retrieve original content and error
            content = hist.GetBinContent(ix, iy)
            error = hist.GetBinError(ix, iy)
            # Scale content and error by the bin's area
            hist.SetBinContent(ix, iy, content / area)
            hist.SetBinError(ix, iy, error / area)
    return hist
def get_bin_edges(hist, axis='x'):
    if axis == 'x':
        axis = hist.GetXaxis()
    else:
        axis = hist.GetYaxis()
    
    bin_edges = [axis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)]  # +2 to include the upper edge of the last bin
    return bin_edges

def numpy_to_hist2D(arr, inp_hist):
    x_bin_edges = get_bin_edges(inp_hist, 'x')
    y_bin_edges = get_bin_edges(inp_hist, 'y')
    if arr.ndim == 3 and arr.shape[1] == 2:
        arr = arr[:, 0, :]
    #print(arr)
    
    # Create histogram with non-uniform binning
    hist = ROOT.TH2D(inp_hist.GetName() + "_copy", inp_hist.GetTitle(),
                len(x_bin_edges) - 1, np.array(x_bin_edges, 'd'),  # 'd' indicates double precision
                len(y_bin_edges) - 1, np.array(y_bin_edges, 'd'))
    
    # Set the bin content from the numpy array
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            hist.SetBinContent(i+1, j+1, arr[i, j])
            
    return hist

def hist2D_to_numpy(hist):
    arr = np.zeros((hist.GetNbinsX(), hist.GetNbinsY()))
    for i in range(hist.GetNbinsX()):
        for j in range(hist.GetNbinsY()):
            arr[i, j] = hist.GetBinContent(i+1, j+1)
    return arr
def compute_jacobian_withBW(hist, total):
    nbins = hist.GetNbinsX()
    jacobian = np.zeros((nbins, nbins))
    
    for i in range(nbins):
        for j in range(nbins):
            bin_width_i = hist.GetXaxis().GetBinWidth(i+1)
            if i == j:
                jacobian[i, j] = 1.0 / (total * bin_width_i)
            else:
                jacobian[i, j] = -hist.GetBinContent(i+1) / (total**2 * bin_width_i)
                
    return jacobian


def compute_jacobian_shape_plus_binwidth(hist, withOF=False):
    """
    Returns the Jacobian for shape normalization *and* dividing
    by bin widths.  That is: y -> (1/bin_width) * (y / Sum(y)).
    """
    # Step 1) shape-only Jacobian
    J = compute_jacobian(hist, withOF)
    
    # Step 2) diagonal for 1/bin_width
    nbins = hist.GetNbinsX()
    D = np.eye(nbins)
    for i in range(nbins):
        bw = hist.GetBinWidth(i+1)
        # Avoid zero or negative widths
        if bw <= 0:
            raise ValueError(f"Bin {i+1} has non-positive width={bw}")
        D[i, i] = 1./bw
    
    # Combined transform M = D * J
    M = D @ J
    return M

def compute_jacobian(hist,ndim=0,withOF=False):
    #if doing ndim:
    #since we want to normalise everything by the XS; 
    #so the total number of events (N) in histos used to compute 
    #the jacobian must be modified such that the terms (N-bin content(i)) 
    #are not incorrect
    
    N = hist.Integral(0, hist.GetNbinsX() + (0 if withOF==False else 1))
    
    if ndim>0:
        N/=ndim
    nbins = hist.GetNbinsX()
    #h = hist.Clone()
    jacobian = np.zeros((nbins, nbins))
        
    for i in range(nbins):
        for j in range(nbins):
            bc = hist.GetBinContent(i+1) if not hist.GetBinContent(i+1)<0 else 0.
            
            if i == j:
                jacobian[i, j] = (N-bc) / (N**2)
            else:
                jacobian[i, j] = - bc / (N**2 )
                
    return jacobian


def compute_jacobian_from_nparray(x):
    """
    Compute the Jacobian matrix for the normalization transformation:
    
        f_i(y) = y_i / sum(y)
    
    for an n-dimensional vector y. The Jacobian is given by:
    
        J_{ij} = δ_ij/s - y_i/s^2,
    
    where s = sum(y) and δ_ij is the Kronecker delta.
    
    Parameters
    ----------
    y : array_like
        A 1D array of n elements.
    
    Returns
    -------
    J : ndarray
        An (n x n) Jacobian matrix.
    
    """
    y = np.atleast_1d(y).astype(float)
    s = np.sum(y)
    if s == 0:
        raise ValueError("The sum of the elements in y is zero; cannot normalize.")
    
    # For diagonal elements: 1/s - y[i]/s^2, and for off-diagonals: - y[i]/s^2.
    J = np.eye(y.size) / s - np.outer(y, np.ones_like(y)) / (s**2)
    return J

def get_normalised_cov( unfHisto, ematrix, ndim=0):
    """
    Generic method to normalise any covariance matrix using a Jacobian
    """
    J = compute_jacobian(unfHisto.Clone(),ndim=ndim)
    
    if isinstance(ematrix, ROOT.TH2):
        cov_abs, _ = th2_to_ndarray(ematrix)
    else:
        cov_abs = ematrix
    
    cov_norm = J @ cov_abs @ J.T
    
    ematrix_norm_th2 = numpy_to_hist2D(cov_norm, ematrix.Clone(ematrix.GetName()+'_jacTrafo'))
    
    return cov_norm, ematrix_norm_th2

def get_normalised_cov_with_binwidth(hist, cov):
    """
    Like get_normalised_cov(), but also divides by bin width.
    So the final distribution is "1/N  * 1/bin_width".
    """
    # Convert TH2 -> numpy if needed
    if isinstance(cov, ROOT.TH2):
        cov_abs, _ = th2_to_ndarray(cov)
    else:
        cov_abs = cov

    # Build the combined Jacobian M = D@J
    M = compute_jacobian_shape_plus_binwidth(hist)

    # Transform V' = M V M^T
    cov_norm = M @ cov_abs @ M.T

    # Convert back to TH2 if desired
    ematrix_norm_th2 = numpy_to_hist2D(cov_norm, cov.Clone(cov.GetName()+'_jacTrafoBW'))

    return cov_norm, ematrix_norm_th2

def get_th1_normedCovErrors(hist, covnorm):
    for i in range( 1, hist.GetNbinsX() + 1):
        unc_tot_norm = ROOT.TMath.Sqrt( covnorm[i-1][i-1] )
        #hist.SetBinContent(i, unc_tot_norm )   
        #print(i,hist.GetBinContent(i),hist.GetBinError(i),unc_tot_norm)
        hist.SetBinError(i, unc_tot_norm )
        

def get_uncertainties_from_cov(cov_matrix):
    # Convert the covariance matrix to a numpy array
    if cov_matrix.ndim == 3 and cov_matrix.shape[1] == 2:
        cov_matrix = cov_matrix[:, 0, :]
    
    #cov_matrix_arr = hist_to_numpy(cov_matrix_hist)

    # Extract the variances from the diagonal and then compute the standard deviations
    uncertainties = np.sqrt(np.diag(cov_matrix))
    
    return uncertainties


def compute_relative_uncertainty(cov_source, cov_total):
    """
    Compute the relative uncertainty contribution per bin from a given source.
    
    Parameters
    ----------
    cov_source : ndarray
        The covariance matrix from a single uncertainty source (n x n).
    cov_total : ndarray
        The total covariance matrix (n x n), assumed to be the sum of all sources.
    
    Returns
    -------
    rel_unc : ndarray
        1D array of length n, where each element is the relative uncertainty 
        contribution in that bin, i.e., sqrt(var_source) / sqrt(var_total).
        If the total variance in a bin is zero, that bin is set to zero.
    """
    # Extract the diagonal elements (variances) for each bin.
    var_source = np.diag(cov_source)
    var_total = np.diag(cov_total)
    
    # Compute standard deviations (uncertainties)
    unc_source = np.sqrt(var_source)
    unc_total = np.sqrt(var_total)
    
    # Avoid division by zero. Where unc_total is zero, set the relative uncertainty to 0.
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_unc = np.divide(unc_source, unc_total, out=np.zeros_like(unc_source), where=(unc_total != 0))
    
    return rel_unc
    

def compute_relative_uncertainty_cov_TH2(cov_source, cov_total, name="relative_unc_cov", 
                                           title="Relative Uncertainty Covariance"):
    """
    Compute the per-bin (diagonal) relative uncertainty contribution from a specific source.
    
    For each bin i the relative uncertainty is defined as:
    
        r_i = sqrt( cov_source(i,i) ) / sqrt( cov_total(i,i) )
    
    Off-diagonal elements of the output TH2 are set to zero.
    
    Parameters
    ----------
    cov_source : ROOT.TH2
        TH2 covariance matrix from a single uncertainty source.
    cov_total : ROOT.TH2
        TH2 total covariance matrix (sum of all sources). Must have the same binning as cov_source.
    name : str, optional
        Name for the output TH2.
    title : str, optional
        Title for the output TH2.
    
    Returns
    -------
    rel_unc_array : np.array
        A numpy arraywhose elements contain the relative uncertainty contributions for a source.
    """
    # Check that cov_total is square
    nbinsX = cov_total.GetNbinsX()
    nbinsY = cov_total.GetNbinsY()
    if nbinsX != nbinsY:
        raise ValueError("cov_total is not square: number of X bins != number of Y bins")
    
    # Clone cov_total to preserve binning, axis labels, etc.
    rel_cov = cov_total.Clone(name)
    rel_cov.SetTitle(title)
    rel_cov.Reset()
    
    for i in range(1, nbinsX+1):
        var_source = cov_source.GetBinContent(i, i)
        var_total  = cov_total.GetBinContent(i, i)
        if var_total > 0:
            ratio = (var_source**0.5) / (var_total**0.5)
        else:
            ratio = 0.0
        rel_cov.SetBinContent(i, i, ratio)
        rel_cov.SetBinError(i, i, 0.0)
    
    return rel_unc_array




def machineEps():
    import sys
    deps = sys.float_info.epsilon
    return deps

########################################################################
###################### Histo modification helpers ######################
########################################################################

import ROOT

def subtract_background_covariance(fVyyData, bkgHistos, bkgScales, scales=None):
    """
    Subtracts background from the data covariance and adds both uncorrelated
    and correlated background uncertainties, exactly as in TUnfoldSys::DoBackgroundSubtraction.
    (https://root.cern.ch/doc/master/classTUnfoldSys.html#a5c1e23aca5ffe7729ce78e6de92a8467)
    The modified data input covariance provided as output is used for the BLT. 
    
    Parameters
    ----------
    covDet_wCorr : TH2
        Input data covariance (fVyyData), including any correlations.
    bkgHistos : dict[str, TH1]
        For each background source name, the TH1 histogram of the background shape.
    bkgScales : dict[str, float]
        For each background source name, the relative scale (/rate) uncertainty.
    scales : dict[str, float], optional
        For each background source name, the normalization factor applied to the
        background histogram before subtraction (default: 1.0 for all).
    Returns
    -------
    TH2
        A modified fVyyData provided as input to unfolding, but with background uncertainties added:
          vyy_ij = vyyData_ij
                  + sum_k [ (Δbkg_k(i))**2 * δᵢⱼ ]
                  + sum_k [ (δscale_k * bkg_k(i))·(δscale_k * bkg_k(j)) ]
        where Δbkg_k(i) = bkghist_k.GetBinError(i+1) and bkg_k(i)=bkghist_k.GetBinContent(i+1),
        δscale_k is the scale/rate unc. on the bkg_k, and δᵢⱼ is the Kronecker delta.
    """
    covDet_wCorr = fVyyData
    nbins = covDet_wCorr.GetNbinsX()
    if scales is None:
        scales = {name: 1.0 for name in bkgHistos}

    #copy input covariance into a TMatrixD; to keep things congruent with TUnfoldSys 
    #and determine bins that are nonzer (used bins in TUnfoldSys)
    vyy = ROOT.TMatrixD(nbins, nbins)
    usedBin = [False]*nbins
    for i in range(nbins):
        for j in range(nbins):
            val = covDet_wCorr.GetBinContent(i+1, j+1)
            vyy[i][j] = val
            if val > 0.0:
                usedBin[i] = True
                usedBin[j] = True

    #add uncorrelated errors:  (scale(=1.) * hist.GetBinError)^2  on the diagonal
    for name, hist in bkgHistos.items():
        scale = scales.get(name, 1.0)
        for i in range(nbins):
            if not usedBin[i]:
                continue
            err = scale * hist.GetBinError(i+1)
            vyy[i][i] += err*err

    #add correlated errors:
    #for each background k, add (δscale_k * bkg_k(i))*(δscale_k * bkg_k(j))
    for name, hist in bkgHistos.items():
        delta = bkgScales[name]
        for i in range(nbins):
            if not usedBin[i]:
                continue
            corr_i = delta * hist.GetBinContent(i+1)
            for j in range(nbins):
                if not usedBin[j]:
                    continue
                corr_j = delta * hist.GetBinContent(j+1)
                vyy[i][j] += corr_i * corr_j

    out = covDet_wCorr.Clone(covDet_wCorr.GetName()+"_withBkgErr")
    out.Reset()
    for i in range(nbins):
        for j in range(nbins):
            out.SetBinContent(i+1, j+1, vyy[i][j])

    return out



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
    
    newHist = ROOT.TH1D(hist.Rebin(len(bins)-1, "%s_rebin"%hist_name, bins))
    #newHist.SetStats(ROOT.kFALSE)
    if isMC==True: newHist.Scale(scale)
    newHist.SetDirectory(0)

    return newHist

def getRebinned_TH1(hist, xrange, scale, isMC=True):
    
    hist.SetDirectory(0)
    ROOT.TH1.StatOverflow(ROOT.kTRUE)
    ROOT.TH1.AddDirectory(ROOT.kFALSE);
    
    bins = array('d', np.array(Bin(xrange)))
    
    newHist = ROOT.TH1D(hist.Rebin(len(bins)-1, "%s_rebin"%hist_name, bins))
    newHist.SetDirectory(0)
    
    return newHist


def rebin_RM_withUF(h_resp,genBin,recoBin):
    #### fancy way to create variable binning TH2D
    tmpHisto = ROOT.TH2D( h_resp.GetName()+"_Rebin", h_resp.GetName()+"_Rebin", len(genBin)-1, array( 'd', genBin), len(recoBin)-1, array( 'd', recoBin) )
    tmpHisto.Sumw2()      
    ROOT.TH2.StatOverflows(ROOT.kTRUE)
    ROOT.TH1.StatOverflows(ROOT.kTRUE)

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
    
    #binnings, such that last bin edge extends beyond last bin of initial binning to capture OF contents in new last bin of histo
    #if isinstance(h_resp, ROOT.TH2):
    #hackey-hack-hack-hackeroo
    nbX, nbY = h_resp.GetNbinsX(), h_resp.GetNbinsY()
    ofX = sum(h_resp.GetBinContent(nbX+1, y) for y in range(1, nbY+1))
    ofY = sum(h_resp.GetBinContent(x, nbY+1) for x in range(1, nbX+1))
    ofXY = h_resp.GetBinContent(nbX+1, nbY+1)
    if ofX or ofY or ofXY:
        print("########################################")
        print(f"WARNING: — TH2 {h_resp.GetName()} has overflow beyond last desired bin edge in rebinning:")
        print(f"  x-overflow total = {ofX:.1f}, y-overflow total = {ofY:.1f}, corner = {ofXY:.1f}")
        print("########################################")
        # — begin overflow handling —

        #nbX = h_resp.GetNbinsX()
        #nbY = h_resp.GetNbinsY()
        #if ('dijet' in tmpHisto.GetName() and '1p5_2' in tmpHisto.GetName()) or ('21' in tmpHisto.GetName() or '32' in tmpHisto.GetName()):
        #    print(f"DEBUG: original {h_resp.GetName()} edge bins:")
        #    for by in range(0, nbY+2):
        #        for bx in range(0, nbX+2):
        #            # skip interior
        #            if 1 <= bx <= nbX and 1 <= by <= nbY:
        #                continue
        #            w = h_resp.GetBinContent(bx, by)
        #            if abs(w) > 1e-9:
        #                print(f"  orig bin ({bx:2d},{by:2d}) = {w:.6f}")
        # 1) x‐overflow for each “normal” y‐bin
        for by in range(1, nbY+1):
            w = h_resp.GetBinContent(nbX+1, by)
            if w == 0: continue
            e = h_resp.GetBinError  (nbX+1, by)
            cy = h_resp.GetYaxis().GetBinCenter(by)
            # find which new reco‐bin this y‐center lands in
            for iY in range(len(recoBin)-1):
                if recoBin[iY] < cy <= recoBin[iY+1]:
                    tgtY = iY+1
                    break
            tgtX = tmpHisto.GetNbinsX()          # last new x‐bin
            ib = tmpHisto.GetBin(tgtX, tgtY)
            tmpHisto.SetBinContent(ib, tmpHisto.GetBinContent(ib) + w)
            tmpHisto.SetBinError  (ib, np.hypot(tmpHisto.GetBinError(ib), e))

        # 2) y‐overflow for each “normal” x‐bin
        for bx in range(1, nbX+1):
            w = h_resp.GetBinContent(bx, nbY+1)
            if w == 0: continue
            e = h_resp.GetBinError  (bx, nbY+1)
            cx = h_resp.GetXaxis().GetBinCenter(bx)
            for iX in range(len(genBin)-1):
                if genBin[iX] < cx <= genBin[iX+1]:
                    tgtX = iX+1
                    break
            tgtY = tmpHisto.GetNbinsY()          # last new y‐bin
            ib = tmpHisto.GetBin(tgtX, tgtY)
            tmpHisto.SetBinContent(ib, tmpHisto.GetBinContent(ib) + w)
            tmpHisto.SetBinError  (ib, np.hypot(tmpHisto.GetBinError(ib), e))

        # 3) corner overflow (both x & y)
        w = h_resp.GetBinContent(nbX+1, nbY+1)
        if w:
            e     = h_resp.GetBinError(nbX+1, nbY+1)
            tgtX  = tmpHisto.GetNbinsX()
            tgtY  = tmpHisto.GetNbinsY()
            ib = tmpHisto.GetBin(tgtX, tgtY)
            tmpHisto.SetBinContent(ib, tmpHisto.GetBinContent(ib) + w)
            tmpHisto.SetBinError  (ib, np.hypot(tmpHisto.GetBinError(ib), e))

        # 4) x‐overflow in the y‐underflow row (bottom‐right “underflow‐overflow”)
        w = h_resp.GetBinContent(nbX+1, 0)
        if w:
            e     = h_resp.GetBinError(nbX+1, 0)
            tgtX  = tmpHisto.GetNbinsX()        # last new x‐bin
            tgtY  = 0                           # the underflow row in tmpHisto
            ib = tmpHisto.GetBin(tgtX, tgtY)
            tmpHisto.SetBinContent(ib, tmpHisto.GetBinContent(ib) + w)
            tmpHisto.SetBinError  (ib, np.hypot(tmpHisto.GetBinError(ib), e))
        # 5) y‐overflow in the x‐underflow column (underflow‐overflow corner)
        w = h_resp.GetBinContent(0, nbY+1)
        if w:
            e     = h_resp.GetBinError(0, nbY+1)
            tgtX  = 0                       # underflow column in tmpHisto
            tgtY  = tmpHisto.GetNbinsY()   # last new y‐bin
            ib    = tmpHisto.GetBin(tgtX, tgtY)
            tmpHisto.SetBinContent(ib, tmpHisto.GetBinContent(ib) + w)
            tmpHisto.SetBinError  (ib, np.hypot(tmpHisto.GetBinError(ib), e))

        # — end overflow handling —
        #tmpHisto.Sumw2()
        tmpHisto.SetDirectory(0)
        #h_resp = copy.deepcopy(tmpHisto.Clone())
        # — DIAGNOSTIC: print every non-zero edge bin of the rebinned tmpHisto —
        if ('dijet' in tmpHisto.GetName() and '1p5_2' in tmpHisto.GetName()) or ('21' in tmpHisto.GetName() or '32' in tmpHisto.GetName()):
            print(f"DEBUG: rebinned {tmpHisto.GetName()} edge bins:")
            for by in range(0, tmpHisto.GetNbinsY()+2):
                for bx in range(0, tmpHisto.GetNbinsX()+2):
                    if 1 <= bx <= tmpHisto.GetNbinsX() and 1 <= by <= tmpHisto.GetNbinsY():
                        continue
                    w = tmpHisto.GetBinContent(bx, by)
                    if abs(w) > 1e-9:
                        print(f"  rebin bin ({bx:2d},{by:2d}) = {w:.6f}")
    orig_total = h_resp.Integral(0, nbX+1, 0, nbY+1)
    reb_total  = tmpHisto.Integral(0, tmpHisto.GetNbinsX()+1,
                                   0, tmpHisto.GetNbinsY()+1)
    if not math.isclose(orig_total, reb_total,rel_tol=1e-6, abs_tol=1e-3):
        print(f"WARNING on L618 of histoHelpers TH2 rebinning for {h_resp.GetName()}: Total weight changed between orig and rebin: {orig_total} vs. {reb_total}!")
    
    return tmpHisto


def rebin_RM_withUFandOF(h_resp, genBin, recoBin):
    #### fancy way to create variable binning TH2D
    

    tmpHisto = ROOT.TH2D(h_resp.GetName()+"_Rebin",
                         h_resp.GetName()+"_Rebin",
                         len(genBin)-1, array('d', genBin),
                         len(recoBin)-1, array('d', recoBin))
    tmpHisto.Sumw2()      
    ROOT.TH2.StatOverflows(ROOT.kTRUE)
    ROOT.TH1.StatOverflows(ROOT.kTRUE)

    nBinsX = h_resp.GetNbinsX()
    nBinsY = h_resp.GetNbinsY()

    for biny in range(0, nBinsY+2):
        by = h_resp.GetYaxis().GetBinCenter(biny)
        for binx in range(0, nBinsX+2):
            bx = h_resp.GetXaxis().GetBinCenter(binx)
            # Central bins: neither under- nor overflow
            if not (binx == 0 or biny == 0 or binx == nBinsX+1 or biny == nBinsY+1):
                for iX in range(len(genBin)-1):
                    for iY in range(len(recoBin)-1):
                        if (bx <= genBin[iX+1] and bx > genBin[iX]) and (by <= recoBin[iY+1] and by > recoBin[iY]):
                            jbin = h_resp.GetBin(binx, biny)
                            newBin = tmpHisto.GetBin(iX+1, iY+1)
                            tmpHisto.SetBinContent(newBin,
                                tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                            tmpHisto.SetBinError(newBin,
                                np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
            else:
                # Special bins: one or both coordinates are in under-/overflow regions.
                jbin = h_resp.GetBin(binx, biny)
                # --- First, treat the explicitly defined corner cases ---
                # Bottom right corner: x overflow (binx==nBinsX+1) and y underflow (biny==0)
                if binx == nBinsX+1 and biny == 0:
                    newX = tmpHisto.GetNbinsX()+1
                    newY = 0
                    newBin = tmpHisto.GetBin(newX, newY)
                    tmpHisto.SetBinContent(newBin,
                        tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                    tmpHisto.SetBinError(newBin,
                        np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
                # (Optionally, you might add a branch for the top‐left corner here.)
                # --- Now treat the remaining cases ---
                # Underflow in at least one coordinate (but not already caught by a corner)
                elif (binx == 0 or biny == 0):
                    # Here we mimic your original underflow treatment for y.
                    # Loop over the new histogram’s X bins to find where bx belongs.
                    j = 0  # new Y index for underflow is 0
                    for i in range(tmpHisto.GetNbinsX()+2):
                        if (bx > tmpHisto.GetXaxis().GetBinLowEdge(i) and 
                            bx <= tmpHisto.GetXaxis().GetBinLowEdge(i+1)):
                            newBin = tmpHisto.GetBin(i, j)
                            tmpHisto.SetBinContent(newBin,
                                tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                            tmpHisto.SetBinError(newBin,
                                np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
                # X overflow only (with central y: 1 <= biny <= nBinsY)
                elif binx == nBinsX+1 and 1 <= biny <= nBinsY:
                    newX = tmpHisto.GetNbinsX()+1
                    # Find the proper new Y bin from the bin center by
                    # looping over the user-defined reco bins.
                    for iY in range(len(recoBin)-1):
                        if (by <= recoBin[iY+1] and by > recoBin[iY]):
                            newY = iY+1
                            newBin = tmpHisto.GetBin(newX, newY)
                            tmpHisto.SetBinContent(newBin,
                                tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                            tmpHisto.SetBinError(newBin,
                                np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
                            break
                # Y overflow only (with central x: 1 <= binx <= nBinsX)
                elif biny == nBinsY+1 and 1 <= binx <= nBinsX:
                    newY = tmpHisto.GetNbinsY()+1
                    for iX in range(len(genBin)-1):
                        if (bx <= genBin[iX+1] and bx > genBin[iX]):
                            newX = iX+1
                            newBin = tmpHisto.GetBin(newX, newY)
                            tmpHisto.SetBinContent(newBin,
                                tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                            tmpHisto.SetBinError(newBin,
                                np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
                            break
                # Top right corner: both x and y are in overflow
                elif binx == nBinsX+1 and biny == nBinsY+1:
                    newX = tmpHisto.GetNbinsX()+1
                    newY = tmpHisto.GetNbinsY()+1
                    newBin = tmpHisto.GetBin(newX, newY)
                    tmpHisto.SetBinContent(newBin,
                        tmpHisto.GetBinContent(newBin) + h_resp.GetBinContent(jbin))
                    tmpHisto.SetBinError(newBin,
                        np.sqrt((tmpHisto.GetBinError(newBin))**2 + (h_resp.GetBinError(jbin))**2))
    tmpHisto.SetDirectory(0)
    return tmpHisto


def rebin_RM_withUFandOF2(h_resp, genBin, recoBin):
    """
    Rebin a 2D response (TH2) histogram into new variable binning (genBin for X, recoBin for Y),
    properly handling all under/overflow bins (including corners).
    """
    

    tmpHisto = ROOT.TH2D(h_resp.GetName()+"_Rebin",
                         h_resp.GetName()+"_Rebin",
                         len(genBin)-1, array('d', genBin),
                         len(recoBin)-1, array('d', recoBin))
    tmpHisto.Sumw2()
    ROOT.TH2.StatOverflows(ROOT.kTRUE)
    ROOT.TH1.StatOverflows(ROOT.kTRUE)

    nBinsX = h_resp.GetNbinsX()
    nBinsY = h_resp.GetNbinsY()

    for biny in range(0, nBinsY+2):
        by = h_resp.GetYaxis().GetBinCenter(biny)
        for binx in range(0, nBinsX+2):
            bx = h_resp.GetXaxis().GetBinCenter(binx)

            # The original bin's content & error
            jbin    = h_resp.GetBin(binx, biny)
            content = h_resp.GetBinContent(jbin)
            err     = h_resp.GetBinError(jbin)
            if content == 0 and err == 0:
                continue  # optional optimization

            #   "Normal" interior bin => find variable bin
            if (1 <= binx <= nBinsX) and (1 <= biny <= nBinsY):
                for iX in range(len(genBin) - 1):
                    if bx >= genBin[iX] and bx < genBin[iX+1]:
                        for iY in range(len(recoBin) - 1):
                            if by >= recoBin[iY] and by < recoBin[iY+1]:
                                newBin = tmpHisto.GetBin(iX+1, iY+1)
                                tmpHisto.SetBinContent(newBin,
                                    tmpHisto.GetBinContent(newBin) + content)
                                tmpHisto.SetBinError(newBin,
                                    np.sqrt(tmpHisto.GetBinError(newBin)**2 + err**2))
                continue

            # ---- Otherwise, at least one coordinate is under-/overflow. ----
            # Handle the four corners first
            if binx == 0 and biny == 0:
                # Bottom-left corner
                newX = 0
                newY = 0
            elif binx == nBinsX+1 and biny == 0:
                # Bottom-right corner
                newX = tmpHisto.GetNbinsX() + 1
                newY = 0
            elif binx == 0 and biny == nBinsY+1:
                # Top-left corner
                newX = 0
                newY = tmpHisto.GetNbinsY() + 1
            elif binx == nBinsX+1 and biny == nBinsY+1:
                # Top-right corner
                newX = tmpHisto.GetNbinsX() + 1
                newY = tmpHisto.GetNbinsY() + 1

            # Next handle pure X‐underflow/overflow with central Y
            elif binx == 0 and (1 <= biny <= nBinsY):
                # x underflow, central y
                newX = 0
                # find correct newY from by
                newY = None
                for iY in range(len(recoBin) - 1):
                    if by >= recoBin[iY] and by < recoBin[iY+1]:
                        newY = iY+1
                        break
                if newY is None:
                    # in principle could happen if by is below or above the user bin range
                    # might want to put it in underflow/overflow too
                    newY = 0 if by <= recoBin[0] else tmpHisto.GetNbinsY()+1

            elif binx == nBinsX+1 and (1 <= biny <= nBinsY):
                # x overflow, central y
                newX = tmpHisto.GetNbinsX() + 1
                newY = None
                for iY in range(len(recoBin) - 1):
                    if by >= recoBin[iY] and by < recoBin[iY+1]:
                        newY = iY+1
                        break
                if newY is None:
                    newY = 0 if by <= recoBin[0] else tmpHisto.GetNbinsY()+1

            # Now handle pure Y‐underflow/overflow with central X
            elif biny == 0 and (1 <= binx <= nBinsX):
                # y underflow, central x
                newY = 0
                newX = None
                for iX in range(len(genBin) - 1):
                    if bx >= genBin[iX] and bx < genBin[iX+1]:
                        newX = iX+1
                        break
                if newX is None:
                    newX = 0 if bx <= genBin[0] else tmpHisto.GetNbinsX()+1

            elif biny == nBinsY+1 and (1 <= binx <= nBinsX):
                # y overflow, central x
                newY = tmpHisto.GetNbinsY() + 1
                newX = None
                for iX in range(len(genBin) - 1):
                    if bx >= genBin[iX] and bx < genBin[iX+1]:
                        newX = iX+1
                        break
                if newX is None:
                    newX = 0 if bx <= genBin[0] else tmpHisto.GetNbinsX()+1

            else:
                # If we get here, it means e.g. x=0,y=0 is already handled above, etc.
                # but if you prefer you could put a "continue" or raise an error.
                continue

            # --- Add the content to the determined bin ---
            theBin = tmpHisto.GetBin(newX, newY)
            tmpHisto.SetBinContent(theBin, tmpHisto.GetBinContent(theBin) + content)
            tmpHisto.SetBinError(theBin, np.sqrt(tmpHisto.GetBinError(theBin)**2 + err**2))

    tmpHisto.SetDirectory(0)
    old_sum = h_resp.Integral(0, -1, 0, -1)    # integral in x from bin0..-1, y from bin0..-1 (includes under/overflow)
    new_sum = tmpHisto.Integral(0, -1, 0, -1)  # likewise for the new histogram
    #print("Old sum of bin contents:", old_sum)
    #print("New sum of bin contents:", new_sum)
    
    return tmpHisto

########################################
# For TH1 histograms
########################################

def print_under_over_flow_th1(h):
    """
    Print the underflow (bin 0) and overflow (bin nBins+1)
    contents (and errors) for a TH1 histogram.
    """
    nb = h.GetNbinsX()
    print("TH1 histogram:", h.GetName())
    print("  Underflow (bin 0): content =", h.GetBinContent(0),
          "error =", h.GetBinError(0))
    
    print("  Overflow (bin %d): content = %f, error = %f" %
          (nb+1, h.GetBinContent(nb+1), h.GetBinError(nb+1)))

def draw_under_over_flow_th1(h):
    """
    Draw the underflow and overflow bins of a TH1 as two separate 1D histograms.
    (One pad shows the underflow bin (bin 0) and the other the overflow (bin nBins+1).)
    """
    c = ROOT.TCanvas("c1", "Underflow/Overflow of "+h.GetName(), 800, 400)
    c.Divide(2, 1)
    
    # Underflow histogram:
    h_uf = ROOT.TH1D(h.GetName()+"_uf", "Underflow bin", 1, 0, 1)
    h_uf.SetBinContent(1, h.GetBinContent(0))
    h_uf.SetBinError(1, h.GetBinError(0))
    
    c.cd(1)
    h_uf.Draw("E")
    
    # Overflow histogram:
    nb = h.GetNbinsX()
    h_of = ROOT.TH1D(h.GetName()+"_of", "Overflow bin", 1, 0, 1)
    h_of.SetBinContent(1, h.GetBinContent(nb+1))
    h_of.SetBinError(1, h.GetBinError(nb+1))
    
    c.cd(2)
    h_of.Draw("E")
    
    c.Update()
    return c

########################################
# For TH2 histograms
########################################

def print_under_over_flow_th2(h):
    """
    Print the underflow/overflow for a TH2 histogram.
    (For each “special” bin, we print the content and error.)
    """
    nbx = h.GetNbinsX()
    nby = h.GetNbinsY()
    print("TH2 histogram:", h.GetName())
    
    print("\nX-axis underflow (x bin 0) for all y bins:")
    for y in range(0, nby+2):
        print("  y bin %2d: content = %f, error = %f" %
              (y, h.GetBinContent(0, y), h.GetBinError(0, y)))
        
    print("\nX-axis overflow (x bin %d) for all y bins:" % (nbx+1))
    for y in range(0, nby+2):
        print("  y bin %2d: content = %f, error = %f" %
              (y, h.GetBinContent(nbx+1, y), h.GetBinError(nbx+1, y)))
        
    print("\nY-axis underflow (y bin 0) for all x bins:")
    for x in range(0, nbx+2):
        print("  x bin %2d: content = %f, error = %f" %
              (x, h.GetBinContent(x, 0), h.GetBinError(x, 0)))
        
    print("\nY-axis overflow (y bin %d) for all x bins:" % (nby+1))
    for x in range(0, nbx+2):
        print("  x bin %2d: content = %f, error = %f" %
              (x, h.GetBinContent(x, nby+1), h.GetBinError(x, nby+1)))
        
def draw_under_over_flow_th2(h):
    """
    Draw the special TH2 bins as 1D histograms.
    Four separate pads show:
       - x underflow (x = 0) vs. y,
       - x overflow (x = nbx+1) vs. y,
       - y underflow (y = 0) vs. x,
       - y overflow (y = nby+1) vs. x.
    """
    nbx = h.GetNbinsX()
    nby = h.GetNbinsY()
    
    c = ROOT.TCanvas("c2", "TH2 Under/Overflows "+h.GetName(), 1200, 800)
    c.Divide(2, 2)
    
    # X underflow (x=0) vs y:
    h_xuf = ROOT.TH1D(h.GetName()+"_xuf", "X underflow (x=0)", nby+2, 0, nby+2)
    for y in range(0, nby+2):
        h_xuf.SetBinContent(y+1, h.GetBinContent(0, y))
        h_xuf.SetBinError(y+1, h.GetBinError(0, y))
    c.cd(1)
    h_xuf.Draw("E")
    
    # X overflow (x = nbx+1) vs y:
    h_xof = ROOT.TH1D(h.GetName()+"_xof", "X overflow (x=nbx+1)", nby+2, 0, nby+2)
    for y in range(0, nby+2):
        h_xof.SetBinContent(y+1, h.GetBinContent(nbx+1, y))
        h_xof.SetBinError(y+1, h.GetBinError(nbx+1, y))
    c.cd(2)
    h_xof.Draw("E")
    
    # Y underflow (y=0) vs x:
    h_yuf = ROOT.TH1D(h.GetName()+"_yuf", "Y underflow (y=0)", nbx+2, 0, nbx+2)
    for x in range(0, nbx+2):
        h_yuf.SetBinContent(x+1, h.GetBinContent(x, 0))
        h_yuf.SetBinError(x+1, h.GetBinError(x, 0))
    c.cd(3)
    h_yuf.Draw("E")
    
    # Y overflow (y = nby+1) vs x:
    h_yof = ROOT.TH1D(h.GetName()+"_yof", "Y overflow (y=nby+1)", nbx+2, 0, nbx+2)
    for x in range(0, nbx+2):
        h_yof.SetBinContent(x+1, h.GetBinContent(x, nby+1))
        h_yof.SetBinError(x+1, h.GetBinError(x, nby+1))
    c.cd(4)
    h_yof.Draw("E")
    
    c.Update()
    return c




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

def get_bin_centers(histogram):
    """
    Extracts the bin centers from a ROOT.TH1 histogram.

    Parameters:
    histogram (ROOT.TH1): A ROOT histogram (e.g., TH1F, TH1D).

    Returns:
    numpy.ndarray: An array of bin centers.
    """
    nbins = histogram.GetNbinsX()
    bin_centers = np.zeros(nbins)
    
    for bin_index in range(1, nbins + 1):  # ROOT bins are 1-indexed
        bin_centers[bin_index - 1] = histogram.GetBinCenter(bin_index)
    
    return bin_centers

def fill_histogram_from_array(array, histogram):
    """
    Fills a ROOT.TH1D histogram with the contents of a given array.

    Parameters:
    array (numpy.ndarray): An array of bin contents.
    histogram (ROOT.TH1D): An empty ROOT.TH1D histogram to be filled.

    Returns:
    None
    """
    if len(array) != histogram.GetNbinsX():
        raise ValueError("The number of bins in the array does not match the histogram bins.")

    for bin_index, value in enumerate(array, start=1):  # ROOT bins are 1-indexed
        histogram.SetBinContent(bin_index, value)
        
    return histogram

def fill_histogramErrors_from_array(array, histogram):
    """
    Fills a ROOT.TH1D histogram with the contents of a given array.

    Parameters:
    array (numpy.ndarray): An array of bin contents.
    histogram (ROOT.TH1D): An empty ROOT.TH1D histogram to be filled.

    Returns:
    None
    """
    if len(array) != histogram.GetNbinsX():
        raise ValueError("The number of bins in the array does not match the histogram bins.")

    for bin_index, value in enumerate(array, start=1):  # ROOT bins are 1-indexed
        histogram.SetBinError(bin_index, value)
    return histogram


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

def fill_histogram_from_array(array, histogram):
    """
    Fills a ROOT.TH1D histogram with the contents of a given array.

    Parameters:
    array (numpy.ndarray): An array of bin contents.
    histogram (ROOT.TH1D): An empty ROOT.TH1D histogram to be filled.

    Returns:
    None
    """
    if len(array) != histogram.GetNbinsX():
        raise ValueError("The number of bins in the array does not match the histogram bins.")

    for bin_index, value in enumerate(array, start=1):  # ROOT bins are 1-indexed
        histogram.SetBinContent(bin_index, value)

def normalise_hist_divide_bin_width(h):
    h_new = h.Clone(h.GetName()+"DivideBinWidth")
    normalise_hist(h_new)
    h_new = hist_divide_bin_width(h_new)
    return h_new

def get_syst_shifted_hist(syst_shift, unfolded=None):
    """Get histogram with systematic shift applied to bin contents
    Can specify starting hist, otherwise assumes unfolded w/no error
    """
    hist_shift = syst_shift.Clone(syst_shift.GetTitle()+'_shiftedbyNom')
    hist_shift.Add(unfolded)  
    return hist_shift


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
            #h_new.SetBinContent(ix, 1)
            print("Warning: convert_syst_shift_to_error_ratio_hist bin", ix, "of h_syst, %s, < 0"%h_syst.GetName(),h_syst.GetBinContent(ix))
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

def convert_shift_hist_to_covariance(hist):
    """Convert 1D shift histogram to 2D covariance matrix, using V = x x^T"""
    xax = hist.GetXaxis()
    bins = array('d', [xax.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX()+2)])
    nbins = len(bins) - 1
    h2d = ROOT.TH2D("covariance_" + hist.GetName(), "Covariance;%s;%s" % (xax.GetTitle(), xax.GetTitle()), nbins, 0, nbins, nbins, 0, nbins)
    values = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)])
    values = values.reshape(len(values), 1)  # turn into column vector
    cov_values = values.dot(values.T)
    for ix in range(nbins):
        for iy in range(nbins):
            h2d.SetBinContent(ix+1, iy+1, cov_values[ix][iy])
            h2d.SetBinError(ix+1, iy+1, 0)
    return h2d

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
def ndarray_to_th1(nd_array, has_oflow_x=False, offset=0., bins=None):
    """Convert numpy ndarray row vector to TH1, with shape (1, nbins)
    Use has_oflow_x to include the under/overflow bins
    """
    nbinsx = nd_array.shape[1]
    #print(len(bins),nbinsx)
    nbins_hist = nbinsx
    if has_oflow_x:
        nbins_hist -= 2

    # need the 0.5 offset to match TUnfold
    if bins == None:
        h = ROOT.TH1D(get_unique_str(), "", nbins_hist, offset, nbins_hist+offset)
    else:
        if len(bins) != nbinsx+1:
            raise IndexError("len(bins) != nbinsx + 1")
        h = ROOT.TH1D(get_unique_str(), "", nbins_hist, bins)

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

def make_diag_cov_hist_from_errors(h1d, do_squaring=True, inverse=False):
    """
    Taken from Robin: https://github.com/raggleton/QGAnalysisPlotting/blob/26bb66e690a4a052b9b1acc328059a372fd25c6b/my_unfolder.py#L1875
    Make diagonal TH2 from errors on TH1.

    Assumes off-diag = 0.

    Can also do inverse by doing 1/(err^2)
    """
    nbins = h1d.GetNbinsX()
    bin_edges = array('d', [h1d.GetBinLowEdge(i) for i in range(1, nbins+2)])
    #print("make_diag_cov_hist_from_errors", bin_edges)
    h = ROOT.TH2D(h1d.GetName()+get_unique_str(), "", nbins, bin_edges, nbins, bin_edges)
    for i in range(1, nbins+1):
        err = h1d.GetBinError(i)
        if do_squaring:
            err *= err
        if inverse:
            if err != 0:
                err = 1/err
            else:
                err = 0
        h.SetBinContent(i, i, err)
    return h
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



def create_presel_string(trig_list):
    conditions = ["("]#
    for trig in trig_list:
        conditions.append(f"(recoEvents[\'selRecoMask_{trig}\'])")
    return (" | ".join(conditions)).replace('( | ', '(')+')'

def create_dict_evtIndices_perTrig(events,year,trig_list):
    #Used when building event matrices for cross-observable correlations and input, combined unfolding covariances in dijet data
    weight_temp = np.ones(len(events))*events['totalRecoWeight_nom'].to_numpy()
    for trig in trig_list:
        prescale = checkDict( 'JetHT', t3_samples_dijets )[year]['triggerList'][trig]

        print(trig,year,prescale)
        
        weight_temp[events[f'selRecoMask_{trig}']] *= prescale #* np.ones(len(events['trigWeights'][events[f'selRecoMask_{trig}']==True]))
    events['trigWeights'] = weight_temp #*events['totalRecoWeight_nom'] #np.ones(len(events))
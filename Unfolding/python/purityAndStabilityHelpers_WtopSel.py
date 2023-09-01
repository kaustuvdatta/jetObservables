import sys, os, gc, glob
import array as array
from array import array
from collections import OrderedDict
import numpy as np
from histoHelpers import *
import copy 

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


from random import gauss
from loadAndRunTUnfold_newRM import *
sys.path.insert(0,'../test/')
#from DrawHistogram_dijetSel import *
sys.path.insert(0,'../python/')
import CMS_lumi as CMS_lumi
#import CMSStyle as new_tdrstyle
import tdrstyle as tdrstyle
#from variables_dijetSel_newSel import nSubVariables as variables



#########################################################################
######################## Dictionary of variables ########################
#########################################################################

nSubVariables = OrderedDict()

# Dict keys correspond to string input required to correctly load histos
# 'bins' (not renamed due to changes that'd be required elsewhere) contains min/max values of the observable
# the max prevent the optimizer from going too crazy-wide with the tail bins' widths 

nSubVariables[ 'Jet'+'_tau_0p25_1' ] = {
    'bins' :   array('d',[0., 1.]),
    'label' : ' AK8 jet #tau_{1}^{(0.25)}',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau_0p25_2' ] = {
    'bins' :  array('d',[0., 0.9]),
    'label' : ' AK8 jet #tau_{2}^{(0.25)}',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau_0p25_3' ] = {
    'bins' :  array('d',[0., 0.8]),
    'label' : ' AK8 jet #tau_{3}^{(0.25)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_0p25_4' ] = {
    'bins' :  array('d',[0.0, 0.7]),
    'label' : ' AK8 jet #tau_{4}^{(0.25)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_0p25_5' ] = {
    'bins' :  array('d',[0.0, 0.7]),
    'label' : ' AK8 jet #tau_{5}^{(0.25)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_0p5_1' ] = {
    'bins' :   array('d',[0., 0.9]),
    'label' : ' AK8 jet #tau_{1}^{(0.5)}',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau_0p5_2' ] = {
    'bins' :  array('d',[0., 0.8]),
    'label' : ' AK8 jet #tau_{2}^{(0.5)}',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau_0p5_3' ] = {
    'bins' :  array('d',[0., 0.7]),
    'label' : ' AK8 jet #tau_{3}^{(0.5)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_0p5_4' ] = {
    'bins' :  array('d',[0.0, 0.6]),
    'label' : ' AK8 jet #tau_{4}^{(0.5)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_0p5_5' ] = {
    'bins' :  array('d',[0.0,0.5]),
    'label' : ' AK8 jet #tau_{5}^{(0.5)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_1_1' ] = {
    'bins' :   array('d',[0., 0.7]),
    'label' : ' AK8 jet #tau_{1}^{(1)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1_2' ] = {
    'bins' : array('d',[0., 0.6]),
    'label' : ' AK8 jet #tau_{2}^{(1)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1_3' ] = {
    'bins' :  array('d',[0., .6]),   
    'label' : ' AK8 jet #tau_{3}^{(1)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1_4' ] = {
    'bins' : array('d',[0., 0.5]),
    'label' : ' AK8 jet #tau_{4}^{(1)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1_5' ] = {
    'bins' :  array('d',[0., 0.4]),
    'label' : ' AK8 jet #tau_{5}^{(1)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_1p5_1' ] = {
    'bins' :   array('d',[0., 0.7]),
    'label' : ' AK8 jet #tau_{1}^{(1.5)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_1p5_2' ] = {
    'bins' :  array('d',[0., 0.6]),
    'label' : ' AK8 jet #tau_{2}^{(1.5)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1p5_3' ] = {
    'bins' :  array('d',[0., 0.6]),
    'label' : ' AK8 jet #tau_{3}^{(1.5)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_1p5_4' ] = {
    'bins' :  array('d',[0.0, 0.5]),
    'label' : ' AK8 jet #tau_{4}^{(1.5)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_1p5_5' ] = {
    'bins' :  array('d',[0.0, 0.4]),
    'label' : ' AK8 jet #tau_{5}^{(1.5)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_2_1' ] = {
    'bins' :  array('d',[0., 0.5]),
    'label' : ' AK8 jet #tau_{1}^{(2)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_2_2' ] = {
    'bins' :  array('d',[0., 0.3]),
    'label' : ' AK8 jet #tau_{2}^{(2)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_2_3' ] = {
    'bins' :  array('d',[0., 0.3]),
    'label' : ' AK8 jet #tau_{3}^{(2)}',
    'alignLeg' : 'right'
    }
nSubVariables[ 'Jet'+'_tau_2_4' ] = {
    #'bins' : np.array( [0., .015, .03, .045, .1 ]  ),
    'bins' :   array('d',[0., 0.18]),
    'label' : ' AK8 jet #tau_{4}^{(2)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau_2_5' ] = {
    #'bins' : np.array( [0., .015, .03, .045, .1 ]  ),
    'bins' :   array('d',[0., 0.14]),
    'label' : ' AK8 jet #tau_{5}^{(2)}',
    'alignLeg' : 'right'
    }

nSubVariables[ 'Jet'+'_tau21' ] = {
    'bins' : array('d',[0., 1.2]),
    'label' : ' AK8 jet #tau_{21}',
    'alignLeg' : 'left'
            }
nSubVariables[ 'Jet'+'_tau32' ] = {
    'bins' :  array('d',[0., 1.2]),
    'label' : ' AK8 jet #tau_{32}',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau21_WTA' ] = {
    'bins' : array('d',[0., 1.1]),
    'label' : ' AK8 jet #tau_{21} (WTA)',
    'alignLeg' : 'left'
            }
nSubVariables[ 'Jet'+'_tau32_WTA' ] = {
    'bins' : array('d',[0., 1.1]),
    'label' : ' AK8 jet #tau_{32} (WTA)',
    'alignLeg' : 'left'
    }

nSubVariables[ 'Jet'+'_tau21_exkT' ] = {
    'bins' :  array('d',[0., 1.4]),
    'label' : ' AK8 jet #tau_{21} (excl. kT)',
    'alignLeg' : 'left'
            }
nSubVariables[ 'Jet'+'_tau32_exkT' ] = {
    'bins' :  array('d',[0., 1.4]),
    'label' : ' AK8 jet #tau_{32} (excl. kT)',
    'alignLeg' : 'left'
    }


#########################################################################
############### Binning scheme Purity/Stability optimizer ###############
#########################################################################


def calculate_optimizedBins_purity_stability(aTH2, var, exp_purity, exp_stability,threshold=0.001,PS=OrderedDict(),lowstandards=False):
    
    #initialise dicts per bin to store purity and stability for plotting

    h_resp = aTH2.Clone()
    h_gen = aTH2.ProjectionX()
    h_reco = aTH2.ProjectionY()
    
    resp_array2d,_ = th2_to_np_arr(h_resp)
    resp_array2d_copy = np.copy(resp_array2d)

    # to get old bin edges, assuming at this point that the input TH2 is square with uniform binning on gen/reco axes
    xaxis = h_resp.GetXaxis() 
    yaxis = h_resp.GetYaxis() 
    old_bin_edges = [xaxis.GetBinLowEdge(i) for i  in range(1,xaxis.GetNbins()+2)]
    old_bin_edgesy = [yaxis.GetBinLowEdge(i) for i  in range(1,yaxis.GetNbins()+2)]
    assert(old_bin_edges==old_bin_edgesy)
    new_gen_bin_edges = np.array(old_bin_edges[:])

    totalhistowidth=new_gen_bin_edges[len(new_gen_bin_edges)-1]-new_gen_bin_edges[0]
    h_midpoint = (round(totalhistowidth/2.,3))
    
    
    # Try to find width of distribution
    #mean = h_gen.GetMean()
    #width = np.sqrt(h_gen.GetRMS())

    # Fit the distribution to actually find the width
    #peak = h_gen.GetBinCenter(h_gen.GetMaximumBin())
    #fit_mod = 1
    #fgaus = ROOT.TF1("fgaus", "gaus", peak-fit_mod*h_gen.GetRMS(), peak+fit_mod*h_gen.GetRMS());
    #fgaus.SetNpx(1000)
    #fgaus.SetParameter('mean', peak)
    #fgaus.SetParameter('sigma', width)
    #fgaus.SetLineColor(ROOT.kRed)
    #fgaus.SetLineWidth(1)
    #fit_res = h_gen.Fit("fgaus", "QSR")
    #if fit_res.Status() == 0:
    #    mean = fit_res.Parameter(1)
    #    width = fit_res.Parameter(2)
        
    first_where_nonZero=h_gen.FindFirstBinAbove(0.)*threshold
    last_where_nonZero=(h_gen.FindLastBinAbove(0.))*threshold#+(10 if not h_gen.FindLastBinAbove(0.)>(h_gen.GetNbinsX()-10) else 0. ))*threshold
    totalhistowidth=last_where_nonZero-first_where_nonZero#(h_gen.FindLastBinAbove(0.))*threshold #real total width    @new_gen_bin_edges[len(new_gen_bin_edges)-1]-new_gen_bin_edges[0]
    print (threshold,f"Width, start and end of distribution: {totalhistowidth,first_where_nonZero,last_where_nonZero}, bins corresponding to the same: {h_gen.FindFirstBinAbove(0.), new_gen_bin_edges[h_gen.FindFirstBinAbove(0.)], h_gen.FindLastBinAbove(0), new_gen_bin_edges[h_gen.FindLastBinAbove(0)]} with initial bin widths: {h_gen.GetBinWidth(1)}")

    print("Mean: %.3f" % mean, "width: %.3f" % width, "[fit]" if fit_res.Status() == 0 else "[raw]")
    
    
    ibin = 0 
    c = 0
    
    bc_tot = h_gen.Integral()
    fractionEvt_in_this_bin = 0.
    diagSum=0.
    
    while ibin < len(resp_array2d_copy)-1 and c<1000001:
        c+=1

        resp_array2d_normX = renorm(resp_array2d_copy, axis=0) # normalise axis to 1, renormed per x/gen bin
        resp_array2d_normY = renorm(resp_array2d_copy, axis=1) # normalise axis to 1, renormed per y/reco bin

        purity = resp_array2d_normY[ibin][ibin] #contains fraction of events in a reco bin that are generated in the same gen bin
        stability = resp_array2d_normX[ibin][ibin] #contains fraction of events in a gen bin that are reconstructed in same reco bin

        flag_even_interval = (new_gen_bin_edges[ibin+1]-new_gen_bin_edges[ibin])/threshold % 2 == 0 #(threshold set by initial, uniform bin widths)
    
        if purity>0.3 or stability>0.3: #compute update only when needed
            fractionEvt_in_this_bin = h_gen.Integral(h_gen.FindBin(new_gen_bin_edges[ibin]),c)/bc_tot
        else: fractionEvt_in_this_bin=None
        current_bin_width=h_gen.GetBinLowEdge(c+1)-new_gen_bin_edges[ibin]
                                                                   
        psFlagAND = True if (purity>=exp_purity and stability>=exp_stability) else False
        psFlagOR = True if psFlagAND else False
        
        if lowstandards: 
            psFlagOR = True if (((purity>=0.4 and stability>=0.4) or (purity>=exp_purity and stability>=exp_stability-0.1))) else False
            #if psFlagOR: print(f"Reverting to psFlagOR={psFlagOR} instead of logical AND between P/S ({purity,stability}), since we have lowered P/S standards={lowstandards} for {var}")
        
        
        #if ((psFlagAND) and not(lowstandards) and purity>=exp_purity and stability>=exp_stability and (fractionEvt_in_this_bin>0.05 or current_bin_width>0.1*totalhistowidth)) or ((psFlagOR) and lowstandards and (current_bin_width>0.1*totalhistowidth)):# or fractionEvt_in_this_bin>0.02)):#fractionEvt_in_this_bin>0.04 or # and lowstandards and (purity>resp_array2d_normY[ibin-1][ibin-1] and stability>resp_array2d_normX[ibin-1][ibin-1])):
        if (psFlagAND) or ((psFlagAND) and (fractionEvt_in_this_bin>0.02)) or ((psFlagAND) and (current_bin_width>0.05*totalhistowidth)) or (psFlagOR and lowstandards and (fractionEvt_in_this_bin>=0.002 or current_bin_width>0.01*totalhistowidth)):# or fractionEvt_in_this_bin>0.005)):
            #either stay at high purity/stability (unless specified otherwise via 'lowstandards' for more difficult observables), 
            #or for lowered standard use the psFLagOR whereby the threshold for either of purity/stability is allowed to be greater than 0.5, 
            #while the other is allowed to float down to 0.45
            #or finally, in cases where p/s>0.5 is easily reached, and bin contents accumulated in a bin are greater than 5% of total events
            #then we move onto the next bin to try and prevent spikey distributions which are also therefore overbinned in the tail-end 
            
            x=new_gen_bin_edges[ibin]
            
            if ( (np.round(new_gen_bin_edges[ibin+1]-x,int(abs(np.log10(threshold)))+1) /2.000 /threshold).is_integer() ): 
                #print(new_gen_bin_edges[ibin+1],x,int(abs(np.log10(threshold)))+1, ((np.round(new_gen_bin_edges[ibin+1]-x,int(abs(np.log10(threshold)))+1) /2.000 /threshold).is_integer() ))
                ibin+=1
            
            else:
                # thresholding to prevent accumulation of bins that are not going to be possible to rebin accurately as per
                # how ROOT complains if you rebin a histo to a set of edges that do not correspond to old ones in it
                # ie, we need to take care of, first, when new bin edges the optimizer settles on do not match existing ones in gen histos
                # we'd have issues rebinning; second, keeping in mind reco level has 2x more bins than gen, 
                # we also want to be able to have bins where midpoints between them correspond to existing bin edges in reco
                # so we ensure diff between gen bins edges is an exact multiple of 0.001 (minimal bin width in orig. histos)
                
                while not ((np.round(new_gen_bin_edges[ibin+1]-x,int(abs(np.log10(threshold)))+1)/2.000/threshold).is_integer()):# and not ibin==len(resp_array2d_copy)-3: 

                    #concatenate a row and column corresponding to ibin to push up bin edge
                    resp_array2d_copy = concat_row(resp_array2d_copy, ibin) # row
                    resp_array2d_copy = concat_row(resp_array2d_copy.T, ibin).T # do transposes to apply same as above to columns
                    #remove concatenated bin to keep track of new binning 
                    new_gen_bin_edges = np.delete(new_gen_bin_edges,ibin+1)
                
                    
                resp_array2d_normX = renorm(resp_array2d_copy, axis=0) # normalise axis to 1, renormed per x/gen bin
                resp_array2d_normY = renorm(resp_array2d_copy, axis=1) # normalise axis to 1, renormed per y/reco bin
                purity = resp_array2d_normY[ibin][ibin] #contains fraction in a reco bin that are actually from the same gen bin
                stability = resp_array2d_normX[ibin][ibin] #contains fraction in a gen bin that are actually from the same reco bin
                flag_even_interval = (new_gen_bin_edges[ibin+1]-new_gen_bin_edges[ibin])/threshold % 2 == 0 #(threshold sets min bin width)
      
                fractionEvt_in_this_bin = h_gen.Integral(h_gen.FindBin(new_gen_bin_edges[ibin]),c+1)/bc_tot


                ibin+=1
            
            continue
        
        else:
            
            #concatenate a row and column corresponding to ibin to push up bin edge
            resp_array2d_copy = concat_row(resp_array2d_copy, ibin) # row
            resp_array2d_copy = concat_row(resp_array2d_copy.T, ibin).T # do transposes to apply same as above to columns
            #remove concatenated bin to keep track of new binning 
            new_gen_bin_edges = np.delete(new_gen_bin_edges, ibin+1)
            
            continue
            
    fractionEvt_in_this_bin = h_gen.Integral(h_gen.FindBin(new_gen_bin_edges[ibin]),c+1)/bc_tot

    # to handle the last bins:
    # now, we're in last bin, so remove penultimate bin edge if the last bin is not pure and/or stable enough
    flag=False #flag to say if bin has to be kept or not

    if (purity<exp_purity or stability<exp_stability):# and (exp_purity==0.5 or exp_stability==0.5):
        
        #print(f"trying to salvage a last bin in case either P or S >{exp_purity}")
        if (((purity>=exp_purity-0.05 and stability>=0.4) or (purity>=0.4 and stability>=exp_stability-0.1)) and lowstandards) or (not lowstandards and ((purity>=0.5 and stability>=0.45) or (purity>=0.45 and stability>=0.5))):# and (exp_purity==0.5 or exp_stability==0.5)): 
            #with this we try to hunt down/salvage 'one last bin' in cases where binning is hard to keep pure or stable in some bins/for > 3 bins (ie, keep one of purity/stability>50% and the other not below 40%) 
            flag=True #save last bin and then resize to make it smaller if it's too long
            print ("flagged!!!")
        
        if not flag:
            new_gen_bin_edges = np.delete(new_gen_bin_edges,-2)
            ibin=ibin-1
            

    binresizeTrue=False 
    #if last bin unnecessarily large (irrespective of how the above processing of the last bin turned out), resizing to something that allows for a final bin edge close to the final bin edge as
    #used from last iteration of unfoldings

    lastbinwidth=new_gen_bin_edges[ibin+1]-new_gen_bin_edges[ibin]
    
    print(f"Found that the last binwidth ({lastbinwidth}) is:{lastbinwidth/totalhistowidth*100.000}\% of the total histo's width: {totalhistowidth}, will resize if last bin extends much further than in previous iterations of results.")
    
    
    if flag==True or ((new_gen_bin_edges[ibin+1]>nSubVariables[var]['bins'][-1] or new_gen_bin_edges[ibin+1]>last_where_nonZero+0.1*totalhistowidth) and lastbinwidth/totalhistowidth>0.1):
        #hardcoding to resize if last bin is >10% of full histowidth and/or extends much beyond past known max bin edge, can probably be optimized further, but works for now
        print(new_gen_bin_edges[ibin+1]>nSubVariables[var]['bins'][-1] , new_gen_bin_edges[ibin+1]>last_where_nonZero+threshold*100. , lastbinwidth/totalhistowidth*100.00>0.1)
        print (f"{var} Last bin edge for bin:({ibin+1},{new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]}) extends too far compared to expectations from last version of results where bins end at: {nSubVariables[var]['bins'][-1]}.")
        
        if nSubVariables[var]['bins'][-1]-new_gen_bin_edges[ibin]<lastbinwidth and lastbinwidth/totalhistowidth>0.1:# and nSubVariables[var]['bins'][-1]-new_gen_bin_edges[ibin]>=0.100: 
            #check if diff between 'new' last bin edge and penultimate bin edge would be smaller than the lastbinwidth and >0.1; then set to former
            
            print (f"Resizing last bin from {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]} to {new_gen_bin_edges[ibin]}-{max(nSubVariables[var]['bins'][-1],new_gen_bin_edges[ibin]+0.1)}")
            new_gen_bin_edges[ibin+1]=max(nSubVariables[var]['bins'][-1],new_gen_bin_edges[ibin]+0.1*totalhistowidth)#+0.2
        
        elif lastbinwidth/totalhistowidth>0.1:#else set to latter--ie, +0.1 added to penultimate bin edge
            print (f"Resizing last bin from {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]} to {new_gen_bin_edges[ibin]}-{max(nSubVariables[var]['bins'][-1],new_gen_bin_edges[ibin]+0.1)}")
            new_gen_bin_edges[ibin+1]=max(new_gen_bin_edges[ibin+1],new_gen_bin_edges[ibin]+0.1*totalhistowidth)#+0.2
            
        #elif lastbinwidth<new_gen_bin_edges[ibin]-new_gen_bin_edges[ibin-1]:
        #    print (f"Bin too small, resizing last bin from {new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]} to {new_gen_bin_edges[ibin]}-{max(nSubVariables[var]['bins'][-1],new_gen_bin_edges[ibin]+(new_gen_bin_edges[ibin]-new_gen_bin_edges[ibin-1]))}")           
        #    new_gen_bin_edges[ibin+1]=max(nSubVariables[var]['bins'][-1],new_gen_bin_edges[ibin]+(new_gen_bin_edges[ibin]-new_gen_bin_edges[ibin-1]))
        binresizeTrue=[True,new_gen_bin_edges[ibin+1]]

                
        
    if not type(binresizeTrue)==bool:
        if binresizeTrue[0]==True: 
            old_bin_edges[-1]=binresizeTrue[1]
    #print(f"{chr(10)}Final, new bin edge pairs: {getBinEdgePairs(new_bins=new_gen_bin_edges,old_bin_edges=old_bin_edges)}")
            
            
    #Calculate dicts of P/S for new binngs
    new_bin_edge_pairs = getBinEdgePairs(new_bins=new_gen_bin_edges,old_bin_edges=old_bin_edges)
    rebinned=make_rebinned_2d_hist(h_resp.Clone(),new_bin_edge_pairs,)#rebinning to gen-level bins
    
    rebinned2=make_rebinned_2d_hist(h_resp.Clone(),new_bin_edge_pairs,True)#rebinning to reco and gen-level bins
    
    arr_rebinned,_ = th2_to_np_arr(rebinned.Clone())
    rebinned_array2d_normX = renorm(arr_rebinned, axis=0) # normalise axis to 1, renormed per x/gen bin
    rebinned_array2d_normY = renorm(arr_rebinned, axis=1) # normalise axis to 1, renormed per y/reco bin
    
    axis = rebinned2.GetYaxis() 
    new_reco_bin_edges = [axis.GetBinLowEdge(i) for i  in range(1,axis.GetNbins()+2)]
    
    p_list=[]
    s_list=[]
    
    for ibin in range(len(new_gen_bin_edges)-1):
        purity = rebinned_array2d_normY[ibin][ibin] #contains fraction in a reco bin that are actually from the same gen bin
        stability = rebinned_array2d_normX[ibin][ibin] #contains fraction in a gen bin that are actually from the same reco bin
        fractionEvt_in_this_bin = h_gen.Integral(h_gen.FindBin(new_gen_bin_edges[ibin]),h_gen.FindBin(new_gen_bin_edges[ibin+1]))/bc_tot
        p_list.append(purity)
        s_list.append(stability)
        print (f"Found p/s/frac. evt. per bin: {purity,stability,fractionEvt_in_this_bin}, in final new binning for bin {ibin}:{new_gen_bin_edges[ibin]}-{new_gen_bin_edges[ibin+1]}")

        
    #build P/S dict entriess for purity and stability for each var
    PS[var]={'purity':     p_list,#array('d',
             'stability':  s_list,
             'nbins':      len(s_list),
             'bin_edge_pairs': new_bin_edge_pairs, 
             'genBins':list( new_gen_bin_edges), 
             'recoBins': list(new_reco_bin_edges), 
            }
    
    print (f"Built purity and stability (vs. ibin) dictionary entry for {var}")
    
    return new_reco_bin_edges,new_bin_edge_pairs,new_gen_bin_edges


########################################################################
################# Helpers for running P/S optimization #################
########################################################################

def loadFiles_forPSOptim(year,inpFolder,sampleLabel, t3_samples): 
    #loading files for Purity/Stability optimization
    
    sigFilesALL={}
    lumi=0.    
    
    if not inpFolder:
        print("############################## WARNING #############################")
        print(" Need to provide input directory for histograms root files; exiting ")
        print("####################################################################")
        return -929. 
    
    listOfYears = ['2016_preVFP','2016', '2017', '2018'] if year.startswith('all') else [ year ] 
    
    for iy in listOfYears:
        lumi = lumi + checkDict( ( 'JetHT' if selection.startswith('_dijet') else 'SingleMuon' ), t3_samples )[iy]['lumi']
        
    print (f"Loading MC corresponding to {lumi} /pb of (UL) Run-II data")
    
    for isam in t3_samples:
        for iy in listOfYears:
            if not isam.startswith(sampleLabel): continue
            if not checkDict( isam, t3_samples )[iy]['skimmerHisto'].endswith('root'): continue
            
            print (isam,iy)
            sigFiles = {}
            dictIn = checkDict( isam, t3_samples )[iy]
            dictIn['XS'] = checkDict( isam, t3_samples )['XS']
            sigFiles[isam.split('_Tune')[0]] = [
                            inpFolder+'/FineBins/'+checkDict( isam, t3_samples )[iy]['skimmerHisto'].replace('.root','_FineBins.root'),#ROOT.TFile.Open(  ),
                            dictIn
                    
                        ]

            sigFilesALL[iy] = sigFiles
            print (f"File dict for year={iy}:")
            import pprint
            pprint.pprint(sigFilesALL[iy])
        
    return sigFilesALL,lumi

def closeFiles(listOfFiles, year='all'):
    
    print("Closing files",listOfFiles)
    for i in list(listOfFiles.keys()):
        print (i,f"Closing {list(listOfFiles.keys()).index(i)}/{len(list(listOfFiles.keys()))-1} of the files you asked me to close, namely {i}")
        listOfFiles[i][0].Close()


# get response matrix per year or for all years in one go
def getResponse( inputFiles, outputFolder, variables, sel, year  ): 
    numBins = 0
    numBinsList = [0]
    dict_respWithMissHistos = OrderedDict()
    signalHistos=OrderedDict()

    if not year.startswith('all'):
        for ivar in variables:
            print (f"Getting input histos for: {ivar.split('Jet_')[1]}")    
            ### Getting input histos
            signalHistos = loadHistograms( inputFiles[year], ivar, sel, sysUnc=[], lumi=checkDict('SingleMuon', t3_samples)[year]['lumi'] , 
                                          year=year, process=process, variables=variables, isMC=True, addGenInfo=False, respOnly=True, noRebin=True )

    else:
        for y in (['2016_preVFP', '2016', '2017', '2018']):
            for ivar in variables:
                #print (y,inputFiles[y])
                print(f"Loading {signalLabel+'_respWithMiss'+ivar+'_nom'+sel} for {y} MC, {ivar}")
                
                sigHistos_perYear = loadHistograms(inputFiles[y], ivar, sel, sysUnc=[], 
                                                   lumi=checkDict('SingleMuon', t3_samples)[y]['lumi'], year=y, process=process, 
                                                   variables=variables, isMC=True, addGenInfo=True, respOnly=True, noRebin=True )
                
                h_respWithMiss=sigHistos_perYear[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone()
                h_respWithMiss.SetSumw2()
                h_respWithMiss.SetDirectory(0)
                
                if y=='2016_preVFP' in y:
                    signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel] =  h_respWithMiss.Clone( h_respWithMiss.GetName() )
                                   
                else: 
                    signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel ].Add(h_respWithMiss.Clone())
          
        for ivar in variables:
            tmpHisto = signalHistos[signalLabel+'_respWithMiss'+ivar+'_nom'+sel].Clone(signalLabel+'_respWithMiss'+ivar+'_nom'+sel)
            tmpHisto.SetDirectory(0)
            dict_respWithMissHistos[ivar]=tmpHisto.Clone()

    return dict_respWithMissHistos
    
def getBinEdgePairs(newBins,oldBinEdges):
    binEdgePairs = [list(x) for x in zip(newBins[:-1], newBins[1:])]

    binEdgePairs[-1][1] = oldBinEdges[-1]
    
    return binEdgePairs

#############################################################
########## Simple purity/stability per bin plotter ##########
#############################################################

colorPallete = [ 0, 2, 4, 8, 12, 28, 30 ]

def makePSplot_simple(purity,stability,var,outputDir,ext='pdf',dictHistos=OrderedDict(),bins=[0.,1.],sel='_WSel'):
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    
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
    p.GetXaxis().SetTitle( nSubVariables[var]['label'] )
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
    canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+version+year+'.'+ext)
    if ext.startswith('pdf'):
        canvas.SaveAs(outputDir+var+'_'+signalLabelBegin+sel+'_Purity'+version+'_'+year+'.png')
    return 1#numBins
    
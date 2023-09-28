import time
import json
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea.util import load
from collections import OrderedDict
import numpy as np
import awkward as ak
import uproot
import os
from copy import deepcopy
from pprint import pprint 
import gc
import hist
from hist import Hist

from datasets_WtopSel_RunIISummer20UL_SampleDictPrep import dictSamples, checkDict

######################### import histogram producer/processor ##########################
from nSubBasis_histoProducer_WtopSel import nSubBasis_unfoldingHistoProd_WtopSel as procWtop
########################################################################################




####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--selection', action='store', default='Wtop', help='Selection: WtopSel, dijetSel' )
    parser.add_argument('-v', '--version', action='store', default='V4', help='Version: V1, V2.' )
    parser.add_argument('-p', '--verbosePrint', action='store_true', default=False, help='Verbosity, ie., with or without printouts' )
    parser.add_argument('--noWtUnc', action='store_true', default=False, help='Produce histograms without wt. uncs in signal sample if True' )
    parser.add_argument('-y', '--year', action='store', default='2017', help='Year: 2016_preVFP, 2016, 2017, 2018, all.' )
    parser.add_argument('-e', '--ext', action='store', default='nomAndWts', help='Suffix to filename: NomAndWts, jesCorr,...' )
    parser.add_argument('-j', '--nWorkers', action='store', default=30, help='nWorkers for uproot concatenation,default=30' )
    parser.add_argument("--outputFolder", action='store', dest="outputFolder", default="/work/kadatta/private/CMSSW_10_6_29/src/nSubHistogramming/V9/", help="Output folder" )

    try: 
        args = parser.parse_args()
        if not(os.path.exists(args.outputFolder)):
            os.makedirs(args.outputFolder)
            os.makedirs(args.outputFolder+'WSel')
            os.makedirs(args.outputFolder+'topSel')
    except:
        parser.print_help()
        sys.exit(0)

    year_list = [y for y in ['2016', '2017', '2018', '2016_preVFP'] if (args.year.lower().endswith('all')) else args.year]

    if args.selection=='Wtop':

        ####################################################################################################################################
        ################## Import dictionary of samples with details per sample of number of events, file locations on T3 ##################
        ####################################################################################################################################
        ####################### this dictionary includes details of samples containing nominal or nominal+genWeights #######################
        ####################################################################################################################################
        with open('Unfolding_WtopDict_FinalNominal.json') as dicto:                                                                        #
            t3_samples=json.load(dicto)                                                                                                    #
        if args.verbosePrint:                                                                                                              #
            pprint (t3_samples)                                                                                                            #
        ####################################################################################################################################
        # include reading dictos for systematics as well, trivial and left for later
        ####################################################################################################################################
        #with open('Unfolding_WtopDict_FinalNominal.json') as dicto:                                                                        #
        #    t3_samples=json.load(dicto)                                                                                                    #
        #if self.verbose:                                                                                                                   #
        #    pprint (t3_samples)                                                                                                            #
        ####################################################################################################################################
        ####################### this dictionary includes details of samples containing nominal or nominal+genWeights #######################
        ####################################################################################################################################
        #with open('Unfolding_WtopDict_FinalNominal.json') as dicto:                                                                        #
        #    t3_samples=json.load(dicto)                                                                                                    #
        #if self.verbose:                                                                                                                   #
        #    pprint (t3_samples)                                                                                                            #
        ####################################################################################################################################
        ####################################################################################################################################

        for i in t3_samples.keys():
            gc.collect()    
            for y in year_list:#
                if args.verbosePrint: print (y)
        
                if ('TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8' in i): 
                
        
                    tstart=time.time()
                    
                    procTop = procWtop( sampleName=i,#'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8', 
                                        sampleDict=t3_samples, selection='_topSel', 
                                        isMC=True if not('SingleMuon' in i) else False, 
                                        isSigMC=True if (i.startswith('TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8')) else False , 
                                        year=y,
                                        wtUnc=True if (i.startswith('TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8') and not(args.noWtUnc)) else False, 
                                        sysUnc=False, verbose=False, saveParquet=False,splitCount='0')
                    
                    procW = procWtop(  sampleName=i,#'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8', 
                                       sampleDict=t3_samples,selection='_WSel', 
                                       isMC=True if not('SingleMuon' in i) else False, 
                                       isSigMC=True if (i.startswith('TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8')) else False, 
                                       year=y,
                                       wtUnc=True if (i.startswith('TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8') and not(args.noWtUnc)) else False, 
                                       sysUnc=False, verbose=False, saveParquet=False,splitCount='0')

                    
                    events = uproot.concatenate(t3_samples[i][y]['t3_dirs'][0].split('/0000/')[0]+'/000*/'+'jetObservables_nanoskim_*.root:Events',  
                                                procTop._branchesToRead, step_size="2048 MB", 
                                                library='ak', num_workers=args.nWorkers, 
                                                #decompression_executor=uproot.ThreadPoolExecutor(num_workers=32),
                                                #interpretation_executor=uproot.ThreadPoolExecutor(num_workers=32),
                                               )
                    elapsed = time.time() - tstart

                    if args.verbosePrint:          
                        print (f'Time for loading events for {i} in {y}:{elapsed}')   
                    
                    tstart=time.time()
                    
                    if args.verbosePrint:          
                        print(i,y,'topSel')
                    
                    
                    ######### processing tops #######
                    histosTop=procTop.process(events)
                    #################################
                    
                    if args.verbosePrint:          
                        print(histosTop['respWithMissJet_tau_0p5_1_nom_topSel'].sum(flow=True),histosTop['respWithMissJet_tau_0p5_1_nom_topSel'].sum(flow=False))
                        print(histosTop['recoJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))
                        print(histosTop['truerecoJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))
                        print(histosTop['fakerecoJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))
                        print(histosTop['genJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))
                        print(histosTop['accepgenJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))
                        print(histosTop['missgenJet_pt_nom_topSel'].sum(flow=True))#,histosTop['respJet_pt_nom_topSel'].sum(flow=False))

                    with uproot.recreate(f'{args.outputFolder}/topSel/jetObservables_histograms_{i.split("_13TeV")[0]}_{y}_topSel.root') as fout:
                        for key in histosTop.keys():
                            fout[key]=histosTop[key]
                        print (f'Done with creating output file:{fout}')
                        fout.close()
                    elapsed = time.time() - tstart

                    print (f'Time for topSel {procTop.sampleName} {procTop.year}:{elapsed}')        
                    #del(events)
                    del(procTop)
                    del(histosTop)
                    
                    gc.collect()
                    #############################################################################################################################################################
                    
                    
                    tstart=time.time()
                   
                    print(i,y,'WSel')
                    
                    ###### processing W's #######
                    histosW=procW.process(events)
                    #############################

                    if args.verbosePrint:
                        print(histosW['respWithMissJet_tau_0p5_1_nom_WSel'].sum(flow=True),histosW['respWithMissJet_tau_0p5_1_nom_WSel'].sum(flow=False))
                        print(histosW['recoJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))
                        print(histosW['truerecoJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))
                        print(histosW['fakerecoJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))
                        print(histosW['genJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))
                        print(histosW['accepgenJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))
                        print(histosW['missgenJet_pt_nom_WSel'].sum(flow=True))#,histosW['respJet_pt_nom_WSel'].sum(flow=False))

                    with uproot.recreate(f'{args.outputFolder}/WSel/jetObservables_histograms_{i.split("_13TeV")[0]}_{y}_WSel.root') as fout:
                        for key in histosW.keys():
                            fout[key]=histosW[key]
                        print (f'Done with creating output file:{fout}')
                        fout.close()
                        
                    elapsed = time.time() - tstart

                    print (f'Time for WSel {procW.sampleName} {procW.year}:{elapsed}')   
                    
                    del(events)
                    del(histosW)
                    del(procW)
                    gc.collect()
                    
                    

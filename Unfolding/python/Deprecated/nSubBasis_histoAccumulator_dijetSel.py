from itertools import islice
#import concurrent
#from concurrent import futures
#from coffea.processor import FuturesExecutor,Runner
from coffea.processor import accumulate
#from uproot import ThreadPoolExecutor

#from coffea.processor import defaultdict_accumulator,dict_accumulator#,IterativeExecutor,FuturesExecutor
import time
import re

########modify for data, potentially just put this in as the postprocessor in the histoproducer module
#{writeAccumulatedORhadd: 0=write accumulated hists in a separate file after saving or not saving constituent file chunks,
#1 = hadd a posteriori the chunks}
def histoMaker_BasicButSplit(myProcessor, sampleIdentifier='qcd_ht', y='2017', sampleDict_PFNano=OrderedDict(),
                             sampleDict_local=OrderedDict(), writeChunks=True, writeAccumulatedORhadd=1, verbose=False, saveParquet=False,
                             isSigMC=True, isMC=True, wtUnc=True, sysUnc=False, onlyUnc='', outputdir='processorTests/', sysSource=[],
                             era='', ext='_nomWts', splitchunks=5, nWorkers=10, stepSize="2048 MB",jetType='Central'):#for jetType: options='Central','Forward' (1 at a time)
    nchunk=copy.deepcopy(splitchunks)
    for sample in sampleDict_local.keys():
        yl=['2016_preVFP', '2016', '2017', '2018'] if y=='all' else [y]
        #if verbose: print (sample,sampleIdentifier,splitchunks)
        
        for year in yl:#['2017','2018','2016','2016_preVFP']:
            gc.collect()
            
            if not( sample.lower().startswith(sampleIdentifier.lower())):# or ('100to200' in sample.lower()) or ('200to300' in sample.lower()) or ('300to500' in sample.lower()) or ('500to700' in sample.lower()): 
                continue
            tstart0=time.time()

            print(f'Histogramming for {sample} in year: {year}')
            if verbose: print(f'using nanoskims from {sampleDict_local[sample][year]["t3_dirs"]}{chr(10)}','\n')
            
            
            if not os.path.exists(f'{outputdir}'): os.makedirs(f'{outputdir}')

            inpdir=sampleDict_local[sample][year]['t3_dirs']
            dirfnames=inpdir[0].split('/0000/')[0]+'/000*/jetObservables_nanoskim_*.root'
            #if verbose: print (dirfnames)
            SC='0'
            if verbose: print (sysSource)
            my_processor=myProcessor(sampleName=sample, sampleDict=sampleDict_PFNano,
                                     isMC=isMC, isSigMC=isSigMC,year=year,saveParquet=saveParquet,sysSource=sysSource,
                                     wtUnc=wtUnc,sysUnc=sysUnc,onlyUnc='' if 'jes' in onlyUnc else onlyUnc,era=era,
                                     verbose=False,splitCount=SC )
            if splitchunks>0:
                fl=[]
                n_tosplit = []
                c=0
                for x in range(len(inpdir)):
                    flist=os.listdir(inpdir[x])
                    flist = sorted([i for i in flist if 'nanoskim' in i], key=lambda s: int(re.search(r'\d+', s).group()))
                    
                    for i in flist:
                        if 'nanoskim' in i:

                            fl.append(inpdir[x]+i)
                            if c%splitchunks==0 and c!=0:
                                n_tosplit.append(splitchunks)
                            c=c+1
                
                if sum(n_tosplit)!=len(fl):
                    if len(fl)%splitchunks!=0: n_tosplit.append(len(fl)%splitchunks)
                    else: n_tosplit.append(splitchunks)

                if verbose: print(f"Splitting input filelist into the following sublists of file chunks:{n_tosplit}")
                ifl=iter(fl)
                ifls=[list(islice(ifl,x)) for x in n_tosplit]
                #if verbose: print(ifls)
                c=0
                hists=None
            else:
                ifls=[dirfnames]
                c=0
            if verbose: print(f"Loading events into arrays")# from the following file chunks:{n_tosplit}")    
            
            
            for i in track(ifls):#
                #print (c,i)
                if splitchunks>0:
                    ns = [int(x.split('jetObservables_nanoskim_')[1].split('.root')[0]) for x in i]
                    s = str(ns).split('[')[1].split(']')[0]
                    stringfnames=i[0].split('jetObservables_nanoskim_')[0]+"jetObservables_nanoskim_{"+f'{s}'+"}.root"
                else: stringfnames=i
                    
                stringfnames=stringfnames.replace(' ','')
                
                if verbose:
                    tstart1 = time.time()
                
                events = uproot.concatenate(stringfnames+':Events', my_processor._branchesToRead, 
                                            step_size=stepSize, library='ak', num_workers=nWorkers,
                                            #interpretation_executor=uproot.ThreadPoolExecutor(num_workers=nWorkers), 
                                            #decompression_executor=uproot.ThreadPoolExecutor(num_workers=nWorkers),
                                           )
                
                if verbose:
                    elapsed1 = time.time()-tstart1
                    print(f"Time taken to load {len(events)} events from {n_tosplit[c]} files in {sample} = {elapsed1}")
                
                if verbose: print(f'nEvents in this file chunk: {len(events)}')
                    
                if verbose: tstart2=time.time()

                if verbose: print("#####Processing events#####")
                    
                processed_events=my_processor.process(events)
                 
                if verbose:
                    elapsed2 = time.time()-tstart2
                    print(f"Time taken to build histos from {len(events)} events = {elapsed2}")
                
                my_processor.splitCount=chr(int(SC)+1)
                
                if c==0: 
                    hists=processed_events
                    if verbose: print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                elif c>0 and writeAccumulatedORhadd==0: 
                    hists=accumulate([hists,processed_events])
                    if verbose: print(f"Current integral of pT nominal hist: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")
                
                if not sysUnc: 
                    if 'pt' in sample.lower(): 
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]+sample.split("_Tune")[0].split("_")[2]}_UL{year}{ext}'
                    else:
                        string=f'{sample.split("_Tune")[0].split("_")[0]+sample.split("_Tune")[0].split("_")[1]}_UL{year}{ext}'

                    fnstem = f'{outputdir}/jetObservables_histograms_{string}'    
                    fn = f'{fnstem}_ForwardJet_{c}.root' if 'Forward' in jetType else f'{fnstem}_CentralJet_{c}.root'
                else: 
                    
                    fnstem = f'{outputdir}/{sampleDict_local[sample][year]["skimmerHisto"].split(".root")[0]}{ext}'
                    if 'jes' in onlyUnc: 
                        fnstem=f'{outputdir}/combinedJES/{fnstem.split(outputdir+"/")[1]}' 
                        if not os.path.exists(f'{outputdir}/combinedJES/'): os.makedirs(f'{outputdir}/combinedJES/')

                    fn = f'{fnstem}_ForwardJet_{c}.root' if 'Forward' in jetType else f'{fnstem}_CentralJet_{c}.root'
                    
                if splitchunks==0 or len(n_tosplit)==1: 
                    fn=f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root'
                if verbose: print (fn,fnstem)

                if writeChunks or writeAccumulatedORhadd==1:# or len(n_tosplit)==1:
                    if verbose: print(f"Current integral of pT nominal hists being saved: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")            
                    

                    #print (processed_events)
                    print ("Events proccesed, now writing output file(s)",fn)
                    
                    with uproot.recreate(fn) as fout:#outputTest_{sample.split("_Tune")[0]}_{year}.root'
                        for key in processed_events.keys():
                            fout[key]=processed_events[key]
                        print (f'Done with creating output file:{fn}')
                        fout.close()
                    
                c=c+1
                del(processed_events)
                del(events)
                gc.collect()
                
            if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):

                print(f"Integral of (reco) nominal hists being saved: {hists['recoJet_pt_nom_dijetSel'].sum(flow=True)}")            
                
                if not os.path.exists(f'{outputdir}/combinedJES/'): os.makedirs(f'{outputdir}/combinedJES/')
                
                print (f'Output file stem for accumulated histos: {fnstem}')
                with uproot.recreate(f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root') as fout_hadd:
                    for key in hists.keys():
                        fout_hadd[key]=hists[key]
                    print (f'Done with creating accumulated output file:{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root')
                    fout_hadd.close()
            #elif writeAccumulatedORhadd==1 and c>0 and splitchunks>0:
            #    if 'jes' in onlyUnc: 
            #        fnstem=f'{outputdir}/combinedJES/{fnstem.split(outputdir+"/")[1]}'
            #    input_haddfilelist = []
            #    for file in os.listdir(fnstem.split('jetObservables_')[0]):
            #        full_filepath=os.path.abspath(fnstem.split('jetObservables_')[0])+file
            #        if fnstem in full_filepath: 
            #            input_haddfilelist.append(file)
            #    print (f'Input files for hadded histos: {input_haddfilelist}')
            #    events = uproot.concatenate(input_haddfilelist)
            #    with uproot.recreate(f'{fnstem}.root') as fout_hadd:
            #        for key in events.fields():#keys():
            #            fout_hadd[key]=events[key]
            #        print (f'Done with creating accumulated output file:{fnstem}.root')
            #        fout_hadd.close()
            #else:
            #    print("separated output files written, the rest is in your hands! :)")
                
            elapsed0 = time.time() - tstart0

            print (f'Time for {sample} {year}:{elapsed0}')        

            del(my_processor)
            #del(events)
            #del(processed_events)
            del(hists)
            gc.collect()
            
            clear_output(wait=True)


    return 1


from itertools import islice
#import concurrent
#from concurrent import futures
from coffea.processor import FuturesExecutor,Runner
from coffea.processor import accumulate

########modify for data, potentially just put this in as the postprocessor in the histoproducer module
#{writeAccumulatedORhadd: 0=write accumulated hists in a separate file after saving or not saving constituent file chunks,
#1 = hadd a posteriori the chunks}
def histoMaker_BasicButSplit(myProcessor, sampleIdentifier='qcd_ht', y='2017', sampleDict_PFNano=OrderedDict(),
sampleDict_local=OrderedDict(), writeChunks=True, writeAccumulatedORhadd=1, 
isSigMC=True, isMC=True, wtUnc=True, sysUnc=False, onlyUnc='', outputdir='processorTests/', era='', ext='_nomWts', splitchunks=5, nWorkers=10, jetType='Central'):#for jetType: options='Central','Forward' (1 at a time)

    for sample in sampleDict_local.keys():
        yl=['2016_preVFP', '2016', '2017', '2018'] if y=='all' else [y]
        print (sample)
        
        for year in yl:#['2017','2018','2016','2016_preVFP']:
            if ('500to700' in sample and 'MLM' in sample and '2018' in year ) or ('100to200' in sample and 'MLM' in sample and '2018' in year):
                splitchunks=250
                
            
            if not( sample.lower().startswith(sampleIdentifier)):# or ('100to200' in sample.lower()) or ('200to300' in sample.lower()) or ('300to500' in sample.lower()) or ('500to700' in sample.lower()): 
                continue
            tstart=time.time()
            #print (sample)

            print(f'{chr(10)} Histogramming for {sample} in year: {year}, using nanoskims from {sampleDict_local[sample][year]["t3_dirs"]}{chr(10)}############################################################################')
            #year='2017'
            
            if not os.path.exists(outputdir): os.makedirs(outputdir)
                
            inpdir=sampleDict_local[sample][year]['t3_dirs']
            dirfnames=inpdir[0].split('/0000/')[0]+'/000*/jetObservables_nanoskim_*.root'
            print (dirfnames)
            my_processor=myProcessor(sampleName=sample, sampleDict=sampleDict_PFNano,
                                     isMC=isMC, isSigMC=isSigMC,year=year,
                                     wtUnc=wtUnc,sysUnc=sysUnc,onlyUnc=onlyUnc,era=era )
            if splitchunks>0:
                fl=[]
                n_tosplit = []
                c=0
                for x in range(len(inpdir)):
                    flist=os.listdir(inpdir[x])
                    flist = sorted([i for i in flist if 'nanoskim' in i], key=lambda s: int(re.search(r'\d+', s).group()))
                    #print (flist)
                    #print (flist)
                    for i in flist:
                        if 'nanoskim' in i:

                            fl.append(inpdir[x]+i)
                            if c%splitchunks==0 and c!=0:
                                n_tosplit.append(splitchunks)
                            c=c+1
                print(fl,n_tosplit)

                if sum(n_tosplit)!=len(fl):
                    print (len(fl), sum(n_tosplit), len(fl)-sum(n_tosplit))
                    if len(fl)%splitchunks!=0: n_tosplit.append(len(fl)%splitchunks)
                    else: n_tosplit.append(splitchunks)

                print(f"Splitting input filelist into the following sublists of file chunks:{n_tosplit}")
                ifl=iter(fl)
                ifls=[list(islice(ifl,x)) for x in n_tosplit]
                #print(ifls)
                c=0
                hists=None
            else:
                ifls=[dirfnames]
                c=0
            print(f"Processing events from the following file chunks:{n_tosplit}")    
            for i in ifls:
                #print (c,i)
                if splitchunks>0:
                    ns = [int(x.split('jetObservables_nanoskim_')[1].split('.root')[0]) for x in i]
                    #print (ns)
                    s = str(ns).split('[')[1].split(']')[0]
                    #print (s)
                    stringfnames=i[0].split('jetObservables_nanoskim_')[0]+"jetObservables_nanoskim_{"+f'{s}'+"}.root"
                    #print (stringfnames)
                else: stringfnames=i
                
                events = uproot.concatenate(stringfnames+':Events', my_processor._branchesToRead, step_size="1024 MB", library='ak', num_workers=nWorkers)
                processed_events=my_processor.process(events)
                
                if c==0: 
                    hists=copy.deepcopy(processed_events)
                elif c>0 and writeAccumulatedORhadd==0: 
                    hists=accumulate([hists,copy.deepcopy(processed_events)])
                
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
                    fn = f'{fnstem}_ForwardJet_{c}.root' if 'Forward' in jetType else f'{fnstem}_CentralJet_{c}.root'
                    
                if splitchunks==0 or len(n_tosplit)==1: 
                    fn=f'{fnstem}_ForwardJet.root' if 'Forward' in jetType else f'{fnstem}_CentralJet.root'
                print (fn,fnstem)

                if writeChunks or writeAccumulatedORhadd==1:# or len(n_tosplit)==1:
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
            if c==0 or (splitchunks>0 and writeAccumulatedORhadd==0):# and not len(n_tosplit)==1:
                if 'jes' in onlyUnc: 
                    fnstem=f'{outputdir}/combinedJES/{fnstem.split(outputdir+"/")[1]}'
                
                print (f'Output file for accumulated histos: {fnstem}.root')
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
            else:
                print("separated output files written, the rest is in your hands! :)")
                
            elapsed = time.time() - tstart

            print (f'Time for {sample} {year}:{elapsed}')        

            del(my_processor)
            #del(events)
            #del(processed_events)
            del(hists)

    return 1


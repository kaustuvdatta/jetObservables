##########################################################################
def loadHistograms( samples, var, sel, sysUnc=[], isMC=True, addGenInfo=True, respOnly=False, lumi=1., variables={}, year='2017', process='data',withRebin=False ):
    """docstring for loadHistograms"""

    if sysUnc==[]: SYSUNC = [ '_nom' ] 
    else: SYSUNC = [ s+u for u in ['Up', 'Down'] for s in sysUnc if not s.startswith(('_model', '_hdamp', '_Tune', '_CR', '_erdON', '_mtop')) ]
    flip = False
    
    allHistos = {}
    for isam in samples:
        if sysUnc!=[]:
            for i in ([ '_jer']):#, '_isrWeight', '_fsrWeight', '_puWeight', '_pdfWeight', '_jesTotal' ]):
                if i in isam:
                    flip=True
                    tmpSYSUNC = [i+u for u in ['Up','Down']]
                    continue
                    
        if not flip: tmpList = [ 'reco'+var+syst+sel for syst in SYSUNC]
        else: tmpList = ['reco'+var+syst+sel for syst in tmpSYSUNC]
        #print (tmpList)
        #if isMC and addGenInfo: tmpList = tmpList + [ 'gen'+var+sel ] + [ 'resp'+var+syst+sel for syst in SYSUNC ]
        if isMC and addGenInfo and not flip: tmpList = tmpList + [ 'gen'+var+sel] + [ 'resp'+var+syst+sel for syst in SYSUNC]
        
        elif isMC and addGenInfo and flip: 
            tmpList = tmpList + [ 'gen'+var+sel]
            for syst in tmpSYSUNC:
                tmpList = tmpList + ['resp'+var+syst+sel]
        
        if respOnly and not flip: 
            tmpList = [ 'resp'+var+syst+sel for syst in SYSUNC] 
            #print (SYSUNC,tmpList)
        elif respOnly and flip: 
            tmpList = ['resp'+var+syst+sel for syst in tmpSYSUNC ] #+ [ 'missgen'+var+sel ]
            #print (tmpSYSUNC,tmpList)
        
        for ih in tmpList:
            #print ('Processing '+isam+' '+ih, flip)
            #print (samples[isam][0])
            if isMC:
                allHistos[isam+'_'+ih] = samples[isam][0].Get( 'jetObservables/'+ih )
                tmpIsam = 'TT' if isam.startswith('data') else isam
                MCScale = samples[isam][1]['XS'] * lumi / samples[isam][1][year][('nGenWeights') ]
                allHistos[isam+'_'+ih].Scale( MCScale )
            else:
                if sel.startswith('_dijet'):# and process.startswith('data'):
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

                    for biny in range( 0, allHistos[isam+'_'+ih].GetNbinsY()+1 ):
                        by = allHistos[isam+'_'+ih].GetYaxis().GetBinCenter( biny )
                        for binx in range( 0, allHistos[isam+'_'+ih].GetNbinsX()+1 ):
                            bx = allHistos[isam+'_'+ih].GetXaxis().GetBinCenter(binx)
                            for iX in range(0, len(genBin)-1 ):
                                for iY in range(0, len(recoBin)-1 ):
                                    if (bx<genBin[iX+1] and bx>genBin[iX]) and (by<recoBin[iY+1] and by>recoBin[iY]):
                                        jbin = allHistos[isam+'_'+ih].GetBin(binx,biny)
                                        tmpArrayContent[iX][iY] = tmpArrayContent[iX][iY] + allHistos[isam+'_'+ih].GetBinContent( jbin )
                                        tmpArrayError[iX][iY] = tmpArrayError[iX][iY] + ROOT.TMath.Power( allHistos[isam+'_'+ih].GetBinError( jbin ), 2 )

                    for biny in range( 0, tmpHisto.GetNbinsY()+1 ):
                        for binx in range( 0, tmpHisto.GetNbinsX()+1 ):
                            #if (binx <= ((biny+1)/2.)+4) and (binx >= ((biny+1)/2.)-4):
                            tmpHisto.SetBinContent( tmpHisto.GetBin(binx,biny), tmpArrayContent[binx-1][biny-1] )
                            tmpHisto.SetBinError( tmpHisto.GetBin(binx,biny), ROOT.TMath.Sqrt(tmpArrayError[binx-1][biny-1] ) )

                    tmpHisto.Sumw2()
                    allHistos[isam+'_'+ih] = tmpHisto

                #if isMC: allHistos[isam+'_'+ih].Scale( MCScale )
                ##### For tests, projections directly from 2D
                genBin = variables[var]['bins']
                
                allHistos[isam+'_accepgen'+var+sel] = allHistos[isam+'_'+ih].ProjectionX()
                allHistos[isam+'_truereco'+var+'_nom'+sel+'_genBin'] = allHistos[isam+'_'+ih].ProjectionY()
                allHistos[isam+'_truereco'+var+'_nom'+sel+'_genBin'].Rebin(len(genBin)-1,'_rebinned', genBin)
                #print (allHistos)
                allHistos[isam+'_missgen'+var+sel] = allHistos[isam+'_gen'+var+sel].Clone()
                allHistos[isam+'_missgen'+var+sel].Add(allHistos[isam+'_accepgen'+var+sel],-1)
                allHistos[isam+'_fakereco'+var+'_nom'+sel+'_genBin'] = allHistos[isam+'_reco'+var+'_nom'+sel+'_genBin'].Clone()
                allHistos[isam+'_fakereco'+var+'_nom'+sel+'_genBin'].Add(allHistos[isam+'_truereco'+var+'_nom'+sel+'_genBin'],-1)

    if sel.startswith('_dijet'):
        tmpHistos = { k:v for (k,v) in allHistos.items() if 'Inf' in k }
        for ih in tmpHistos:
            for jh in allHistos:
                if (jh.endswith('0'+ih.split('Inf')[1])) and not ('Inf' in jh ):
                    tmpHistos[ih].Add( allHistos[jh].Clone() )
        if len(tmpHistos)>0: allHistos = tmpHistos
    '''
    if len(samples)<2:
        
        #print (samples[isam][1][year]['skimmerHisto'])
        outFileName = samples[isam][1][year]['skimmerHisto'].split('.root')[0]+'_loadedHistos.root'
        print ("Written loaded/rebinned histos to:", outFileName)
        outFile = ROOT.TFile.Open ( '../test/Samples/'+outFileName ," RECREATE ")
        outHistFile.cd()
        for h in allHistos:
            h.Write()
        outFile.Close()
    '''
    return allHistos



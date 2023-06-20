import ROOT
import math, os, sys
import numpy as np
from collections import OrderedDict
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class triggerEfficiencies(Module):

    def __init__(self, year='2016'):
        self.writeHistFile=True
        self.year = year
        #this is historical
        self.triggers = {
                'AK8PFJet80' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 14504.85901 ],
                        'v2': [ 297557, 299329, 82434.76397 ],
                        },
                    },
                'AK8PFJet140' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 154.3970862 ],
                        'v2': [ 297557, 299329, 1761.200357 ],
                        },
                    },
                'AK8PFJet200' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 10.42278681 ],
                        'v2': [ 297557, 299329, 355.1960767 ],
                        },
                    },
                'AK8PFJet260' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 52.55889822 ],
                        'v2': [ 297557, 299329, 60.42317971 ],
                        },
                    },
                'AK8PFJet320' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 20.95316709 ],
                        'v2': [ 297557, 299329, 23.83882777 ],
                        },
                    },
                'AK8PFJet400' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 1.001746325 ],
                        'v2': [ 297557, 299329, 1.001043199 ],
                        },
                    },
                'AK8PFJet450' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 1.039043707 ],
                        'v2': [ 297557, 299329, 1.154728156 ],
                        },
                    },
                'AK8PFJet500' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 1. ],
                        'v2': [ 297557, 299329, 1. ],
                        },
                    },
                'AK8PFJet550' : {
                    '2017' : {
                        'v1': [ 297046, 297505, 1. ],
                        'v2': [ 297557, 299329, 1. ],
                        },
                    }
                }

        
    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        ### Booking histograms
        self.addObject( ROOT.TH1F('cutflow',   ';cutflow',   10, 0, 10) )
        #### general selection
        self.addP4Hists( 'AK8Jet', '_check' )
        for t in [ x+'_baseline' for x in self.triggers ] + self.triggers.keys():
            self.addObject( ROOT.TH1F( 'AK8Jet1Pt_'+t, ';p_{T} (GeV)',   200, 0, 2000) )
            self.addObject( ROOT.TH1F( 'AK8Jet2Pt_'+t, ';p_{T} (GeV)',   200, 0, 2000) )


    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -4.0, 4.0 ) )
        self.addObject( ROOT.TH1F(s+'_y'+t, s+';y', 100, -4.0, 4.0 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14159) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 1000) )


    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''

        self.isMC = event.run == 1

        AK8jets = list(Collection(event, 'FatJet' ))
        trig = list(Collection(event, 'TrigObj' ))
        listOfTriggers = self.triggers.keys()

        #### Basic AK8 jet selection
        #Shouldn't affect that this is different from measurement PS, since this is done with the leading jet -ale
        ak8jets = [ x for x in AK8jets if abs(x.eta) < 2.5 ]
        AK8HT = sum( [ x.pt for x in ak8jets ] )
        ak8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        for ijet in ak8jets:
            getattr( self, 'AK8Jet_pt_check' ).Fill( ijet.pt )
            getattr( self, 'AK8Jet_eta_check' ).Fill( ijet.eta )
            getattr( self, 'AK8Jet_y_check' ).Fill( ijet.p4().Rapidity()) 
            getattr( self, 'AK8Jet_phi_check' ).Fill( ijet.phi )
            getattr( self, 'AK8Jet_mass_check' ).Fill( ijet.msoftdrop )

        if ( len(ak8jets)>0 ):
        #if ( getattr( event, 'HLT_'+self.baseTriggers[self.year][0] )==1 ) and ( len(ak8jets)>0 ):
            #print( getattr(event, 'TrigObj_id'  )[1] )
            '''historical 130-133'''
            triggerVersion = ''
            for v in self.triggers['AK8PFJet80'][self.year]:
                if self.triggers['AK8PFJet80'][self.year][v][0] <= event.run <= self.triggers['AK8PFJet80'][self.year][v][1]:
                    triggerVersion=v

            getattr( self, 'cutflow' ).Fill( 0 )
            #getattr( self, 'AK8Jet1Pt_baseTrigger' ).Fill( ak8jets[0].pt )
            #if ( len(ak8jets)>1 ): getattr( self, 'AK8Jet2Pt_baseTrigger' ).Fill( ak8jets[1].pt )

            #### Loop over triggers
            for it in range(len(listOfTriggers)-1):
                getattr( self, 'cutflow' ).Fill( it )

                if ( getattr( event, 'HLT_'+listOfTriggers[it] )== 1 ): 
                    # meant to check if given trigger fires FOR THE LEADING JET
                    getattr( self, 'AK8Jet1Pt_'+listOfTriggers[it+1]+'_baseline' ).Fill( ak8jets[0].pt, self.triggers[ listOfTriggers[it] ][self.year][triggerVersion][2] )

                    if ( len(ak8jets)>1 ): 
                        getattr( self, 'AK8Jet2Pt_'+listOfTriggers[it+1]+'_baseline' ).Fill( ak8jets[1].pt, self.triggers[ listOfTriggers[it] ][self.year][triggerVersion][2]  )

                    if ( getattr( event, 'HLT_'+listOfTriggers[it+1] )==1 ): #CHECK TRIGGER ABOVE THE GIVEN ONE IN LOOP
                        # CHECK IF below triggers both given and given+1 trigger
                        getattr( self, 'AK8Jet1Pt_'+listOfTriggers[it+1] ).Fill( ak8jets[0].pt, self.triggers[ listOfTriggers[it+1] ][self.year][triggerVersion][2]  )
                    
                        if ( len(ak8jets)>1 ): 
                            getattr( self, 'AK8Jet2Pt_'+listOfTriggers[it+1] ).Fill( ak8jets[1].pt, self.triggers[ listOfTriggers[it+1] ][self.year][triggerVersion][2]  )

#            for i, it in enumerate(self.triggers):
#                if ( getattr( event, 'HLT_'+it )==1 ):
#                    for x in trig:
#                        print( x.id, x.filterBits, x.pt )
#                        if x.id==6: print( x.id, x.pt )
#                    print( event.run, event.luminosityBlock, int(event.event), 'pass', it )
#                    getattr( self, 'cutflow' ).Fill( i+1 )
#                    getattr( self, 'AK8Jet1Pt_'+it ).Fill( ak8jets[0].pt )
#                    if ( len(ak8jets)>1 ): getattr( self, 'AK8Jet2Pt_'+it ).Fill( ak8jets[1].pt )
#                    #continue
#                #else: print( event.event, 'fail', it )

        return True

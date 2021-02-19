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
        self.baseTriggers = {
            '2016' : [ 'IsoMu24', 'IsoTkMu24'],
            '2017' : [ 'AK8PFJet80' ],  #'IsoMu27' ],
            '2018' : [ 'IsoMu24' ],
        }
        self.triggers = [ 'AK8PFJet80', 'AK8PFJet140','AK8PFJet200', 'AK8PFJet260', 'AK8PFJet320', 'AK8PFJet400', 'AK8PFJet450', 'AK8PFJet500', 'AK8PFJet550']


    #############################################################################
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        ### Booking histograms
        self.addObject( ROOT.TH1F('cutflow',   ';cutflow',   10, 0, 10) )
        #### general selection
        self.addP4Hists( 'AK8Jet', '_check' )
        for t in ['baseTrigger'] + self.triggers:
            self.addObject( ROOT.TH1F( 'AK8Jet1Pt_'+t, ';p_{T} (GeV)',   200, 0, 2000) )
            self.addObject( ROOT.TH1F( 'AK8Jet2Pt_'+t, ';p_{T} (GeV)',   200, 0, 2000) )


    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -4.0, 4.0 ) )
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

        #### Basic AK8 jet selection
        ak8jets = [ x for x in AK8jets if abs(x.eta) < 2.5 ]
        AK8HT = sum( [ x.pt for x in ak8jets ] )
        ak8jets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################


        for ijet in ak8jets:
            getattr( self, 'AK8Jet_pt_check' ).Fill( ijet.pt )
            getattr( self, 'AK8Jet_eta_check' ).Fill( ijet.eta )
            getattr( self, 'AK8Jet_phi_check' ).Fill( ijet.phi )
            getattr( self, 'AK8Jet_mass_check' ).Fill( ijet.msoftdrop )

        if ( getattr( event, 'HLT_'+self.baseTriggers[self.year][0] )==1 ) and ( len(ak8jets)>0 ):

            getattr( self, 'cutflow' ).Fill( 0 )
            getattr( self, 'AK8Jet1Pt_baseTrigger' ).Fill( ak8jets[0].pt )
            if ( len(ak8jets)>1 ): getattr( self, 'AK8Jet2Pt_baseTrigger' ).Fill( ak8jets[1].pt )

            for i, it in enumerate(self.triggers):
                if ( getattr( event, 'HLT_'+it )==1 ):
                    getattr( self, 'cutflow' ).Fill( i+1 )
                    getattr( self, 'AK8Jet1Pt_'+it ).Fill( ak8jets[0].pt )
                    if ( len(ak8jets)>1 ): getattr( self, 'AK8Jet2Pt_'+it ).Fill( ak8jets[1].pt )
                    continue

        return True

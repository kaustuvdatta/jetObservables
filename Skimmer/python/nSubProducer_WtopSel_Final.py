#next check if b-tagging the hadronic hem b changes things compared to the below, particluarly for the top selection


import ROOT
import math, os, sys
import numpy as np
from collections import OrderedDict
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nSubProd(Module):

    def __init__(self, sysSource=[], leptonSF={}, year='2017', isMC=True, topreweight=False, onlyUnc='_nom', onlyTrees=False, evtSelection='Wtop'):
        self.writeHistFile=True
        self.leptonSFhelper = leptonSF
        print(self.leptonSFhelper)
        self.year = year
        self.isMC = isMC
        self.onlyUnc = onlyUnc
        self.onlyTrees = onlyTrees
        self.evtSelection = evtSelection
        self.recoLevel=0
        self.fakes=0
        self.miss=0
        self.genLevel=0
        self.response=0
        self.ufo=0
        self.ufoResponse=0
        self.ufoFake=0
        self.ufoMiss=0
        self.recoLevelW=0
        self.fakesW=0
        self.missW=0
        self.genLevelW=0
        self.responseW=0
        self.ufoW=0
        self.recoLeveltop=0
        self.fakestop=0
        self.misstop=0
        self.genLeveltop=0
        self.responsetop=0
        self.ufotop=0

        ### Kinematics Cuts AK8Jets ###
        self.minAK8JetPt = 170  ### this is the basic minimum, not the final
        self.maxJetAK8Eta = 1.5
        self.maxJetAK8Rap = 1.7
        ### Cuts for selections
        self.minLeadAK8JetPtW = 200.
        self.minSDMassW = 65.#60.#   ### looser to pick low bins
        self.maxSDMassW = 115. #120.#  ### looser to pick higher bins
        self.minLeadAK8JetPtTop= 350.
        self.minSDMassTop = 140.
        self.maxSDMassTop = 220.
        self.METCutWtop = 50.
        self.minLeptonicWPt = 150.


        ### Kinematics Cuts Jets ###
        self.minJetPt = 30.
        self.maxJetEta = 2.4
        self.minBDisc = 0.3040   ### L: 0.0532, M: 0.3040, T: 0.7476, for DeepJet (ie, DeepFlavB)

        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 55.
        self.maxMuonEta = 2.4
        self.minMuonMETPt = 50.

        ### Kinenatic Cuts Electrons ###
        self.minLooseElectronPt = 40.
        self.minTightElectronPt = 120.
        self.minElectronMETPt = 80.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.4]

        self.totalRecoWeight = 1.
        self.totalGenWeight = 1.
        self.btaggingWeight = 1.
        self.topreweight = topreweight 
        self.topweight = 1.
        self.leptonWeight = 1.
        self.puWeight = 1.

        self.pdfWeightUp = 1.
        self.puWeightUp = 1.
        self.isrWeightUp = 1.
        self.fsrWeightUp = 1.
        self.pdfWeightDown = 1.
        self.puWeightDown = 1.
        self.isrWeightDown = 1.
        self.fsrWeightDown = 1.
        
        self.eventCategory = -1
        self.dummy = 0

        ### Defining nsubjetiness basis
        self.maxTau = 5
        self.nSub_labels = {
                        "_tau_0p5_1": [0., 0.9, 900  ],
                        "_tau_0p5_2": [0., 0.9, 900  ],
                        "_tau_0p5_3": [0., 0.8, 800  ],
                        "_tau_0p5_4": [0., 0.8, 800  ],
                        "_tau_0p5_5": [0., 0.8, 800  ],
                        "_tau_1_1": [ 0., 0.9, 900  ],
                        "_tau_1_2": [ 0., 0.7, 700  ],
                        "_tau_1_3": [ 0., 0.5, 500  ],
                        "_tau_1_4": [ 0., 0.5, 500  ],
                        "_tau_1_5": [ 0., 0.5, 500  ],
                        "_tau_2_1": [ 0., 0.9, 900  ],
                        "_tau_2_2": [ 0., 0.5, 500  ],
                        "_tau_2_3": [ 0., 0.5, 500  ],
                        "_tau_2_4": [ 0., 0.5, 500  ],
                        "_tau_2_5": [ 0., 0.5, 500  ]
                }
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )
        self.nSub1_WTA_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6) ##### needed for genjet tau21 or tau32 # NOW USING opm for the CMS comparisons

        ### Softdrop quantities
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, self.minAK8JetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")

        ### Helpers
        self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
        self.nJet = [ 'Jet']#, 'sdJet' ]

        ### Uncertainties
        self.sysSource = ['_nom'] +[ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        if onlyUnc: self.sysSource = [ onlyUnc+i for i in [ 'Up', 'Down' ] ]
        if onlyUnc.startswith('_jes'): self.sysSource = [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        if onlyUnc.startswith('_btag'): self.sysSource = [isys for isys in sysSource if not isys.endswith('nom') ]
        print ("Sys sources:", self.sysSource)
        self.sysWeightList = ( '_pu', '_pdf', '_ps', '_isr', '_fsr', '_btag')
        #print ("Sys weight list", self.sysWeightList)
        
        self.FLAGS = [
        "isPrompt",
        "isDecayedLeptonHadron",
        "isTauDecayProduct",
        "isPromptTauDecayProduct",
        "isDirectTauDecayProduct",
        "isDirectPromptTauDecayProduct",
        "isDirectHadronDecayProduct",
        "isHardProcess",
        "fromHardProcess",
        "isHardProcessTauDecayProduct",
        "isDirectHardProcessTauDecayProduct",
        "fromHardProcessBeforeFSR",
        "isFirstCopy",
        "isLastCopy",
        "isLastCopyBeforeFSR",
        ]
        """bit-packed statusFlags interpretations.  Use `GenParticle.hasFlags` to query"""
        #############################################################################


    def hasFlags(self, topcand, *flags):
        """Check if one or more status flags are set
        Parameters
        ----------
            flags : str or list
                A list of flags that are required to be set true. If the first argument
                is a list, it is expanded and subsequent arguments ignored.
                Possible flags are enumerated in the `FLAGS` attribute
        Returns a boolean array
        """
        if not len(flags):
            raise ValueError("No flags specified")
        elif isinstance(flags[0], list):
            flags = flags[0]
        mask = 0
        for flag in flags:
            mask |= 1 << self.FLAGS.index(flag)
        return (topcand.statusFlags & mask) == mask

    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        #tauBins = 950
        #if not self.onlyTrees:
        ### Booking histograms
        if self.evtSelection.startswith('Wtop'): 
            selList = [ '_WSel', '_topSel' ]
        elif self.evtSelection.startswith('WSel') and not 'top' in self.evtSelection: 
            selList = [ '_WSel' ]
        elif self.evtSelection.startswith('topSel'): 
            selList = [ '_topSel']
        else: 
            print ("Selection %s cannot be handled with this skimmer; what funny business are you even up to?"%self.evtSelection)

        #print ("################################# Working on the %s selection"%self.evtSelection+" #################################")

        if not self.onlyUnc:
        
            self.addObject( ROOT.TH1F('cutflow_test',   ';Categories',   25, 0, 25) )
            self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
            self.addObject( ROOT.TH1F('Lepweight',   ';LepWeight',   20, 0, 2) )
            self.addObject( ROOT.TH1F('Btagweight',   ';BtagWeight',   25, 0, 2) )
            self.addObject( ROOT.TH1F('Topweight',   ';Topweight',   25, 0, 2) )
        
        #### general selection
        for isel in [ '_noSelnoWeight', '_noSel' ] + selList:
            self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
            self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   40, 0, 20) )
            self.addP4Hists( 'muons', isel )
            self.addP4Hists( 'eles', isel )
            self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   40, 0, 20) )
            self.addP4Hists( 'AK8jets', isel )
            #self.addObject( ROOT.TH1F('AK8jets_rap'+isel,   ';rapidity of AK8 jets',   100, -2., 2.) )
            self.addObject( ROOT.TH1F('nAK4jets'+isel,   ';number of AK4 jets',   40, 0, 20) )
            self.addObject( ROOT.TH1F('nAK4bjets'+isel,   ';number of AK4 jets',   40, 0, 20) )
            self.addObject( ROOT.TH1F('nAK4bhadjets'+isel,   ';number of AK4 jets',   40, 0, 20) )
            self.addObject( ROOT.TH1F('nAK4blepjets'+isel,   ';number of AK4 jets',   40, 0, 20) )
            self.addP4Hists( 'AK4jets', isel )
            self.addP4Hists( 'AK4btaggedjets', isel )
            self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )
            self.addObject( ROOT.TH1F('leptonicWMT'+isel,   ';leptonic W m_T(GeV)',   400, 0, 4000) )
            self.addObject( ROOT.TH1F('Mtt'+isel,   '; m_{t#bar{t}}(GeV)',   400, 0, 4000) )
            self.addP4Hists( 'leptonicW', isel)#,   ';Leptonic W p_T (GeV)',   200, 0, 2000) )
            self.addP4Hists( 'leptonicTop', isel)#,   ';Leptonic  p_T (GeV)',   200, 0, 2000) )
            self.addP4Hists( 'lepHemB', isel)#,   ';Leptonic  p_T (GeV)',   200, 0, 2000) )
            self.addP4Hists( 'hadHemB', isel)#,   ';Leptonic  p_T (GeV)',   200, 0, 2000) )
            self.addObject( ROOT.TH1F('HT'+isel,   ';HT (GeV)',   200, 0, 2000) )
            #self.addObject( ROOT.TH1F('pTrel'+isel,   '; p_{T,rel}', 50, 0, 200.) )
            #self.addObject( ROOT.TH1F('mindR'+isel,   '; #Delta R(j,#mu)', 50, 0, 5) )
            
            self.addObject( ROOT.TH1F('leadAK8JetMatched'+isel, ';AK8 reco jet SD mass matched'+isel+' [GeV]', 500, 0, 500) )

        if self.isMC:
            for isel in [ '_noSel' ] + selList :
                self.addObject( ROOT.TH1F('ngenleps'+isel,   ';number of gen leptons',   40, 0, 20) )
                self.addP4Hists( 'genmuons', isel )
                self.addP4Hists( 'geneles', isel )
                self.addObject( ROOT.TH1F('ngenAK8jets'+isel,   ';number of AK8 genjets',   20, 0, 20) )
                self.addP4Hists( 'AK8genjets', isel )
                self.addObject( ROOT.TH1F('ngenAK4jets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                self.addObject( ROOT.TH1F('ngenAK4bjets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                self.addObject( ROOT.TH1F('ngenAK4bhadjets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                self.addObject( ROOT.TH1F('ngenAK4blepjets'+isel,   ';number of AK4 genjets',   20, 0, 20) )
                self.addP4Hists( 'AK4genjets', isel )
                self.addP4Hists( 'AK4btaggedgenjets', isel )
                self.addObject( ROOT.TH1F('genMETPt'+isel,   ';gen MET (GeV)',   200, 0, 2000) )
                self.addObject( ROOT.TH1F('genHT'+isel,   ';genHT (GeV)',   200, 0, 2000) )
                self.addObject( ROOT.TH1F('genleptonicWMT'+isel,   ';gen leptonic W m_T(GeV)',   400, 0, 4000) )
                self.addObject( ROOT.TH1F('genMtt'+isel,   '; gen m_{t#bar{t}}(GeV)',   400, 0, 4000) )
                self.addP4Hists( 'genleptonicW', isel)#,   '; gen leptonic W p_T (GeV)',   200, 0, 2000) )
                self.addP4Hists( 'genleptonicTop', isel)#,   '; gen leptonic W p_T (GeV)',   200, 0, 2000) )
                self.addP4Hists( 'genlepHemB', isel)#,   '; gen leptonic W p_T (GeV)',   200, 0, 2000) )
                self.addP4Hists( 'genhadHemB', isel)#,   '; gen leptonic W p_T (GeV)',   200, 0, 2000) )
                #self.addObject( ROOT.TH1F('genpTrel'+isel,   '; p_{T,rel}', 50, 0, 200.) )
                #self.addObject( ROOT.TH1F('genmindR'+isel,   '; #Delta R(j,#mu)', 50, 0, 5) )
                

        for isel in selList:
            for iJ in self.nJet:
                #self.addP4Hists( 'uforeco'+iJ, isel )
                #self.addP4Hists( 'ufogen'+iJ, isel )
                
                for itype in ( [ 'gen',  'accepgen', 'missgen', 'reco', 'fakereco', 'truereco' ] if self.isMC else [ 'reco' ] ): #'missgen',

                    if itype.endswith(('reco','gen')):

                        for sysUnc in self.sysSource:
                            self.addP4Hists( itype+iJ, sysUnc+isel )

                            binning = np.array([(i/2000) for i in np.arange(0*2000, 2.*2000)])
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{21}', len(binning)-1, binning) )
                            self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+sysUnc+isel, ';AK8 '+itype+' jet #tau_{32}', len(binning)-1, binning) )


                            if itype.startswith('truereco') and self.isMC:
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau21'+sysUnc+isel, ';AK8 accepgen jet #tau_{21};AK8 truereco jet #tau_{21}', len(binning)-1, binning, len(binning)-1, binning) )
                                self.addObject( ROOT.TH2F('resp'+iJ+'_tau32'+sysUnc+isel, ';AK8 accepgen jet #tau_{32};AK8 truereco jet #tau_{32}', len(binning)-1, binning, len(binning)-1, binning) )
                                
                                if sysUnc.endswith('nom'):
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_pt'+isel, ';AK8 reco/gen jet pt', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_sdmass'+isel, ';AK8 reco/gen jet SD mass', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau21'+isel, ';AK8 reco/gen jet #tau_{21}', 500, 0, 5) )
                                    self.addObject( ROOT.TH1F('resol'+iJ+'_tau32'+isel, ';AK8 reco/gen jet #tau_{32}', 500, 0, 5) )

                            for x, y in self.nSub_labels.items():
                                binning = np.array([(i/1000) for i in np.arange(y[0]*1000, y[1]*1000)])
                                self.addObject( ROOT.TH1F(itype+iJ+x+sysUnc+isel, ';AK8 '+itype+' jet #tau', len(binning)-1, binning ) )
                                if itype.startswith('reco') and self.isMC:
                                    self.addObject( ROOT.TH2F('resp'+iJ+x+sysUnc+isel, ';AK8 gen jet '+x+';AK8 reco jet '+x, len(binning)-1, binning, len(binning)-1, binning ) )
                                    if sysUnc.endswith('nom'): self.addObject( ROOT.TH1F('resol'+iJ+x+isel, ';AK8 reco/gen jet '+x, 500, 0, 5) )

                    else:
                        #print ('Creating gen objects:', itype, iJ, isel)
                        self.addP4Hists( itype+iJ, isel )
                        binning = np.array([(i/2000) for i in np.arange(0*2000, 2.*2000)])
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau21'+isel, ';AK8 '+itype+' jet #tau_{21}', len(binning)-1, binning) )
                        self.addObject( ROOT.TH1F(itype+iJ+'_tau32'+isel, ';AK8 '+itype+' jet #tau_{32}', len(binning)-1, binning) )

                        for x, y in self.nSub_labels.items():
                            binning = np.array([(i/1000) for i in np.arange(y[0]*1000, y[1]*1000)])
                            self.addObject( ROOT.TH1F(itype+iJ+x+isel, ';AK8 '+itype+' jet '+x, len(binning)-1, binning ) )


    #############################################################################
    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_rap'+t, s+';y', 100, -2.5, 2.5 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14259) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 500) )


    #############################################################################
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #print (wrappedOutputTree, outputFile, inputFile, inputTree)
        #if self.onlyTrees:
        self.out = wrappedOutputTree
        #self.out.branch('triggerWeight',  "F")
        self.out.branch('eventCategory',  "I")
        self.out.branch( 'totalRecoWeight', "F" )
        self.out.branch( 'totalGenWeight', "F" )
        self.out.branch('leptonWeight',  "F")
        self.out.branch('btagWeight',  "F")
        self.out.branch('puWeight',  "F")
        self.out.branch('top_pTreweight_SF',  "F")
        
        if self.onlyUnc:
            if 'pdf' in self.onlyUnc:
                self.out.branch( 'pdfWeightUp', "F" )
                self.out.branch( 'pdfWeightDown', "F" )

            if 'isr' in self.onlyUnc:
                self.out.branch( 'isrWeightUp', "F" )
                self.out.branch( 'isrWeightDown', "F" )

            if 'fsr' in self.onlyUnc:
                self.out.branch( 'fsrWeightUp', "F" )
                self.out.branch( 'fsrWeightDown', "F" )

            if 'puWeight' in self.onlyUnc:
                self.out.branch( 'puWeightUp', "F" )
                self.out.branch( 'puWeightDown', "F" )

        tmplist = [ 'selRecoJets'+sys for sys in self.sysSource ]
        #print ('selrecojetlist:', tmplist)

        if self.isMC: 
            X = ['trueRecoJets'+sys for sys in self.sysSource]
            

            #print ("appended MC reco/genjets to nano output file")
        
            if not self.onlyUnc: 
                tmplist.append(X[0])
                tmplist.append( 'selGenJets_nom' )
                tmplist.append( 'accepGenJets_nom' )
            else:
                tmpGenList = [ 'selGenJets'+sys for sys in self.sysSource ]
                tmpAccepGenList = [ 'accepGenJets'+sys for sys in self.sysSource ] 
                for x in X:
                    tmplist.append(x)
                for x in tmpGenList:
                    tmplist.append(x)
                for x in tmpAccepGenList:
                    tmplist.append(x)
            
        print ("Stored tree branches:", tmplist)
        for iJ in tmplist:
            #print (iJ)
            self.out.branch('n'+iJ,  'I')  ### dummy for nanoAOD Tools
            self.out.branch(iJ+'_pt',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_eta',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_phi',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_rap',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_mass',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_msoftdrop',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_Tau21',  'F', lenVar='n'+iJ)
            self.out.branch(iJ+'_Tau32',  'F', lenVar='n'+iJ)
            #self.out.branch(iJ+'_SD_Tau21',  'F', lenVar='n'+iJ)
            #self.out.branch(iJ+'_SD_Tau32',  'F', lenVar='n'+iJ)

            for x in self.nSub_labels:
                self.out.branch(iJ+x, 'F', lenVar='n'+iJ )
                #self.out.branch(iJ+'SD'+x, 'F', lenVar='n'+iJ )
        pass

    #############################################################################
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

        #if not self.onlyTrees:
        if self.isMC and not self.onlyUnc:

            self.genLevel = self.response+self.miss
            self.genLevelW = self.responseW+self.missW
            self.genLeveltop = self.responsetop+self.misstop

            getattr( self, 'cutflow_test' ).SetBinContent( 1, self.recoLevel )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 1, "allReco" )
            getattr( self, 'cutflow_test' ).SetBinContent( 2, self.response )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 2, "allTrueRecoAccepGen" )
            getattr( self, 'cutflow_test' ).SetBinContent( 3, self.fakes )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 3, "allFakes" )
            getattr( self, 'cutflow_test' ).SetBinContent( 4, self.miss )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 4, "allMiss" )
            getattr( self, 'cutflow_test' ).SetBinContent( 5, self.genLevel )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 5, "allGen" )

            getattr( self, 'cutflow_test' ).SetBinContent( 7, self.recoLevelW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 7, "WReco" )
            getattr( self, 'cutflow_test' ).SetBinContent( 8, self.responseW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 8, "WtrueRecoAccepGen" )
            getattr( self, 'cutflow_test' ).SetBinContent( 9, self.fakesW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 9, "WFakes" )
            getattr( self, 'cutflow_test' ).SetBinContent( 10, self.missW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 10, "WMiss" )
            getattr( self, 'cutflow_test' ).SetBinContent( 11, self.genLevelW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 11, "WGen" )

            getattr( self, 'cutflow_test' ).SetBinContent( 13, self.recoLeveltop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 13, "topReco" )
            getattr( self, 'cutflow_test' ).SetBinContent( 14, self.responsetop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 14, "toptrueRecoAccepGen" )
            getattr( self, 'cutflow_test' ).SetBinContent( 15, self.fakestop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 15, "topFakes" )
            getattr( self, 'cutflow_test' ).SetBinContent( 16, self.misstop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 16, "topMiss" )
            getattr( self, 'cutflow_test' ).SetBinContent( 17, self.genLeveltop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 17, "topGen" )

            getattr( self, 'cutflow_test' ).SetBinContent( 19, self.ufo )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 19, "allUfo" )
            getattr( self, 'cutflow_test' ).SetBinContent( 20, self.ufoResponse )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 20, "ufoResponse" )
            getattr( self, 'cutflow_test' ).SetBinContent( 21, self.ufoFake )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 21, "ufoFake" )
            getattr( self, 'cutflow_test' ).SetBinContent( 22, self.ufoMiss )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 22, "ufoMiss" )
            getattr( self, 'cutflow_test' ).SetBinContent( 23, self.ufoW )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 23, "ufoW" )
            getattr( self, 'cutflow_test' ).SetBinContent( 24, self.ufotop )
            getattr( self, 'cutflow_test' ).GetXaxis().SetBinLabel( 24, "ufotop" )

    #############################################################################
    def leptonSF(self, lepton, leptonP4 ):

        leptonP4eta = abs(leptonP4.eta)
        leptonP = ROOT.TMath.Sqrt(leptonP4.p4().Px()**2 + leptonP4.p4().Py()**2 + leptonP4.p4().Pz()**2)

        SFFileTrigger = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['Trigger'][0] )
        histoSFTrigger = SFFileTrigger.Get( self.leptonSFhelper[lepton]['Trigger'][1] )
        SFTrigger = histoSFTrigger.GetBinContent( histoSFTrigger.GetXaxis().FindBin( leptonP4eta ), histoSFTrigger.GetYaxis().FindBin( leptonP4.pt ) )

        SFFileID = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ID'][0] )
        histoSFID = SFFileID.Get( self.leptonSFhelper[lepton]['ID'][1] )
        histoSFID_X = histoSFID.GetXaxis().FindBin( leptonP4eta)
        histoSFID_Y = histoSFID.GetYaxis().FindBin( leptonP4.pt )
        SFID = histoSFID.GetBinContent( histoSFID_X, histoSFID_Y )
        SFID = SFID if SFID>0 else 1

        SFFileISO = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['ISO'][0] )
        histoSFISO = SFFileISO.Get( self.leptonSFhelper[lepton]['ISO'][1] )
        histoSFISO_X = histoSFISO.GetXaxis().FindBin( leptonP4eta )
        histoSFISO_Y = histoSFISO.GetYaxis().FindBin( leptonP4.pt )
        SFISO = histoSFISO.GetBinContent( histoSFISO_X, histoSFISO_Y )
        SFISO = SFISO if SFISO>0 else 1
        
        SFFileRecoEff = ROOT.TFile( os.environ['CMSSW_BASE']+"/src/jetObservables/Skimmer/data/leptonSF/"+self.leptonSFhelper[lepton]['RecoEff'][0] )
        histoSFRecoEff = SFFileRecoEff.Get( self.leptonSFhelper[lepton]['RecoEff'][1] )
        histoSFRecoEff_X = histoSFRecoEff.GetXaxis().FindBin( leptonP4eta )
        histoSFRecoEff_Y = histoSFRecoEff.GetYaxis().FindBin( leptonP )
        SFRecoEff = histoSFRecoEff.GetBinContent( histoSFRecoEff_X, histoSFRecoEff_Y )
        SFRecoEff = SFRecoEff if SFRecoEff>0 else 1

        #print (SFTrigger * SFID * SFISO), SFTrigger , SFID , SFISO, leptonP4.pt, leptonP4.eta
        return [SFTrigger , SFID , SFISO, SFRecoEff]

     #############################################################################



    def mindRandIndex(self, recoAK4Coll, muon):
        mindR=999
        c=0
        for x in recoAK4Coll:
            if x.p4().DeltaR(muon)<mindR:
                mindR=x.p4().DeltaR(muon)
                c+=1
        return mindR, c-1

    def getpTRel(self, recoAK4, muon): 

        p_mu = muon.Vect()#ROOT.TVector3(muon.Px(),muon.Py(),muon.Pz())
        p_j = recoAK4.Vect()#ROOT.TVector3(recoAK4.Px(),recoAK4.Py(),recoAK4.Pz())
        
        pTrel2 = p_mu.Dot(p_mu) - (p_mu.Dot(p_j)/p_j.Mag())**2.
        
        if pTrel2<0.: return 0.
        else: return ROOT.TMath.Sqrt(pTrel2)

    #############################################################################
    '''
    def getBTagWeight(self, nBTagged=0, jet_SFs=[0]): # Now redundant to allow for more than just two btagged jets per event
        bTagWeight=0
        if len(jet_SFs)>2 or nBTagged>2:
            print "Error, only leading and subleading AK4 jets are considered: # of btagged jets cannot exceed 2"
        elif nBTagged>len(jet_SFs):
            print "Number of b-tagged jets cannot be greater than number of them for which SFs are provided!"
            return 0
        elif nBTagged==0 and len(jet_SFs)==0: return 1

        if len(jet_SFs)==1:
            SF = jet_SFs[0]

            for i in range(0,2):
                if i!=nBTagged: continue
                bTagWeight+=pow(SF,i)*pow(1-SF,1-i)

        elif len(jet_SFs)==2:
            SF1, SF2 = jet_SFs[0], jet_SFs[1]
            for i in range(0,2):
                for j in range(0,2):
                    if (i+j)!=nBTagged: continue
                    bTagWeight+=pow(SF1,i)*pow(1-SF1,1-i)*pow(SF2,j)*pow(1-SF2,1-j)

        return bTagWeight

    '''
    #############################################################################
    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''
        self.isMC = event.run == 1

        if not self.isMC:
            passGenSel=False
            iGenSel=None
        else:
            passGenSel, iGenSel, selGenMuons, selGenElectrons, selGenAK4bjets, selGenJets, selGenMET = self.genSelection(event)
        passRecoSel, iRecoSel, selRecoMuons, selRecoElectrons, selRecoAK4bjets, selRecoJets, selRecoMET, selbtagweights = self.recoSelection( event )
        #print (event.genWeight)
        if not self.isMC and not passRecoSel['_nom']: return False
        
        
        #if self.isMC and not passGenSel and not passRecoSel['_nom']: return False

        for sys in self.sysSource:

            if sys.startswith(self.sysWeightList):
    
                selRecoJets[sys] = selRecoJets['_nom']
                iRecoSel[sys] = iRecoSel['_nom']
                passRecoSel[sys] = passRecoSel['_nom']

                if sys.startswith('_pdfWeight'):
                    arr =  getattr( event, 'LHEPdfWeight' )
                    size = arr.GetSize() # force loading of array contents for this entry
                    buf = arr.GetAddress() # still not sure this is always equivalent to `&arr[0]`
                    buf.SetSize(size)
                    npLHEPdfWeight = np.frombuffer(buf, count=size, dtype='float32')
                    pdfWeight = np.std(npLHEPdfWeight)

                    if sys.endswith('Up'):
                        self.pdfWeightUp = 1 + pdfWeight
                        self.totalGenWeight = event.genWeight * self.pdfWeightUp
                        WEIGHT = self.totalRecoWeight * self.pdfWeightUp 
                    else:
                        self.pdfWeightDown = 1 - pdfWeight
                        self.totalGenWeight = event.genWeight * self.pdfWeightDown
                        WEIGHT = self.totalRecoWeight * self.pdfWeightDown 

                elif sys.endswith('isrWeightDown'):
                    self.isrWeightDown = getattr( event, 'PSWeight' )[0]
                    self.totalGenWeight = event.genWeight * self.isrWeightDown
                    WEIGHT = self.totalRecoWeight * self.isrWeightDown 

                elif sys.endswith('fsrWeightDown'):
                    self.fsrWeightDown = getattr( event, 'PSWeight' )[1]
                    self.totalGenWeight = event.genWeight * self.fsrWeightDown
                    WEIGHT = self.totalRecoWeight * self.fsrWeightDown 

                elif sys.endswith('isrWeightUp'):
                    self.isrWeightUp = getattr( event, 'PSWeight' )[2]
                    self.totalGenWeight = event.genWeight * self.isrWeightUp
                    WEIGHT = self.totalRecoWeight * self.isrWeightUp 
                elif sys.endswith('fsrWeightUp'):
                    self.fsrWeightUp = getattr( event, 'PSWeight' )[3]
                    self.totalGenWeight = event.genWeight * self.fsrWeightUp
                    WEIGHT = self.totalRecoWeight * self.fsrWeightUp 

                elif sys.endswith('puWeightUp'):
                    self.puWeightUp = event.puWeightUp
                    self.totalGenWeight = event.genWeight 
                    WEIGHT = self.totalRecoWeight * self.puWeightUp / self.puWeight #divide out puweight nominal in this case, since that's multiplied into the weight in recoSelection() for the nominal case and other non-PU variations
                elif sys.endswith('puWeightDown'):
                    self.puWeightDown = event.puWeightDown
                    self.totalGenWeight = event.genWeight 
                    WEIGHT = self.totalRecoWeight * self.puWeightDown / self.puWeight
                elif sys.startswith('_btag'):
                    self.totalGenWeight = event.genWeight
                    WEIGHT = self.totalRecoWeight # the up down modulation of btagweights is controlled in the WtopSel function
            else: 

                WEIGHT =  self.totalRecoWeight*(selbtagweights[sys]/selbtagweights['_nom']) 
                # the funky btag weight ratio is to handle small variations if same # of ak4b's aren't tagged in case of AK8 jet sys variations
                # this is because of how we decide to id the hadronic hemispheric AK4s via a deltaR selection vs the leading AK8 leading to this dependency
                if self.isMC: self.totalGenWeight = event.genWeight

            
            genJet = OrderedDict()
            recoJet = OrderedDict()
            tmpRecoJets = OrderedDict()
            
            if self.isMC and not passGenSel and not passRecoSel[sys]: return False

            if self.isMC: 
                deltaRmatch = False



            if passRecoSel[sys]:  #### Detector level dist.
            
                #print ("All reco weights: self.totalRecoWeight, WEIGHT, event.puWeight , event.genWeight , self.btaggingWeight, self.leptonWeight , self.topreweight")
                #print ("All reco weights:", self.totalRecoWeight, WEIGHT, event.puWeight , event.genWeight , self.btaggingWeight, self.leptonWeight , self.topreweight) 

                tmpRecoJets[sys] = {}

                if sys.endswith('nom'):
                    tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )
                else:
                    if '_nom' in tmpRecoJets.keys():
                        tmpRecoJets[sys][0] = tmpRecoJets['_nom'][0]
                    else:
                        tmpRecoJets[sys][0] = self.createNsubBasis( selRecoJets[sys][0], event, 'PFCands' )

                recoJet['Jet'] = tmpRecoJets[sys][0] #self.createNsubBasis( tmpRecoJets[sys][0], event, 'PFCands' )#tmpRecoJet[sys][0]
                self.recoLevel = self.recoLevel+1       #### counting ALL the recoLevel
                if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                    self.recoLevelW = self.recoLevelW+1
                elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                    self.recoLeveltop = self.recoLeveltop+1
                else: print('Weird reco', iRecoSel[sys])

                self.eventCategory=1

                
                
                for iRJ,ireco in recoJet.items():

                    if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                        minLeadAK8JetPt = self.minLeadAK8JetPtW
                    elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                        minLeadAK8JetPt = self.minLeadAK8JetPtTop

                    if ( getattr(ireco['jet'], 'pt'+( '_nom' if sys.startswith(self.sysWeightList) else sys ) ) < minLeadAK8JetPt ): continue

                    if sys.startswith(self.sysWeightList):
                        
                        getattr( self, 'reco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                        if iRecoSel[sys].startswith('_W'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                        elif iRecoSel[sys].startswith('_top'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                        
                    else:
                        getattr( self, 'reco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                        if iRecoSel[sys].startswith('_W'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                        elif iRecoSel[sys].startswith('_top'): getattr( self, 'reco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop'+sys )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT) #?? regarding the puppi corr division

                    getattr( self, 'reco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_rap'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'rapidity'), WEIGHT )
                    getattr( self, 'reco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                    getattr( self, 'reco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                    
                    for tauN in range(1, self.maxTau+1):
                        getattr( self, 'reco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                        getattr( self, 'reco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                        
                self.fillBranches( 'selRecoJets'+sys, recoJet )
                
                if self.isMC:
                    if not passGenSel:
                        # if evt does not pass gensel but passes recosel[sys] add to fakes
                        self.fakes = self.fakes+1

                        if iRecoSel[sys].startswith('_W'): self.fakesW = self.fakesW+1
                        elif iRecoSel[sys].startswith('_top'): self.fakestop = self.fakestop+1
                        else: self.ufoFake = self.ufoFake+1
                        
                        self.eventCategory = 2
                        
                        for iRJ,ireco in recoJet.items():
                            if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                                minLeadAK8JetPt = self.minLeadAK8JetPtW
                            elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                                minLeadAK8JetPt = self.minLeadAK8JetPtTop

                            if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < minLeadAK8JetPt ): continue

                            if sys.startswith(self.sysWeightList):
                                getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                if iRecoSel[sys].startswith('_W'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                elif iRecoSel[sys].startswith('_top'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                            else:
                                getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_rap'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'rapidity'), WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                            getattr( self, 'fakereco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                                
                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                    else:
                        #print (sys, self.totalRecoWeight, self.totalGenWeight, WEIGHT, self.responseW, self.responsetop)

                        genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )

                        for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):
                            if iGenSel.startswith('_W'):           #### counting recoLevelW
                                minLeadAK8JetPt = self.minLeadAK8JetPtW
                            elif iGenSel.startswith('_top'):       #### counting recoLeveltop
                                minLeadAK8JetPt = self.minLeadAK8JetPtTop
                            if getattr(igen['jet'], 'pt') < minLeadAK8JetPt: 
                                continue
                            #if sys.endswith('nom'):
                            getattr( self, 'gen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_rap'+sys+iGenSel ).Fill( igen['jet'].rapidity, self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                          
                            for tauN in range(1, self.maxTau+1):
                                getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                             
                        self.fillBranches( 'selGenJets'+sys, genJet )

                    
                        if ( genJet['Jet']['jet'].p4().DeltaR( recoJet['Jet']['jet'].p4()) < 0.8 ): # we should tighten/loosen/remove this?
                            deltaRmatch = True
                            self.response = self.response+1
                            if iGenSel.startswith('_W'): self.responseW = self.responseW+1
                            elif iGenSel.startswith('_top'): self.responsetop = self.responsetop+1
                            else: print('Weird response', iRecoSel[sys])
                            self.eventCategory = 4

                        #if deltaR matched, fill hists for accepted gen jet and truereco jets, as well as the response matrices, of course!

                        if  deltaRmatch:

                            for (iGJ,igen), (iRJ,ireco) in zip( genJet.items(), recoJet.items() ):

                                if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                                    minLeadAK8JetPt = self.minLeadAK8JetPtW
                                elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                                    minLeadAK8JetPt = self.minLeadAK8JetPtTop

                                if getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < minLeadAK8JetPt or getattr(igen['jet'], 'pt') < minLeadAK8JetPt: 
                                    continue 

                                #if sys.endswith('nom'):
                                #    print ("yayx4!", iGenSel, iRecoSel[sys] , sys, passGenSel, iGJ, igen['jet'], iRJ, ireco )

                                getattr( self, 'accepgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_rap'+sys+iGenSel ).Fill( igen['jet'].rapidity, self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                                getattr( self, 'accepgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                                
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'accepgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'accepgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                    
                                

                                    
                                if sys.startswith(self.sysWeightList):
                                    getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                    if iGenSel.startswith('_W'): getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                    elif iGenSel.startswith('_top'): getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                                    
                                else:
                                    getattr( self, 'truereco'+iRJ+'_pt'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_mass'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                
                                getattr( self, 'truereco'+iRJ+'_eta'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_phi'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_rap'+sys+iGenSel ).Fill( getattr(ireco['jet'], 'rapidity'), WEIGHT )

                                getattr( self, 'truereco'+iRJ+'_tau21'+sys+iGenSel ).Fill( ireco['tau21'], WEIGHT )
                                getattr( self, 'truereco'+iRJ+'_tau32'+sys+iGenSel ).Fill( ireco['tau32'], WEIGHT )
                                
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'truereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'truereco'+iRJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                

                                #### Filling response matrices and resolution plots
                                if sys.endswith('nom'):
                                    getattr( self, 'resol'+iGJ+'_pt'+iGenSel ).Fill( getattr(ireco['jet'], 'pt')/getattr(igen['jet'], 'pt'))#, WEIGHT )
                                    getattr( self, 'resol'+iGJ+'_sdmass'+iGenSel ).Fill( getattr(ireco['jet'], 'msoftdrop')/getattr(igen['jet'], 'mass'))#, WEIGHT )
                                    getattr( self, 'resol'+iGJ+'_tau21'+iGenSel ).Fill( ireco['tau21']/igen['tau21'])#, WEIGHT )
                                    getattr( self, 'resol'+iGJ+'_tau32'+iGenSel ).Fill( ireco['tau32']/igen['tau32'])#, WEIGHT )
                                    
                                getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], ireco['tau21'], WEIGHT )
                                getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], ireco['tau32'], WEIGHT )
                                getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight - WEIGHT )
                                getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight - WEIGHT )
                                
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], ireco['2'+str(tauN)], WEIGHT )
                                    #counter weights below
                                    getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight - WEIGHT  )
                                    getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight - WEIGHT  )
                                    getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight - WEIGHT  )
                                    
                                    if sys.endswith('nom'):
                                        getattr( self, 'resol'+iGJ+'_tau_0p5_'+str(tauN)+iGenSel ).Fill( ( ireco['0p5'+str(tauN)]/igen['0p5'+str(tauN)] if igen['0p5'+str(tauN)]!=0 else -999 ))#, WEIGHT )
                                        getattr( self, 'resol'+iGJ+'_tau_1_'+str(tauN)+iGenSel ).Fill( ( ireco['1'+str(tauN)]/igen['1'+str(tauN)] if igen['1'+str(tauN)]!=0 else -999 ))#, WEIGHT )
                                        getattr( self, 'resol'+iGJ+'_tau_2_'+str(tauN)+iGenSel ).Fill( ( ireco['2'+str(tauN)]/igen['2'+str(tauN)] if igen['2'+str(tauN)]!=0 else -999 ))#, WEIGHT )
                            self.fillBranches( 'accepGenJets'+sys, genJet )    
                            self.fillBranches( 'trueRecoJets'+sys, recoJet )    


                        #if not deltaR matched, fill hists for misreconstructed gen jet and fakereco jets
                        else:# and not self.onlyUnc.startswith(self.sysWeightList):


                            for iGJ,igen in genJet.items():
                                
                                if iGenSel.startswith('_W'):           #### counting recoLevelW
                                    minLeadAK8JetPt = self.minLeadAK8JetPtW
                                elif iGenSel.startswith('_top'):       #### counting recoLeveltop
                                    minLeadAK8JetPt = self.minLeadAK8JetPtTop
                                else: print ("something funky at missgen filling, found event from %s sel!?"%iGenSel )

                                if getattr(igen['jet'], 'pt') < minLeadAK8JetPt: continue

                                #if sys.endswith('nom'):

                                getattr( self, 'missgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_rap'+sys+iGenSel ).Fill( igen['jet'].rapidity, self.totalGenWeight )

                                getattr( self, 'missgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                                getattr( self, 'missgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                                
                                getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                                getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )


                                for tauN in range(1, self.maxTau+1):
                                    #if sys.endswith('nom'):
                                    getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                    getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight  )
                                    getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight  )
                                    getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight  )

                                #print (self.totalGenWeight - WEIGHT, self.totalGenWeight, WEIGHT, self.responseW, self.responsetop)
                            self.fakes = self.fakes+1

                            if iRecoSel[sys].startswith('_W'): self.fakesW = self.fakesW+1
                            elif iRecoSel[sys].startswith('_top'): self.fakestop = self.fakestop+1
                            else: self.ufoFake = self.ufoFake+1
                            
                            self.eventCategory = 2
                            
                            for iRJ,ireco in recoJet.items():
                                if iRecoSel[sys].startswith('_W'):           #### counting recoLevelW
                                    minLeadAK8JetPt = self.minLeadAK8JetPtW
                                elif iRecoSel[sys].startswith('_top'):       #### counting recoLeveltop
                                    minLeadAK8JetPt = self.minLeadAK8JetPtTop

                                if ( getattr(ireco['jet'], 'pt'+(sys if not sys.startswith(self.sysWeightList) else '_nom')) < minLeadAK8JetPt ): continue

                                if sys.startswith(self.sysWeightList):
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt_nom' ), WEIGHT )
                                    if iRecoSel[sys].startswith('_W'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop_nom' ), WEIGHT )
                                    elif iRecoSel[sys].startswith('_top'): getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( (getattr(ireco['jet'], 'msoftdrop_nom' )/getattr(ireco['jet'], 'msoftdrop_corr_PUPPI' )), WEIGHT)
                                else:
                                    getattr( self, 'fakereco'+iRJ+'_pt'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'pt'+sys ), WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_mass'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'msoftdrop'+sys ), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_eta'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'eta'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_phi'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'phi'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_rap'+sys+iRecoSel[sys] ).Fill( getattr(ireco['jet'], 'rapidity'), WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau21'+sys+iRecoSel[sys] ).Fill( ireco['tau21'], WEIGHT )
                                getattr( self, 'fakereco'+iRJ+'_tau32'+sys+iRecoSel[sys] ).Fill( ireco['tau32'], WEIGHT )
                                    
                                for tauN in range(1, self.maxTau+1):
                                    getattr( self, 'fakereco'+iRJ+'_tau_0p5_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['0p5'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_1_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['1'+str(tauN)], WEIGHT )
                                    getattr( self, 'fakereco'+iRJ+'_tau_2_'+str(tauN)+sys+iRecoSel[sys] ).Fill( ireco['2'+str(tauN)], WEIGHT )
                                    
            if self.isMC and not passRecoSel[sys]: # handle misreconstructed events here, ie, the ones that do not pass the reco sel but do pass the gen sel
                if passGenSel:               
                    genJet['Jet'] = self.createNsubBasis( selGenJets[0], event, 'GenCands' )

                    for (iGJ,igen) in genJet.items():
                        if iGenSel.startswith('_W'):           #### counting recoLevelW
                            minLeadAK8JetPt = self.minLeadAK8JetPtW
                        elif iGenSel.startswith('_top'):       #### counting recoLeveltop
                            minLeadAK8JetPt = self.minLeadAK8JetPtTop
                        if getattr(igen['jet'], 'pt') < minLeadAK8JetPt: 
                            continue
                        #if sys.endswith('nom'):
                        getattr( self, 'gen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_rap'+sys+iGenSel ).Fill( igen['jet'].rapidity, self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                        getattr( self, 'gen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                            
                        getattr( self, 'missgen'+iGJ+'_pt'+sys+iGenSel ).Fill( igen['jet'].pt, self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_eta'+sys+iGenSel ).Fill( igen['jet'].eta, self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_mass'+sys+iGenSel ).Fill( igen['jet'].mass, self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_phi'+sys+iGenSel ).Fill( igen['jet'].phi, self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_rap'+sys+iGenSel ).Fill( igen['jet'].rapidity, self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], self.totalGenWeight )
                        getattr( self, 'missgen'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], self.totalGenWeight )
                        
                        getattr( self, 'resp'+iGJ+'_tau21'+sys+iGenSel ).Fill( igen['tau21'], -1., self.totalGenWeight )
                        getattr( self, 'resp'+iGJ+'_tau32'+sys+iGenSel ).Fill( igen['tau32'], -1., self.totalGenWeight )


                        for tauN in range(1, self.maxTau+1):
                            #if sys.endswith('nom'):
                            getattr( self, 'gen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'gen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )
                                         
                            getattr( self, 'missgen'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'missgen'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], self.totalGenWeight )
                            getattr( self, 'missgen'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], self.totalGenWeight )

                            getattr( self, 'resp'+iGJ+'_tau_0p5_'+str(tauN)+sys+iGenSel ).Fill( igen['0p5'+str(tauN)], -1., self.totalGenWeight  )
                            getattr( self, 'resp'+iGJ+'_tau_1_'+str(tauN)+sys+iGenSel ).Fill( igen['1'+str(tauN)], -1., self.totalGenWeight  )
                            getattr( self, 'resp'+iGJ+'_tau_2_'+str(tauN)+sys+iGenSel ).Fill( igen['2'+str(tauN)], -1., self.totalGenWeight  )
                    
                    self.fillBranches( 'selGenJets'+sys, genJet )

            

    
        self.out.fillBranch( 'eventCategory', self.eventCategory )
        self.out.fillBranch( 'totalRecoWeight', self.totalRecoWeight if self.isMC else 1)
        self.out.fillBranch( 'totalGenWeight', self.totalGenWeight if self.isMC else -1)
        self.out.fillBranch( 'leptonWeight', self.leptonWeight if self.isMC else 1)
        self.out.fillBranch( 'btagWeight', self.btaggingWeight if self.isMC else 1)
        if not ('puWeight' in self.onlyUnc):
           self.out.fillBranch( 'puWeight', self.puWeight if self.isMC else 1)
        self.out.fillBranch( 'top_pTreweight_SF', self.topweight if self.isMC else 1)

        if self.onlyUnc:
            if 'pdf' in self.onlyUnc:
                self.out.fillBranch( 'pdfWeightUp', self.pdfWeightUp if self.isMC else 1)
                self.out.fillBranch( 'pdfWeightDown', self.pdfWeightDown if self.isMC else 1)
            if 'isr' in self.onlyUnc:
                self.out.fillBranch( 'isrWeightUp', self.isrWeightUp if self.isMC else 1)
                self.out.fillBranch( 'isrWeightDown', self.isrWeightDown if self.isMC else 1)
            if 'fsr' in self.onlyUnc:
                self.out.fillBranch( 'fsrWeightUp', self.fsrWeightUp if self.isMC else 1)
                self.out.fillBranch( 'fsrWeightDown', self.fsrWeightDown if self.isMC else 1)
            if 'puWeight' in self.onlyUnc:
                self.out.fillBranch( 'puWeightUp', event.puWeightUp if self.isMC else 1 )
                self.out.fillBranch( 'puWeightDown', event.puWeightDown if self.isMC else 1 )
            
        return True


    #############################################################################
    def recoSelection( self, event, sysUnc=[]  ):
        '''Analyzing reco information'''

        AK8jets = list(Collection(event, 'FatJet' ))

        for ijets in AK8jets: ijets.rapidity = self.etaToRapidity(ijets)



        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))
        met = Object(event, 'MET')
        if self.isMC: genParticles = Collection(event, 'GenPart')
        
        
        ########### Lepton selection ###############
        recoElectrons  = [x for x in electrons if x.pt > self.minLooseElectronPt and x.cutBased >= 2 and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #only loose selection for e
        recoMuons = [x for x in muons if x.pt > self.minTightMuonPt and abs(x.p4().Eta()) < self.maxMuonEta and x.tightId and abs(x.dxy)<0.2 and abs(x.dz)<0.5 and x.isGlobal and x.highPtId and x.tkRelIso<0.3] #  applying tight selection on muons already here since we only veto for loose muons and loose electrons in event

        recoElectrons.sort(key=lambda x:x.pt, reverse=True)
        recoMuons.sort(key=lambda x:x.pt,reverse=True)

        nleptons = len(recoMuons)+len(recoElectrons)

        ############################################

        ################### MET  #######################
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiM(met.pt, 0., met.phi, 0.)#met.sumEt)
        ################################################


        ############ Basic AK4 selection ###############
        recoAK4jets = [ x for x in jets if x.pt > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta and (x.jetId >= 2)]# and x.btagDeepFlavB > self.minBDisc]
        recoAK4jets.sort(key=lambda x:x.pt,reverse=True)
        ################################################

        recoAK8jets = {}
        btagweights = {}
        passSel = {}
        iSel = {}
        flagTopWt=False
        recoAK4bjets = []
        recoAK4_lepTop_jets = []
        recoAK4_hadrTop_jets = [] 
        recottbarp4 = []
        recoleptW = []
        recoleptonicTop = []
        for sys in self.sysSource:
             
            if sys.startswith(self.sysWeightList): sys = '_nom'

            ################### Basic AK8 jet selection ############################
            recoAK8jets[sys] = [ x for x in AK8jets if getattr( x, 'pt'+sys ) > self.minAK8JetPt and abs(x.rapidity) < self.maxJetAK8Rap and (x.jetId >= 2) ]
            recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)
            AK8HT = sum( [ getattr( x, 'pt'+sys ) for x in recoAK8jets[sys] ] )
            ###################################################################

           
            ########################## Selections on jets ############################## 
            if len(recoMuons)!=0: 

                # isolate leptonic hemisphere AK4's via deltaR(mu,AK4jet), require a medium b-tag on AK4's in this hemisphere
                recoAK4_lepTop_jets = [x for x in recoAK4jets if (x.p4().DeltaR(recoMuons[0].p4())<1.6 and x.p4().DeltaR( recoMuons[0].p4() )>0.4)]# and x.btagDeepFlavB > self.minBDisc] 
                recoAK4_lepTop_jets.sort(key=lambda x:x.pt,reverse=True) 
                #######################################################

                
                # leptonic W definition
                recoleptW = [MET+recoMuons[0].p4()]

                #######################################################

                # leptonic top definition
                if len(recoAK4_lepTop_jets)!=0: 
                    #leptonic top reco'd using leading leptonic hemisphere AK4, and leptonically decaying W 
                    recoleptonicTop = [x+recoAK4_lepTop_jets[0].p4() for x in recoleptW] 


                    # Ensure separation of the AK8s from the leptonic hemisphere by deltaPhi(mu,AK8) and deltaR(AK8,lept hem AK4)
                    recoAK8jets[sys] = [ x for x in recoAK8jets[sys] if abs(recoMuons[0].p4().DeltaPhi(x.p4()))>2. and recoAK4_lepTop_jets[0].p4().DeltaR(x.p4())>1.6]
                    recoAK8jets[sys].sort(key=lambda x:getattr( x, 'pt'+sys ), reverse=True)
                    ####################################################### (consider putting the below in this above if to check)


                    # Isolate hadronic hemisphere AK4s via deltaR(leading AK8, AK4jet) 
                    # change angular selection  based on whether we're looking for W's or tops to make the next steps more efficient
                    # require at least a medium b-tagged jet to ensure we're catching the t->W(->qq~)b rather than a light quark from the W->qq~

                    #removing deltaR(W,b) for check with and without this cut
                    if len(recoAK8jets[sys])>0 and self.evtSelection.startswith('WSel'): 
                        recoAK4_hadrTop_jets = [x for x in recoAK4jets if x.p4().DeltaR(recoAK8jets[sys][0].p4())>0.8 and x.p4().DeltaR(recoAK8jets[sys][0].p4())<1.6]# and x.btagDeepFlavB > self.minBDisc] #
                        recoAK4_hadrTop_jets.sort(key=lambda x:x.pt,reverse=True)
                    elif len(recoAK8jets[sys])>0 and self.evtSelection.startswith('topSel'):  
                        recoAK4_hadrTop_jets = [x for x in recoAK4jets if x.p4().DeltaR(recoAK8jets[sys][0].p4())<0.8]# and x.btagDeepFlavB > self.minBDisc] 
                        recoAK4_hadrTop_jets.sort(key=lambda x:x.pt,reverse=True)
                    #######################################################

                    #Not a cut, just so we can plot m_tt when we can reconstruct a ttbar pair
                    if len(recoAK4_hadrTop_jets)!=0: #len(recoAK4_lepTop_jets)!=0 and 

                        
                        if self.evtSelection.startswith('WSel'): 
                            # Considering the leading b-tagged AK4 in the hadronic hemisphere, and outside the W->qq jet
                            # to reconstruct the hadronically decaying top
                            recottbarp4 = [(x+(recoAK4_hadrTop_jets[0].p4()+recoAK8jets[sys][0].p4())) for x in recoleptonicTop]
                        
                        elif self.evtSelection.startswith('topSel'): 
                            # Here, the hadronically decaying top is reconstructed in it's entirety as an AK8 jet
                            # the leading b-tagged hadronic hemisphere AK4 is required to be within DeltaR(AK8,AK4)<0.8 already above
                            recottbarp4 = [(x+recoAK8jets[sys][0].p4()) for x in recoleptonicTop]
                        
                    recoAK4bjets = recoAK4_lepTop_jets+recoAK4_hadrTop_jets#[x for x in [recoAK4_lepTop_jets[0],recoAK4_hadrTop_jets[0]] if x.btagDeepFlavB > self.minBDisc]
                    recoAK4bjets = [x for x in recoAK4bjets if x.btagDeepFlavB > self.minBDisc]
                    recoAK4bjets.sort(key=lambda x:x.pt,reverse=True)
                    #######################################################


            ############################################################################



            ################################################### Applying selection #########################################################
            
            passSel[sys], iSel[sys] = self.WtopSelection( False, event, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets[sys], MET, sys)
            
            ################################################################################################################################
            
            # Since number of b-tagged jets may vary based on the systematics that modify the FatJet pt, calculate them dynamically here in the loop
            #### Calculate btagging weights on the fly for MC reco
            if self.isMC:  
                
                if self.onlyUnc.startswith('_btagWeightUp'): bTagSFs = [x.btagSF_deepjet_M_up for x in recoAK4bjets]
                elif self.onlyUnc.startswith('_btagWeightDown'): bTagSFs = [x.btagSF_deepjet_M_down for x in recoAK4bjets]
                else: bTagSFs = [x.btagSF_deepjet_M for x in recoAK4bjets]
                w=1.
                for i in bTagSFs:
                    w *= i  
                self.btaggingWeight = w 
                btagweights[sys] = w
            ########################################################
           


        ######################### Weights ############################
        weight=1.
        
        #### b-tagging Weights #####
        #if self.isMC: weight = weight#self.btaggingWeight 
        ##################################################
        

        #### Lepton Weights ####
        if self.isMC:
            if len(recoMuons)>0: leptonWeights = self.leptonSF( "muon", recoMuons[0] )
            else: leptonWeights = [0, 0, 0, 0]
        else: leptonWeights = [1, 1, 1, 1.]
        self.leptonWeight = np.prod(leptonWeights)
        ##################################################

        if self.isMC:
            
            self.puWeight = event.puWeight
            #### Applying ALL remaining object-related weights ####
            weight = event.puWeight * event.genWeight * self.leptonWeight * self.topweight * btagweights['_nom'] # last term includes btag weight updated in WtopSel function            
        
        else:
            weight = 1

        self.totalRecoWeight = weight
        
        ##############################################################

        # filling test histos
        if not self.onlyUnc:

            #### Checking no selection without weights
            getattr( self, 'nPVs_noSelnoWeight' ).Fill( getattr( event, 'PV_npvsGood') )
            
            getattr( self, 'nleps_noSelnoWeight' ).Fill( nleptons )
            
            for imuon in recoMuons:
                getattr( self, 'muons_pt_noSelnoWeight' ).Fill( imuon.pt )
                getattr( self, 'muons_eta_noSelnoWeight' ).Fill( imuon.eta )
                getattr( self, 'muons_phi_noSelnoWeight' ).Fill( imuon.phi )

            for iele in recoElectrons:
                getattr( self, 'eles_pt_noSelnoWeight' ).Fill( iele.pt )
                getattr( self, 'eles_eta_noSelnoWeight' ).Fill( iele.eta )
                getattr( self, 'eles_phi_noSelnoWeight' ).Fill( iele.phi )

            getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recoAK8jets['_nom']) )
            getattr( self, 'HT_noSelnoWeight' ).Fill( AK8HT )

            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK8jets_rap_noSelnoWeight' ).Fill( ijet.rapidity )
                getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop_nom if self.evtSelection.startswith('W') else ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI )
            
            getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4jets) )
            getattr( self, 'nAK4bjets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
            getattr( self, 'nAK4bhadjets_noSelnoWeight' ).Fill( len(recoAK4_hadrTop_jets) )
            getattr( self, 'nAK4blepjets_noSelnoWeight' ).Fill( len(recoAK4_lepTop_jets) )

            for ijet in recoAK4jets:
                getattr( self, 'AK4jets_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'AK4jets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK4jets_phi_noSelnoWeight' ).Fill( ijet.phi )

            for ijet in recoAK4bjets:
                getattr( self, 'AK4btaggedjets_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'AK4btaggedjets_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'AK4btaggedjets_phi_noSelnoWeight' ).Fill( ijet.phi )

            getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )
            for lW in recoleptW:
                getattr( self, 'leptonicW_pt_noSelnoWeight' ).Fill( lW.Pt() )
                getattr( self, 'leptonicW_eta_noSelnoWeight' ).Fill( lW.Eta() )
                getattr( self, 'leptonicW_phi_noSelnoWeight' ).Fill( lW.Phi() )
                getattr( self, 'leptonicW_mass_noSelnoWeight' ).Fill( lW.M() )
                getattr( self, 'leptonicWMT_noSelnoWeight').Fill( lW.Mt() )

            for ttbar in recottbarp4:
                getattr( self, 'Mtt_noSelnoWeight').Fill( ttbar.M())

            for lT in recoleptonicTop:
                getattr( self, 'leptonicTop_pt_noSelnoWeight' ).Fill( lT.Pt() )
                getattr( self, 'leptonicTop_eta_noSelnoWeight' ).Fill( lT.Eta() )
                getattr( self, 'leptonicTop_phi_noSelnoWeight' ).Fill( lT.Phi() )
                getattr( self, 'leptonicTop_mass_noSelnoWeight' ).Fill( lT.M() )

            for ijet in recoAK4_lepTop_jets:
                getattr( self, 'lepHemB_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'lepHemB_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'lepHemB_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'lepHemB_mass_noSelnoWeight' ).Fill( ijet.mass )
            for ijet in recoAK4_hadrTop_jets:
                getattr( self, 'hadHemB_pt_noSelnoWeight' ).Fill( ijet.pt )
                getattr( self, 'hadHemB_eta_noSelnoWeight' ).Fill( ijet.eta )
                getattr( self, 'hadHemB_phi_noSelnoWeight' ).Fill( ijet.phi )
                getattr( self, 'hadHemB_mass_noSelnoWeight' ).Fill( ijet.mass )

            #### Checking no selection with weights
            getattr( self, 'nPVs_noSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
            
            getattr( self, 'nleps_noSel' ).Fill( nleptons , weight)
            
            for imuon in recoMuons:
                getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt , weight)
                getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta , weight)
                getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi , weight)

            for iele in recoElectrons:
                getattr( self, 'eles_pt_noSel' ).Fill( iele.pt , weight)
                getattr( self, 'eles_eta_noSel' ).Fill( iele.eta , weight)
                getattr( self, 'eles_phi_noSel' ).Fill( iele.phi , weight)

            getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets['_nom']) , weight)
            getattr( self, 'HT_noSel' ).Fill( AK8HT , weight)

            for ijet in recoAK8jets['_nom']:
                getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt , weight)
                getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta , weight)
                getattr( self, 'AK8jets_rap_noSel' ).Fill( ijet.rapidity , weight)
                getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi , weight)
                getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop_nom if self.evtSelection.startswith('W') else ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI , weight)
            
            getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4jets) , weight )
            getattr( self, 'nAK4bjets_noSelnoWeight' ).Fill( len(recoAK4bjets)  , weight)
            getattr( self, 'nAK4bhadjets_noSelnoWeight' ).Fill( len(recoAK4_hadrTop_jets) , weight )
            getattr( self, 'nAK4blepjets_noSelnoWeight' ).Fill( len(recoAK4_lepTop_jets) , weight )

            for ijet in recoAK4jets:
                getattr( self, 'AK4jets_pt_noSel' ).Fill( ijet.pt , weight)
                getattr( self, 'AK4jets_eta_noSel' ).Fill( ijet.eta , weight)
                getattr( self, 'AK4jets_phi_noSel' ).Fill( ijet.phi , weight)
                
            for ijet in recoAK4bjets:
                getattr( self, 'AK4btaggedjets_pt_noSel' ).Fill( ijet.pt , weight)
                getattr( self, 'AK4btaggedjets_eta_noSel' ).Fill( ijet.eta , weight)
                getattr( self, 'AK4btaggedjets_phi_noSel' ).Fill( ijet.phi , weight)

            getattr( self, 'METPt_noSel' ).Fill( MET.Pt() , weight)
            
            for lW in recoleptW:
                getattr( self, 'leptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                getattr( self, 'leptonicW_eta_noSel' ).Fill( lW.Eta() , weight)
                getattr( self, 'leptonicW_phi_noSel' ).Fill( lW.Phi() , weight)
                getattr( self, 'leptonicW_mass_noSel' ).Fill( lW.M() , weight)
                getattr( self, 'leptonicWMT_noSel').Fill( lW.Mt() , weight)
            
            for ttbar in recottbarp4:
                getattr( self, 'Mtt_noSel').Fill( ttbar.M() , weight)

            for lT in recoleptonicTop:
                getattr( self, 'leptonicTop_pt_noSel' ).Fill( lT.Pt() , weight)
                getattr( self, 'leptonicTop_eta_noSel' ).Fill( lT.Eta() , weight)
                getattr( self, 'leptonicTop_phi_noSel' ).Fill( lT.Phi() , weight)
                getattr( self, 'leptonicTop_mass_noSel' ).Fill( lT.M() , weight)

            
            

            for ijet in recoAK4_lepTop_jets:
                getattr( self, 'lepHemB_pt_noSel' ).Fill( ijet.pt, weight)
                getattr( self, 'lepHemB_eta_noSel' ).Fill( ijet.eta, weight)
                getattr( self, 'lepHemB_phi_noSel' ).Fill( ijet.phi, weight)
                getattr( self, 'lepHemB_mass_noSel' ).Fill( ijet.mass, weight)
            for ijet in recoAK4_hadrTop_jets:
                getattr( self, 'hadHemB_pt_noSel' ).Fill( ijet.pt, weight)
                getattr( self, 'hadHemB_eta_noSel' ).Fill( ijet.eta, weight)
                getattr( self, 'hadHemB_phi_noSel' ).Fill( ijet.phi, weight)
                getattr( self, 'hadHemB_mass_noSel' ).Fill( ijet.mass, weight)
                


            reweight = self.totalRecoWeight 
            #### Checking nominal selection with weights
            if ('_nom' in iSel.keys()) and passSel['_nom'] and iSel['_nom']:

                # basic reco histos
                getattr( self, 'nPVs'+iSel['_nom'] ).Fill( getattr( event, 'PV_npvsGood'), reweight )
            
                getattr( self, 'nleps'+iSel['_nom'] ).Fill( nleptons , reweight)
                
                for imuon in recoMuons:
                    getattr( self, 'muons_pt'+iSel['_nom'] ).Fill( imuon.pt , reweight)
                    getattr( self, 'muons_eta'+iSel['_nom'] ).Fill( imuon.eta , reweight)
                    getattr( self, 'muons_phi'+iSel['_nom'] ).Fill( imuon.phi , reweight)

                for iele in recoElectrons:
                    getattr( self, 'eles_pt'+iSel['_nom'] ).Fill( iele.pt , reweight)
                    getattr( self, 'eles_eta'+iSel['_nom'] ).Fill( iele.eta , reweight)
                    getattr( self, 'eles_phi'+iSel['_nom'] ).Fill( iele.phi , reweight)

                getattr( self, 'nAK8jets'+iSel['_nom'] ).Fill( len(recoAK8jets['_nom']) , reweight)
                getattr( self, 'HT'+iSel['_nom'] ).Fill( AK8HT , reweight)

                for ijet in recoAK8jets['_nom']:
                    getattr( self, 'AK8jets_pt'+iSel['_nom'] ).Fill( ijet.pt , reweight)
                    getattr( self, 'AK8jets_eta'+iSel['_nom'] ).Fill( ijet.eta , reweight)
                    getattr( self, 'AK8jets_rap'+iSel['_nom'] ).Fill( ijet.rapidity , reweight)
                    getattr( self, 'AK8jets_phi'+iSel['_nom'] ).Fill( ijet.phi , reweight)
                    getattr( self, 'AK8jets_mass'+iSel['_nom'] ).Fill( ijet.msoftdrop_nom if self.evtSelection.startswith('W') else ijet.msoftdrop_nom/ijet.msoftdrop_corr_PUPPI , reweight)
                
                getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4jets) , reweight )
                getattr( self, 'nAK4bjets_noSelnoWeight' ).Fill( len(recoAK4bjets)  , reweight)
                getattr( self, 'nAK4bhadjets_noSelnoWeight' ).Fill( len(recoAK4_hadrTop_jets) , reweight )
                getattr( self, 'nAK4blepjets_noSelnoWeight' ).Fill( len(recoAK4_lepTop_jets) , reweight )

                for ijet in recoAK4jets:
                    getattr( self, 'AK4jets_pt'+iSel['_nom'] ).Fill( ijet.pt , reweight)
                    getattr( self, 'AK4jets_eta'+iSel['_nom'] ).Fill( ijet.eta , reweight)
                    getattr( self, 'AK4jets_phi'+iSel['_nom'] ).Fill( ijet.phi , reweight)
                    
                for ijet in recoAK4bjets:
                    getattr( self, 'AK4btaggedjets_pt'+iSel['_nom'] ).Fill( ijet.pt , reweight)
                    getattr( self, 'AK4btaggedjets_eta'+iSel['_nom'] ).Fill( ijet.eta , reweight)
                    getattr( self, 'AK4btaggedjets_phi'+iSel['_nom'] ).Fill( ijet.phi , reweight)

                getattr( self, 'METPt'+iSel['_nom'] ).Fill( MET.Pt() , reweight)
                
                for lW in recoleptW:
                    getattr( self, 'leptonicW_pt'+iSel['_nom'] ).Fill( lW.Pt(), reweight )
                    getattr( self, 'leptonicW_eta'+iSel['_nom'] ).Fill( lW.Eta() , reweight)
                    getattr( self, 'leptonicW_phi'+iSel['_nom'] ).Fill( lW.Phi() , reweight)
                    getattr( self, 'leptonicW_mass'+iSel['_nom'] ).Fill( lW.M() , reweight)
                    getattr( self, 'leptonicWMT'+iSel['_nom']).Fill( lW.Mt() , reweight)

                for ttbar in recottbarp4:
                    getattr( self, 'Mtt'+iSel['_nom']).Fill( ttbar.M() , reweight)                

                for lT in recoleptonicTop:
                    getattr( self, 'leptonicTop_pt'+iSel['_nom'] ).Fill( lT.Pt() , reweight)
                    getattr( self, 'leptonicTop_eta'+iSel['_nom'] ).Fill( lT.Eta() , reweight)
                    getattr( self, 'leptonicTop_phi'+iSel['_nom'] ).Fill( lT.Phi() , reweight)
                    getattr( self, 'leptonicTop_mass'+iSel['_nom'] ).Fill( lT.M() , reweight)

                for ijet in recoAK4_lepTop_jets:
                    getattr( self, 'lepHemB_pt'+iSel['_nom'] ).Fill( ijet.pt, reweight)
                    getattr( self, 'lepHemB_eta'+iSel['_nom'] ).Fill( ijet.eta, reweight)
                    getattr( self, 'lepHemB_phi'+iSel['_nom'] ).Fill( ijet.phi, reweight)
                    getattr( self, 'lepHemB_mass'+iSel['_nom'] ).Fill( ijet.mass, reweight)
                for ijet in recoAK4_hadrTop_jets:
                    getattr( self, 'hadHemB_pt'+iSel['_nom'] ).Fill( ijet.pt, reweight)
                    getattr( self, 'hadHemB_eta'+iSel['_nom'] ).Fill( ijet.eta, reweight)
                    getattr( self, 'hadHemB_phi'+iSel['_nom'] ).Fill( ijet.phi, reweight)
                    getattr( self, 'hadHemB_mass'+iSel['_nom'] ).Fill( ijet.mass, reweight)
            
                #print (met.pt, MET.Pt())

        return passSel, iSel, recoMuons, recoElectrons, recoAK4bjets, recoAK8jets, MET, btagweights

    #############################################################################
    def genSelection( self, event ):
        '''Analyzing gen information'''

        genJetsAK8 = list(Collection( event, 'GenJetAK8' ))
        for ijets in genJetsAK8: ijets.rapidity = self.etaToRapidity(ijets)

        genLeptons = list(Collection( event, 'GenDressedLepton' ))
        genJets = list(Collection( event, 'GenJet'))
        genParticles = Collection(event, 'GenPart')
        genmet = Object( event, 'GenMET')
        
        tops =  [x for x in genParticles if x.pdgId==6 and x.statusFlags==14]
        tops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)
        antiTops =  [x for x in genParticles if x.pdgId==-6 and x.statusFlags==14]
        antiTops.sort(key=lambda x:getattr( x, 'pt' ), reverse=True)

        #### Lepton selection
        genElectrons  = [x for x in genLeptons if abs(x.pdgId)==11 and x.pt>self.minLooseElectronPt and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))] #and abs(x.eta)<self.maxMuonEta ]
        genElectrons.sort(key=lambda x:x.pt, reverse=True)

        genMuons = [ x for x in genLeptons if abs(x.pdgId)==13 and  x.pt > self.minTightMuonPt and abs(x.eta) < self.maxMuonEta ]
        genMuons.sort(key=lambda x:x.pt, reverse=True)
        
        ngenLeptons = len(genMuons)+len(genElectrons)
        ##################################################

        #### MET (not sure if needed)
        genMET = ROOT.TLorentzVector()
        genMET.SetPtEtaPhiM(genmet.pt, 0., genmet.phi, 0)
        ##################################################

        #### Basic AK4 jet selection
        genAK4jets = [ x for x in genJets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta]# and abs(x.hadronFlavour)==5 ] 
        ##################################################

        #### Basic AK8 jet selection 
        genAK8jets = [ x for x in genJetsAK8 if x.pt > self.minAK8JetPt and abs(x.rapidity) < self.maxJetAK8Rap]
        genAK8HT = sum( [ x.pt for x in genAK8jets ] )
        ##################################################    
        
        genAK4bjets = []
        genAK4_lepTop_jets = []
        genAK4_hadrTop_jets = []
        genleptonicTop = []#ROOT.TLorentzVector(0.,0.,0.,0.)
        genleptW = []#ROOT.TLorentzVector(0.,0.,0.,0.) 
        genttbarp4 = []

        ################## selections on jets ####################
        if len(genMuons)!=0:

            # isolate leptonic hemisphere AK4's via 0.4<deltaR(mu,AK4jet)<1.6, require |hadronFlavour|==5 matched AK4's in this hemisphere
            genAK4_lepTop_jets = [x for x in genAK4jets if (x.p4().DeltaR(genMuons[0].p4())<1.6 and x.p4().DeltaR( genMuons[0].p4() )>0.4)]# and abs(x.hadronFlavour)==5 ] 
            genAK4_lepTop_jets.sort(key=lambda x:x.pt,reverse=True) 
            #######################################################

            # leptonic W definition
            genleptW = [genMET+genMuons[0].p4()]

            # leptonic top definition
            if len(genAK4_lepTop_jets)!=0: 
                genleptonicTop = [x+genAK4_lepTop_jets[0].p4() for x in genleptW]

                # Ensure separation of the AK8s from the leptonic hemisphere by deltaPhi(mu,AK8) and deltaR(AK8,lept hem AK4)
                genAK8jets = [ x for x in genAK8jets if abs(genMuons[0].p4().DeltaPhi(x.p4()))>2. and genAK4_lepTop_jets[0].p4().DeltaR(x.p4())>1.6]
                genAK8jets.sort(key=lambda x:x.pt,reverse=True)

                # Isolate hadronic hemisphere AK4 bjets via deltaR(leading AK8, AK4bjet) 
                # change angular selection  based on whether we're looking for W's or tops to make the next steps more efficient
                # require hadronFlavour matched b jets to ensure we're catching the t->W(->qq~)b rather than a light quark from the W->qq~
                
                # v2 select tops only if no further AK4's in x.p4().DeltaR(genAK8jets[0].p4())>0.8 and x.p4().DeltaR(genAK8jets[0].p4())<1.6
                # do in Wtop sel function

                #genAK8jets = [ x for x in genAK8jets if abs(genMuons[0].p4().DeltaPhi(x.p4()))>2. and genAK4_lepTop_jets[0].p4().DeltaR(x.p4())>1.6]
                #genAK8jets.sort(key=lambda x:x.pt,reverse=True)

                if len(genAK8jets)>0 and self.evtSelection.startswith('WSel'): 
                    genAK4_hadrTop_jets = [x for x in genAK4jets if x.p4().DeltaR(genAK8jets[0].p4())>0.8 and x.p4().DeltaR(genAK8jets[0].p4())<1.6]# and abs(x.hadronFlavour)==5] #
                    genAK4_hadrTop_jets.sort(key=lambda x:x.pt,reverse=True)
                elif len(genAK8jets)>0 and self.evtSelection.startswith('topSel'):  
                    genAK4_hadrTop_jets = [x for x in genAK4jets if x.p4().DeltaR(genAK8jets[0].p4())<0.8]# and abs(x.hadronFlavour)==5] 
                    genAK4_hadrTop_jets.sort(key=lambda x:x.pt,reverse=True)


                # this is not a cut, used just to reconstruct ttbar pair, if possible, for m_tt measurement 
                if len(genAK4_hadrTop_jets)!=0: #len(genAK4_lepTop_jets)!=0 and 

                    
                    if self.evtSelection.startswith('WSel'):
                        # Considering the leading AK4 b in the hadronic hemisphere, and outside the W->qq jet
                        # to reconstruct the hadronically decaying top 
                        genttbarp4 = [(x+(genAK4_hadrTop_jets[0].p4()+genAK8jets[0].p4())) for x in genleptonicTop] 
                    
                    elif self.evtSelection.startswith('topSel'): 
                        # Here, the hadronically decaying top is reconstructed in its entirety as an AK8 jet
                        # the leading hadronic hemisphere AK4 b jet is required to be within DeltaR(AK8,AK4)<0.8 already above
                        genttbarp4 = [(x+genAK8jets[0].p4()) for x in genleptonicTop]
                    
        
                genAK4bjets = genAK4_lepTop_jets+genAK4_hadrTop_jets#[x for x in [genAK4_lepTop_jets[0],genAK4_hadrTop_jets[0]] if abs(x.hadronFlavour)==5]
                genAK4bjets = [x for x in genAK4bjets if abs(x.hadronFlavour)==5]
                genAK4bjets.sort(key=lambda x:x.pt,reverse=True)
               
        #####################################################################

        
        
        ##### Applying selection
        #### Creating Nsub basis, filling histos and creating branches IF passSel
        passSel, iSel = self.WtopSelection( True, event, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET,'' )
        #### Weight
        weight = event.genWeight
        #self.totalGenWeight = weight

        if not self.onlyUnc: 
            #### Checking no selection
            getattr( self, 'ngenleps_noSel' ).Fill( len(genMuons)+len(genElectrons), weight )
            
            for imuon in genMuons:
                getattr( self, 'genmuons_pt_noSel' ).Fill( imuon.pt, weight )
                getattr( self, 'genmuons_eta_noSel' ).Fill( imuon.eta, weight )
                getattr( self, 'genmuons_phi_noSel' ).Fill( imuon.phi, weight )
            
            for iele in genElectrons:
                getattr( self, 'geneles_pt_noSel' ).Fill( iele.pt, weight )
                getattr( self, 'geneles_eta_noSel' ).Fill( iele.eta, weight )
                getattr( self, 'geneles_phi_noSel' ).Fill( iele.phi, weight )
            
            getattr( self, 'ngenAK8jets_noSel' ).Fill( len(genAK8jets), weight )
            getattr( self, 'genHT_noSel' ).Fill( genAK8HT, weight )
            
            for ijet in genAK8jets:
                getattr( self, 'AK8genjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK8genjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK8genjets_rap_noSel' ).Fill( ijet.rapidity, weight )
                getattr( self, 'AK8genjets_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'AK8genjets_mass_noSel' ).Fill( ijet.mass, weight )
            
            getattr( self, 'ngenAK4jets_noSel' ).Fill( len(genAK4jets), weight )
            getattr( self, 'ngenAK4bjets_noSel' ).Fill( len(genAK4bjets), weight )
            getattr( self, 'ngenAK4bhadjets_noSel' ).Fill( len(genAK4_hadrTop_jets), weight )
            getattr( self, 'ngenAK4blepjets_noSel' ).Fill( len(genAK4_lepTop_jets), weight )
            
            for ijet in genAK4jets:
                getattr( self, 'AK4genjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK4genjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK4genjets_phi_noSel' ).Fill( ijet.phi, weight )
            for ijet in genAK4bjets:
                getattr( self, 'AK4btaggedgenjets_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'AK4btaggedgenjets_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'AK4btaggedgenjets_phi_noSel' ).Fill( ijet.phi, weight )
            
            getattr( self, 'genMETPt_noSel' ).Fill( genMET.Pt(), weight )

            for lW in genleptW:
                getattr( self, 'genleptonicW_pt_noSel' ).Fill( lW.Pt(), weight )
                getattr( self, 'genleptonicW_eta_noSel' ).Fill( lW.Eta(), weight )
                getattr( self, 'genleptonicW_phi_noSel' ).Fill( lW.Phi(), weight )
                getattr( self, 'genleptonicW_mass_noSel' ).Fill( lW.M(), weight )
                getattr( self, 'genleptonicWMT_noSel' ).Fill( lW.Mt(), weight )

            for ttbar in genttbarp4:
                getattr( self, 'genMtt_noSel' ).Fill( ttbar.M(), weight )

            for lT in genleptonicTop:
                getattr( self, 'genleptonicTop_pt_noSel' ).Fill( lT.Pt(), weight )
                getattr( self, 'genleptonicTop_eta_noSel' ).Fill( lT.Eta(), weight )
                getattr( self, 'genleptonicTop_phi_noSel' ).Fill( lT.Phi(), weight )
                getattr( self, 'genleptonicTop_mass_noSel' ).Fill( lT.M(), weight )

            for ijet in genAK4_lepTop_jets:
                getattr( self, 'genlepHemB_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'genlepHemB_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'genlepHemB_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'genlepHemB_mass_noSel' ).Fill( ijet.mass, weight )
            for ijet in genAK4_hadrTop_jets:
                getattr( self, 'genhadHemB_pt_noSel' ).Fill( ijet.pt, weight )
                getattr( self, 'genhadHemB_eta_noSel' ).Fill( ijet.eta, weight )
                getattr( self, 'genhadHemB_phi_noSel' ).Fill( ijet.phi, weight )
                getattr( self, 'genhadHemB_mass_noSel' ).Fill( ijet.mass, weight )

            ##### Filling histograms
            if passSel and iSel:
                #print ("Filling gen histos")
                #### Checking nominal selection with weights
                getattr( self, 'ngenleps'+iSel ).Fill( len(genMuons)+len(genElectrons), weight )
                for imuon in genMuons:
                    getattr( self, 'genmuons_pt'+iSel ).Fill( imuon.pt, weight )
                    getattr( self, 'genmuons_eta'+iSel ).Fill( imuon.eta, weight )
                    getattr( self, 'genmuons_phi'+iSel ).Fill( imuon.phi, weight )
                for iele in genElectrons:
                    getattr( self, 'geneles_pt'+iSel ).Fill( iele.pt, weight )
                    getattr( self, 'geneles_eta'+iSel ).Fill( iele.eta, weight )
                    getattr( self, 'geneles_phi'+iSel ).Fill( iele.phi, weight )
                getattr( self, 'ngenAK8jets'+iSel ).Fill( len(genAK8jets), weight )
                getattr( self, 'genHT'+iSel ).Fill( genAK8HT, weight )
                for ijet in genAK8jets:
                    getattr( self, 'AK8genjets_pt'+iSel ).Fill( ijet.pt, weight )
                    getattr( self, 'AK8genjets_eta'+iSel ).Fill( ijet.eta, weight )
                    getattr( self, 'AK8genjets_rap'+iSel ).Fill( ijet.rapidity, weight )
                    getattr( self, 'AK8genjets_phi'+iSel ).Fill( ijet.phi, weight )
                    getattr( self, 'AK8genjets_mass'+iSel ).Fill( ijet.mass, weight )

                getattr( self, 'ngenAK4jets_noSel' ).Fill( len(genAK4jets), weight )
                getattr( self, 'ngenAK4bjets_noSel' ).Fill( len(genAK4bjets), weight )
                getattr( self, 'ngenAK4bhadjets_noSel' ).Fill( len(genAK4_hadrTop_jets), weight )
                getattr( self, 'ngenAK4blepjets_noSel' ).Fill( len(genAK4_lepTop_jets), weight )

                for ijet in genAK4jets:
                    getattr( self, 'AK4genjets_pt'+iSel  ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4genjets_eta'+iSel  ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4genjets_phi'+iSel  ).Fill( ijet.phi, weight )
                for ijet in genAK4bjets:
                    getattr( self, 'AK4btaggedgenjets_pt'+iSel  ).Fill( ijet.pt, weight )
                    getattr( self, 'AK4btaggedgenjets_eta'+iSel  ).Fill( ijet.eta, weight )
                    getattr( self, 'AK4btaggedgenjets_phi'+iSel ).Fill( ijet.phi, weight )
                getattr( self, 'genMETPt'+iSel ).Fill( genMET.Pt(), weight )
                
                for lW in genleptW:
                    getattr( self, 'genleptonicW_pt'+iSel  ).Fill( lW.Pt(), weight )
                    getattr( self, 'genleptonicW_eta'+iSel  ).Fill( lW.Eta(), weight )
                    getattr( self, 'genleptonicW_phi'+iSel  ).Fill( lW.Phi(), weight )
                    getattr( self, 'genleptonicW_mass'+iSel  ).Fill( lW.M(), weight )
                    getattr( self, 'genleptonicWMT'+iSel ).Fill( lW.Mt(), weight )

                for ttbar in genttbarp4:
                    getattr( self, 'genMtt'+iSel).Fill( ttbar.M(), weight )

                for lT in genleptonicTop:
                    getattr( self, 'genleptonicTop_pt'+iSel  ).Fill( lT.Pt(), weight )
                    getattr( self, 'genleptonicTop_eta'+iSel  ).Fill( lT.Eta(), weight )
                    getattr( self, 'genleptonicTop_phi'+iSel  ).Fill( lT.Phi(), weight )
                    getattr( self, 'genleptonicTop_mass'+iSel  ).Fill( lT.M(), weight )

                for ijet in genAK4_lepTop_jets:
                    getattr( self, 'genlepHemB_pt'+iSel  ).Fill( ijet.pt, weight )
                    getattr( self, 'genlepHemB_eta'+iSel  ).Fill( ijet.eta, weight )
                    getattr( self, 'genlepHemB_phi'+iSel  ).Fill( ijet.phi, weight )
                    getattr( self, 'genlepHemB_mass'+iSel  ).Fill( ijet.mass, weight )
                for ijet in genAK4_hadrTop_jets:
                    getattr( self, 'genhadHemB_pt'+iSel  ).Fill( ijet.pt, weight )
                    getattr( self, 'genhadHemB_eta'+iSel  ).Fill( ijet.eta, weight )
                    getattr( self, 'genhadHemB_phi'+iSel  ).Fill( ijet.phi, weight )
                    getattr( self, 'genhadHemB_mass'+iSel  ).Fill( ijet.mass, weight )
                #print (genmet.pt, genMET.Pt())


        return passSel, iSel, genMuons, genElectrons, genAK4bjets, genAK8jets, genMET

    #############################################################################
    def WtopSelection( self, isGen, event, muons, electrons, AK4bjets, AK8jets, MET, ptLabel):
    
        if (len(muons)==1) and (len(electrons)==0) and (len(AK8jets)>0) and (len(AK4bjets)>= 1):                
            
            leadJetpT = getattr( AK8jets[0], 'pt'+ptLabel )
            
            leptWpT=muons[0].p4()+MET
            
            if leptWpT.Pt()>self.minLeptonicWPt and (MET.Pt()>self.METCutWtop) and leadJetpT>self.minLeadAK8JetPtW:
                
                ################ Cuts on AK8 to finally consider whether the fatjet is a W/top candidate #################                  
                jetMass = AK8jets[0].mass if isGen else AK8jets[0].msoftdrop_nom
                if 'WSel' in self.evtSelection: 
                    if (((jetMass<self.maxSDMassW) and (jetMass>=self.minSDMassW)) and (leadJetpT>self.minLeadAK8JetPtW)): 
                        #print (leptWpT.Mt(),isGen)                           
                        return True, '_WSel' 
                    else: return False, None

                elif 'topSel' in self.evtSelection:
                    if ((jetMass/ (1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) >= self.minSDMassTop) and (jetMass/(1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop) and (leadJetpT> self.minLeadAK8JetPtTop)): # 
                    #if ((jetMass/ (1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) > 150. ) and (leadJetpT> 400.)):#self.minLeadAK8JetPtTop)): # >= self.minSDMassTop) and (jetMass/(1 if isGen else AK8jets[0].msoftdrop_corr_PUPPI) < self.maxSDMassTop
                        #print (leptWpT.Mt(),isGen)                           
                        return True, '_topSel'
                    else: return False, None    
                ###########################################################################################################

            else: return False, None 
        else: return False, None
        
        
    #############################################################################
    def createNsubBasis(self, AK8jet, event, PFCollection ):
        '''Generic, taking a AK8 jet and computing Nsub basis from PFCollection'''

        pfCands = list(Collection(event, PFCollection ))
        ak8jet = {}          ### Storing good jet as list for later use

        ##### Computing quantities
        ak8jet['jet'] = AK8jet

        #### Run calculations of NSub bases and store for ungroomed AK8jets (default in CMS)

        #### Applying PUPPI weights to the PF candidates
        constituents = ROOT.vector("TLorentzVector")()
        CandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()
        for p in pfCands :
            tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
            try: tp = tp * p.puppiWeight
            except RuntimeError: tp = tp    ### for genjets
            CandsPUPPIweightedVec.push_back(tp)

        #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
        for x in CandsPUPPIweightedVec:
            if abs(AK8jet.p4().DeltaR( x )) < 0.8: constituents.push_back(x) #Shouldn't this be DrRapidityPhi?

        #### Computing n-subjetiness basis from PF PUPPI constituents
        nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
        nsub1 = self.nSub1.getTau( self.maxTau, constituents )
        nsub2 = self.nSub2.getTau( self.maxTau, constituents )
        nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, constituents )    ### needed for genjet tau21 tau32


        ### default in CMS WTA_kT https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/python/nJettinessAdder_cfi.py
        try: ak8jet['tau21'] = nsub1_WTA_kT[1]/nsub1_WTA_kT[0]
        except ZeroDivisionError: ak8jet['tau21'] = -1.
        try: ak8jet['tau32'] = nsub1_WTA_kT[2]/nsub1_WTA_kT[1]
        except ZeroDivisionError: ak8jet['tau32'] = -1.


        #### filling histos and branches with nsub basis
        for tauN in range(self.maxTau):
            ak8jet['0p5'+str(tauN+1)] = nsub0p5[tauN]
            ak8jet['1'+str(tauN+1)] = nsub1[tauN]
            ak8jet['2'+str(tauN+1)] = nsub2[tauN]
        '''
        #### Computing Softdrop jets
            sdAK8jets = self.sd.result( constituents ) #CandsPUPPIweightedVec )

            ### Storing good jet as list for later use
            if len(sdAK8jets)>0:

                ak8jet['sdjet'] = sdAK8jets[0]

                # Cluster only the particles near the appropriate jet to save time
                sd_constituents =  ROOT.vector("TLorentzVector")()

                for x in ak8jet['sdjet'].constituents():
                    sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
                    if abs(ak8jet['sdjet'].delta_R( x )) < 0.8:
                        sd_constituents.push_back(sd_constits)
                sd_nsub0p5 = self.nSub0p5.getTau( self.maxTau, sd_constituents )
                sd_nsub1 = self.nSub1.getTau( self.maxTau, sd_constituents )
                sd_nsub2 = self.nSub2.getTau( self.maxTau, sd_constituents )
                sd_nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( 3, sd_constituents )

                try: ak8jet['sdtau21'] = sd_nsub1_WTA_kT[1]/sd_nsub1_WTA_kT[0]
                except ZeroDivisionError: ak8jet['sdtau21'] = -1
                try: ak8jet['sdtau32'] = sd_nsub1_WTA_kT[2]/sd_nsub1_WTA_kT[1]
                except ZeroDivisionError: ak8jet['sdtau32'] = -1

                for tauN in range(self.maxTau):
                    ak8jet['sd0p5'+str(tauN+1)] = sd_nsub0p5[tauN]
                    ak8jet['sd1'+str(tauN+1)] = sd_nsub1[tauN]
                    ak8jet['sd2'+str(tauN+1)] = sd_nsub2[tauN]
            else:
                ak8jet['sdjet'] = ROOT.TLorentzVector( )
                ak8jet['sdtau21'] = -1
                ak8jet['sdtau32'] = -1
                for tauN in range(self.maxTau):
                    ak8jet['sd0p5'+str(tauN+1)] = -1
                    ak8jet['sd1'+str(tauN+1)] = -1
                    ak8jet['sd2'+str(tauN+1)] = -1
        '''
        return ak8jet
        
    #############################################################################
    '''
    def dijetSelection( self, event, muons, electrons, AK8jets, ptAsym, deltaPhi ):

        if (len(muons)+len(electrons)==0) and (len(AK8jets)>1) and (ptAsym<0.3) and (deltaPhi>2) and (AK8jets[1].pt>self.minLeadAK8JetPtDijet):         ### PT REQUIREMENT NEEDS TO BE REVISIT
            return True
        else: return False

    '''
    #############################################################################
    def fillBranches( self, jetLabel, jetInfo ): #update later if you ever decide to go back to trees

        #if 'true' in jetLabel : print (jetLabel) 
        #### Filling branch with passAK8jet info after selection
        self.out.fillBranch( 'n'+jetLabel, len(jetInfo) )
        self.out.fillBranch(jetLabel+"_mass", [ iJ['jet'].mass for i,iJ in jetInfo.items() ] )
        if not 'gen' in jetLabel.lower(): 
            if not 'top'in self.evtSelection: self.out.fillBranch(jetLabel+"_msoftdrop", [ iJ['jet'].msoftdrop_nom for i,iJ in jetInfo.items() ] )
            else: self.out.fillBranch(jetLabel+"_msoftdrop", [ iJ['jet'].msoftdrop_nom/iJ['jet'].msoftdrop_corr_PUPPI for i,iJ in jetInfo.items() ] )
        self.out.fillBranch(jetLabel+"_pt", [ iJ['jet'].pt for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_eta", [ iJ['jet'].eta for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_phi", [ iJ['jet'].phi for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau21", [ iJ['tau21'] for i,iJ in jetInfo.items()  ] )
        self.out.fillBranch(jetLabel+"_Tau32", [ iJ['tau32'] for i,iJ in jetInfo.items()  ] )
        #self.out.fillBranch(jetLabel+"_SD_Tau21", [ iJ['sdtau21'] for i,iJ in jetInfo.items()  ] )
        #self.out.fillBranch(jetLabel+"_SD_Tau32", [ iJ['sdtau32'] for i,iJ in jetInfo.items()  ] )
        
        for tauN in range(1, self.maxTau+1):
            self.out.fillBranch(jetLabel+"_tau_0p5_"+str(tauN),  [ iJ['0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_1_"+str(tauN),  [ iJ['1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            self.out.fillBranch(jetLabel+"_tau_2_"+str(tauN),  [ iJ['2'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_0p5_"+str(tauN),  [ iJ['sd0p5'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_1_"+str(tauN),  [ iJ['sd1'+str(tauN)] for i,iJ in jetInfo.items() ] )
            #self.out.fillBranch(jetLabel+"SD_tau_2_"+str(tauN),  [ iJ['sd2'+str(tauN)] for i,iJ in jetInfo.items() ] )



    #############################################################################
    def matchRecoGenParticle( self, event, recoJet ):

        genParticles = Collection(event, "GenPart")

        quarksFromW = GenQuarkFromW( genParticles )
        bquarksFromTop = GenBquarkFromTop( genParticles )
        allQuarksFromWtop = quarksFromW + bquarksFromTop

        listMatched = []
        for q in allQuarksFromWtop:
            #print event.event, q.pdgId, genParticles[q.genPartIdxMother].pdgId, q.p4().Pt()
            if recoJet.p4().DeltaR( q.p4() )<0.3: listMatched.append( True )
            else: listMatched.append( False )

        if (len(listMatched)==4) and all(listMatched[:2]) and not all(listMatched[2:]): boosted = 2  ## only boostedW
        elif (len(listMatched)==4) and all(listMatched[:2]) and any(listMatched[2:]): boosted = 4      ## boosted Top
        else: boosted = 0
        #print listMatched, boosted

        return boosted

    ##################### Helpful functions
    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        elif hasattr( c, "Perp"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        else:
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.perp(), c.eta(), c.phi(), c.m() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )

    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret
    
    #############################################################################
    def etaToRapidity( self, ijet ):
        nom = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt * np.cosh(ijet.eta),2) ) + ijet.pt * np.sinh(ijet.eta)
        den = np.sqrt( np.power(ijet.mass,2) + np.power(ijet.pt,2) )
        return np.log(nom/den)


    ###############
    

###################### Different gen functions
def getDaughters(GenParticle,gp):
    ret = []
    tmpListGenParticles = list(GenParticle)
    for part in GenParticle:
        if part != gp:
            if part.genPartIdxMother == tmpListGenParticles.index(gp):
                ret.append(part)
    return ret

def GenBquarkFromTop(GenParticle):
    ret = []
    for gP in GenParticle:
        if (abs(gP.pdgId) == 5) and gP.genPartIdxMother>0:
            mom = GenParticle[gP.genPartIdxMother]
            if abs(mom.pdgId)  == 6:
                dauids = [abs(dau.pdgId) for dau in getDaughters(GenParticle,mom)]
                if 6 not in dauids:
                    ret.append(gP)
    return ret

def GenQuarkFromW(GenParticle):
    ret = []
    for gP in GenParticle:
        if gP.genPartIdxMother>0:
            if (abs(gP.pdgId) <= 5) and (abs(GenParticle[gP.genPartIdxMother].pdgId) == 24):
                ret.append(gP)
        if (abs(gP.pdgId) <= 5) and gP.genPartIdxMother>0:
            mom = GenParticle[gP.genPartIdxMother]
            if abs(mom.pdgId)  == 24:
                dauids = [abs(dau.pdgId) for dau in getDaughters(GenParticle,mom)]
                if 24 not in dauids:
                    ret.append(gP)



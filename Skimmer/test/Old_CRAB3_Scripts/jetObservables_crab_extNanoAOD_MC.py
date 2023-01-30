#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor


from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puAutoWeight_2016 #puWeight_2016
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncertAll #jetRecalib2016BCD, jetRecalib2016BCDAK8Puppi

#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *  

# our module
#from Analysis.jetObservables.nSubProducer_axscan import nsubjettinessProducer #for groomed and ungroomed gen- and reco-jets in semileptonic ttbar with axisdef={0,3,6}  
from Analysis.jetObservables.nSubProducer_boostedW_axdefScan import nsubjettinessProducer
#from Analysis.jetObservables.nSubProducer import nsubjettinessProducer #vanilla, for ungroomed recojets in semileptonic ttbar 
#from Analysis.jetObservables.nSubProducer_gen_reco import nsubjettinessProducer #for ungroomed gen- and reco-jets in semileptonic ttbar with Escheme+excl. kT 
#from Analysis.jetObservables.nSubProducer_QCD import nsubjettinessProducer
#from Analysis.jetObservables.nSubProducer_matched_TTSemilept import nsubjettinessProducer
#from Analysis.jetObservables.nSubProducer_matched_bbveto_TTSemilept import nsubjettinessProducer
# Output file name, and decide if this is MC or data


haddname = "jetObservables_nanoskim.root"

#if any( ['JetHT' in i for i in inputFiles() ] ):
#    isData = True
#else :
#    isData = False

print '---------------------------------------------------'
print 'Input files:'
print inputFiles()

jetmetCorrector = createJMECorrector(isMC=True, dataYear=2016, jesUncert="All", redojec=True)
fatJetCorrector = createJMECorrector(isMC=True, dataYear=2016, jesUncert="All", redojec=True, jetType = "AK8PFPuppi") 

modulesToRun = [puWeight_2016(),
        jetmetCorrector(),
        fatJetCorrector(),
        btagSF2016(),
        nsubjettinessProducer()
        ]
p1=PostProcessor( ".", inputFiles(), "", "keep_and_drop.txt",
                    modules=modulesToRun,
                    provenance=True,
                    fwkJobReport=True,
                    #jsonInput=runsAndLumis(),
                    #prefetch=True,
                    #maxEntries=20000,
                    haddFileName=haddname )

p1.run()
print "DONE"
os.system("ls -lR")


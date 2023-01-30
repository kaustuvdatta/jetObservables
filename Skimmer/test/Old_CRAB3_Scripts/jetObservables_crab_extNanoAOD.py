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

#from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016 #puAutoWeight_2016 #puWeight_2016     
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016         
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib import jetRecalib2016BCDAK8Puppi, jetRecalib2016EFAK8Puppi, jetRecalib2016GHAK8Puppi

# our module

#from Analysis.jetObservables.nSubProducer import nsubjettinessProducer
from Analysis.jetObservables.nSubProducer_axscan_data import nsubjettinessProducer
#from Analysis.jetObservables.nSubproducer_WP import nsubjettinessProducer     
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


#if any( ['2017B' in i for i in inputFiles() ] ):
#    is2017B = True
#else :
#    is2017B = False
#
#if any( ['2016' in i for i in inputFiles() ] ):
#    is2016 = True
#else :
#    is2016 = False
#
#ptcuts = "(nFatJet >=2 && FatJet_pt[0]>170 && FatJet_pt[1] > 170) "
#if not is2016:
#    trig60  = "( HLT_AK8PFJet60 ==1 && FatJet_pt[0] >  70 )"
#    trig80  = "( HLT_AK8PFJet80 ==1 && FatJet_pt[0] >  90 )"
#    trig140 = "( HLT_AK8PFJet140==1 && FatJet_pt[0] > 150 )"
#    trig200 = "( HLT_AK8PFJet200==1 && FatJet_pt[0] > 210 )"
#    trig260 = "( HLT_AK8PFJet260==1 && FatJet_pt[0] > 280 )"
#    trig320 = "( HLT_AK8PFJet320==1 && FatJet_pt[0] > 350 )"
#    trig400 = "( HLT_AK8PFJet400==1 && FatJet_pt[0] > 430 )"
#    trig450 = "( HLT_AK8PFJet450==1 && FatJet_pt[0] > 480 )"
#    trig500 = "( HLT_AK8PFJet500==1 && FatJet_pt[0] > 550 )"
#    if not is2017B:
#        trigHT = "(HLT_PFHT1050 == 1 || HLT_AK8PFHT850_TrimMass50 == 1 || HLT_AK8PFHT900_TrimMass50 == 1)"
#    else:
#        trigHT = "(HLT_PFHT1050 == 1)"
#    triglist = [ trig60, trig80, trig140, trig200, trig260, trig320, trig400, trig450, trig500, trigHT]
#else :
#    trig60  = "( HLT_PFJet60 ==1 && FatJet_pt[0] > 100 )"
#    trig80  = "( HLT_PFJet80 ==1 && FatJet_pt[0] > 160 )"
#    trig140 = "( HLT_PFJet140==1 && FatJet_pt[0] > 220 )"
#    trig200 = "( HLT_PFJet200==1 && FatJet_pt[0] > 310 )"
#    trig260 = "( HLT_PFJet260==1 && FatJet_pt[0] > 420 )"
#    trig320 = "( HLT_PFJet320==1 && FatJet_pt[0] > 510 )"
#    trig400 = "( HLT_PFJet400==1 && FatJet_pt[0] > 610 )"
#    trig450 = "( HLT_PFJet450==1 && FatJet_pt[0] > 720 )"
#    trigHT = "(HLT_PFHT900 == 1 || HLT_AK8PFHT650_TrimR0p1PT0p03Mass50 == 1 || HLT_AK8PFHT700_TrimR0p1PT0p03Mass50 == 1)"
#    triglist = [ trig60, trig80, trig140, trig200, trig260, trig320, trig400, trig450, trigHT]
#
#
#
#
#if isData  :
#    trigcuts = "(" +  ' || '.join( triglist ) + ")"
#    cuts = " && ".join( [ptcuts,trigcuts] )
#else :
#    cuts = ptcuts
#
#print 'Applying cuts : ' + cuts
#jetRecalib2016BCDAK8Puppi(),
#puWeight_2016(), btagSF2016()


jetmetCorrector = createJMECorrector(isMC=False, dataYear=2016, jesUncert="All", redojec=True)
fatJetCorrector = createJMECorrector(isMC=False, dataYear=2016, jesUncert="All", redojec=True, jetType = "AK8PFPuppi") 

modulesToRun = [
        jetmetCorrector(),
        fatJetCorrector(),
        btagSF2016(),
        nsubjettinessProducer()
        ]

p1=PostProcessor( ".", inputFiles(), "", "keep_and_drop.txt",
                    modules=modulesToRun,#modules=[nsubjettinessProducer()], 
                    provenance=True,
                    fwkJobReport=True,
                    #jsonInput=runsAndLumis(),
                    haddFileName=haddname )

p1.run()
print "DONE"
os.system("ls -lR")

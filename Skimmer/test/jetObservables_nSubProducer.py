#!/usr/bin/env python
import os, sys
import psutil
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puAutoWeight_2016, puWeight_2017, puAutoWeight_2017
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016, btagSF2017
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

import argparse

parser = argparse.ArgumentParser(description='Runs MEAnalysis')
parser.add_argument(
    '--sample',
    action="store",
    help="Sample to process",
    default='ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8'
)
parser.add_argument(
    '--numEvents',
    action="store",
    type=int,
    help="Number of events to process",
    default=1000000000000,
)
parser.add_argument(
    '--iFile',
    action="store",
    help="Input file (for condor)",
    default=""
)
parser.add_argument(
    '--oFile',
    action="store",
    help="Output file (for condor)",
    default=""
)
parser.add_argument(
    '--local',
    action="store_true",
    help="Run local or condor/crab"
)
parser.add_argument(
    '--year',
    action="store",
    help="year of data",
    choices=["2016", "2017", "2018"],
    default="2017",
    required=False
)
parser.add_argument(
    '--selection',
    action="store",
    help="Event selection",
    choices=["dijet", "Wtop"],
    default="Wtop",
    required=False
)
args = parser.parse_args(sys.argv[1:])
if args.sample.startswith(('/EGamma', '/Single', 'EGamma', 'Single', 'UL16_Single', '/UL16_Single', 'UL17_Single', '/UL17_Single', 'UL18_Single', '/UL18_Single', '/JetHT', 'JetHT', '/UL17_Jet', 'UL17_Jet' )) or ('EGamma' in args.iFile or 'SingleMuon' in args.iFile or ('JetHT' in args.iFile)):
    isMC = False
    print "sample is data"
else: isMC = True

### General selections:
PV = "(PV_npvsGood>0)"
METFilters = "( (Flag_goodVertices==1) && (Flag_globalSuperTightHalo2016Filter==1) && (Flag_HBHENoiseFilter==1) && (Flag_HBHENoiseIsoFilter==1) && (Flag_EcalDeadCellTriggerPrimitiveFilter==1) && (Flag_BadPFMuonFilter==1) )"
if not isMC: METFilters = METFilters + ' && (Flag_eeBadScFilter==1)'

if args.selection.startswith('dijet'):
    Triggers =  '( (HLT_PFJet140==1) || (HLT_PFJet200==1) || (HLT_PFJet260==1) || (HLT_PFJet320==1) || (HLT_PFJet400==1) || (HLT_PFJet450==1) || (HLT_PFJet500==1) || (HLT_PFJet550==1) )'
else:
    Triggers = '(HLT_Mu50==1)'
#if args.year.startswith('2016'): Triggers = ...
#elif args.year.startswith('2017'): Triggers =  ...
#elif args.year.startswith('2018'): Triggers = ...

cuts = PV + " && " + METFilters + " && " + Triggers

systSources = [ '_jesTotal', '_jer', '_pu' ] if isMC else []

### Lepton scale factors
LeptonSF = {
    '2016' : {
        'muon' : {
            'Trigger' : [ "EfficienciesAndSF_RunBtoF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio" ],
            'ID' : [ "MuonID_2016_RunBCDEF_SF_ID.root", "NUM_TightID_DEN_genTracks_eta_pt", False ],       ### True: X:pt Y:eta
            'ISO' : [ "MuonID_2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt", True ],
        },
        'electron' : {
            'Trigger' : [ "TriggerSF_Run2016All_v1.root", "Ele27_WPTight_Gsf" ],
            'ID' : [ "2016LegacyReReco_ElectronTight_Fall17V2.root", "EGamma_SF2D", False ],
            'ISO' : [ "EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root", "EGamma_SF2D", False ],
        },
    },
    '2017' : {
        'muon' : {
            'Trigger' : [ "EfficienciesAndSF_RunBtoF_Nov17Nov2017.root", "Mu50_PtEtaBins/pt_abseta_ratio" ],
            'ID' : [ "Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt", False ],
            'ISO' : [ "Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt", False ],
        },
        'electron' : {
            'Trigger' : [ "SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5.0.histo.root", "SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF" ],
            'ID' : [ "egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "EGamma_SF2D", False ],
            'ISO' : [ "2017_ElectronTight.root", "EGamma_SF2D", False ],
        },
    },
    '2018' : {
        'muon' : {
            'Trigger' : [ "EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root", "IsoMu24_PtEtaBins/pt_abseta_ratio" ],
            'ID' : [ "MuonID_2018_RunABCD_SF_ID.root", "NUM_TightID_DEN_TrackerMuons_pt_abseta", True ],
            'ISO' : [ "MuonID_2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", True ],
        },
        'electron' : {
            'Trigger' : [ "SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v5.0.root", "SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF" ],
            'ID' : [ "egammaEffi.txt_EGM2D_updatedAll.root", "EGamma_SF2D", False ],
            'ISO' : [ "2018_ElectronTight.root", "EGamma_SF2D", False ],
        },
    },
}


#### Modules to run
if not args.selection.startswith('dijet'): jetmetCorrector = createJMECorrector(isMC=isMC, applySmearing=False, dataYear='UL'+args.year, jesUncert="All")
fatJetCorrector = createJMECorrector(isMC=isMC, applySmearing=False, dataYear='UL'+args.year, jesUncert="All", jetType = "AK8PFPuppi") #redojec=True,??? also do I need to fill the runPeriod var with something beyond the default value, currently it is always set to the default of run 'B' irrespective of year?

modulesToRun = []
if isMC:
    if args.year=='2018':
        modulesToRun.append( puWeight_2018() )
        print "###Running with btag SF calc.###"
        if not args.selection.startswith('dijet'): modulesToRun.append( btagSF2018() )
    if args.year=='2017':
        modulesToRun.append( puAutoWeight_2017() )
        print "###Running with btag SF calc.###"
        if not args.selection.startswith('dijet'): modulesToRun.append( btagSF2017() )
    if args.year=='2016':
        modulesToRun.append( puWeight_2016() )
        print "Running with btag SF calc."
        if not args.selection.startswith('dijet'): modulesToRun.append( btagSF2016() )
modulesToRun.append( fatJetCorrector() )
if not args.selection.startswith('dijet'): modulesToRun.append( jetmetCorrector() )

# our module
if args.selection.startswith('dijet'):
    if args.local: triggerFile = "../../TriggerEfficiencies/test/Rootfiles/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017ALL.root"
    #if args.local: triggerFile = "/eos/home-a/algomez/tmpFiles/jetObservables/triggerPrescales/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017B.pkl"
    else: triggerFile = 'triggerEfficiencies_histograms_MiniAOD_JetHTRun2017B.pkl'
    from jetObservables.Skimmer.nSubProducer_dijetSel import nSubProd
    modulesToRun.append( nSubProd( sysSource=systSources, leptonSF=LeptonSF[args.year], isMC=isMC, triggerFile=triggerFile ) )
else:
    from jetObservables.Skimmer.nSubProducer_WtopSel import nSubProd
    modulesToRun.append( nSubProd( sysSource=systSources, leptonSF=LeptonSF[args.year], isMC=isMC ) )


#### Make it run
p1=PostProcessor(
        '.', (inputFiles() if not args.iFile else args.iFile.split(',')),
        cut          = cuts,
        outputbranchsel   = "keep_and_drop.txt",
        modules      = modulesToRun,
        provenance   = True,
        #jsonInput   = runsAndLumis(),
        maxEntries   = args.numEvents,
        prefetch     = args.local,
        longTermCache= args.local,
        fwkJobReport = True,
        haddFileName = "jetObservables_"+args.selection+"_nanoskim.root" if args.local else 'jetObservables_nanoskim.root',
        histFileName = "jetObservables_"+args.selection+"_histograms.root" if args.local else 'jetObservables_histograms.root',
        histDirName  = 'jetObservables',
        )
p1.run()
print "DONE"
#if not args.local: os.system("ls -lR")

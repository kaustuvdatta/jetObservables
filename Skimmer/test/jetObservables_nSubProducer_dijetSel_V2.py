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

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer, puAutoWeight_UL2016, puAutoWeight_UL2017, puAutoWeight_UL2018
#from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016, btagSF2017
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
    '--runEra',
    action="store",
    help="Run era for data",
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
    choices=["2016_preVFP", "2016", "2017", "2018"],
    default="2017",
    required=False
)
parser.add_argument(
    '--selection',
    action="store",
    help="Event selection",
    choices=["dijet"], #, "Wtop"],
    default="dijet",
    required=False
)
parser.add_argument(
    '--onlyUnc',
    action="store",
    default='',
    help="Run only specific uncertainty variations"
)
parser.add_argument(
    '--onlyTrees',
    action="store_true",
    help="Do not save histograms, only trees"
)
parser.add_argument(
    '--isSigMC',
    action="store_true",
    #default='',
    help="Save branches with sys+puWeights for signal/main MC else store basic reco/gen branches"
)


args = parser.parse_args(sys.argv[1:])
print(args)
if args.sample.startswith(('/EGamma', '/Single', 'EGamma', 'Single', 'UL16_Single', '/UL16_Single', 'UL17_Single', '/UL17_Single', 'UL18_Single', '/UL18_Single', '/JetHT', 'JetHT', '/UL17_Jet', 'UL17_Jet' )) or ('EGamma' in args.iFile or 'SingleMuon' in args.iFile or ('JetHT' in args.iFile)):
    isMC = False
    print "sample is data"
else: isMC = True

### General selections:
PV = "(PV_npvsGood>0)"
if not('2016' in args.year): METFilters = "( (Flag_goodVertices==1) && (Flag_globalSuperTightHalo2016Filter==1) && (Flag_HBHENoiseFilter==1) && (Flag_HBHENoiseIsoFilter==1) && (Flag_EcalDeadCellTriggerPrimitiveFilter==1) && (Flag_BadPFMuonFilter==1) && (Flag_BadPFMuonDzFilter==1) )"
else: METFilters = "( (Flag_goodVertices==1) && (Flag_globalSuperTightHalo2016Filter==1) && (Flag_HBHENoiseFilter==1) && (Flag_HBHENoiseIsoFilter==1) && (Flag_EcalDeadCellTriggerPrimitiveFilter==1) && (Flag_BadPFMuonFilter==1)  )"

if not isMC: METFilters = METFilters + ' && (Flag_eeBadScFilter==1)'

if not('2016' in args.year): 
    Triggers =  '( (HLT_AK8PFJet80==1) || (HLT_AK8PFJet140==1) || (HLT_AK8PFJet200==1) || (HLT_AK8PFJet260==1) || (HLT_AK8PFJet320==1) || (HLT_AK8PFJet400==1) || (HLT_AK8PFJet450==1) || (HLT_AK8PFJet500==1 ) || (HLT_AK8PFJet550==1) )' 
else:
    Triggers =  '( (HLT_AK8PFJet80==1) || (HLT_AK8PFJet140==1) || (HLT_AK8PFJet200==1) || (HLT_AK8PFJet260==1) || (HLT_AK8PFJet320==1) || (HLT_AK8PFJet400==1) || (HLT_AK8PFJet450==1) || (HLT_AK8PFJet500==1) )'

cuts = PV + " && " + METFilters + " && " + Triggers

#systSources = ['_jesTotal', '_jer', '_puWeight'] if isMC else []   ######### NEEDS TO BE REVIEWED FOR WTOP
#if args.selection.startswith('dijet'):
    #systSources = [ '_puWeight', '_isrWeight', '_fsrWeight', '_pdfWeight' ] if args.sample.startswith('QCD_HT') else []
systSources = [] #['_jesTotal', '_jer', '_puWeight', '_isrWeight', '_fsrWeight', '_pdfWeight' ] if args.sample.startswith('QCD_HT') else []
if args.onlyUnc: 
    if '_jesCorr' in args.onlyUnc:
        systSources = [ '_jesAbsolute', '_jesBBEC1', '_jesEC2', '_jesFlavorQCD', '_jesHF', '_jesRelativeBal']
    elif '_jesUncorr' in args.onlyUnc:
        if '2018' in args.year: systSources=systSources+['_jesAbsolute_2018', '_jesBBEC1_2018', '_jesEC2_2018', '_jesHF_2018', '_jesRelativeSample_2018']
        elif '2017' in args.year: systSources=systSources+['_jesAbsolute_2017', '_jesBBEC1_2017', '_jesEC2_2017', '_jesHF_2017', '_jesRelativeSample_2017']
        elif '2016' in args.year: systSources=systSources+['_jesAbsolute_2016', '_jesBBEC1_2016', '_jesEC2_2016', '_jesHF_2016', '_jesRelativeSample_2016']

    else: systSources = [ args.onlyUnc ]


### Lepton scale factors
LeptonSF = {

    '2016_preVFP' : {
        'muon' : {
            'Trigger' : [ "MuonSF_ULRunII.root", "UL16_preVFP_Trigger", False ],
            'ID' : [ "MuonSF_ULRunII.root", "UL16_preVFP_ID", False ],
            'ISO' : [ "MuonSF_ULRunII.root", "UL16_preVFP_ISO", False ],
            'RecoEff' : [ "MuonSF_ULRunII.root", "UL16_preVFP_Reco", False ],
        },

    },

    '2016' : {
        'muon' : {
            'Trigger' : [ "MuonSF_ULRunII.root", "UL16_postVFP_Trigger", False ],
            'ID' : [ "MuonSF_ULRunII.root", "UL16_postVFP_ID", False ],
            'ISO' : [ "MuonSF_ULRunII.root", "UL16_postVFP_ISO", False ],
            'RecoEff' : [ "MuonSF_ULRunII.root", "UL16_postVFP_Reco", False ],
        },
    },

    '2017' : {
        'muon' : {
            'Trigger' : [ "MuonSF_ULRunII.root", "UL17_Trigger", False ],
            'ID' : [ "MuonSF_ULRunII.root", "UL17_ID", False ],
            'ISO' : [ "MuonSF_ULRunII.root", "UL17_ISO", False ],
            'RecoEff' : [ "MuonSF_ULRunII.root", "UL17_Reco", False ],
        },

    },
    '2018' : {
        'muon' : {
            'Trigger' : [ "MuonSF_ULRunII.root", "UL18_Trigger", False ],
            'ID' : [ "MuonSF_ULRunII.root", "UL18_ID", False ],
            'ISO' : [ "MuonSF_ULRunII.root", "UL18_ISO", False ],
            'RecoEff' : [ "MuonSF_ULRunII.root", "UL18_Reco", False ],
        },

    },
}


#### Modules to run
jesUncert = 'Merged' if any('_jes' in iunc for iunc in systSources) else 'Total'
print(jesUncert, systSources)
fatJetCorrector = createJMECorrector(isMC=isMC, applySmearing=False, dataYear='UL'+args.year, jesUncert=jesUncert, jetType = "AK8PFPuppi", runPeriod=args.runEra, applyHEMfix=(True if args.year.startswith('2018') else False))

modulesToRun = []
modulesToRun.append( fatJetCorrector() )

if isMC:
    if args.year=='2018':
        modulesToRun.append( puAutoWeight_UL2018() )
    if args.year=='2017':
        modulesToRun.append( puAutoWeight_UL2017() )
    if '2016' in args.year:
        modulesToRun.append( puAutoWeight_UL2016() )
        
# our module
if args.selection.startswith('dijet'):
    from jetObservables.Skimmer.nSubProducer_dijetSel_RunIISummer20UL_CustomTrees_V2 import nSubProd
    modulesToRun.append( nSubProd( sysSource=systSources, isMC=isMC, year=args.year, onlyUnc=args.onlyUnc, onlyTrees=args.onlyTrees, isSigMC=True if args.isSigMC else False ) )

#### Make it run
p1=PostProcessor(
        '.', (inputFiles() if not args.iFile else args.iFile.split(',')),
        cut          = cuts,
        outputbranchsel   = "keep_and_drop_dijet.txt" if args.selection.startswith('dijet') else "keep_and_drop.txt",
        modules      = modulesToRun,
        provenance   = True,
        #jsonInput   = runsAndLumis(),
        maxEntries   = args.numEvents,
        prefetch     = args.local,
        longTermCache= args.local,
        fwkJobReport = True,
        haddFileName = "jetObservables_"+args.selection+args.year+args.onlyUnc+("sigMC" if args.isSigMC else "")+("isMC" if isMC else "isData")+"_nanoskim_V2.root" if args.local else 'jetObservables_nanoskim.root',
        histFileName = "jetObservables_"+args.selection+args.year+args.onlyUnc+("sigMC" if args.isSigMC else "")+("isMC" if isMC else "isData")+"_histograms_V2.root" if args.local else 'jetObservables_histograms.root',
        histDirName  = 'jetObservables',
        )
p1.run()
print "DONE"
#if not args.local: os.system("ls -lR")

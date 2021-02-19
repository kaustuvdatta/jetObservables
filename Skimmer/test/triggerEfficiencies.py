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


# our module
from jetObservables.Skimmer.triggerEff_dijet import triggerEfficiencies

import argparse

parser = argparse.ArgumentParser(description='Runs TriggerEfficiencies')
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

args = parser.parse_args(sys.argv[1:])
if args.sample.startswith(('/EGamma', '/Single', 'EGamma', 'Single', 'UL16_Single', '/UL16_Single', 'UL17_Single', '/UL17_Single', 'UL18_Single', '/UL18_Single', '/JetHT', 'JetHT' )) or ('EGamma' in args.iFile or 'SingleMuon' in args.iFile or ('JetHT' in args.iFile)):
    isMC = False
    print "sample is data"
else: isMC = True

modulesToRun = []
modulesToRun.append( triggerEfficiencies( year=args.year ) )

#### Make it run
p1=PostProcessor(
        '.', (inputFiles() if not args.iFile else [args.iFile]),
        outputbranchsel   = "keep_and_drop.txt",
        modules      = modulesToRun,
        provenance   = True,
        jsonInput   = runsAndLumis(),
        maxEntries   = args.numEvents,
        prefetch     = args.local,
        longTermCache= args.local,
        fwkJobReport = True,
        haddFileName = 'triggerEfficiencies.root',
        histFileName = "triggerEfficiencies_histograms.root",  #if args.local else 'triggerEfficiencies_histograms.root',
        histDirName  = 'triggerEfficiencies',
        )
p1.run()
print "DONE"
#if not args.local: os.system("ls ")

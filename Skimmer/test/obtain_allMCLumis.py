import os, sys

#from six.moves.urllib.parse import urlparse, urlencode
from FWCore.PythonUtilities.LumiList import LumiList
#from CRABClient.ClientUtilities import DBSURLS, LOGLEVEL_MUTE, colors
#from CRABClient.ClientUtilities import execute_command, getUserProxy
#from CRABClient.ClientExceptions import ClientException, UsernameException
from crab_UserUtils_hacked import config, getLumiListInValidFiles
#from WMCore.Configuration import Configuration

import argparse

parser = argparse.ArgumentParser(description='Runs calculation of missing lumis in processed PFNano files')
parser.add_argument(
    '--inputPFNano', '-i',
    action="store",
    help="Input sample to check for list of possible input lumis",
    default='/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER'
)

args = parser.parse_args(sys.argv[1:])

print (args)

def retrieve_availableLumis(inputDatasetTag=''):
    #data file name example:/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD
    # config = config()
    inputLumis = getLumiListInValidFiles(dataset=inputDatasetTag,dbsurl='phys03')

    outputJSON = "%s_allValidLumis.json"%(inputDatasetTag.split('/')[1]+inputDatasetTag.split('/')[2].split('-')[1])
    print ([x for x in inputLumis])
    inputLumis.writeJSON(outputJSON)

    return outputJSON

if (not(args.inputPFNano.startswith('/')) or not('PF' in args.inputPFNano)):
    print ("Not a valid PFNano sample; check and try again")


out_JSON_name=retrieve_availableLumis(inputDatasetTag=args.inputPFNano)

print ("Lumis to be run over for %s are written to %s"%(args.inputPFNano.split('/')[1],out_JSON_name))

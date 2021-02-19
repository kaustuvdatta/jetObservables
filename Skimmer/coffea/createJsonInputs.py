#####################################
##### This script sum the genWeights from nanoAOD
##### To run: python computeGenWeights.py -d NAMEOFDATASET
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import numpy as np
import json
from dbs.apis.dbsClient import DbsApi  ## talk to DBS to get list of files in this dataset
dbsGlobal = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
dbsPhys03 = DbsApi('https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader')

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--datasets", action='store', dest="datasets", default="all", help="Name of dataset to process" )
parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
parser.add_argument("-y", "--year", action='store', choices=[ '2016', '2017', '2018' ],  default="2017", help="Version" )

try: args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

dictSamples = {}
dictSamples['SingleMuon_Run2017B'] = '/SingleMuon/Run2017B-UL2017_MiniAODv1_NanoAODv2-v1/NANOAOD'
dictSamples['SingleMuon_Run2017C'] = '/SingleMuon/Run2017C-UL2017_02Dec2019-v1/NANOAOD'
dictSamples['SingleMuon_Run2017D'] = '/SingleMuon/Run2017D-UL2017_MiniAODv1_NanoAODv2-v1/NANOAOD'
dictSamples['SingleMuon_Run2017E'] = '/SingleMuon/Run2017E-UL2017_MiniAODv1_NanoAODv2-v1/NANOAOD'
dictSamples['SingleMuon_Run2017F'] = '/SingleMuon/Run2017F-UL2017_MiniAODv1_NanoAODv2-v2/NANOAOD'


### To choose dataset from dictSamples
processingSamples = {}
for sam in dictSamples:
    if 'all' in args.datasets:
        processingSamples[ sam ] = dictSamples[ sam ]
    elif sam.startswith( args.datasets ):
        processingSamples[ sam ] = dictSamples[ sam ]
if len(processingSamples)==0: print('No sample found. \n Have a nice day :)')

print(processingSamples)
allfiles = {}
for isample, jsample  in processingSamples.items():

    fileDictList = ( dbsPhys03 if jsample.endswith('USER') else dbsGlobal).listFiles(dataset=jsample,validFileOnly=1)
    # DBS client returns a list of dictionaries, but we want a list of Logical File Names
    #allfiles = [ "root://cms-xrd-global.cern.ch/"+dic['logical_file_name'] for dic in fileDictList ]
    allfiles[isample] = [ "root://xrootd-cms.infn.it/"+dic['logical_file_name'] for dic in fileDictList ]
    print ("dataset %s has %d files" % (jsample, len(allfiles)))

outputFile='samples_'+args.year+'.json'
with open( outputFile, 'w') as fp:
    json.dump(allfiles, fp, indent=4, sort_keys=True)
print('Output file:', outputFile)

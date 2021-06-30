#####################################
##### This script sum the genWeights from nanoAOD
##### To run: python computeGenWeights.py -d NAMEOFDATASET
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from datasets import dictSamples, checkDict
import numpy as np
import ROOT
import uproot
import concurrent.futures
#from root_numpy import root2array, tree2array
from dbs.apis.dbsClient import DbsApi  ## talk to DBS to get list of files in this dataset
dbsGlobal = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
dbsPhys03 = DbsApi('https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader')


#####################################
def computeGenWeights( inputFiles, genEventSumw = 0, genEventSumw2 = 0, genEventCount = 0 ):
    """docstring for computeGenWeights"""

    #executor = concurrent.futures.ThreadPoolExecutor()
    notProcessedFiles = []
    for i, fi in enumerate(inputFiles):
        try:
            #bl = uproot.open(fi, executor=executor)['Runs']
            bl = uproot.open(fi)['Runs']
            genEventSumw += bl.array("genEventSumw")[0]
            genEventSumw2 += bl.array("genEventSumw2")[0]
            genEventCount += bl.array("genEventCount")[0]
            print("file",i,genEventSumw)
        except (OSError, IOError) as e:
            notProcessedFiles.append(fi)
            print('file not processed',fi)
            continue

    print ('Total number of genEventCount in sample: ', genEventCount)
    print ('Total number of genEventSumw in sample: ', genEventSumw)
    print ('Total number of genEventSumw2 in sample: ', genEventSumw2)
    print ('List of not processed files\n', notProcessedFiles)



#####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datasets", action='store', dest="datasets", default="ttHTobb", help="Name of dataset to process" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument("-l", "--local", action='store_true', help="Run local or through xrootd" )
    parser.add_argument("-y", "--year", action='store', choices=[ '2016', '2017', '2018' ],  default="2017", help="Version" )
    parser.add_argument( '--selection', action="store", help="Event selection", choices=["dijet", "Wtop"], default="dijet" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    ### To choose dataset from dictSamples
    processingSamples = {}
    for sam in dictSamples:
        if sam.startswith( args.datasets ) | args.datasets.startswith('all'):
            if checkDict( sam, dictSamples )['selection'] != args.selection: continue
            processingSamples[ sam ] = checkDict( sam, dictSamples )[args.year]['skimmer'][0]

    if len(processingSamples)==0: print ('No sample found. \n Have a nice day :)')

    print(processingSamples)
    for isample, jsample  in processingSamples.items():
        if not jsample:
            print('|--------> Sample '+isample+' does not have skimmer sample in datasets.py')
            continue

        prefix = '/pnfs/psi.ch/cms/trivcat' if args.local else "root://xrootd-cms.infn.it/"

        ### Create a list from the dataset
        if isinstance( jsample, list ):
            allfiles = []
            for jsam in jsample:
                fileDictList = ( dbsPhys03 if jsam.endswith('USER') else dbsGlobal).listFiles(dataset=jsam,validFileOnly=1)
                tmpfiles = [ prefix+dic['logical_file_name'] for dic in fileDictList ]
                #print tmpfiles
                allfiles = allfiles + tmpfiles
            #print allfiles
        else:
            fileDictList = ( dbsPhys03 if jsample.endswith('USER') else dbsGlobal).listFiles(dataset=jsample,validFileOnly=1)
            # DBS client returns a list of dictionaries, but we want a list of Logical File Names
            #allfiles = [ {"root://cms-xrd-global.cern.ch/"+dic['logical_file_name']: "Runs"} for dic in fileDictList ]
            allfiles = [ prefix+dic['logical_file_name'] for dic in fileDictList  ]
            #allfiles = [ prefix+dic['logical_file_name'] for dic in fileDictList if 'nanoskim_1' in dic['logical_file_name'] ]
        print ("dataset %s has %d files" % (jsample, len(allfiles)))

        computeGenWeights( allfiles )


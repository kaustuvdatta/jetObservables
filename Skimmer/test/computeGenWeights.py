#####################################
##### This script sum the genWeights from nanoAOD
##### To run: python computeGenWeights.py -d NAMEOFDATASET
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import numpy as np
import ROOT
import uproot
from root_numpy import root2array, tree2array
from dbs.apis.dbsClient import DbsApi  ## talk to DBS to get list of files in this dataset
dbsGlobal = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
dbsPhys03 = DbsApi('https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader')


#####################################
def computeGenWeights( inputFiles, outputName ):
    """docstring for computeGenWeights"""

#    ### Open root files
#    intree = ROOT.TChain("Runs")
#    intree2 = ROOT.TChain("Events")
#    if isinstance(inputFiles, list):
#        for iFile in inputFiles:
#            intree.Add(iFile)
#            intree2.Add(iFile)
#    else:
#        intree.Add(inputFiles)
#        intree2.Add(inputFiles)
#    print intree.GetEntries()
#
#    ### Convert root tree to numpy array, applying cuts
#    arrays = tree2array( intree,
#                            branches=['genEventCount_', 'genEventSumw_', 'genEventSumw2_'],
#                            selection='',
#                            #stop=1000,  #### to test only, run only 1000 events
#                            )
#    arrays2 = tree2array( intree2,
#                            branches=['genWeight'],
#                            selection='',
#                            #stop=1000,  #### to test only, run only 1000 events
#                            )
#
#    #tmp = arrays['genWeight']/arrays['genWeight'][0]
#    #print 'Total number of events in sample: ', intree.GetEntries()
#    #print 'Event weights per file: ', arrays['genEventSumw']
#    print 'Total number of genEventCount in sample: ', sum(arrays['genEventCount_'])
#    print 'Total number of genEventSumw in sample: ', sum(arrays['genEventSumw_'])
#    print 'Total number of genEventSumw2 in sample: ', sum(arrays['genEventSumw2_'])
#    print 'Total sum of genWeights in sample: ', sum(arrays2['genWeight'])
#    print 'Total sum of sign(genWeights) in sample: ', sum(np.sign(arrays2['genWeight']))

    genEventSumw = 0
    genEventSumw2 = 0
    genEventCount = 0
    genWeight = 0
    for fi in inputFiles:
        ff = uproot.open(fi)
        bl = ff.get("Runs")
        genEventSumw += bl.array("genEventSumw").sum()
        genEventSumw2 += bl.array("genEventSumw2").sum()
        genEventCount += bl.array("genEventCount").sum()
        bl = ff.get("Events")
        genWeight += bl.array("genWeight").sum()

    print 'Total number of genEventCount in sample: ', genEventCount
    print 'Total number of genEventSumw in sample: ', genEventSumw
    print 'Total number of genEventSumw2 in sample: ', genEventSumw2
    print 'Total sum of genWeights in sample: ', genWeight


#####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datasets", action='store', dest="datasets", default="ttHTobb", help="Name of dataset to process" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument("-y", "--year", action='store', choices=[ '2016', '2017', '2018' ],  default="2017", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dictSamples = {}
    dictSamples['QCDPt170to300'] = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER'
    dictSamples['ST_s-channel'] = '/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER'

    ### To choose dataset from dictSamples
    processingSamples = {}
    for sam in dictSamples:
        if 'all' in args.datasets:
            processingSamples[ sam ] = dictSamples[ sam ]
        elif sam.startswith( args.datasets ):
            processingSamples[ sam ] = dictSamples[ sam ]
    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    print(processingSamples)
    for isample, jsample  in processingSamples.items():

        ### Create a list from the dataset
        if isinstance( jsample, list ):
            allfiles = []
            for jsam in jsample:
                fileDictList = ( dbsPhys03 if jsam.endswith('USER') else dbsGlobal).listFiles(dataset=jsam,validFileOnly=1)
                tmpfiles = [ "root://xrootd-cms.infn.it/"+dic['logical_file_name'] for dic in fileDictList ]
                #print tmpfiles
                allfiles = allfiles + tmpfiles
            #print allfiles
        else:
            fileDictList = ( dbsPhys03 if jsample.endswith('USER') else dbsGlobal).listFiles(dataset=jsample,validFileOnly=1)
            # DBS client returns a list of dictionaries, but we want a list of Logical File Names
            #allfiles = [ "root://cms-xrd-global.cern.ch/"+dic['logical_file_name'] for dic in fileDictList ]
            allfiles = [ "root://xrootd-cms.infn.it/"+dic['logical_file_name'] for dic in fileDictList ]
        print ("dataset %s has %d files" % (jsample, len(allfiles)))

        computeGenWeights( allfiles, isample )

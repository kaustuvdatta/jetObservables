#####################################
##### This script computes cross section
##### based on this script: https://github.com/cms-sw/genproductions/blob/master/Utilities/calculateXSectionAndFilterEfficiency/compute_cross_section.py
##### To run: python computeXSection.py -d NAMEOFDATASET -y YEAR
##### The datasets are stored in datasets.py
##### and needs to clone this package first: https://github.com/cms-sw/genproductions
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from datasets import dictSamples, checkDict
from optparse import OptionParser
import commands
import re
import datetime
from time import sleep

def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False

def compute_cross_section():
    """docstring for comput_cross_section"""

    debug = str_to_bool(str(args.debug))
    mcm = str_to_bool(str(args.mcm))
    skipexisting = str_to_bool(str(args.skipexisting))
    if debug: print 'args.mcm',args.mcm,'mcm',mcm,'debug',debug
    if debug: print 'debug is True!'
    if debug and mcm: print 'mcm is True!'

    if debug:
        print
        print 'RUNNING PARAMS: '
        print '                debug                 = ' + str(debug)
        print '                dataset               = ' + args.inputdataset
        print '                MC campaign           = ' + args.campaign
        print '                Datatier              = ' + args.datatier
        print '                number of events      = ' + str(args.events)
        print '                use McM prepID        = ' + str(mcm)
        print '                skipexisting          = ' + str(skipexisting)
        print

    das_cmd = "/cvmfs/cms.cern.ch/common/dasgoclient"

    # if mcm is specified, retrieve dataset name from prepID:
    if mcm:
        if "/" in str(args.inputdataset):
            print "not a McM prepID format, please check"
            sys.exit(1)
        # load McM
        sys.path.append('/afs/cern.ch/cms/PPD/PdmV/tools/McM/')
        from rest import McM
        mcm = McM()
        # retrieve request with given prepid
        temp = sys.stdout
        f = open('/dev/null', 'w')
        sys.stdout = f
        request = mcm.get('requests', str(args.inputdataset))
        sys.stdout = temp
        if debug: print 'request prepid',request['prepid']
        # search dataset name as returned by mcm
        dataset_used = str(request['output_dataset'][0])
        primary_dataset_name = dataset_used.split('/')[1]
    else:
        # search dataset name as name + campaign + datatier
        primary_dataset_name = args.inputdataset.split('/')[1]
        command=das_cmd+" --limit=0 --query=\"dataset dataset="+args.inputdataset+'"'
        dataset_used = commands.getstatusoutput(command)[1].split("\n")
        if debug: print 'command',command,'\n'
        dataset_used = [x.strip() for x in dataset_used][0]

    if skipexisting and os.path.isfile("xsec_"+primary_dataset_name+".log"):
        print "xsec_"+primary_dataset_name+".log existing and NO skipexisting asked, skipping"
        # sys.exit(0)
    else:
        if debug: print 'dataset_used',dataset_used
        if debug: print 'primary_dataset_name',primary_dataset_name,'\n'
        # pick up only the first dataset of the list
        if debug: print 'dataset_used',dataset_used
        # retrieve filelist
        command=das_cmd+" --limit=50 --query=\"file dataset="+dataset_used+"\" "
        if debug: print 'command',command
        filelist_used = "/store"+commands.getstatusoutput(command)[1].replace("\n",",").split("/store",1)[1]
        if debug:
            print 'filelist_used',filelist_used.split(',')[0]
            filelist_used = filelist_used.split(',')[0]
        # compute cross section
        command = 'cmsRun calculateXSectionAndFilterEfficiency/genXsec_cfg.py inputFiles=\"'+filelist_used+'\" maxEvents='+str(args.events)+" 2>&1 | tee calculateXSectionAndFilterEfficiency/xsec_"+primary_dataset_name+".log"
        print command
        os.system(command)

#####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--datasets", action='store', dest="datasets", default="ttHTobb", help="Name of dataset to process" )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )
    parser.add_argument("-y", "--year", action='store', choices=[ '2016', '2017', '2018' ],  default="2017", help="Version" )
    parser.add_argument('-n', '--events'        , dest="events",        default=int(1e6),         help='number of events to calculate the cross section')
    parser.add_argument(      '--mcm'           , dest="mcm",           default=False,            help='use McM prepID instead of dataset name')
    parser.add_argument('-s', '--skipexisting'  , dest="skipexisting",  default=False,            help='skipexisting existing output files containing xsec results')
    parser.add_argument(      '--debug'         , dest="debug",         default=False,            help='use debug options (debug couts...)')

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)


    ### To choose dataset from dictSamples
    processingSamples = {}
    for sam in dictSamples:
        if sam.startswith( args.datasets ) | args.datasets.startswith('all'):
            processingSamples[ sam ] = checkDict( sam, dictSamples )[args.year]['miniAOD'][0]

    if len(processingSamples)==0: print 'No sample found. \n Have a nice day :)'

    print(processingSamples)
    for isample, jsample  in processingSamples.items():
        if not jsample:
            print('|--------> Sample '+isample+' does not have miniAOD sample in datasets.py')
            continue
        args.inputdataset = jsample
        args.datatier = 'MINIAODSIM'
        compute_cross_section()


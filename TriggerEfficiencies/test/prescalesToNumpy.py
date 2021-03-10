#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import numpy as np
import ROOT
import uproot
from root_numpy import root2array, tree2array


#####################################
def computeGenWeights( inputFiles, outputName ):
    """docstring for computeGenWeights"""

    ff = uproot.open(inputFiles)
    bl = ff.get("TriggerEfficiencies/TriggerObjects")
    data = bl.pandas.df(['run', 'lumi', 'event'])
    data = data.assign( AK8PFJet80 = bl.array('triggerPrescales')[:,0] )
    data = data.assign( AK8PFJet140 = bl.array('triggerPrescales')[:,1] )
    data = data.assign( AK8PFJet200 = bl.array('triggerPrescales')[:,2] )
    data = data.assign( AK8PFJet260 = bl.array('triggerPrescales')[:,3] )
    data = data.assign( AK8PFJet320 = bl.array('triggerPrescales')[:,4] )
    data = data.assign( AK8PFJet400 = bl.array('triggerPrescales')[:,5] )
    data = data.assign( AK8PFJet450 = bl.array('triggerPrescales')[:,6] )
    data = data.assign( AK8PFJet500 = bl.array('triggerPrescales')[:,7] )
    data = data.assign( AK8PFJet550 = bl.array('triggerPrescales')[:,8] )
    #print data[ (data.run==297050) & ( data.lumi==68 ) ]
    data.to_pickle( outputName )



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

    allfiles = 'Rootfiles/triggerEfficiencies_histograms_MiniAOD_JetHTRun2017E.root'
    outputFile = allfiles.split('/')[1].replace('root','pkl')
    computeGenWeights( allfiles, outputFile )

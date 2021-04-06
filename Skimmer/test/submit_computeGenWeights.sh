#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc820
export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Skimmer/test/
eval `scramv1 runtime -sh`

echo "Running: python computeGenWeights.py -d ${1} -y ${2}"
time python computeGenWeights.py -d ${1} -y ${2}


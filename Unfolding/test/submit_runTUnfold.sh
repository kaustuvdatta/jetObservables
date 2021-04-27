#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc820
export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/
eval `scramv1 runtime -sh`
env

echo "Running: python runTUnfold.py -v v02 --only ${1}}"
time python runTUnfold.py -v v02 --only ${1} --QCDHT -p MC
time python runTUnfold.py -v v02 --only ${1} --QCDHT -p MC --noUnc
time python runTUnfold.py -v v02 --only ${1} --QCDHT -p MC --selfClosure
time python runTUnfold.py -v v02 --only ${1} --QCDHT -p MC --selfClosure --noUnc
time python runTUnfold.py -v v02 --only ${1} --QCDHT
time python runTUnfold.py -v v02 --only ${1} --QCDHT --noUnc


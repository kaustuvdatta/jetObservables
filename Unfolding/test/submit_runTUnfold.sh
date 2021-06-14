#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc820
export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/
eval `scramv1 runtime -sh`
env
echo "Running: python runTUnfold.py -v v03 --only ${1} -y ${2} --outputFolder ${3}"
time python runTUnfold.py -v v03 --only ${1} -p MCSelfClosure -y ${2} --outputFolder ${3} --main Ptbin
time python runTUnfold.py -v v03 --only ${1} -p MCSelfClosure -y ${2} --outputFolder ${3} --main HTbin
time python runTUnfold.py -v v03 --only ${1} -p MCSelfClosure -y ${2} --outputFolder ${3} --main herwig
time python runTUnfold.py -v v03 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} --main Ptbin --alt herwig
time python runTUnfold.py -v v03 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} --main Ptbin --alt HTbin
time python runTUnfold.py -v v03 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} --main HTbin --alt herwig
time python runTUnfold.py -v v03 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} --main HTbin --alt Ptbin
time python runTUnfold.py -v v03 --only ${1} -p data -y ${2} --outputFolder ${3} --main Ptbin --alt herwig
time python runTUnfold.py -v v03 --only ${1} -p data -y ${2} --outputFolder ${3} --main HTbin --alt herwig

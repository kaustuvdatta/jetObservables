#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc820
export X509_USER_PROXY=/afs/cern.ch/user/k/kadatta/tmp/x509up_u113232

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/
eval `scramv1 runtime -sh`
env
echo "Running: python runTUnfold.py -v v1 --only ${1} -y ${2} --outputFolder ${3}"
time python runTUnfold.py -v v1 --only ${1} -p MCSelfClosure -y ${2} --outputFolder ${3} -s '_WSel' 
time python runTUnfold.py -v v1 --only ${1} -p MCSelfClosure -y ${2} --outputFolder ${3} -s '_topSel' 
time python runTUnfold.py -v v1 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} -s '_WSel' 
time python runTUnfold.py -v v1 --only ${1} -p MCClosure -y ${2} --outputFolder ${3} -s '_topSel' 
time python runTUnfold.py -v v1 --only ${1} -p data -y ${2} --outputFolder ${3} -s '_WSel' 
time python runTUnfold.py -v v1 --only ${1} -p data -y ${2} --outputFolder ${3} -s '_topSel' 


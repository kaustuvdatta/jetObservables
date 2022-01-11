#!/bin/bash
export SCRAM_ARCH=slc7_amd64_gcc820
export X509_USER_PROXY=/afs/cern.ch/user/a/algomez/x509up_u15148

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/
eval `scramv1 runtime -sh`
env
ivar=${1}
folder=${2}
echo "Running: python runTUnfold.py -v v04 --only ${ivar} --outputFolder ${folder}"
for year in 2017 2018;
do
    time python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y ${year} --outputFolder ${folder} --main Ptbin
    time python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y ${year} --outputFolder ${folder} --main HTbin
    time python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y ${year} --outputFolder ${folder} --main herwig
    time python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y ${year} --outputFolder ${folder} --main Ptbin --alt HTbin --noUnc
    time python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y ${year} --outputFolder ${folder} --main HTbin --alt herwig
    time python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y ${year} --outputFolder ${folder} --main HTbin --alt Ptbin
    time python runTUnfold.py -v v04 --only ${ivar} -p data -y ${year} --outputFolder ${folder} --main HTbin --alt herwig
done

python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main Ptbin
python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main HTbin
python runTUnfold.py -v v04 --only ${ivar} -p MCSelfClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main herwig
python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main Ptbin --alt HTbin --noUnc
python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main HTbin --alt Ptbin
python runTUnfold.py -v v04 --only ${ivar} -p MCClosure -y all --inputFolder ${folder} --outputFolder ${folder} --main HTbin --alt herwig
python runTUnfold.py -v v04 --only ${ivar} -p MCCrossClosure -y all --main HTbin --alt herwig --inputFolder ${folder} --outputFolder ${folder} --noUnc
python runTUnfold.py -v v04 --only ${ivar} -p data -y all --inputFolder ${folder} --outputFolder ${folder} --main HTbin --alt herwig
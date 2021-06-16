cp $CMSSW_BASE/src/Analysis/jetObservables/CMS_2021_PAS_SMP_21_XXX.cc $CMSSW_BASE/src/Rivet/SMP/src/
cp $CMSSW_BASE/src/Analysis/jetObservables/CMS_2021_PAS_SMP_21_XXX.plot $CMSSW_BASE/src/Rivet/SMP/data/
cp $CMSSW_BASE/src/Analysis/jetObservables/CMS_2021_PAS_SMP_21_XXX.info $CMSSW_BASE/src/Rivet/SMP/data/
cp $CMSSW_BASE/src/Analysis/jetObservables/CMS_2021_PAS_SMP_21_XXX.yoda $CMSSW_BASE/src/Rivet/SMP/data/

tmpDir=${PWD}
source $CMSSW_BASE/src/Rivet/rivetSetup.sh
cd $CMSSW_BASE/src/
scram b -j 8
cmsenv
cd ${tmpDir}
source $CMSSW_BASE/src/Rivet/rivetSetup.sh


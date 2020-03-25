#jssMeasurements

This analysis uses the extended nanoAOD format (info available [here] (https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1))

[based on fork from https://github.com/alefisico/jetObservables.git 106X branch] 

This repo includes, in addition to the above, Jupyter notebooks for unfolding the nSub bases and ratios, using TUnfold

##Instructions for running

###0) Basic setup:
To set up the code:
```bash
cmsrel  CMSSW_10_6_5
cd  CMSSW_10_6_5/src
cmsenv
git cms-addpkg GeneratorInterface/Core
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone https://github.com/kaustuvdatta/jssMeasurements.git -b 106X jssMeasurements/
git clone ssh://git@gitlab.cern.ch:7999/asparker/QJetMass.git jssMeasurements/QJetMass    ### This is just as example
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py jssMeasurements/Skimmer/test/
scram b -j 6
```

###1) Change input file paths/formats as required in the Unfolding/notebooks/nSubExtractor.py script; should also be easily modifiable to run over files via xrootd

###2) To run notebooks and tree extractor script, ensure following packages are vailable (either via conda env or cmsenv):
   ROOT 6.18/00 +
   root_numpy

### Further information

This package contains two folders: 
1. Skimmer: where the trees and histograms are created for the step 2. More information in the following [README](Skimmer/README.md)
2. Unfolding: where takes the input from step 1 and uses datasets on the T3 to do the unfolding procedure with TUnfold in jupyter notebooks.
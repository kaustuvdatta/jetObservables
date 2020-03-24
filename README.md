# jetObservables


This analysis uses extended NanoAOD format. More info about this format [here](https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1)

## Instructions

### Basic setup 
To set up the code:
```bash
cmsrel  CMSSW_10_6_5
cd  CMSSW_10_6_5/src
cmsenv
git cms-addpkg GeneratorInterface/Core
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone https://github.com/alefisico/jetObservables.git -b 106X jetObservables/
git clone ssh://git@gitlab.cern.ch:7999/asparker/QJetMass.git jetObservables/QJetMass    ### This is just as example
git clone git@github.com:kaustuvdatta/JSS-Unfolding.git jetObservables/Unfolding/notebooks/
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py jetObservables/Skimmer/test/
scram b -j 6
```

### Further information

This package contains two folders: 
1. Skimmer: where the trees and histograms are created for the step 2. More information in the following [README](Skimmer/README.md)
2. Unfolding: where takes the input from step 1 and uses combine to do the unfolding procedure. More information in the following [README](Unfolding/README.md)


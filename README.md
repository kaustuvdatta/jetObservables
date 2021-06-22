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
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone https://github.com/alefisico/jetObservables.git -b 106X jetObservables/
git clone ssh://git@gitlab.cern.ch:7999/asparker/QJetMass.git jetObservables/QJetMass    ### This is just as example
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py jetObservables/Skimmer/test/
scram b -j 6
```

### Further information

This package contains two folders: 
1. Skimmer: where the trees and histograms are created for the step 2. More information in the following [README](Skimmer/README.md)
2. Unfolding: where takes the input from step 1 and uses combine to do the unfolding procedure. More information in the following [README](Unfolding/README.md)
=======
# RIVET routine for jetObservables

This branch is for the rivet routines only.

## How to install

We will use the `CMSSW` infrastructure and follow this [tutorial](https://indico.cern.ch/event/962610/contributions/4049790/attachments/2131081/3588988/rivet_tutorial_mseidel.pdf). Dont forget to fork the main [cms-gen/Rivet](https://gitlab.cern.ch/cms-gen/Rivet) repository. Then:
```
cmsrel CMSSW_11_1_0
cd CMSSW_11_1_0/src
cmsenv
git clone ssh://git@gitlab.cern.ch:7999/${USER}/Rivet.git
cd Rivet
# add main repository and ensure sync
git remote add cms-gen ssh://git@gitlab.cern.ch:7999/cms-gen/Rivet.git
git pull cms-gen master
git push origin master # for those of you who forked long time ago ;)# set some paths for Rivet and LaTeX
source Rivet/rivetSetup.sh
scram b -j8
wget https://raw.githubusercontent.com/cms-sw/genproductions/22cdc0aba20c34087929cef168da715dad25581a/python/rivet_customize.py  -P Configuration/GenProduction/python/
git clone git@github.com:alefisico/jetObservables.git -b Rivet Analysis/jetObservables
source Analysis/jetObservables/copyFiles.sh
```

Although the files will be located in the folder `Analysis/jetObservables/`, the files must be located inside the `Rivet` folder from CMSSW. In order to copy them to the correct folder, the script [copyFiles.sh](copyFiles.sh) does this for you. Dont forget to run it _everytime_ that something changes in the files needed.


## Files needed

To create this routines one needs the following files:
 * [XXX.yoda](CMS_2021_PAS_SMP_21_XXX.yoda) Unfolded histograms from analysis. 
 * [XXX.cc](CMS_2021_PAS_SMP_21_XXX.cc) c++ code with RIVET syntax
 * [XXX.info](CMS_2021_PAS_SMP_21_XXX.info) describing the analysis
 * [XXX.plot](CMS_2021_PAS_SMP_21_XXX.plot) 

## Yoda file

Unfolded histograms needs to be saved in `yoda` format. For this in your unfolding script you can do something like:
``````
import yoda 

histToYoda = [  yoda.root.to_yoda( UnfoldHistogram ) ]
yoda.writeYODA( histToYoda, "outputName.yoda" )
``````
Dont forget to put all the unfold histograms together in a single `yoda` file. 

## Run routine

To run the routine within CMSSW, one can run:
```
cmsRun RIVET_QCD_Pt-15To7000_TuneCUETP8M1_Flat_13TeV-pythia8_cff.py
```
which will create a yoda file with the outputs defined in the `XXX.c` code. To create simple plots one can run:
```
rivet-mkhtml file.yoda
```
After running this command, a folder called `rivet-plots` contains a `index.html` file with all the plots.

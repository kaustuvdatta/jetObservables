# jetObservables


This analysis uses extended NanoAOD format. More info about this format [here](https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1)

## Instructions

### Basic setup 
To set up the code:
```bash
cmsrel  CMSSW_10_6_30
cd  CMSSW_10_6_30/src
cmsenv
git cms-addpkg GeneratorInterface/Core
git clone https://github.com/kaustuvdatta/nanoAOD-tools.git -b nSubMeasurements PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone https://github.com/kaustuvdatta/jetObservables.git -b RunIISummer20UL jetObservables/
git clone https://github.com/cms-sw/RecoLuminosity-LumiDB/tree/master 
cp jetObservables/TriggerEfficiencies/scripts/*.sh RecoLuminosity-LumiDB/scripts/ 
cd RecoLuminosity-LumiDB/scripts/ && chmod u+x {calc,PU}*.sh && cd -
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py jetObservables/Skimmer/test/
scram b -j 6
cmsenv
```

### Further information
(----> V1 of rep worked on CMSSW_10_6_14)
This package contains four folders: 
1. TriggerEfficiencies: where one can skim MiniAOD data samples to create histograms necessary for calculating the 99% efficiency (in pT) point for prescaled triggers used in the dijet event selection. Notebooks to carry the necessary fits are included therein. (calculation of prescale values is done after computing the luminosities accepted per trigger per given run period in RecoLuminosity-LumiDB/scripts/calc_prescale_lumi.sh, where for 2016 (non-)HIPM one must add up the luminosities also written out to the (2016-)2016HIPM-ext' csv files since some runs were not correctly labelled as pre/post VFP runs as per https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessingUL2016).
2. Skimmer: where the trees and histograms are created for the step 2. More information in the following [README](Skimmer/README.md).
3. SaturationROC: which notebooks showing how to pre-process parquet datasets and save in hdf5 format for training simple feedforward deep neural networks on M-body bases of N-subjettiness observables (at gen/reco level); this is to obtain ROC curves to understand for which value of 'M' discrimination power saturates between light flavour jets, as per the dijet selections, and boosted hadronic W/top decays. 
4. Unfolding: where takes the input from step 1 and uses combine to do the unfolding procedure. More information in the following [README](Unfolding/README.md).
5. Note: Most/all .py files in the python directories in the Unfolding and SaturationROC directories require python3 to be run; this will lead scram to throw errors, just move them out of the directories and replace after compiling CMSSSW
========================================================

# RIVET routine for jetObservables

This branch is for the rivet routines only.

## How to install
(to be updated?)
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

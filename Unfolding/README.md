## Unfolding step

### Unfolding with TUnfold

Several notebooks with the unfolding procedure are in  `jetObservables/Unfolding/notebooks/`. More information in each notebook.

In addition, there is a python script to compute the unfolding in [runTUnfold.py](test/runTUnfold.py).


### Unfolding with combine

A work-in-progress script is in [makeDatacards.py](test/makeDatacards.py)


### How to compute cross sections

We are using a modified version of [these instructions](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer#Automated_scripts_to_compute_the). To set up the code, if you did not include it before:

```bash
cd $CMSSW_BASE/src/
git cms-addpkg GeneratorInterface/Core
scram b -j 8
```

To run:
```bash
cd jetObservables/Skimmer/test/
git clone https://github.com/cms-sw/genproductions.git
mv genproductions/test/calculateXSectionAndFilterEfficiency/ .
rm -rf genproductions
```

To run:
```bash
cd $CMSSW_BASE/src/Analysis/jetObservables/test/calculateXSectionAndFilterEfficiency/
./calculateXSectionAndFilterEfficiency.sh -f datasets.txt -c Moriond17 -d MINIAODSIM -n 1000000  ## run using list of dataset names mode
```
where the input parameters are:
 * -f wants the input file containing the list of dataset names (default) or McM prepID (requires -m)
 * -c specifies the campaign, i.e. the string to be used to search for the secondary dataset name /.../*Moriond17*/*
 * -d specifies the datatier to be used, i.e.  /.../*/MINIAODSIM
 * -n number of events to be used for each dataset to compute the cross section
 * -m use the McM prepID instead of the dataset names


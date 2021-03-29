## Skimmer step


This step creates skim files for different selections: dijet, boosted W and top. The main files are [nSubProducer_dijet.py](python/nSubProducer_dijet.py), [nSubProducer_Wtop.py](python/nSubProducer_Wtop.py), [jetObservables_nSubProducer.py](test/jetObservables_nSubProducer.py) and [multicrab_nSubProducer.py](test/multicrab_nSubProducer.py)

### How to run skimmer

To run local:
```bash
cd jetObservables/Skimmer/test/
python jetObservables_nSubProducer.py --sample TTJets --local --selection dijet
```
The local option is to run in your local machine. Selection can be `dijet` or `Wtop`.


### Run multicrab jobs and storing datasets

To run in crab:
```bash
cd jetObservables/Skimmer/test/
python multicrab_nSubProducer.py --datasets TTJets -v 106X_v01 --selection dijet
```

while `-v` and `--selection` are self explanatory options, the `--datasets` option is linked to the [datasets.py](test/datasets.py) file.
This is an important file containing a list of datasets (miniAOD, nanoAOD, Skimmer), cross sections, number of events, etc. The `multicrab_nSubproducer.py` and the scripts below *require* to have that information in the `datasets.py` file.


### How to compute the total number of genWeights (events)

To compute the total number of processed events from genWeights, there is a script called [computeGenWeights.py](test/computeGenWeights.py). There, one needs to include the name of the dataset as written in the `datasets.py` script:

```bash
python computeGenWeights.py -d QCD
```

The output of this command looks like this:
```
{'QCDPt170to300': '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER'}
dataset /QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER has 28 files
Total number of genEventCount in sample:  6875125
Total number of genEventSumw in sample:  6875125.0
Total number of genEventSumw2 in sample:  6875125.0
Total sum of genWeights in sample:  27095.0
```
where the number needed for normalization of the samples is the one `Total number of genEventSumw in sample`.


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
mv genproductions/Utilities/calculateXSectionAndFilterEfficiency/ .
rm -rf genproductions
```

A modified version of [this script](https://github.com/cms-sw/genproductions/blob/master/Utilities/calculateXSectionAndFilterEfficiency/compute_cross_section.py) will compute the cross sections:
```bash
cd $CMSSW_BASE/src/jetObservables/Skimmer/test/
python computeXSection.py -y 2017 -d QCD_HT300
```
where the input parameters are `-d` for the name of the dataset listed in `datasets.py`. 

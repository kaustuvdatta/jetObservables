## Skimmer step


This step creates skim files for different selections: dijet, boosted W and top. The main files are [nSubProducer_withAllSel.py](python/nSubProducer_withAllSel.py) and [multicrab_nSubProducer.py](test/multicrab_nSubProducer.py)

### How to run skimmer

To run local:
```bash
cd jetObservables/Skimmer/test/
python jetObservables_nSubProducer.py --sample TTJets --local --selection dijet
```
The local option is to run in your local machine. Selection can be `dijet` or `Wtop`.

To run in crab:
```bash
cd jetObservables/Skimmer/test/
python multicrab_nSubProducer.py --datasets TTJets -v 106X_v01 --selection dijet
```

### How to compute the total number of genWeights (events)

To compute the total number of processed events from genWeights, there is a script called [computeGenWeights.py](test/computeGenWeights.py). There, one needs to include the name of the dataset processed in the skimmer step inside the `dictSamples` dictionary and then, if the key of the sample starts with `QCD`, one can run:

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

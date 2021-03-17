import awkward as ak
import coffea
from coffea import processor, hist
import numpy as np


class Analyzer(processor.ProcessorABC):

    def __init__(self, ismc, year):
        self.ismc = ismc
        self.year = year
        self.baseTriggers = {
            '2016' : [ 'IsoMu24', 'IsoTkMu24'],
            '2017' : [ 'IsoMu27' ],
            '2018' : [ 'IsoMu24' ],
        }
        self.triggers = [ 'AK8PFJet140','AK8PFJet200', 'AK8PFJet260', 'AK8PFJet320', 'AK8PFJet400', 'AK8PFJet450', 'AK8PFJet500']

        self._columns = [ 'FatJet_pt', 'FatJet_eta', 'FatJet_msoftdrop', 'FatJet_jetID' ] + [ 'HLT_'+t for t in self.triggers ]
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),

            "jets": hist.Hist(
                "entries",
                hist.Cat("dataset", "Dataset"),
                hist.Cat("selection", "Selection"),
                hist.Bin("pt", r"Leading jet $p^{T}$ [GeV]", 100, 0, 1000),
                hist.Bin("eta", r"Leading jet $\eta$", 10, -5, 5),
                hist.Bin("mass", r"Leading jet softdrop mass [GeV]", 51, 0, 500),
            ),
            'cutflow': hist.Hist(
                'entries',
                hist.Cat('dataset', 'Dataset'),
                hist.Bin('cut', 'Cut index', 11, 0, 11),
            ),

        })

    @property
    def accumulator(self):
        return self._accumulator

#    @property
#    def columns(self):
#        return self._columns


    def process(self, events): #, parameters=parameters, samples_info={}, \
                #lumimask=None, cat=False, boosted=False, uncertainty=None, \
                #uncertaintyName=None, parametersName=None, extraCorrection=None):

        dataset = events.metadata['dataset']
        #isRealData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        weights = processor.Weights(len(events))
        output = self.accumulator.identity()

        nEvents = len(events.event)
        run = events.run
        luminosityBlock = events.luminosityBlock
        #print("Processing %d %s events" % (nEvents, dataset))

        basetrigger = np.zeros(nEvents, dtype=bool)
        for t in self.baseTriggers[self.year]:
            basetrigger = basetrigger | events.HLT[t]
        selection.add('basetrigger', np.array(basetrigger))

        dummyTriggerDict = {}
        for it in self.triggers:
            dummyTriggerDict[ it ] = np.zeros(nEvents, dtype='bool')
            dummyTriggerDict[ it ] = dummyTriggerDict[ it ] | events.HLT[it]
            selection.add('trigger'+it, np.array(dummyTriggerDict[ it ]))

        fatjets = ak.zip({
            'pt' : events.FatJet.pt,
            'eta' : events.FatJet.eta,
            'phi' : events.FatJet.phi,
            'mass' : events.FatJet.mass,
            'msoftdrop' : events.FatJet.msoftdrop,
            'isTight' : events.FatJet.isTight,
            }, with_name='PtEtaPhiMCandidate')

        candidatejet = fatjets[
                (fatjets.pt > 200)
                & (abs(fatjets.eta) < 2.5)
                & (fatjets.msoftdrop >= 40)
                & fatjets.isTight
                ]
        canJetSel = ( ak.num(candidatejet) > 0 )
        selection.add( 'njetsBasic', np.array(canJetSel ) )

        baseCuts = [ 'njetsBasic', 'basetrigger'  ]
        triggerCuts = [ 'trigger'+i for i in self.triggers ]
        cuts = baseCuts + triggerCuts

        allcuts = set()
        output['cutflow'].fill(dataset=dataset, cut=0 )
        for i, cut in enumerate(cuts):
            allcuts.add(cut)
            cut = selection.all(*allcuts)
            output['cutflow'].fill(dataset=dataset, cut=i)

        allcuts = set()
        for icut in cuts:
            allcuts.add(icut)
            jcut = selection.all(*allcuts)
            FinalJet = candidatejet[ jcut ][:,0]
            output['jets'].fill(
                dataset=dataset,
                selection=icut,
                pt= FinalJet.pt,
                eta= FinalJet.eta,
                mass= FinalJet.msoftdrop,
#                 weight=weight,
            )


        return output

    def postprocess(self, accumulator):
        return accumulator

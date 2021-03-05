import os
import sys
import json
import argparse
import time

import numpy as np

import uproot as uproot
import awkward as ak
import coffea
from coffea import processor, hist
from coffea.nanoevents import NanoEventsFactory
from coffea.util import load, save
from coffea import processor

##############################################################################

class Analyzer(processor.ProcessorABC):

    def __init__(self, ismc, year):
        self.ismc = ismc
        self.year = year
        self.baseTriggers = {
            '2016' : [ 'IsoMu24', 'IsoTkMu24'],
            '2017' : [ 'AK8PFJet80' ], #'IsoMu27' ],
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
                (fatjets.pt > 170)
                & (abs(fatjets.eta) < 2.5)
                & (fatjets.msoftdrop >= 40)
                & fatjets.isTight
                ]
        canJetSel = ( ak.num(candidatejet) > 1 )
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

##############################################################################


def validate(file):
    try:
        fin = uproot.open(file)
        return fin['Events'].num_entries
    except:
        print("Corrupted file: {}".format(file))
        return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor coffea files')
    # Inputs
    parser.add_argument('-o', '--output', default=r'hists.coffea', help='Output histogram filename (default: %(default)s)')
    parser.add_argument('--samples', '--json', dest='samplejson', default='dummy_samples.json', help='JSON file containing dataset and file locations (default: %(default)s)')

    # Scale out
    parser.add_argument('--executor', choices=['iterative', 'futures', 'parsl/slurm', 'dask/condor', 'dask/slurm'], default='futures', help='The type of executor to use (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers (cores/threads) to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    parser.add_argument('-s', '--scaleout', type=int, default=6, help='Number of nodes to scale out to if using slurm/condor. Total number of concurrent threads is ``workers x scaleout`` (default: %(default)s)')
    parser.add_argument('--voms', default=None, type=str, help='Path to voms proxy, accsessible to worker nodes. By default a copy will be made to $HOME.')

    # Debugging
    parser.add_argument('--validate', action='store_true', help='Do not process, just check all files are accessible')
    parser.add_argument('--skipbadfiles', action='store_true', help='Skip bad files.')
    parser.add_argument('--only', type=str, default=None, help='Only process specific dataset or file')
    parser.add_argument('--limit', type=int, default=None, metavar='N', help='Limit to the first N files of each dataset in sample JSON')
    parser.add_argument('--chunk', type=int, default=500000, metavar='N', help='Number of events per process chunk')
    parser.add_argument('--max', type=int, default=None, metavar='N', help='Max number of chunks to run in total')

    args = parser.parse_args()
    if args.output == parser.get_default('output'):
        args.output = f'hists_{(args.samplejson).rstrip(".json")}.coffea'
        if args.only: args.output = f'hists_{(args.samplejson).rstrip(".json")}_{args.only}.coffea'


    # load dataset
    with open(args.samplejson) as f:
        sample_dict = json.load(f)

    for key in sample_dict.keys():
        sample_dict[key] = sample_dict[key][:args.limit]

    # For debugging
    if args.only is not None:
        if args.only in sample_dict.keys():  # is dataset
            sample_dict = dict([(args.only, sample_dict[args.only])])
        if "*" in args.only: # wildcard for datasets
            _new_dict = {}
            print("Will only proces the following datasets:")
            for k, v in sample_dict.items():
                if k.lstrip("/").startswith(args.only.rstrip("*")):
                    print("    ", k)
                    _new_dict[k] = v
            sample_dict = _new_dict
        else:  # is file
            for key in sample_dict.keys():
                if args.only in sample_dict[key]:
                    sample_dict = dict([(key, [args.only])])


    # Scan if files can be opened
    if args.validate:
        start = time.time()
        from p_tqdm import p_map
        all_invalid = []
        for sample in sample_dict.keys():
            _rmap = p_map(validate, sample_dict[sample], num_cpus=args.workers,
                      desc=f'Validating {sample[:20]}...')
            _results = list(_rmap)
            counts = np.sum([r for r in _results if np.isreal(r)])
            all_invalid += [r for r in _results if type(r) == str]
            print("Events:", np.sum(counts))
        print("Bad files:")
        for fi in all_invalid:
            print(f"  {fi}")
        end = time.time()
        print("TIME:", time.strftime("%H:%M:%S", time.gmtime(end-start)))
        if input("Remove bad files? (y/n)") == "y":
            print("Removing:")
            for fi in all_invalid:
                print(f"Removing: {fi}")
                os.system(f'rm {fi}')
        sys.exit(0)

    processor_instance = Analyzer( ismc=False, year='2017' )

    if args.executor not in ['futures', 'iterative']:
        # dask/parsl needs to export x509 to read over xrootd
        if args.voms is not None:
            _x509_path = args.voms
        else:
            _x509_localpath = [l for l in os.popen('voms-proxy-info').read().split("\n") if l.startswith('path')][0].split(":")[-1].strip()
            _x509_path = os.environ['HOME'] + f'/.{_x509_localpath.split("/")[-1]}'
            os.system(f'cp {_x509_localpath} {_x509_path}')

        env_extra = [
            'export XRD_RUNFORKHANDLER=1',
            f'export X509_USER_PROXY={_x509_path}',
            f'export X509_CERT_DIR={os.environ["X509_CERT_DIR"]}',
            'ulimit -u 32768',
        ]

    #########
    # Execute
    if args.executor in ['futures', 'iterative']:
        if args.executor == 'iterative':
            _exec = processor.iterative_executor
        else:
            _exec = processor.futures_executor
        output = processor.run_uproot_job(sample_dict,
                                    treename='Events',
                                    processor_instance=processor_instance,
                                    executor=_exec,
                                    executor_args={
                                        'skipbadfiles':args.skipbadfiles,
                                        'schema': processor.NanoAODSchema,
                                        'workers': args.workers},
                                    chunksize=args.chunk, maxchunks=args.max
                                    )
    elif args.executor == 'parsl/slurm':
        import parsl
        from parsl.providers import LocalProvider, CondorProvider, SlurmProvider
        from parsl.channels import LocalChannel
        from parsl.config import Config
        from parsl.executors import HighThroughputExecutor
        from parsl.launchers import SrunLauncher
        from parsl.addresses import address_by_hostname

        slurm_htex = Config(
            executors=[
                HighThroughputExecutor(
                    label="coffea_parsl_slurm",
                    address=address_by_hostname(),
                    prefetch_capacity=0,
                    provider=SlurmProvider(
                        channel=LocalChannel(script_dir='logs_parsl'),
                        launcher=SrunLauncher(),
                        max_blocks=(args.scaleout)+10,
                        init_blocks=args.scaleout,
                        partition='all',
                        worker_init="\n".join(env_extra) + "\nexport PYTHONPATH=$PYTHONPATH:$PWD",
                        walltime='00:120:00'
                    ),
                )
            ],
            retries=20,
        )
        dfk = parsl.load(slurm_htex)

        output = processor.run_uproot_job(sample_dict,
                                    treename='Events',
                                    processor_instance=processor_instance,
                                    executor=processor.parsl_executor,
                                    executor_args={
                                        'skipbadfiles':True,
                                        'schema': processor.NanoAODSchema,
                                        'config': None,
                                    },
                                    chunksize=args.chunk, maxchunks=args.max
                                    )

    elif 'dask' in args.executor:
        from dask_jobqueue import SLURMCluster, HTCondorCluster
        from distributed import Client
        from dask.distributed import performance_report
        import socket

        if 'slurm' in args.executor:
            cluster = SLURMCluster(
                queue='all',
                cores=args.workers,
                processes=args.workers,
                memory="200 GB",
                retries=10,
                walltime='00:30:00',
                env_extra=env_extra,
            )
        elif 'condor' in args.executor:
            cluster = HTCondorCluster(
                 cores=1, #args.workers,
                 memory='2GB',
                 disk='20GB',
                 env_extra=env_extra,
                 nanny=True,
                 extra=[ '--worker-port 8786', '--nanny-port 8785'],
                 scheduler_options={
                     'port' : 8786,
                     'host' : socket.gethostname()
                     },
                 #death_timeout='3200',
                 job_extra={
                    'log': 'dask_job_output.log',
                    'output': 'dask_job_output.out',
                    'error': 'dask_job_output.err',
                    '+JobFlavour' : "workday",
                    'should_transfer_files': 'Yes',
                    'when_to_transfer_output': 'ON_EXIT'
                }
            )
        cluster.scale(jobs=args.scaleout)

        client = Client(cluster)
        #with performance_report(filename="dask-report.html"):
        output = processor.run_uproot_job(sample_dict,
                                        treename='Events',
                                        processor_instance=processor_instance,
                                        executor=processor.dask_executor,
                                        executor_args={
                                            'client': client,
                                            'skipbadfiles':args.skipbadfiles,
                                            'schema': processor.NanoAODSchema,
                                        },
                                        #chunksize=args.chunk, maxchunks=args.max
                            )


    save(output, args.output)

    print(output)
    print(f"Saving output to {args.output}")

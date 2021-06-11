#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {
        ####### Wtop (xsecs not available on xsdb and/or not matching between xsdb and xseccalc are cross-checked with an approved analysis a la CMS AN-2018/129 )
    'SingleMuon' : {
        '2017' :  {
            'miniAOD' : {
                'B' : '/SingleMuon/Run2017B-UL2017_MiniAODv2-v1/MINIAOD',
                'C' : '/SingleMuon/Run2017C-UL2017_MiniAODv2-v1/MINIAOD',
                'D' : '/SingleMuon/Run2017D-UL2017_MiniAODv2-v1/MINIAOD',
                'E' : '/SingleMuon/Run2017E-UL2017_MiniAODv2-v1/MINIAOD',
                'F' : '/SingleMuon/Run2017F-UL2017_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'B': '/SingleMuon/algomez-Run2017B-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'C': '/SingleMuon/kadatta-Run2017C-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'D': '/SingleMuon/kadatta-Run2017D-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'E': '/SingleMuon/kadatta-Run2017E-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'F': '/SingleMuon/kadatta-Run2017F-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                },
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/SingleMuonAll_2017UL.root',
            'lumi' : 39412.23,
            },
        '2018' :  {
            'miniAOD' : {
                'A' : '/SingleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD',
                'B' : '/SingleMuon/Run2018B-UL2018_MiniAODv2-v1/MINIAOD',
                'C' : '/SingleMuon/Run2018C-UL2018_MiniAODv2-v1/MINIAOD',
                'D' : '/SingleMuon/Run2018D-UL2018_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'A' : '/SingleMuon/kadatta-Run2018A-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'B' : '/SingleMuon/kadatta-Run2018B-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'C' : '/SingleMuon/kadatta-Run2018C-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'D' : '/SingleMuon/kadatta-Run2018D-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                },
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/SingleMuonAll_2018UL.root',
            'lumi' : 56385.82159,
            },
        'color': ROOT.kWhite
    },

    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTJets_amcatnloFXFX-pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 34320777426.760,
            },
        '2018' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTJets_amcatnloFXFX-pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 34320777426.760,
            },
        'XS' : 722.8, #some notes say 831.76?
        'label' : 't#bar{t} amcatnloFXFX' ,
        'color': ROOT.kBlue
    },

    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 33674904582.836,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 33674904582.836,
            },
        'XS' : 365.24,
        'label' : 't#bar{t} powheg' ,
        'color': ROOT.kMagenta
    },

    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTTo2L2Nu_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 4778168320.403,
            },
        '2018' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTTo2L2Nu_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 4778168320.403,
            },
        'XS' : 88.29,
        'label' : 't#bar{t} DL' ,
        'color': ROOT.kCyan+1
    },

    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WJetsToLNu_madgraphMLM_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 10501241224508.998,
            },
        '2018' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WJetsToLNu_madgraphMLM_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 10501241224508.998,
            },
        'XS' : 5.368e+04,
        'label' : 'WJets' ,
        'color': ROOT.kOrange-2
    },
    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_s-channel_amcatnlo_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 70465639.296,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_s-channel_amcatnlo_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 70465639.296,
            },
        'XS' : 3.549e+00,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_t-channel_top_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 655197632.514,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_t-channel_top_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 655197632.514,
            },
        'XS' : 1.197e+02,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_t-channel_antitop_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 265124234.844,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_t-channel_antitop_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 265124234.844,
            },
        'XS' : 7.174e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_tW_top_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 219634313.571,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_tW_top_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 219634313.571,
            },
        'XS' : 3.245e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_tW_antitop_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 277735003.1208,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ST_tW_antitop_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 277735003.1208,
            },
        'XS' : 3.251e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'WW_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WW_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 7876265.259,
            },
        '2018' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WW_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 7876265.259,
            },
        'XS' : 7.577e+01,
        'label' : 'Diboson' ,
        'color': ROOT.kBlue
    },
    'ZZ_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ZZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 7876265.259,
            },
        '2018' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/ZZ_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 1981800.000000,
            },
        'XS' : 2.748e+00,
        'label' : 'Diboson' ,
        'color': ROOT.kBlue
    },
    'WZ_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 3970000.00,
            },
        '2018' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WZ_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 3970000.00,
            },
        'XS' : 1.21e+00,
        'label' : 'Dibosons' ,
        'color': ROOT.kBlue
    },
    'QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_170to300_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 36027673.00,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 7055.0,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_300to470_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 494796.07288,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 619.3,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_470to600_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 517383.713837,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 59.24,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_600to800_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 17318812.571925,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 3970000.00,
            },
        'XS' : 18.21,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_800to1000_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.0,
            },
        'XS' : 7055.0,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/mmarcheg-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/QCD_Muen_Pt_1000toInf_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 14642553.00,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/WZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.00,
            },
        'XS' : 1.078,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    ## MC variations
    'TTToSemileptonic_TuneCP5_erdON_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_erdON_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_erdON_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 365.24,
        'label' : 'CR Model (erdON)' ,
        'color': ROOT.kBlue-3
    },

    'TTToSemileptonic_TuneCP5down_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_TuneCP5down_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_TuneCP5down_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 365.24,
        'label' : 'Tune' ,
        'color': ROOT.kBlue+7
    },


    'TTToSemileptonic_TuneCP5up_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_TuneCP5up_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_TuneCP5up_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 365.24,
        'label' : 'Tune' ,
        'color': ROOT.kBlue+7
    },

    'TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_hdampDOWN_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_hdampDOWN_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 365.24,
        'label' : '#h_{damp}' ,
        'color': ROOT.kBlue+7
    },


    'TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_hdampUP_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-754f8ab81f6f0298c5fa7c45094d30e4/USER' ],
            'skimmerHisto' : '/eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Samples/TTToSemileptonic_hdampUP_powheg_pythia8_2018UL.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 365.24,
        'label' : '#h_{damp}' ,
        'color': ROOT.kBlue+7
    },

}

def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


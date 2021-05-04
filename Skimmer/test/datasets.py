#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {

    'JetHT' : {
        '2017' :  {
            'miniAOD' : {
                'B' : '/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD',
                'C' : '/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD',
                'D' : '/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD',
                'E' : '/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD',
                'F' : '/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'B': '/JetHT/algomez-Run2017B-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'C': '/JetHT/algomez-Run2017C_UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'D': '/JetHT/algomez-Run2017D-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'E': '/JetHT/algomez-Run2017E-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'F': '/JetHT/algomez-Run2017F-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                },
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_JetHT2017ALL.root',
            'lumi' : 41473.09,
            'triggerList' : {
                'AK8PFJet80' : 30000.00, #30000,  #16419.91,
                'AK8PFJet140' : 3060.15, #1000., #1560.15,
                'AK8PFJet200' : 249.70, #80., #219.70,
                'AK8PFJet260' : 320, #110., #88.43,
                'AK8PFJet320' : 80, #45.,#33.83,
                'AK8PFJet400' : 5.40, #10,  #5.40,
                'AK8PFJet450' : 4.4 , #6.,  # 4.299,
                'AK8PFJet500' : 1.,
                }
            },
        '2018' :  {
            'miniAOD' : {
                'A' : '/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD',
                'B' : '/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD',
                'C' : '/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD',
                'D' : '/JetHT/Run2018D-UL2018_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'A' : '/JetHT/algomez-Run2018A-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'B' : '/JetHT/algomez-Run2018B-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'C' : '/JetHT/algomez-Run2018C-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                'D' : '/JetHT/algomez-Run2018D-UL2018_MiniAODv2-v1_PFNanoAOD-690f1d1c4f2456484a2616002f2d2b6d/USER',
                },
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_JetHT2017ALL.root',
            'lumi' : 58918.895,
            'triggerList' : {
                'AK8PFJet80' : 30000.00, #30000,  #16419.91,
                'AK8PFJet140' : 3060.15, #1000., #1560.15,
                'AK8PFJet200' : 249.70, #80., #219.70,
                'AK8PFJet260' : 320, #110., #88.43,
                'AK8PFJet320' : 80, #45.,#33.83,
                'AK8PFJet400' : 5.40, #10,  #5.40,
                'AK8PFJet450' : 4.4 , #6.,  # 4.299,
                'AK8PFJet500' : 1.,
                }
            },
        'color': ROOT.kWhite
    },

    'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-jetObservables_Skimmer_UL17_QCDHerwig_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHerwig_UL2017.root',
            'nevents' : 19880000.,
            'nGenWeights' :397097.833,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.32E+09,
        'label' : 'QCD Herwig7',
        'color': ROOT.kBlue
    },

    'QCD_Pt_170to300_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt170to300_UL2017.root',
            'nevents' : 28768700.,
            'nGenWeights' : 28768700.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1_ext1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.035e+05,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_300to470_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt300_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt300to470_UL2017.root',
            'nevents' : 33257500.,
            'nGenWeights' : 33257505.738,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 6.760e+03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_470to600_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt470to600_UL2017.root',
            'nevents' : 25797300.,
            'nGenWeights' : 25797383.655,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 5.516e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_600to800_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt600to800_UL2017.root',
            'nevents' : 58340300.,
            'nGenWeights' : 58340301.431,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.564e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_800to1000_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_ext1-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt800to1000_UL2017.root',
            'nevents' : 32708500.,
            'nGenWeights' : 32708500.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v02-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 2.624e+01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt1000to1400_UL2017.root',
            'nevents' : 18644400.,
            'nGenWeights' : 18644400.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 7.477e+00,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt1400to1800.root',
            'nevents' : 5193100.,
            'nGenWeights' : 5193100.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 6.423e-01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1.,
            'nGenWeights' : 1.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt2400_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt2400to3200_UL2017.root',
            'nevents' : 1916800.,
            'nGenWeights' : 1916800.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 5.233e-03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt3200_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDPt3200toInf_UL2017.root',
            'nevents' : 800000.,
            'nGenWeights' : 800000.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.351e-04,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT100to200_UL2017.root',
            'nevents' : 76745700.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT100to200_UL2017_Closure.root',
            'nGenWeightsClosure' : 9106631.,
            'nGenWeights' : 76745700.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 2.362e+07,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT200to300_UL2017.root',
            'nevents' : 56892662.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT200to300_UL2017_Closure.root',
            'nGenWeightsClosure' : 9153696.0,
            'nGenWeights' : 56892662.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.555e+06,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT300to500_UL2017.root',
            'nevents' : 52225204.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT300to500_UL2017_Closure.root',
            'nGenWeightsClosure' : 9110458.0,
            'nGenWeights' : 52225204.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 3.243e+05,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT500to700_UL2017.root',
            'nevents' : 8928812.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT500to700_UL2017_Closure.root',
            'nGenWeightsClosure' : 2290261.0,
            'nGenWeights' : 8928812.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 3.048e+04,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT700to1000_UL2017.root',
            'nevents' : 41753837.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT700to1000_UL2017_Closure.root',
            'nGenWeightsClosure' : 8532247.0,
            'nGenWeights' : 41753837.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 6.433e+03,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT1000to1500_UL2017.root',
            'nevents' : 11964428.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT1000to1500_UL2017_Closure.root',
            'nGenWeightsClosure' : 5634302.0,
            'nGenWeights' : 11964428.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.116e+03,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT1500to2000_UL2017.root',
            'nevents' : 10496669.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT1500to2000_UL2017_Closure.root',
            'nGenWeightsClosure' : 3881216.0,
            'nGenWeights' : 10496669.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.081e+02,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT2000toInf_UL2017.root',
            'nevents' : 5366803.,
            'skimmerHistoClosure' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHT2000toInf_UL2017_Closure.root',
            'nGenWeightsClosure' : 732412.0,
            'nGenWeights' : 5366803.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 2.193e+01,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kMagenta
    },

    ####### Wtop
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
                'B': '',
                'C': '',
                'D': '',
                'E': '',
                'F': '',
                },
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/SingleMuonB_2017UL.root',
            'lumi' : 4803.36,
            },
        '2018' :  {
            'miniAOD' : {
                'A' : '/SingleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD',
                'B' : '/SingleMuon/Run2018B-UL2018_MiniAODv2-v1/MINIAOD',
                'C' : '/SingleMuon/Run2018C-UL2018_MiniAODv2-v1/MINIAOD',
                'D' : '/SingleMuon/Run2018D-UL2018_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'A' : '',
                'B' : '',
                'C' : '',
                'D' : '',
                },
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'lumi' : 1.,
            },
        'color': ROOT.kWhite
    },

    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/TTJets_amcatnloFXFX-pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 34320777426.760,
            },
        'XS' : 722.8,
        'label' : 't#bar{t} amcatnloFXFX' ,
        'color': ROOT.kBlue
    },

    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/TTToSemileptonic_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 33674904582.836,
            },
        'XS' : 365.34,
        'label' : 't#bar{t} powheg' ,
        'color': ROOT.kMagenta
    },

    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/TTTo2L2Nu_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 4778168320.403,
            },
        'XS' : 88.29,
        'label' : 't#bar{t} DL' ,
        'color': ROOT.kMagenta
    },

    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/WJetsToLNu_madgraphMLM_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 10501241224508.998,
            },
        'XS' : 5.368e+04,
        'label' : 'WJets' ,
        'color': ROOT.kOrange
    },
    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/ST_s-channel_amcatnlo_pythia8_2017UL.root',
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
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/ST_t-channel_top_powheg_pythia8_2017UL.root',
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
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/ST_t-channel_antitop_powheg_pythia8_2017UL.root',
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
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/ST_tW_top_powheg_pythia8_2017UL.root',
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
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/ST_tW_antitop_powheg_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 277735003.1208,
            },
        'XS' : 3.251e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'WW_TuneCP5_DoubleScattering_13TeV-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/WW_TuneCP5_DoubleScattering_13TeV-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/WW_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 7876265.259,
            },
        'XS' : 7.577e+01,
        'label' : 'Diboson' ,
        'color': ROOT.kBlue
    },
    'WZ' : {
        '2017' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/WZ_pythia8_2017UL.root',
            'nevents' : 1.,
            'nGenWeights' : 3970000.00,
            },
        'XS' : 1.21e+00,
        'label' : 'Dibosons' ,
        'color': ROOT.kBlue
    },
#    '' : {
#        '2017' :  {
#            'miniAOD' : [ '' ],
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
#            'nevents' : 1.,
#            'nGenWeights' : 33674904582.836,
#            },
#        'XS' : 365.34,
#        'label' : 't#bar{t} powheg' ,
#        'color': ROOT.kBlue
#    },

}

def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


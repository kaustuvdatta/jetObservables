#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {

    'JetHT' : {
        'selection' : 'dijet',
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
                'C': '/JetHT/algomez-Run2017C-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'D': '/JetHT/algomez-Run2017D-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'E': '/JetHT/algomez-Run2017E-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                'F': '/JetHT/algomez-Run2017F-UL2017_MiniAODv2-v1_PFNanoAOD-aeb179a6519d272d4336b9926381baad/USER',
                },
            'skimmerHisto' : 'jetObservables_histograms_JetHT2017ALL.root',
            'lumi' : 41527.74,
            'triggerList' : {
                'AK8PFJet80' : 16419.91,
                'AK8PFJet140' : 1560.15,
                'AK8PFJet200' : 219.70,
                'AK8PFJet260' : 88.43,
                'AK8PFJet320' : 33.827,
                'AK8PFJet400' : 5.404,
                'AK8PFJet450' : 4.299,
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
            'skimmerHisto' : 'jetObservables_histograms_JetHT2018ALL.root',
            'lumi' : 58918.895,
            'triggerList' : {
                'AK8PFJet80' : 27585.68,
                'AK8PFJet140' : 1268.74,
                'AK8PFJet200' : 295.24,
                'AK8PFJet260' : 127.99,
                'AK8PFJet320' : 48.17,
                'AK8PFJet400' : 16.08,
                'AK8PFJet450' : 8.09,
                'AK8PFJet500' : 1.,
                }
            },
        'color': ROOT.kWhite
    },

#    'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7' : {
#        'selection' : 'dijet',
#        '2017' :  {
#            'miniAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
#            'skimmer' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03-7a1edb72314467730e458def0bc98536/USER' ],
#            'skimmerHisto' : 'jetObservables_histograms_QCDPt15to7000_herwig_UL2017.root',
#            'nevents' : 12718100.,
#            'nGenWeights' : 53765.806,
#            },
#        '2018' :  {
#            'miniAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
#            'skimmer' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p2-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
#            'skimmerHisto' : 'jetObservables_histograms_QCDPt15to7000_herwig_UL2018.root',
#            'nevents' : 18744500.,
#            'nGenWeights' : 375606.55,
#            },
#        'XS' : 7.417e+08,
#        'label' : 'QCD Herwig7',
#        'color': ROOT.kOrange
#    },


    'QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt150to3000_herwig_UL2017.root',
            'nevents' : 19880000.,
            'nGenWeights' :397097.833,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt-150to3000_TuneCH3_FlatPower7_13TeV-herwig7/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt150to3000_herwig_UL2018.root',
            'nevents' : 6649600.,
            'nGenWeights' : 28765.923,
            },
        'selection' : 'dijet',
        'XS' : 8.637e+03, ## for 2018 1.086e+04
        'label' : 'QCD Herwig7',
        'color': ROOT.kPink
    },

    'QCD_Pt_170to300_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2017.root',
            'nevents' : 27963000.0,
            'nGenWeights' : 27963000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2018.root',
            'nevents' : 28425000.0,
            'nGenWeights' : 28425000.0,
            },
        'selection' : 'dijet',
        'XS' : 1.035e+05,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_300to470_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2017.root',
            'nevents' : 53017000.,
            'nGenWeights' : 53017008.47,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2018.root',
            'nevents' : 55315000.,
            'nGenWeights' : 55315008.680,
            },
        'selection' : 'dijet',
        'XS' : 6.760e+03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_470to600_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2017.root',
            'nevents' : 46970000.,
            'nGenWeights' : 46970154.474,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2018.root',
            'nevents' : 51191000.,
            'nGenWeights' : 51191166.543,
            },
        'selection' : 'dijet',
        'XS' : 5.516e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_600to800_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt600to800_UL2017.root',
            'nevents' : 63458000.,
            'nGenWeights' : 63458001.544,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' :  'jetObservables_histograms_QCDPt600to800_UL2018.root',
            'nevents' : 65300000.,
            'nGenWeights' : 65300001.492,
            },
        'selection' : 'dijet',
        'XS' : 1.564e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_800to1000_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2017.root',
            'nevents' : 35696000.0,
            'nGenWeights' : 35696000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2018.root',
            'nevents' : 36056000.0,
            'nGenWeights' : 36056000.0,
            },
        'selection' : 'dijet',
        'XS' : 2.624e+01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2017.root',
            'nevents' : 18653000.0,
            'nGenWeights' : 18653000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2018.root',
            'nevents' : 19106000.0,
            'nGenWeights' : 19106000.0,
            },
        'selection' : 'dijet',
        'XS' : 7.477e+00,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2017.root',
            'nevents' : 10358000.0,
            'nGenWeights' : 10358000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2018.root',
            'nevents' : 10550000.0,
            'nGenWeights' : 10550000.0,
            },
        'selection' : 'dijet',
        'XS' : 6.423e-01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '',
            'nevents' : 5191000.0,
            'nGenWeights' : 5191000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1800to2400_UL2018.root',
            'nevents' : 5152000.0,
            'nGenWeights' : 5152000.0,
            },
        'selection' : 'dijet',
        'XS' : 8.746e-02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2017.root',
            'nevents' : 2757000.0,
            'nGenWeights' : 2757000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2018.root',
            'nevents' : 2901000.0,
            'nGenWeights' : 2901000.0,
            },
        'selection' : 'dijet',
        'XS' : 5.233e-03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v03p1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2017.root',
            'nevents' : 1000000.0,
            'nGenWeights' : 1000000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2018.root',
            'nevents' : 934000.0,
            'nGenWeights' : 934000.0,
            },
        'selection' : 'dijet',
        'XS' : 1.351e-04,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    },

    'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2017.root',
            'nevents' : 76745700.0,
            'nGenWeights' : 76745700.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2018.root',
            'nevents' : 70101398.0,
            'nGenWeights' : 70101398.0,
            },
        'selection' : 'dijet',
        'XS' : 2.362e+07,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2017.root',
            'nevents' : 56892662.,
            'nGenWeights' : 56892662.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2018.root',
            'nevents' : 22793054.0,
            'nGenWeights' : 22793054.0,
            },
        'selection' : 'dijet',
        'XS' : 1.555e+06,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2017.root',
            'nevents' : 52225204.,
            'nGenWeights' : 52225204.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2018.root',
            'nevents' : 55102146.0,
            'nGenWeights' : 55102146.0,
            },
        'selection' : 'dijet',
        'XS' : 3.243e+05,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2017.root',
            'nevents' : 8928812.0,
            'nGenWeights' : 8928812.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2018.root',
            'nevents' : 51866729.0,
            'nGenWeights' : 51866729.0,
            },
        'selection' : 'dijet',
        'XS' : 3.048e+04,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2017.root',
            'nevents' : 41753837.0,
            'nGenWeights' : 41753837.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2018.root',
            'nevents' : 47028925.0,
            'nGenWeights' : 47028925.0,
            },
        'selection' : 'dijet',
        'XS' : 6.433e+03,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2017.root',
            'nevents' : 12643258.0,
            'nGenWeights' : 12643258.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2018.root',
            'nevents' : 15589003.0,
            'nGenWeights' : 15589003.0,
            },
        'selection' : 'dijet',
        'XS' : 1.116e+03,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2017.root',
            'nevents' : 10496669.0,
            'nGenWeights' : 10496669.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2018.root',
            'nevents' : 9825850.0,
            'nGenWeights' : 9825850.0,
            },
        'selection' : 'dijet',
        'XS' : 1.081e+02,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v03p2-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2017.root',
            'nevents' : 5366803.0,
            'nGenWeights' : 5366803.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL18PFNanoAOD_jetObservables_Skimmer_v03p1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2018.root',
            'nevents' : 4484780.0,
            'nGenWeights' : 4484780.0,
            },
        'selection' : 'dijet',
        'XS' : 2.193e+01,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    },

    ####### Wtop (xsecs not available on xsdb and/or not matching between xsdb and xseccalc are cross-checked with an approved analysis a la CMS AN-2018/129 )
    'SingleMuon' : {
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
        'XS' : 1.078,
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    # MC variations
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
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
        'selection' : 'Wtop',
        'XS' : 365.24,
        'label' : '#h_{damp}' ,
        'color': ROOT.kBlue+7
    },

}

def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


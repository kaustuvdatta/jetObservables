#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {

    'JetHT' : {
        '2017' :  {
            'miniAOD' : [ '' ],
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
                'AK8PFJet80' : 16419.91, #30000,  #16419.91,
                'AK8PFJet140' : 1560.15, #1000., #1560.15,
                'AK8PFJet200' : 219.70, #80., #219.70,
                'AK8PFJet260' : 200, #110., #88.43,
                'AK8PFJet320' : 50, #45.,#33.83,
                'AK8PFJet400' : 5.40, #10,  #5.40,
                'AK8PFJet450' : 4.299, #6.,  # 4.299,
                'AK8PFJet500' : 1.,
                }
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'color': ROOT.kWhite
    },

    'QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7' : {
        '2016' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt-15to7000_TuneCH3_Flat_13TeV_herwig7/algomez-jetObservables_Skimmer_UL17_QCDHerwig_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/jetObservables_histograms_QCDHerwig2017.root',
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
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1.,
            'nGenWeights' : 1.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt300_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 33257500.,
            'nGenWeights' : 33257505.738,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 25797300.,
            'nGenWeights' : 25797383.655,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1.,
            'nGenWeights' : 1.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 32708500.,
            'nGenWeights' : 32708500.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 18644400.,
            'nGenWeights' : 18644400.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1.,
            'nGenWeights' : 1.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt2400_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1916800.,
            'nGenWeights' : 1916800.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-jetObservables_Skimmer_UL17_QCDPt3200_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 800000.,
            'nGenWeights' : 800000.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
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
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 76745700.,
            'nGenWeights' : 76745700.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 2.362e+07,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 56892662.,
            'nGenWeights' : 56892662.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.555e+06,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 52225204.,
            'nGenWeights' : 52225204.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 3.243e+05,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 8928812.,
            'nGenWeights' : 8928812.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 3.048e+04,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 41753837.,
            'nGenWeights' : 41753837.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 6.433e+03,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 1.,
            'nGenWeights' : 1.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.116e+03,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 10496669.,
            'nGenWeights' : 10496669.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 1.081e+02,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

    'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer19UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/algomez-RunIISummer19UL17PFNanoAOD_jetObservables_Skimmer_v02-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : '/afs/cern.ch/work/a/algomez/jetObservables/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Rootfiles/',
            'nevents' : 5366803.,
            'nGenWeights' : 5366803.,
            },
        '2018' :  {
            'miniAOD' : [ '' ],
            'nanoAOD' : [ '' ],
            'skimmer' : [ '' ],
            'nevents' : 00,
            'nGenWeights' : 00,
            },
        'XS' : 2.193e+01,
        'label' : 'QCD madgraphMLM' ,
        'color': ROOT.kBlue
    },

#    '' : {
#        '2016' :  {
#            'miniAOD' : [ '' ],
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2017' :  {
#            'miniAOD' : [ '' ],
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'miniAOD' : [ '' ],
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 0,
#    },
#    'TT' : {
#        '2016' :  {
#            'miniAOD' : [ '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/algomez-TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'skimmer' : [ '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 76915549.,
#            'nGenWeights' : 76489216.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 7.306e+02, # +- 5.572e-01 pb
#    },
#    'TTJets' : {
#        '2016' :  {
#            'miniAOD' : [ '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv3-PUMoriond1794XmcRun2-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'skimmer' : [ '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/algomez-jetObservables_Skimmer_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 10199051.,
#            'nGenWeights' : 7662481.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 5.093e+02, # +- 4.456e-01 pb
#    },
#    'ST_s-channel_4f_InclusiveDecays' : {
#        '2016' :  {
#            'miniAOD' : [ '/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/algomez-STs-channel4fInclusiveDecays13TeV-amcatnlo-pythia8RunIISummer16MiniAODv3-PUMoriond17-dafc15ff64439ee3efd0c8e48ce3e57e/USER ' ],
#            'skimmer' : [ '/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 2989199.,
#            'nGenWeights' : 30289975.20,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 1.012e+01, # +- 1.334e-02 pb
#    },
#    'ST_t-channel_antitop' : {
#        '2016' :  {
#            'miniAOD' : [ '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/algomez-STt-channelantitop4finclusiveDecays13TeV-powhegV2-madspin-pythia8TuneCUETP8M1-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'skimmer' : [ '/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 38811017.,
#            'nGenWeights' : 37877067.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 7.441e+01, # +- 4.171e-01 pb
#    },
#    'ST_t-channel_top' : {
#        '2016' :  {
#            'miniAOD' : [ '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ "/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/algomez-STt-channeltop4finclusiveDecays13TeV-powhegV2-madspin-pythia8TuneCUETP8M1-dafc15ff64439ee3efd0c8e48ce3e57e/USER" ],
#            'skimmer' : [ '/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 67105876.,
#            'nGenWeights' : 66826980.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 1.233e+02, # +- 7.721e-01 pb
#    },
#    'ST_tW_antitop' : {
#        '2016' :  {
#            'miniAOD' : [ "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM" ],
#            'nanoAOD' : [ "/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/algomez-STtWantitop5finclusiveDecays13TeV-powheg-pythia8TuneCUETP8M2T4-dafc15ff64439ee3efd0c8e48ce3e57e/USER" ],
#            'skimmer' : [ '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 998276.,
#            'nGenWeights' : 998276.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 3.806e+01, # +- 3.055e-02 pb
#    },
#    'ST_tW_top' : {
#        '2016' :  {
#            'miniAOD' : [ '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM' ],
#            'nanoAOD' : [ "/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/algomez-STtWtop5finclusiveDecays13TeV-powheg-pythia8TuneCUETP8M2T4-dafc15ff64439ee3efd0c8e48ce3e57e/USER" ],
#            'skimmer' : [ '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 992024.,
#            'nGenWeights' : 992024.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 3.809e+01, # +- 3.050e-02 pb
#    },
#    'WJetsToLNu' : {
#        '2016' :  {
#            'miniAOD' : [ '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM', '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM' ],
#            'nanoAOD' : [ '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/algomez-WJetsToLNuTuneCUETP8M113TeV-amcatnloFXFX-pythia8RunIISummer16MiniAODv3-PUMoriond1794X-dafc15ff64439ee3efd0c8e48ce3e57e/USER', '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/algomez-WJetsToLNuTuneCUETP8M113TeV-amcatnloFXFX-pythia8RunIISummer16MiniAODv3-PUMoriond1794X_ext2-v1-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'skimmer' : [ '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/algomez-jetObservables_Skimmer_v02-dafc15ff64439ee3efd0c8e48ce3e57e/USER' ],
#            'nevents' : 24120319.,
#            'nGenWeights' : 611509190272.00,
#            },
#        '2017' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        '2018' :  {
#            'nanoAOD' : [ '' ],
#            'skimmer' : [ '' ],
#            'nevents' : 00,
#            'nGenWeights' : 00,
#            },
#        'XS' : 6.038e+04, # +- 1.238e+02 pb
#    },

}

def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


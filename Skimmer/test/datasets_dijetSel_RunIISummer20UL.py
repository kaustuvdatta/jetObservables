#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {

	'JetHT' : 
		{
        'selection' : 'dijet',
        '2016_preVFP' :  {
            'miniAOD' : {
                'B': '/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'C': '/JetHT/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'D': '/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'E': '/JetHT/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'F': '/JetHT/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                },
            'nanoAOD' : {
                'B': '/JetHT/algomez-Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'C': '/JetHT/algomez-Run2016C-HIPM_UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'D': '/JetHT/algomez-Run2016D-HIPM_UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'E': '/JetHT/algomez-Run2016E-HIPM_UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'F': '/JetHT/algomez-Run2016F-HIPM_UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                },
            'skimmerHisto' : 'jetObservables_histograms_JetHT2016ALL_HIPM.root',
            'lumi' : 16707.068,
            'triggerList' : {
                'AK8PFJet80' : 27661.64,
                'AK8PFJet140' : 2689.77,
                'AK8PFJet200' : 278.62,
                'AK8PFJet260' : 56.78,
                'AK8PFJet320' : 19.59,
                'AK8PFJet400' : 6.41,
                'AK8PFJet450' : 1.,
                'AK8PFJet500' : 1.,
                }
            },
        '2016' :  {
            'miniAOD' : {
                'F': '/JetHT/Run2016F-UL2016_MiniAODv2-v2/MINIAOD',
                'G': '/JetHT/Run2016G-UL2016_MiniAODv2-v2/MINIAOD',
                'H': '/JetHT/Run2016H-UL2016_MiniAODv2-v2/MINIAOD',
                },
            'nanoAOD' : {
                'F': '/JetHT/algomez-Run2016F-UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'G': '/JetHT/algomez-Run2016G-UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                'H': '/JetHT/algomez-Run2016H-UL2016_MiniAODv2-v2_PFNanoAOD-f51f8a19a3ab12c4f8e76ddb12529cf0/USER',
                },
            'skimmerHisto' : 'jetObservables_histograms_JetHT2016ALL.root',
            'lumi' : 16795.659,
            'triggerList' : {
                'AK8PFJet80' : 41995.50,
                'AK8PFJet140' : 4319.21,
                'AK8PFJet200' : 652.64,
                'AK8PFJet260' : 75.17,
                'AK8PFJet320' : 25.00,
                'AK8PFJet400' : 8.47,
                'AK8PFJet450' : 1.,
                'AK8PFJet500' : 1.,
                }
            },
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

    'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 

    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT100to200_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 2.362e+07,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT200to300_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 1.555e+06,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT300to500_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 3.243e+05,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT500to700_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 3.048e+04,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
	    '2017' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT700to1000_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 6.433e+03,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
   		},

    'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1000to1500_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
    
        'selection' : 'dijet',
        'XS' : 1.116e+03,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT1500to2000_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
    	          },
        'selection' : 'dijet',
        'XS' : 1.081e+02,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2017.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDHT2000toInf_UL2018.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        'selection' : 'dijet',
        'XS' : 2.193e+01,
        'label' : 'QCD MG5+Pythia8' ,
        'color': ROOT.kMagenta
    	},

    'QCD_Pt_170to300_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2017.root',
            'nevents' : 27963000.0,
            'nGenWeights' : 27963000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt170to300_UL2018.root',
            'nevents' : 28425000.0,
            'nGenWeights' : 28425000.0,
            },
        'selection' : 'dijet',
        'XS' : 1.035e+05,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_300to470_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2017.root',
            'nevents' : 53017000.,
            'nGenWeights' : 53017008.47,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt300to470_UL2018.root',
            'nevents' : 55315000.,
            'nGenWeights' : 55315008.680,
            },

        'selection' : 'dijet',
        'XS' : 6.760e+03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_470to600_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2017.root',
            'nevents' : 46970000.,
            'nGenWeights' : 46970154.474,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt470to600_UL2018.root',
            'nevents' : 51191000.,
            'nGenWeights' : 51191166.543,
            },
        'selection' : 'dijet',
        'XS' : 5.516e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_600to800_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt600to800_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' :  'jetObservables_histograms_QCDPt600to800_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt600to800_UL2017.root',
            'nevents' : 63458000.,
            'nGenWeights' : 63458001.544,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' :  'jetObservables_histograms_QCDPt600to800_UL2018.root',
            'nevents' : 65300000.,
            'nGenWeights' : 65300001.492,
            },
        'selection' : 'dijet',
        'XS' : 1.564e+02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_800to1000_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2017.root',
            'nevents' : 35696000.0,
            'nGenWeights' : 35696000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt800to1000_UL2018.root',
            'nevents' : 36056000.0,
            'nGenWeights' : 36056000.0,
            },
        'selection' : 'dijet',
        'XS' : 2.624e+01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2017.root',
            'nevents' : 18653000.0,
            'nGenWeights' : 18653000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1000to1400_UL2018.root',
            'nevents' : 19106000.0,
            'nGenWeights' : 19106000.0,
            },
        'selection' : 'dijet',
        'XS' : 7.477e+00,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2017.root',
            'nevents' : 10358000.0,
            'nGenWeights' : 10358000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1400to1800_UL2018.root',
            'nevents' : 10550000.0,
            'nGenWeights' : 10550000.0,
            },
        'selection' : 'dijet',
        'XS' : 6.423e-01,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1800to2400_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1800to2400_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1800to2400_UL2017.root',
            'nevents' : 5191000.0,
            'nGenWeights' : 5191000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt1800to2400_UL2018.root',
            'nevents' : 5152000.0,
            'nGenWeights' : 5152000.0,
            },
        'selection' : 'dijet',
        'XS' : 8.746e-02,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2017.root',
            'nevents' : 2757000.0,
            'nGenWeights' : 2757000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt2400to3200_UL2018.root',
            'nevents' : 2901000.0,
            'nGenWeights' : 2901000.0,
            },
        'selection' : 'dijet',
        'XS' : 5.233e-03,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},

    'QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8' : 
    	{
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2016_preVFP.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2016.root',
            'nevents' : 0.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2017.root',
            'nevents' : 1000000.0,
            'nGenWeights' : 1000000.0,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmer' : [ '/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER' ],
            'skimmerHisto' : 'jetObservables_histograms_QCDPt3200toInf_UL2018.root',
            'nevents' : 934000.0,
            'nGenWeights' : 934000.0,
            },
        'selection' : 'dijet',
        'XS' : 1.351e-04,
        'label' : 'QCD Pythia8',
        'color': ROOT.kBlue
    	},    
    

}


def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


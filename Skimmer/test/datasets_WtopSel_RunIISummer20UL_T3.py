#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt

#cross checked XS vs AN2022_118_v11 (https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=TOP-22-014&tp=an&id=2632&ancode=TOP-22-014)
import ROOT

dictSamples = {

'SingleMuon' : {
        
        '2016_preVFP' :  {
            'miniAOD' : {
                'B' : '/SingleMuon/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'C' : '/SingleMuon/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'D' : '/SingleMuon/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'E' : '/SingleMuon/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                'F' : '/SingleMuon/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD',
                },
            'nanoAOD' : {
                'B': '/SingleMuon/kadatta-Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2_PFNanov2pt3-b43a2f9e6d8b6e76e19138f1f2a38dac/USER',
                'C': '/SingleMuon/kadatta-Run2016C-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt3-b43a2f9e6d8b6e76e19138f1f2a38dac/USER',
                'D': '/SingleMuon/kadatta-Run2016D-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt3-b43a2f9e6d8b6e76e19138f1f2a38dac/USER',
                'E': '/SingleMuon/kadatta-Run2016E-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt3-b43a2f9e6d8b6e76e19138f1f2a38dac/USER',
                'F': '/SingleMuon/kadatta-Run2016F-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt3-b43a2f9e6d8b6e76e19138f1f2a38dac/USER',
                },
            'skimmerHisto' : 'SingleMuonAll_2016UL_preVFP.root',
            'lumi' : 19495.076,
            },
        '2016' :  {
            'miniAOD' : {
                'F' : '/SingleMuon/Run2016F-UL2016_MiniAODv2-v2/MINIAOD',
                'G' : '/SingleMuon/Run2016G-UL2016_MiniAODv2-v2/MINIAOD',
                'H' : '/SingleMuon/Run2016H-UL2016_MiniAODv2-v2/MINIAOD',
                
                },
            'nanoAOD' : {
                'F' : '/SingleMuon/kadatta-Run2016F-UL2016_MiniAODv2-v2_PFNanov2pt3-a30de466bdcb4fc34de6668b02d7ad33/USER',
                'G' : '/SingleMuon/kadatta-Run2016G-UL2016_MiniAODv2-v2_PFNanov2pt3-a30de466bdcb4fc34de6668b02d7ad33/USER',
                'H' : '/SingleMuon/kadatta-Run2016H-UL2016_MiniAODv2-v2_PFNanov2pt3-a30de466bdcb4fc34de6668b02d7ad33/USER',
                },
            'skimmerHisto' : 'SingleMuonAll_2016UL.root',
            'lumi' : 16810.813,
            },
        '2017' :  {
            'miniAOD' : {
                'B' : '/SingleMuon/Run2017B-UL2017_MiniAODv2-v1/MINIAOD',
                'C' : '/SingleMuon/Run2017C-UL2017_MiniAODv2-v1/MINIAOD',
                'D' : '/SingleMuon/Run2017D-UL2017_MiniAODv2-v1/MINIAOD',
                'E' : '/SingleMuon/Run2017E-UL2017_MiniAODv2-v1/MINIAOD',
                'F' : '/SingleMuon/Run2017F-UL2017_MiniAODv2-v1/MINIAOD',
                },
            'nanoAOD' : {
                'B': '/SingleMuon/kadatta-Run2017B-UL2017_MiniAODv2-v1_PFNanov2pt3-80b78937c3a3b7c462c5f2c6e6c26583/USER',
                'C': '/SingleMuon/kadatta-Run2017C-UL2017_MiniAODv2-v1_PFNanov2pt3-80b78937c3a3b7c462c5f2c6e6c26583/USER',
                'D': '/SingleMuon/kadatta-Run2017D-UL2017_MiniAODv2-v1_PFNanov2pt3-80b78937c3a3b7c462c5f2c6e6c26583/USER',
                'E': '/SingleMuon/kadatta-Run2017E-UL2017_MiniAODv2-v1_PFNanov2pt3-80b78937c3a3b7c462c5f2c6e6c26583/USER',
                'F': '/SingleMuon/kadatta-Run2017F-UL2017_MiniAODv2-v1_PFNanov2pt3-80b78937c3a3b7c462c5f2c6e6c26583/USER',
                },
            'skimmerHisto' : 'SingleMuonAll_2017UL.root',
            'lumi' : 41475.262,
            },
        '2018' :  {
            'miniAOD' : {
                'A' : '/SingleMuon/Run2018A-UL2018_MiniAODv2_GT36-v1/MINIAOD',
                'B' : '/SingleMuon/Run2018B-UL2018_MiniAODv2_GT36-v1/MINIAOD',
                'C' : '/SingleMuon/Run2018C-UL2018_MiniAODv2_GT36-v2/MINIAOD',
                'D' : '/SingleMuon/Run2018D-UL2018_MiniAODv2_GT36-v1/MINIAOD',
                },
            'nanoAOD' : {
                'A' : '/SingleMuon/kadatta-Run2018A-UL2018_MiniAODv2_GT36-v1_PFNano-05315415ba22ce0eb9676f5f753b2b2a/USER',
                'B' : '/SingleMuon/kadatta-Run2018B-UL2018_MiniAODv2_GT36-v1_PFNano-05315415ba22ce0eb9676f5f753b2b2a/USER',
                'C' : '/SingleMuon/kadatta-Run2018C-UL2018_MiniAODv2_GT36-v2_PFNano-05315415ba22ce0eb9676f5f753b2b2a/USER',
                'D' : '/SingleMuon/kadatta-Run2018D-UL2018_MiniAODv2_GT36-v1_PFNanov2pt2-05315415ba22ce0eb9676f5f753b2b2a/USER',
                },
            'skimmerHisto' : 'SingleMuonAll_2018UL.root',
            'lumi' : 59673.938,
            },
        'color': ROOT.kWhite,
        'selection' : 'Wtop'
    },

    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 833.9,#7.568e+02,
        'selection' : 'Wtop',
        'label' : 't#bar{t} amcatnloFXFX' ,
        'color': ROOT.kBlue
    },


    'TTJets_TuneCP5_13TeV-madgraphMLM-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTJets_madgraphMLM-pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v2-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v2-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTJets_madgraphMLM-pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTJets_madgraphMLM-pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTJets_madgraphMLM-pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 4.761e+02,
        'selection' : 'Wtop',
        'label' : 't#bar{t} amcatnloFXFX' ,
        'color': ROOT.kBlue
    },

    
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt4-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt4-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        #'XS' : 364.328, #43.8% of 831.8 a la PDG https://pdg.lbl.gov/2021/reviews/rpp2021-rev-top-quark.pdf
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 't#bar{t} powheg' ,
        'color': ROOT.kMagenta
    },

    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTToHadronic_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToHadronic_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToHadronic_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ], #pfnano has 1 extra file not counted in DAS filelist (stupid crab publishing issues c.a. 08/2022), ensure number of events run on ends up matching nevents in skimmed pfnano for a full successful run over the entire dataset
            'nanoAOD' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToHadronic_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'XS' : 381.0923,#45.7% of 833.9 a la PDG
        'selection' : 'Wtop',
        'label' : 't#bar{t} powheg' ,
        'color': ROOT.kMagenta
    },

    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' :  0.,
            },
        '2016' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' :  0.,
            },
        '2018' :  {
            'miniAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'selection' : 'Wtop',
        'XS' : 87.5595,#10.5% of 833.9  a la PDG,
        'label' : 't#bar{t} DL' ,
        'color': ROOT.kCyan+1
    },

    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 9949891586212.455078,
            },
        '2016' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 10484986870276.357422,
            },
        '2017' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 9949891586212.455078,
            },
        '2018' :  {
            'miniAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 10484986870276.357422,
            },
        'selection' : 'Wtop',
        'XS' : 5.359e+04, #xsec tool
        'label' : 'WJets' ,
        'color': ROOT.kOrange-2
    },

    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_leptonDecays_amcatnlo_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 70439807.475972,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_leptonDecays_amcatnlo_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 64014403.924644,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_leptonDecays_amcatnlo_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 70439807.475972,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_leptonDecays_amcatnlo_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 64014403.924644,
            },
        'selection' : 'Wtop',
        'XS' : 3.549e+00,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    
    'ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 656310119.274000,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 656310119.274000,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'selection' : 'Wtop',
        'XS' : 134.2,#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef #1.197e+02, #xsec tool
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 0.,
            },
        'selection' : 'Wtop',
        'XS' : 80.0,#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef #7.175e+01,xsec tool
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 331239851.967700,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 328321574.253000,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 331239851.967700,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 328321574.253000,
            },
        'selection' : 'Wtop',
        'XS' : 21.609, # 79.3 for t+tbar from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef, multiplied by (43.8+10.5)% a la PDG BRs,# 3.245e+01, #xsec tool
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 303960792.435600,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 296655487.557600,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 303960792.435600,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 296655487.557600,
            },
        'selection' : 'Wtop',
        'XS' : 21.609, # 79.3 for t+tbar from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef, multiplied by (43.8+10.5)% a la PDG BRs,# 3.245e+01: xsec tool
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'WW_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'WW_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 15634116.199514,
            },
        '2016' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'WW_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7959266.140092,
            },
        '2017' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'WW_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7959266.140092,
            },
        '2018' :  {
            'miniAOD' : [ '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/WW_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'WW_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7959266.140092,
            },
        'XS' : 7.592e+01, #xsec tool
        'selection' : 'Wtop',
        'label' : 'Diboson' ,
        'color': ROOT.kPink-1
    },
    'ZZ_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'ZZ_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 2706000.000000,
            },
        '2016' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ZZ_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 2000000.000000,
            },
        '2017' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'ZZ_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 2706000.000000,
            },
        '2018' :  {
            'miniAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ZZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ZZ_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 2000000.000000,
            },
        'selection' : 'Wtop',
        'XS' : 1.214e+01,#xsec tool
        'label' : 'Diboson' ,
        'color': ROOT.kPink-1
    },
    'WZ_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'WZ_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7889000.000000,
            },
        '2016' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'WZ_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 3976000.000000,
            },
        '2017' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER' ],
            'skimmerHisto' : 'WZ_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7889000.000000,
            },
        '2018' :  {
            'miniAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/WZ_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'WZ_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 3976000.000000,
            },
        'XS' : 2.758e+01, #xsec tool
        'selection' : 'Wtop',
        'label' : 'Dibosons' ,
        'color': ROOT.kPink-1
    },
    'QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 36027673.00,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER'],
            'skimmer' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER'],
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 35868369.000000,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 36027673.00,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER'],
            'skimmer' : [ '/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER'],
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 35868369.000000,
            },
    'XS' : 7.020e+03,#all QCD from XS tool 
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 494796.07288,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 49232711.259666,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 494796.07288,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d733' ],
            'skimmer' : [ '/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d733' ],
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 49232711.259666,
            },
        'XS' : 6.197e+02,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 517383.713837,
        },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 20585907.002748,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 517383.713837,
        },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 20585907.002748,
            },
        'XS' : 5.908e+01,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 17318812.571925,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 16983542.512253,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 17318812.571925,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 16983542.512253,
            },
        'XS' : 1.820e+01,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 16962615.000000,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_UL2016.root',
            'nevents' : 1.,
           'nGenWeights' : 17178227.000000,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 16962615.000000,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_UL2018.root',
            'nevents' : 1.,
           'nGenWeights' : 17178227.000000,
            },
       'XS' : 3.274e+00,
       'selection' : 'Wtop',
       'label' : 'QCD' ,
       'color': ROOT.kRed-7
    },
    'QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 14642553.00,
            },
        '2016' :  {
            'miniAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 14556168.000000,
            },
        '2017' :  {
            'miniAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 14642553.00,
            },
        '2018' :  {
            'miniAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 14556168.000000,
            },
        'XS' : 1.078e+00,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    
    # MC variations
    'varTTToSemileptonic_TuneCP5_erdON_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 8032868875.856007,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 8066810085.685998,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 8032868875.856007,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5_erdON_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 8066810085.685998,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 'CR (erdON)' ,
        'color': ROOT.kBlue+3
    },

    'varTTToSemileptonic_TuneCP5Down_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7685702947.076004,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7783384223.075996,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7685702947.076004,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5down_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7783384223.075996,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 'Tune CP5 down' ,
        'color': ROOT.kBlue+7
    },


    'varTTToSemileptonic_TuneCP5Up_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7926747573.580008,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 8024917255.960003,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7926747573.580008,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 8024917255.960003,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 'TuneCP5 up' ,
        'color': ROOT.kBlue+7
    },

    'varTTToSemileptonic_hdampDown_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7560511524.540000,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7407695568.779999,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7560511524.540000,
            },
        '2018' :  {#could re-up this dataset... missing nearly 15 million events (~7.8%) vs. MiniAOD
            'miniAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7407695568.779999,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : '#h_{damp} down' ,
        'color': ROOT.kBlue+5
    },


    'varTTToSemileptonic_hdampUp_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7947260377.072005,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7743168497.928005,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7947260377.072005,
            },
        '2018' :  { #could also be topped up, missing ~10M events (~5.3%) vs. MiniAOD
            'miniAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7743168497.928005,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : '#h_{damp} up' ,
        'color': ROOT.kBlue+5
    },
   
     'varTTToSemileptonic_TuneCP5CR1_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 8105362943.832006,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 8196473485.551998,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 8105362943.832006,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER'],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 8196473485.551998,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 'CR1' ,
        'color': ROOT.kBlue+3
    },
    
     'varTTToSemileptonic_TuneCP5CR2_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7685702947.076004,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7834556470.664001,
            },
        '2017' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7685702947.076004,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt4-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_TuneCP5CR2_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt4-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7834556470.664001,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : 'CR2' ,
        'color': ROOT.kBlue+3
    },
    
    'varTTToSemileptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7519104783.423996,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7343389413.151999,
            },
        '2017' :  {#sh/could be topped up, is missing ~9M events (~9.62%) vs. MiniAOD
            'miniAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7519104783.423996,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7343389413.151999,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : '#m_{top}=171.5' ,
        'color': ROOT.kBlue-3
    },
    
    'varTTToSemileptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 7259862511.980000,
            },
        '2016' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 7159151017.746006,
            },
        '2017' :  {#sh/could be topped up, mising 11M events (~10%) vs. MiniAOD
            'miniAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL17PFNanov2pt2-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 7259862511.980000,
            },
        '2018' :  {
            'miniAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' ],
            'nanoAOD' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 7159151017.746006,
            },
        'XS' : 365.2482, #43.8% of 833.9 a la PDG '22/'23 https://pdg.lbl.gov/2023/reviews/rpp2023-rev-top-quark.pdf
        'selection' : 'Wtop',
        'label' : '#m_{top}=173.5' ,
        'color': ROOT.kBlue-3
    },
   

}

def checkDict( string, dictio ):
    return next(dictio[k] for k in list(dictio.keys()) if k.startswith(string))

"""    
    
    'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8' : {
        '2016_preVFP' :  {
            'miniAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v2-5b62f87422b4b3633317ac4b6822e17c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_hadronicDecays_amcatnlo_pythia8_UL2016_preVFP.root',
            'nevents' : 1.,
            'nGenWeights' : 70439807.475972,
            },
        '2016' :  {
            'miniAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v2-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v2-36773dbd6a9d41d00afbee31f82bc069/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_hadronicDecays_amcatnlo_pythia8_UL2016.root',
            'nevents' : 1.,
            'nGenWeights' : 64014403.924644,
            },
        '2017' :  {
            'miniAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v2-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_hadronicDecays_amcatnlo_pythia8_UL2017.root',
            'nevents' : 1.,
            'nGenWeights' : 70439807.475972,
            },
        '2018' :  {
            'miniAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM' ],
            'nanoAOD' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmer' : [ '/ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER' ],
            'skimmerHisto' : 'ST_s-channel_4f_hadronicDecays_amcatnlo_pythia8_UL2018.root',
            'nevents' : 1.,
            'nGenWeights' : 64014403.924644,
            },
        'selection' : 'Wtop',
        'XS' :  7.104e+00,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    }, 
"""
#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples = {
    
    'SingleMuon' : {
        
        '2017' :  {
            'skimmerHisto' : 'SingleMuonAll_2017UL.root',
            'lumi' : 41472.052,
            },
        '2018' :  {
           'skimmerHisto' : 'SingleMuonAll_2018UL.root',
            'lumi' : 59816.0#56385.82159,
            },
        'all' :  {
           'skimmerHisto' : 'SingleMuonAll_UL17-18.root',
           'lumi' : 59816.0#56385.82159,
            },
        'color': ROOT.kWhite,
        'selection' : 'Wtop'
    },

    'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_2017UL.root',
            'nGenWeights' : 286629789811.7677,
            },
        '2018' :  {
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_2018UL.root',
            'nGenWeights' : 297902725633.38257,
            },
        
        'all' :  {
            'skimmerHisto' : 'TTJets_amcatnloFXFX-pythia8_UL17-18.root',
            'nGenWeights' : 584532515445.1503,
            },
    
        'XS' : 7.535e+02, #some notes say 722.8 and a few say 831.76?
        'selection' : 'Wtop',
        'label' : 't#bar{t} amcatnloFXFX' ,
        'color': ROOT.kBlue
    },

    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_NoTopReWt.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_NoTopReWt.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_nom_UL17-18_NoTopReWt.root',
            'nGenWeights' : 63770862770.80198,
              },
        'XS' : 301.64,#364.0,#5.34,
        'selection' : 'Wtop',
        'label' : 't#bar{t} powheg' ,
        'color': ROOT.kMagenta
    },

    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_2017UL.root',
            'nGenWeights' :  4759777104.469413,
            },
        '2018' :  {
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_2018UL.root',
            'nGenWeights' : 5063754431.9676,
            },
        'all' :  {
            'skimmerHisto' : 'TTTo2L2Nu_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 9823531536.437012,
            },
        'selection' : 'Wtop',
        'XS' : 72.83,#26,
        'label' : 't#bar{t} DL' ,
        'color': ROOT.kCyan+1
    },

    'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_2017UL.root',
            'nGenWeights' : 9949891586212.455078,
            },
        '2018' :  {
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_2018UL.root',
            'nGenWeights' : 10484986870276.357422,
            },
        'all' :  {
            'skimmerHisto' : 'WJetsToLNu_madgraphMLM_pythia8_UL17-18.root',
            'nGenWeights' : 20434878456488.812,
            },
        'selection' : 'Wtop',
        'XS' : 5.371e+04,
        'label' : 'WJets' ,
        'color': ROOT.kOrange-2
    },
    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ST_s-channel_amcatnlo_pythia8_2017UL.root',
            'nGenWeights' : 70439807.475972,
            },
        '2018' :  {
            'skimmerHisto' : 'ST_s-channel_amcatnlo_pythia8_2018UL.root',
            'nGenWeights' : 64014403.924644,
            },
        'all' :  {
            'skimmerHisto' : 'ST_s-channel_amcatnlo_pythia8_UL17-18.root',
            'nGenWeights' : 134454211.400616,
            },
        'selection' : 'Wtop',
        'XS' : 3.549,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_2017UL.root',
            'nGenWeights' : 656310119.274000,
            },
        '2018' :  {
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_2018UL.root',
            'nGenWeights' : 683747124.990000,
            },
        'all' :  {
            'skimmerHisto' : 'ST_t-channel_top_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 1340057244.264,
            },
        'selection' : 'Wtop',
        'XS' : 1.197e+02,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_2017UL.root',
            'nGenWeights' : 303960792.435600,
            },
        '2018' :  {
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_2018UL.root',
            'nGenWeights' : 271219283.652000,
            },
        'all' :  {
            'skimmerHisto' : 'ST_t-channel_antitop_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 575180076.0876,
            },
        'selection' : 'Wtop',
        'XS' : 7.174e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_2017UL.root',
            'nGenWeights' : 331239851.967700,
            },
        '2018' :  {
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_2018UL.root',
            'nGenWeights' : 328321574.253000,
            },
        'all' :  {
            'skimmerHisto' : 'ST_tW_top_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 659561426.2207,
            },
        'selection' : 'Wtop',
        'XS' : 3.245e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_2017UL.root',
            'nGenWeights' : 303960792.435600,
            },
        '2018' :  {
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_2018UL.root',
            'nGenWeights' : 296655487.557600,
            },
        'all' :  {
            'skimmerHisto' : 'ST_tW_antitop_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 600616279.9932001,
            },
        'selection' : 'Wtop',
        'XS' : 3.251e+01,
        'label' : 'Single top' ,
        'color': ROOT.kGreen+2,
    },
    'WW_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'WW_pythia8_2017UL.root',
            'nGenWeights' : 15634116.199514,
            },
        '2018' :  {
            'skimmerHisto' : 'WW_pythia8_2018UL.root',
            'nGenWeights' : 7959266.140092,
            },
        'all' :  {
            'skimmerHisto' : 'WW_pythia8_UL17-18.root',
            'nGenWeights' : 23593382.339606002,
            },
        'selection' : 'Wtop',
        'XS': 7.56e+01,
        'label' : 'Diboson' ,
        'color': ROOT.kPink-1
    },
    'ZZ_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'ZZ_pythia8_2017UL.root',
            'nGenWeights' : 2706000.000000,
            },
        '2018' :  {
            'skimmerHisto' : 'ZZ_pythia8_2018UL.root',
            'nGenWeights' : 2000000.000000,
            },
        'all' :  {
            'skimmerHisto' : 'ZZ_pythia8_UL17-18.root',
            'nGenWeights' : 4706000.000000,
            },
        'selection' : 'Wtop',
        'XS' : 1.210e+01,
        'label' : 'Diboson' ,
        'color': ROOT.kPink-1
    },
    'WZ_TuneCP5_13TeV-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'WZ_pythia8_2017UL.root',
            'nGenWeights' : 7889000.000000,
            },
        '2018' :  {
            'skimmerHisto' : 'WZ_pythia8_2018UL.root',
            'nGenWeights' : 3976000.000000,
            },
        'all' :  {
            'skimmerHisto' : 'WZ_pythia8_UL17-18.root',
            'nGenWeights' : 11865000.0,
            },
        'XS' : 2.751e+01,
        'selection' : 'Wtop',
        'label' : 'Dibosons' ,
        'color': ROOT.kPink-1
    },
    'QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_2017UL.root',
            'nGenWeights' : 36027673.00,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_2018UL.root',
            'nGenWeights' : 35868369.000000,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_170to300_pythia8_UL17-18.root',
            'nGenWeights' : 71896042.0,
            },
        'XS' : 7.014e+03,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_2017UL.root',
            'nGenWeights' : 494796.07288,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_2018UL.root',
            'nGenWeights' : 49232711.259666,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_300to470_pythia8_UL17-18.root',
            'nGenWeights' : 49727507.332546,
            },
        'XS' : 6.205e+02,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_2017UL.root',
            'nGenWeights' : 517383.713837,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_2018UL.root',
            'nGenWeights' : 20585907.002748,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_470to600_pythia8_UL17-18.root',
            'nGenWeights' : 21103290.716585003,
            },
        'XS' : 5.910e+01,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_2017UL.root',
            'nGenWeights' : 17318812.571925,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_2018UL.root',
            'nGenWeights' : 16983542.512253,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_600to800_pythia8_UL17-18.root',
            'nGenWeights' : 34302355.0842,
            },
        'XS' : 1.825e+01,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_2017UL.root',
            'nGenWeights' : 16962615.000000,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_2018UL.root',
            'nGenWeights' : 17178227.000000,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_800to1000_pythia8_UL17-18.root',
            'nGenWeights' : 34140842.0,
            },
        'XS' : 3.287e+00,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },
    'QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_2017UL.root',
            'nGenWeights' : 14642553.00,
            },
        '2018' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_2018UL.root',
            'nGenWeights' : 14556168.000000,
            },
        'all' :  {
            'skimmerHisto' : 'QCD_Muen_Pt_1000toInf_pythia8_UL17-18.root',
            'nGenWeights' : 29198721,
            },
        'XS' : 1.087e+00,
        'selection' : 'Wtop',
        'label' : 'QCD' ,
        'color': ROOT.kRed-7
    },    
    
    # MC variations
    'varTTToSemileptonic_TuneCP5_erdON_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_2017UL.root',
            'nGenWeights' : 8032868875.856008,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_2018UL.root',
            'nGenWeights' : 8066810085.685997,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_erdON_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 16099678961.542004,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'CR Model (erdON)' ,
        'color': ROOT.kBlue-3
    },

    'varTTToSemileptonic_TuneCP5Down_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7685702947.076002,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7783384223.075996,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5down_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 15469087170.151999,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'Tune' ,
        'color': ROOT.kBlue+7
    },


    'varTTToSemileptonic_TuneCP5Up_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7926747573.580006,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_2018UL.root',
            'nGenWeights' : 8024917255.960002,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_TuneCP5up_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 15951664829.540009,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'Tune' ,
        'color': ROOT.kBlue+7
    },

    'varTTToSemileptonic_hdampDOWN_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
             'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7819513446.420001,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7422136644.780002,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_hdampDOWN_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 15241650091.200003,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : '#h_{damp}' ,
        'color': ROOT.kBlue+7
    },


    'varTTToSemileptonic_hdampUP_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7998353023.320009,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7786479144.880011,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_hdampUP_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 15784832168.20002,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : '#h_{damp}' ,
        'color': ROOT.kBlue+7
    },
    
    #'varTTToSemileptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8' : {
    #    '2017' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop169p5_powheg_pythia8_2017UL.root',
    #        'nGenWeights' : 6297999091.747996,
    #        },
    #    '2018' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop169p5_powheg_pythia8_2018UL.root',
    #        'nGenWeights' : 6307297092.055999,
    #        },
    #    'all' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop169p5_powheg_pythia8_UL17-18.root',
    #        'nGenWeights' : 12605296183.803995,
    #        },
    #    'XS' : 301.64,
    #    'selection' : 'Wtop',
    #    'label' : '#m_{top}=169.5' ,
    #     'color': ROOT.kBlue+7
    #},
    
    'varTTToSemileptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7519104783.423997,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7343389413.151999,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop171p5_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 14862494196.575996,
            },
        'XS' : 309.977,
        'selection' : 'Wtop',
        'label' : '#m_{top}=171.5' ,
        'color': ROOT.kBlue+7
    },
    
    'varTTToSemileptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_2017UL.root',
            'nGenWeights' : 7259862511.980000,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7159151017.746006,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_mtop173p5_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 14419013529.726006,
            },
        'XS' : 293.5154,
        'selection' : 'Wtop',
        'label' : '#m_{top}=173.5' ,
        'color': ROOT.kBlue+7
    },
    
    #'varTTToSemileptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8' : {
    #    '2017' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop175p5_powheg_pythia8_2017UL.root',
    #        'nGenWeights' : 5129149295.824001,
    #        },
    #    '2018' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop175p5_powheg_pythia8_2018UL.root',
    #        'nGenWeights' : 5418843197.792000,
    #        },
    #    'all' :  {
    #        'skimmerHisto' : 'TTToSemiLeptonic_mtop175p5_powheg_pythia8_UL17-18.root',
    #         'nGenWeights' : 10547992493.616001,
    #          },
    #    'XS' : 301.64,
    #    'selection' : 'Wtop',
    #    'label' : '#m_{top}=175.5' ,
    #    'color': ROOT.kBlue+7
    #},
    
     'varTTToSemileptonic_TuneCP5CR1_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_2017UL.root',
            'nGenWeights' : 8105362943.832007,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_2018UL.root',
            'nGenWeights' : 8196473485.551998,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR1_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 16301836429.384007,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'CR1' ,
        'color': ROOT.kBlue+7
    },
    
     'varTTToSemileptonic_TuneCP5CR2_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_2017UL.root',
            'nGenWeights' : 8180582986.796001,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_2018UL.root',
            'nGenWeights' : 7834556470.664003,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemiLeptonic_TuneCP5CR2_powheg_pythia8_UL17-18.root',
            'nGenWeights' : 16015139457.460005,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'CR2' ,
        'color': ROOT.kBlue+7
    },
    
    'sysTTToSemileptonic_puWeight_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_puWeight.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_puWeight.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_puWeight_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'puWeight' ,
        'color': ROOT.kMagenta
    },
    
    'sysTTToSemileptonic_jer_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_jer.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_jer.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_jer_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'JER' ,
        'color': ROOT.kMagenta
    },
    
    'sysTTToSemileptonic_jesTotal_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_jesTotal.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_jesTotal.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_jesTotal_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'JESTotal' ,
        'color': ROOT.kMagenta
    },
    
    'sysTTToSemileptonic_isrWeight_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_isrWeight.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_isrWeight.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_isrWeight_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'isrWeight' ,
        'color': ROOT.kMagenta
    },
    
    'sysTTToSemileptonic_fsrWeight_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_fsrWeight.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_fsrWeight.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_fsrWeight_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'fsrWeight' ,
        'color': ROOT.kMagenta
    },
    
    'sysTTToSemileptonic_pdfWeight_TuneCP5_13TeV-powheg-pythia8' : {
        '2017' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2017UL_fsrWeight.root',
            'nGenWeights' : 32392787477.907993,
            },
        '2018' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_2018UL_fsrWeight.root',
            'nGenWeights' : 31378075292.893986,
            },
        'all' :  {
            'skimmerHisto' : 'TTToSemileptonic_powheg_pythia8_fsrWeight_UL17-18.root',
            'nGenWeights' : 63770862770.80198,
            },
        'XS' : 301.64,
        'selection' : 'Wtop',
        'label' : 'fsrWeight' ,
        'color': ROOT.kMagenta
    },
    
}


def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


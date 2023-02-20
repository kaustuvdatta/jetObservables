#!/usr/bin/env python
### For xsec calculation: grep miniAOD datasets.py | awk '{ print $4 }' | sed "s/'//g" | sed  's/"//g' | sed 's/\,//g' > calculateXSectionAndFilterEfficiency/datasets.txt
import ROOT

dictSamples={
"JetHT": 
	{
"selection": "dijet", 
	"2016_preVFP": 
			
{"miniAOD": 
					{"B": "/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD", 
			 		 "C": "/JetHT/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD", 
					 "D": "/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD", 
					 "E": "/JetHT/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD", 
					 "F": "/JetHT/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD"}, 
			 
"nanoAOD": 
			 		{"B": "/JetHT/kadatta-Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2_PFNanov2pt2-b43a2f9e6d8b6e76e19138f1f2a38dac/USER", 
			 		 "C": "/JetHT/kadatta-Run2016C-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt2-b43a2f9e6d8b6e76e19138f1f2a38dac/USER", 
			 		 "D": "/JetHT/kadatta-Run2016D-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt2-b43a2f9e6d8b6e76e19138f1f2a38dac/USER", 
			 		 "E": "/JetHT/kadatta-Run2016E-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt2-b43a2f9e6d8b6e76e19138f1f2a38dac/USER", 
			 		 "F": "/JetHT/kadatta-Run2016F-HIPM_UL2016_MiniAODv2-v2_PFNanov2pt2-b43a2f9e6d8b6e76e19138f1f2a38dac/USER"}, 
			 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_JetHT2016_HIPM_All.root", 
			 "lumi": 19497.905, 
			 "triggerList": 
			 	{"AK8PFJet80": 27661.64, 
			 	 "AK8PFJet140": 2689.77, 
			 	 "AK8PFJet200": 278.62, 
			 	 "AK8PFJet260": 56.78, 
			 	 "AK8PFJet320": 19.59, 
			 	 "AK8PFJet400": 6.41, 
			 	 "AK8PFJet450": 1.0, 
			 	 "AK8PFJet500": 1.0}, 
			 
"t3_dirs": 
			  	{"B": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160438/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160438/0001/"], 
			 	 "C": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016C-HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160327/0000/"], 
			 	 "D": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016D-HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160108/0000/"], 
			 	 "E": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016E-HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_155936/0000/"], 
			 	 "F": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016F-HIPM_UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160215/0000/"]}}, 
	
"2016": 
		
{"miniAOD": 
				{"F": "/JetHT/Run2016F-UL2016_MiniAODv2-v2/MINIAOD", 
			 	 "G": "/JetHT/Run2016G-UL2016_MiniAODv2-v2/MINIAOD", 
			 	 "H": "/JetHT/Run2016H-UL2016_MiniAODv2-v2/MINIAOD"}, 
		 
"nanoAOD": 
		 		{"F": "/JetHT/kadatta-Run2016F-UL2016_MiniAODv2-v2_PFNanov2pt2-a30de466bdcb4fc34de6668b02d7ad33/USER", 
		 	 	 "G": "/JetHT/kadatta-Run2016G-UL2016_MiniAODv2-v2_PFNanov2pt2-a30de466bdcb4fc34de6668b02d7ad33/USER", 
		 	 	 "H": "/JetHT/kadatta-Run2016H-UL2016_MiniAODv2-v2_PFNanov2pt2-a30de466bdcb4fc34de6668b02d7ad33/USER"}, 
		 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_JetHT2016_All.root", 
		 "lumi": 16810.813, 
		 "triggerList": 
		 		{"AK8PFJet80": 41995.5, 
		 	 	 "AK8PFJet140": 4319.21, 
		 		 "AK8PFJet200": 652.64,
		 		 "AK8PFJet260": 75.17, 
		 		 "AK8PFJet320": 25.0, 
		 		 "AK8PFJet400": 8.47, 
		 		 "AK8PFJet450": 1.0, 
		 		 "AK8PFJet500": 1.0}, 
		 
"t3_dirs": 
		 	 	{"F": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016F-UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_155937/0000/"], 
		  	 	 "G": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016G-UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160110/0000/"], 
		 	 	 "H": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2016H-UL2016_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160215/0000/"]}}, 
	
"2017": 
		
{"miniAOD": 
				{"B": "/JetHT/Run2017B-UL2017_MiniAODv2-v1/MINIAOD", 
		  		 "C": "/JetHT/Run2017C-UL2017_MiniAODv2-v1/MINIAOD", 
		  		 "D": "/JetHT/Run2017D-UL2017_MiniAODv2-v1/MINIAOD", 
		  		 "E": "/JetHT/Run2017E-UL2017_MiniAODv2-v1/MINIAOD", 
		  		 "F": "/JetHT/Run2017F-UL2017_MiniAODv2-v1/MINIAOD"}, 
		  
"nanoAOD": 
		  		{"B": "/JetHT/kadatta-Run2017B-UL2017_MiniAODv2-v1_PFNanov2pt2-80b78937c3a3b7c462c5f2c6e6c26583/USER", 
		  		 "C": "/JetHT/kadatta-Run2017C-UL2017_MiniAODv2-v1_PFNanov2pt2-80b78937c3a3b7c462c5f2c6e6c26583/USER", 
		  		 "D": "/JetHT/kadatta-Run2017D-UL2017_MiniAODv2-v1_PFNanov2pt2-80b78937c3a3b7c462c5f2c6e6c26583/USER", 
		  		 "E": "/JetHT/kadatta-Run2017E-UL2017_MiniAODv2-v1_PFNanov2pt2-80b78937c3a3b7c462c5f2c6e6c26583/USER", 
		  		 "F": "/JetHT/kadatta-Run2017F-UL2017_MiniAODv2-v1_PFNanov2pt2-80b78937c3a3b7c462c5f2c6e6c26583/USER"}, 
		  

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_JetHT2017_All.root", 
		  "lumi": 41473.12, 
		  "triggerList": 
		  		{"AK8PFJet80": 16419.91,
		  		 "AK8PFJet140": 1560.15, 
		  		 "AK8PFJet200": 219.7, 
		  		 "AK8PFJet260": 88.43, 
		  		 "AK8PFJet320": 33.827, 
		  		 "AK8PFJet400": 5.404, 
		  		 "AK8PFJet450": 4.299, 
		  		 "AK8PFJet500": 1.0}, 
		  
"t3_dirs": 
		  		{"B": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017B-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160110/0000/"], 
		  		 "C": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017C-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_155935/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017C-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_155935/0001/"], 
		  		 "D": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017D-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160435/0000/"], 
		  		 "E": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017E-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160327/0000/"], 
		  		 "F": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017F-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160216/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2017F-UL2017_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160216/0001/"]}}, 
	
"2018": 
		
{"miniAOD": 
				{"A": "/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD", 
				 "B": "/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD", 
				 "C": "/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD", 
				 "D": "/JetHT/Run2018D-UL2018_MiniAODv2-v1/MINIAOD"}, 
		 
"nanoAOD": 
		 		{"A": "/JetHT/kadatta-Run2018A-UL2018_MiniAODv2-v1_PFNanov2pt2-384370c5705060bbc71e621a12dca622/USER", 
		 		 "B": "/JetHT/kadatta-Run2018B-UL2018_MiniAODv2-v1_PFNanov2pt2-384370c5705060bbc71e621a12dca622/USER", 
		 		 "C": "/JetHT/kadatta-Run2018C-UL2018_MiniAODv2-v1_PFNanov2pt2-384370c5705060bbc71e621a12dca622/USER", 
		 		 "D": "/JetHT/kadatta-Run2018D-UL2018_MiniAODv2-v2_PFNanov2pt2-384370c5705060bbc71e621a12dca622/USER"}, 
		 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_JetHT2018_All.root", 
		 "lumi": 59817.407, 
		 "triggerList": 
		 		{"AK8PFJet80": 27585.68, 
		 		 "AK8PFJet140": 1268.74, 
		 		 "AK8PFJet200": 295.24, 
		 		 "AK8PFJet260": 127.99, 
		 		 "AK8PFJet320": 48.17, 
		 		 "AK8PFJet400": 16.08, 
		 		 "AK8PFJet450": 8.09, 
		 		 "AK8PFJet500": 1.0}, 
		 
"t3_dirs": 
		 		{"A": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018A-UL2018_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160216/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018A-UL2018_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160216/0001/"], 
		 		 "B": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018B-UL2018_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_155936/0000/"], 
		 		 "C": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018C-UL2018_MiniAODv2-v1_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160111/0000/"], 
		 		 "D": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160327/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160327/0001/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/JetHT/Run2018D-UL2018_MiniAODv2-v2_PFNano_jetObsSkim_Summer20_nomTreesV3/221114_160327/0002/"]}},
	 "color": 0}, 
"QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
	{"2016_preVFP"	: 
		{"miniAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
		 "nanoAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
		 "skimmer": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
		 
"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT100to200_nom_UL2016_preVFP.root", 
		 "nevents": 0.0, 
	     "nGenWeights": 28066870, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_141028/0000/"], 
"nevents": 28066870.0, 
"MCScaling": 16422.582005047232}, 
"2016": 
{"miniAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT100to200_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 23519301, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_141028/0000/"], 
"nevents": 23519301.0, 
"MCScaling": 16897.084625091535}, 
"2017": 
{"miniAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT100to200_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 54760426, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212054/0000/"], 
"nevents": 54760426.0, 
"MCScaling": 17903.88841752254}, 
"2018": 
{"miniAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT100to200_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 84083628, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212054/0000/"], 
"nevents": 84083628.0, 
"MCScaling": 16817.58429215257}, 
"selection": "dijet", 
"XS": 23640000.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT200to300_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 18273591, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_141439/0000/"], 
"nevents": 18273591.0, 
"MCScaling": 1654.9155913033185}, 
"2016": 
{"miniAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT200to300_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 17569141, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_141437/0000/"], 
"nevents": 17569141.0, 
"MCScaling": 1484.0549667738449}, 
"2017": 
{"miniAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanov2pt3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT200to300_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 42714435, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212302/0000/"], 
"nevents": 42714435.0, 
"MCScaling": 1505.9267228982428}, 
"2018": 
{"miniAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt2-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT200to300_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 57336623, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221106_212302/0000/"], 
"nevents": 57336623.0, 
"MCScaling": 1618.1071259986832}, 
"selection": "dijet", 
"XS": 1551000.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT300to500_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 15341307, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_140833/0000/"], 
"nevents": 15341307.0, 
"MCScaling": 411.65797865201444}, 
"2016": 
{"miniAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT300to500_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 16747056, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_140833/0000/"], 
"nevents": 16747056.0, 
"MCScaling": 325.13310582468944}, 
"2017": 
{"miniAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT300to500_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 43589739, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211952/0000/"], 
"nevents": 43589739.0, 
"MCScaling": 308.17214959694985}, 
"2018": 
{"miniAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanov2pt3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT300to500_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 61705174, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_211952/0000/"], 
"nevents": 61705174.0, 
"MCScaling": 313.99081910538}, 
"selection": "dijet", 
"XS": 323900.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT500to700_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 15775001, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_141247/0000/"], 
"nevents": 15775001.0, 
"MCScaling": 37.512607241673074}, 
"2016": 
{"miniAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT500to700_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 15222746, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_141246/0000/"], 
"nevents": 15222746.0, 
"MCScaling": 33.5161720855094}, 
"2017": 
{"miniAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT500to700_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 36194860, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212157/0000/"], 
"nevents": 36194860.0, 
"MCScaling": 34.775909949644785}, 
"2018": 
{"miniAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT500to700_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 49184771, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212157/0000/"], 
"nevents": 49184771.0, 
"MCScaling": 36.910984142835595}, 
"selection": "dijet", 
"XS": 30350.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT700to1000_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 15808790, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_141855/0000/"], 
"nevents": 15808790.0, 
"MCScaling": 7.930494943003227}, 
"2016": 
{"miniAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT700to1000_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 13905714, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_141855/0000/"], 
"nevents": 13905714.0, 
"MCScaling": 7.7733173276827054}, 
"2017": 
{"miniAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT700to1000_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 34051754, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212537/0000/"], 
"nevents": 34051754.0, 
"MCScaling": 7.83137813106485}, 
"2018": 
{"miniAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM", 
"nanoAOD": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v2-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT700to1000_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 48506751, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212538/0000/"], 
"nevents": 48506751.0, 
"MCScaling": 7.929327755017028}, 
"selection": "dijet", 
"XS": 6430.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1000to1500_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 4773503, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_140640/0000/"], 
"nevents": 4773503.0, 
"MCScaling": 4.578849433005488}, 
"2016": 
{"miniAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1000to1500_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 4365993, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_140640/0000/"], 
"nevents": 4365993.0, 
"MCScaling": 4.316296744635184}, 
"2017": 
{"miniAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1000to1500_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 10256089, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_211852/0000/"], 
"nevents": 10256089.0, 
"MCScaling": 4.5330503196686385}, 
"2018": 
{"miniAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1000to1500_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 14527915, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_211852/0000/"], 
"nevents": 14527915.0, 
"MCScaling": 4.615618500452405}, 
"selection": "dijet", 
"XS": 1121.0, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1500to2000_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 3503675, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_141640/0000/"], 
"nevents": 3503675.0, 
"MCScaling": 0.6015750691773637}, 
"2016": 
{"miniAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1500to2000_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 3217830, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_141641/0000/"], 
"nevents": 3217830.0, 
"MCScaling": 0.5647435959326627}, 
"2017": 
{"miniAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1500to2000_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 7326541, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212437/0000/"], 
"nevents": 7326541.0, 
"MCScaling": 0.6119182670239612}, 
"2018": 
{"miniAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT1500to2000_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 10871473, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221106_212439/0000/"], 
"nevents": 10871473.0, 
"MCScaling": 0.5947916806397807}, 
"selection": "dijet", 
"XS": 108.1, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT2000toInf_nom_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 1629000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_142049/0000/"], 
"nevents": 1629000.0, 
"MCScaling": 0.2632037636279926}, 
"2016": 
{"miniAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT2000toInf_nom_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 1847781, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_142049/0000/"], 
"nevents": 1847781.0, 
"MCScaling": 0.20006146717062245}, 
"2017": 
{"miniAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 
"skimmer": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL17PFNanoV3-106X_mc2017_realistic_v9-v1-3a26c9a4394b9b50b8bf5dda6bb9a62c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT2000toInf_nom_UL2017.root", 
"nevents": 0.0, 
"nGenWeights": 4112573, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212639/0000/"], 
"nevents": 4112573.0, 
"MCScaling": 0.22175750042613224}, 
"2018": 
{"miniAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 
"skimmer": "/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/kadatta-RunIISummer20UL18PFNanoV3-106X_upgrade2018_realistic_v16_L1v1-v1-654ac7a30df3f7a4474385d73390277c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDHT2000toInf_nom_UL2018.root", 
"nevents": 0.0, 
"nGenWeights": 5374711, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212640/0000/"], 
"nevents": 5374711.0, 
"MCScaling": 0.24473590857815425}, 
"selection": "dijet", 
"XS": 21.99, 
"label": "QCD MG5+Pythia8", "color": 616}, "QCD_Pt_170to300_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt170to300_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 27953000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_143142/0000/"], 
"nevents": 27953000.0, 
"MCScaling": 72.26354802704539}, 
"2016": 
{"miniAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt170to300_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 25809000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_143143/0000/"], 
"nevents": 25809000.0, 
"MCScaling": 67.48034510442093}, 
"2017": 
{"miniAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt170to300_UL2017.root", 
"nevents": 27963000.0, 
"nGenWeights": 29223000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211849/0000/"], 
"nevents": 29223000.0, 
"MCScaling": 147.02854710330902}, 
"2018": 
{"miniAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt170to300_UL2018.root", 
"nevents": 28425000.0, 
"nGenWeights": 29628000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211856/0000/"], 
"nevents": 29628000.0, 
"MCScaling": 209.16306754421493}, 
"selection": "dijet", 
"XS": 103600.0, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_300to470_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt300to470_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 49674000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_144456/0000/"], 
"nevents": 49674008.02757648, 
"MCScaling": 2.6828558335346457}, 
"2016": 
{"miniAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt300to470_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 55264000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_144457/0000/"], 
"nevents": 55264009.34843346, 
"MCScaling": 2.0791456799182106}, 
"2017": 
{"miniAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt300to470_UL2017.root", 
"nevents": 53017000.0, 
"nGenWeights": 55438000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212537/0000/"], 
"nevents": 55438008.8400189, 
"MCScaling": 5.113257606695769}, 
"2018": 
{"miniAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt300to470_UL2018.root", 
"nevents": 55315000.0, 
"nGenWeights": 57670000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212538/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212538/0001/"], 
"nevents": 57670009.131807365, 
"MCScaling": 7.089508875411825}, 
"selection": "dijet", 
"XS": 6835.0, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_470to600_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt470to600_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 50563000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_144732/0000/"], 
"nevents": 50563169.25699422, 
"MCScaling": 0.21282150524098647}, 
"2016": 
{"miniAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt470to600_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 52408000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_144734/0000/"], 
"nevents": 52408168.58703068, 
"MCScaling": 0.1770318977007327}, 
"2017": 
{"miniAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt470to600_UL2017.root", 
"nevents": 46970000.0, 
"nGenWeights": 49709000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212639/0000/"], 
"nevents": 49709163.34203214, 
"MCScaling": 0.46046017678891143}, 
"2018": 
{"miniAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt470to600_UL2018.root", 
"nevents": 51191000.0, 
"nGenWeights": 52391000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212639/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212639/0001/"], 
"nevents": 52391170.42485976, 
"MCScaling": 0.6301316432841518}, 
"selection": "dijet", 
"XS": 551.9, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_600to800_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt600to800_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 52147000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_143336/0000/"], 
"nevents": 52147001.17667397, 
"MCScaling": 0.05851577526032178}, 
"2016": 
{"miniAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt600to800_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 64633000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_143335/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_143335/0001/"], 
"nevents": 64633001.48389972, 
"MCScaling": 0.040705092359939965}, 
"2017": 
{"miniAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt600to800_UL2017.root", 
"nevents": 63458000.0, 
"nGenWeights": 67043000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211952/0000/"], 
"nevents": 67043001.60200574, 
"MCScaling": 0.09681164745014394}, 
"2018": 
{"miniAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt600to800_UL2018.root", 
"nevents": 65300000.0, 
"nGenWeights": 67460000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211953/0000/", "/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_211953/0001/"], 
"nevents": 67460001.54427473, 
"MCScaling": 0.1387699999332938}, 
"selection": "dijet", 
"XS": 156.5, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_800to1000_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt800to1000_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 35063000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_145241/0000/"], 
"nevents": 35063000.0, 
"MCScaling": 0.01459159305250549}, 
"2016": 
{"miniAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt800to1000_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 37714000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_145241/0000/"], 
"nevents": 37714000.0, 
"MCScaling": 0.011696339108023543}, 
"2017": 
{"miniAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt800to1000_UL2017.root", 
"nevents": 35696000.0, 
"nGenWeights": 36890000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212838/0000/"], 
"nevents": 36890000.0, 
"MCScaling": 0.029499991021957172}, 
"2018": 
{"miniAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt800to1000_UL2018.root", 
"nevents": 36056000.0, 
"nGenWeights": 37016000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212838/0000/"], 
"nevents": 37016000.0, 
"MCScaling": 0.04240352171169224}, 
"selection": "dijet", 
"XS": 26.24, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1000to1400_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 19077000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_143800/0000/"], 
"nevents": 19077000.0, 
"MCScaling": 0.007638902446401427}, 
"2016": 
{"miniAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1000to1400_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 19813000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_143802/0000/"], 
"nevents": 19813000.0, 
"MCScaling": 0.006341493784989653}, 
"2017": 
{"miniAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1000to1400_UL2017.root", 
"nevents": 18653000.0, 
"nGenWeights": 18989000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212157/0000/"], 
"nevents": 18989000.0, 
"MCScaling": 0.016323666274158724}, 
"2018": 
{"miniAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1000to1400_UL2018.root", 
"nevents": 19106000.0, 
"nGenWeights": 19538000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212157/0000/"], 
"nevents": 19538000.0, 
"MCScaling": 0.02288234721660354}, 
"selection": "dijet", 
"XS": 7.474, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1400to1800_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 3945000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_145015/0000/"], 
"nevents": 3945000.0, 
"MCScaling": 0.0032036861903675536}, 
"2016": 
{"miniAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1400to1800_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 10722000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_145015/0000/"], 
"nevents": 10722000.0, 
"MCScaling": 0.0010163000360567057}, 
"2017": 
{"miniAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1400to1800_UL2017.root", 
"nevents": 10358000.0, 
"nGenWeights": 10994000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212739/0000/"], 
"nevents": 10994000.0, 
"MCScaling": 0.002445231615790431}, 
"2018": 
{"miniAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1400to1800_UL2018.root", 
"nevents": 10550000.0, 
"nGenWeights": 10886000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212739/0000/"], 
"nevents": 10886000.0, 
"MCScaling": 0.003561789749898953}, 
"selection": "dijet", 
"XS": 0.6482, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1800to2400_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 5262000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_144014/0000/"], 
"nevents": 5262000.0, 
"MCScaling": 0.0003240387290478905}, 
"2016": 
{"miniAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1800to2400_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 5236000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2_jetObsSkim_Summer20_nomTreesV3/221105_144014/0000/"], 
"nevents": 5236000.0, 
"MCScaling": 0.00028076883056722684}, 
"2017": 
{"miniAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1800to2400_UL2017.root", 
"nevents": 5191000.0, 
"nGenWeights": 5488000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212301/0000/"], 
"nevents": 5488000.0, 
"MCScaling": 0.0006608644941690962}, 
"2018": 
{"miniAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt1800to2400_UL2018.root", 
"nevents": 5152000.0, 
"nGenWeights": 5461000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212302/0000/"], 
"nevents": 5461000.0, 
"MCScaling": 0.0009578890756546419}, 
"selection": "dijet", 
"XS": 0.08745, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 
"skimmer": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/kadatta-RunIISummer20UL16PFNanov2pt3APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-a11d409e36cbafe26c98db696a6b5c66/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt2400to3200_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 2999000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3APVv2_jetObsSkim_Summer20_nomTreesV3/221105_143550/0000/"], 
"nevents": 2999000.0, 
"MCScaling": 3.404169075691897e-05}, 
"2016": 
{"miniAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt2400to3200_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 2836000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_143550/0000/"], 
"nevents": 2836000.0, 
"MCScaling": 3.103717096897038e-05}, 
"2017": 
{"miniAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt2400to3200_UL2017.root", 
"nevents": 2757000.0, 
"nGenWeights": 2997000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212054/0000/"], 
"nevents": 2997000.0, 
"MCScaling": 7.245687564898232e-05}, 
"2018": 
{"miniAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt2400to3200_UL2018.root", 
"nevents": 2901000.0, 
"nGenWeights": 2997000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212054/0000/"], 
"nevents": 2997000.0, 
"MCScaling": 0.0001045058201708375}, 
"selection": "dijet", 
"XS": 0.005236, 
"label": "QCD Pythia8", "color": 600}, "QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8": 
{"2016_preVFP"	: 
{"miniAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 
"skimmer": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt2APVv2-106X_mcRun2_asymptotic_preVFP_v11-v1-5b62f87422b4b3633317ac4b6822e17c/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt3200toInf_UL2016_preVFP.root", 
"nevents": 0.0, 
"nGenWeights": 1000000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt2APVv2_jetObsSkim_Summer20_nomTreesV3/221105_144230/0000/"], 
"nevents": 1000000.0, 
"MCScaling": 2.636116756e-06}, 
"2016": 
{"miniAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 
"skimmer": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/spigazzi-RunIISummer20UL16PFNanov2pt3-106X_mcRun2_asymptotic_v17-v1-36773dbd6a9d41d00afbee31f82bc069/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt3200toInf_UL2016.root", 
"nevents": 0.0, 
"nGenWeights": 996000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL16PFNanov2pt3_jetObsSkim_Summer20_nomTreesV3/221105_144231/0000/"], 
"nevents": 996000.0, 
"MCScaling": 2.2819497164658635e-06}, 
"2017": 
{"miniAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD-106X_mc2017_realistic_v9-v1-7a1edb72314467730e458def0bc98536/USER", 
"skimmer": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL17PFNanoAOD_jetObservables_Skimmer_v04p3-7a1edb72314467730e458def0bc98536/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt3200toInf_UL2017.root", 
"nevents": 1000000.0, 
"nGenWeights": 1000000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL17PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212435/0000/"], 
"nevents": 1000000.0, 
"MCScaling": 5.607165824000001e-06}, 
"2018": 
{"miniAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM", 
"nanoAOD": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD-106X_upgrade2018_realistic_v16_L1v1-v1-fff189d3e67d18da8f7301eb2c0e2940/USER", 
"skimmer": "/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/algomez-RunIISummer20UL18PFNanoAOD_jetObservables_Skimmer_v04p3-fff189d3e67d18da8f7301eb2c0e2940/USER", 

"localNano": "UL17and18_nano_tests/ROOT/jetObservables_nanoskim_QCDPt3200toInf_UL2018.root", 
"nevents": 934000.0, 
"nGenWeights": 1000000, 
"t3_dirs": ["/pnfs/psi.ch/cms/trivcat/store/user/kadatta/jetObservables/QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8/RunIISummer20UL18PFNanoAOD_jetObsSkim_Summer20_nomTreesV3/221106_212436/0000/"], 
"nevents": 1000000.0, 
"MCScaling": 8.0873134264e-06}, 
"selection": "dijet", 
"XS": 0.0001352, 
"label": "QCD Pythia8", "color": 600}}


def checkDict( string, dictio ):
    return next(v for k,v in dictio.items() if string.startswith(k))


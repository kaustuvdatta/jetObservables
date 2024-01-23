 #time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/QCD_MuEn470to600_UL18_1.root --local --selection '_WtopSel' --year '2018' --numEvents 10000000 --isBkgMC
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile  'root://cms-xrd-global.cern.ch///store/user/kadatta/PFNano/106x_v02/SingleMuon/Run2017C-UL2017_MiniAODv2-v1_PFNanov2pt3/221116_104120/0000/nano_data2017_1.root' --local --selection '_WtopSel' --year '2017' --numEvents 50000 --runEra 'C'

python jetObservables_nSubProducer_WtopSel_V3.py --iFile nano_mc2017_250.root  --local --selection '_WtopSel' --year '2017' --numEvents 1000000 --effMap
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/TTToSemiLeptonic_UL18_testPFNano.root --local --selection '_WSel' --year '2018' --numEvents 10000000 --isSigMC
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/TTToSemiLeptonic_UL18_testPFNano.root --local --selection '_WtopSel' --year '2018' --numEvents 10000000 --isSigMC
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile  nano_mc2017_250.root --local --selection '_WtopSel' --year '2017' --numEvents 250000 --isSigMC #nano_mc2017_250
#--onlyUnc '_jer' --effMap
#isSigMC 
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/TTToSemiLeptonic_UL18_testPFNano.root --local --selection '_WtopSel' --year '2018' --numEvents 1000000 --onlyUnc '_jer'
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/TTToSemiLeptonic_UL18_testPFNano.root --local --selection '_WtopSel' --year '2018' --numEvents 1000000 --onlyUnc '_jesCorr'
#time python jetObservables_nSubProducer_WtopSel_V3.py --iFile skimTest_nano/TTToSemiLeptonic_UL18_testPFNano.root --local --selection '_WtopSel' --year '2018' --numEvents 1000000 --onlyUnc '_jesUncorr'

#TTToSemiLeptonic_UL18_testPFNano.root
#python jetObservables_nSubProducer_V3.py --iFile ttbar_testfiles/ttsemileptonic_test_nano_mc2017.root --local --selection 'topSel' --year '2017'  --numEvents 1000000

#python jetObservables_nSubProducer_V3.py --iFile 'root://cms-xrd-global.cern.ch///store/user/kadatta/PFNano/SingleMuon/Run2017E-UL2017_MiniAODv2-v1_PFNanoAOD/210427_094611/0000/nano_data2017_787.root' --local --selection '_WtopSel' --year '2017' --numEvents 10000 --isSigMC



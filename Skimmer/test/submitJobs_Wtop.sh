echo 'Please input version for skimmer run'
read version
echo 'Chosen versioning:' $version

#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v nom$version --year '2016' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v nom$version --year '2016_preVFP' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v nom$version --year '2017' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v nom$version --year '2018' --selection '_WtopSel' --isSigMC

#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016' --selection '_WtopSel' --onlyUnc '_jer' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016_preVFP' --selection '_WtopSel' --onlyUnc '_jer' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2017' --selection '_WtopSel' --onlyUnc '_jer' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2018' --selection '_WtopSel' --onlyUnc '_jer' 

#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2018' --selection '_WtopSel' --onlyUnc '_jesCorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2017' --selection '_WtopSel' --onlyUnc '_jesCorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016' --selection '_WtopSel' --onlyUnc '_jesCorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016_preVFP' --selection '_WtopSel' --onlyUnc '_jesCorr' 

#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016' --selection '_WtopSel' --onlyUnc '_jesUncorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2016_preVFP' --selection '_WtopSel' --onlyUnc '_jesUncorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2017' --selection '_WtopSel' --onlyUnc '_jesUncorr' 
#python multicrab_nSubProducer_WtopSel_V2.py --datasets TTToSemiLeptonic_TuneCP5_13TeV -v sys$version --year '2018' --selection '_WtopSel' --onlyUnc '_jesUncorr' 

#python multicrab_nSubProducer_WtopSel_V2.py --datasets varTT -v nom$version --year '2016' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets varTT -v nom$version --year '2016_preVFP' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets varTT -v nom$version --year '2017' --selection '_WtopSel' --isSigMC
#python multicrab_nSubProducer_WtopSel_V2.py --datasets varTT -v nom$version --year '2018' --selection '_WtopSel' --isSigMC

#python multicrab_nSubProducer_WtopSel_V2.py --datasets TT -v nom$version --year '2016' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets TT -v nom$version --year '2016_preVFP' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets TT -v nom$version --year '2017' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets TT -v nom$version --year '2018' --selection '_WtopSel' 

python multicrab_nSubProducer_WtopSel_V2.py --datasets ST -v nom$version --year '2016' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets ST -v nom$version --year '2016_preVFP' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets ST -v nom$version --year '2017' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets ST -v nom$version --year '2018' --selection '_WtopSel' 

python multicrab_nSubProducer_WtopSel_V2.py --datasets W -v nom$version --year '2016' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets W -v nom$version --year '2016_preVFP' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets W -v nom$version --year '2017' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets W -v nom$version --year '2018' --selection '_WtopSel' 

python multicrab_nSubProducer_WtopSel_V2.py --datasets Z -v nom$version --year '2016' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets Z -v nom$version --year '2016_preVFP' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets Z -v nom$version --year '2017' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets Z -v nom$version --year '2018' --selection '_WtopSel' 

python multicrab_nSubProducer_WtopSel_V2.py --datasets QCD -v nom$version --year '2016' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets QCD -v nom$version --year '2016_preVFP' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets QCD -v nom$version --year '2017' --selection '_WtopSel' 
python multicrab_nSubProducer_WtopSel_V2.py --datasets QCD -v nom$version --year '2018' --selection '_WtopSel'


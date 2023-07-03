#!/usr/bin/bash
for year in "2016HIPM" "2016HIPM-ext1" "2016HIPM-ext2" "2016" "2016-ext1" "2016-ext2" "2017" "2018"
do
    # "Mu50" "TkMu50" 
    for k in "AK8PFJet60" "AK8PFJet80" "AK8PFJet140" "AK8PFJet200" "AK8PFJet260" "AK8PFJet320" "AK8PFJet400" "AK8PFJet450" "AK8PFJet500" "AK8PFJet550";
    do
        CERT=" "
        TRIG="$k"
        echo "$TRIG"
        echo "$year"
        puLatest=" "
        runStart=0
		runStop=0
	    if  [[ "$year" == "2016HIPM" ]];
	        then
	            CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
	            puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
	            runStart=272760
			    runStop=278761
		
		elif  [[ "$year" == "2016HIPM-ext1" ]];
	        then
	            CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
	            puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
	            runStart=278770
			    runStop=278770

		elif  [[ "$year" == "2016HIPM-ext2" ]];
	        then
	            CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
	            puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
	            runStart=278806
			    runStop=278807

		elif  [[ "$year" == "2016" ]]; 
	        then
			    CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
			    puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
			    runStart=278808
			    runStop=284044

		elif  [[ "$year" == "2016-ext1" ]]; 
	        then
			    CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
			    puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
			    runStart=278769
			    runStop=278769

		elif  [[ "$year" == "2016-ext2" ]]; 
	        then
			    CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
			    puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
			    runStart=278801
			    runStop=278805	

		elif [[ "$year" == "2017" ]] ; 
			then 
		        CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
			    puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
			    runStart=297047
			    runStop=306460

		elif [[ "$year" == "2018" ]]; 
			then
			    CERT="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
			    puLatest="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/UltraLegacy/pileup_latest.txt"
				runStart=315257
			    runStop=325175	        

	    else
	    	echo Year not available;
			fi
	
		echo "$CERT"              
		echo "$puLatest"
		echo "$TRIG"
        echo "$year"
	 	echo "Run start: $runStart, run end: $runStop"
		brilcalc lumi --byls --begin $runStart --end $runStop --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i "$CERT" --hltpath HLT_${TRIG}_v* -o output_${TRIG}_${year}.csv   
		
    
    done
done



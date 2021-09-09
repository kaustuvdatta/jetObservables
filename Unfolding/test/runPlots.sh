for iyear in '2017' '2018' 'all';
do
    python DrawHistogram.py -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

    python DrawHistogram.py -L -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -L -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

    python DrawHistogram.py -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

    python DrawHistogram.py -L -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -L -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext png --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s WSel
    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext pdf --outputFolder /eos/home-k/kadatta/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/20210908/ -s topSel

done

for iyear in '2017' '2018' 'all';
do
    #python DrawHistogram.py -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    #python DrawHistogram.py -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

    #python DrawHistogram.py -L -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    #python DrawHistogram.py -L -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

    #python DrawHistogram.py -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    #python DrawHistogram.py -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

    #python DrawHistogram.py -L -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    #python DrawHistogram.py -L -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext png --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s WSel
    python DrawHistogram.py -p resol -y ${iyear} -v v1 --ext pdf --outputFolder /home/kaustuv1993/cernbox/PhD_Projects/CMSSW_10_6_14/src/jetObservables/Unfolding/test/Results/ -s topSel

done

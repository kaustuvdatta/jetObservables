#for iyear in 2017 2018 all;
#do
#    python DrawHistogram.py -y ${iyear} -v v04 -L --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20211011/
#    python DrawHistogram.py -y ${iyear} -v v04 --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20211011/
#    python DrawHistogram.py -p resol -y ${iyear} -v v04 --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20211011/
#    python DrawHistogram.py -p tauComp -y ${iyear} -v v04 --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20211011/
#done

for iyear in 2017 2018;
do
    python Purity.py -v v04 -e pdf -y ${iyear} --main Ptbin --outputFolder ~/cernbox/JetObservables/Updates/20211011/
    python Purity.py -v v04 -e pdf -y ${iyear} --main HTbin --outputFolder ~/cernbox/JetObservables/Updates/20211011/
    python Purity.py -v v04 -e pdf -y ${iyear} --main herwig --outputFolder ~/cernbox/JetObservables/Updates/20211011/
done
for iyear in 2017 2018 all;
do
    python DrawHistogram.py -y ${iyear} -v v02 -L --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20210527/
    python DrawHistogram.py -y ${iyear} -v v02 --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20210527/
    python DrawHistogram.py -p resol -y ${iyear} -v v02 --ext pdf --outputFolder ~/cernbox/JetObservables/Updates/20210527/
done

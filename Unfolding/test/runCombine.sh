variableList="tau21 tau32 tau_0p5_1 tau_0p5_2 tau_0p5_3 tau_0p5_4 tau_1_1 tau_1_2 tau_1_3 tau_1_4 tau_2_1 tau_2_2 tau_2_3 tau_2_4"
jettype="Jet1 Jet2"

workingDir=$PWD
unfoldDir=$1
year=$2
unfold=$3
typeQCD=$4

for jet in $jettype; do
    for tau in $variableList; do
        echo "|========> Runnning maxLikelihood unfold for ${year}/${jet}_${tau}/${unfold}"
        cd ${workingDir}/${unfoldDir}/${year}/${jet}_${tau}/${unfold}/
        source runCombine${4}.sh
        cd ${workingDir}
    done
done

#!/bin/bash
read -p 'What dataset do you want to run on? ' dataset
read -p 'Which run to start on? ' runStart
read -p 'Which run to end on? ' runStop

echo 'Running following type of command over the given run range: dasgoclient -query=file run=273403 dataset=/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD >>' run${runStart}to${runStop}_MiniAOD.txt

for (( c=$runStart; c<=$runStop; c++))
do 
    echo 'Running the following command for the given run and dataset: dasgoclient -query=file run='${c}' dataset='${dataset}
    echo $c
    dasgoclient -query="file run=${c} dataset=${dataset}">>run${runStart}to${runStop}_MiniAOD.txt
done

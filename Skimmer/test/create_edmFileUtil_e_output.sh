#!/bin/bash

echo "Run number?"
read -r run 
dasgoclient -query="file run=${run} dataset=/JetHT/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD" > run${run}_MiniAODs.txt
input=run${run}_MiniAODs.txt
echo "Input file containing listed data files: ${input}"

while IFS= read -r line
do
  echo "$line"
  edmFileUtil -e "root://cms-xrd-global.cern.ch//${line}" > edmFU_out_run${run}.out 
done < "$input"

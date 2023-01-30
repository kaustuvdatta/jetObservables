#/bin/bash
# example usage:
# ./createInvalidationShellScript.sh -f QCD_Summer19UL_datasetsToInvalidate.txt -o QCD_DAS_Invalidation_Summer19.sh

#$usage() { echo "Usage: $0 [-s <45|90>] [-p <string>]" 1>&2; exit 1; }

FILES="QCD_Summer19UL_datasetsToInvalidate.txt"
oFile="QCD_DAS_Invalidation_Summer19.sh"

while getopts ":f:o:" option; do

    case $option in
	       f) FILES=$OPTARG;;
	       o) oFile=$OPTARG;;
    esac
done

echo "Writing dasgoclient invalidation output lines to '${oFile}'  "
echo "#!/bin/bash">>${oFile}

while read -r dataset; do

    name="$dataset"
    echo "Name of dataset being configured for deletion - $name"
    echo -e "python \$DBS3_CLIENT_ROOT/examples/DataOpsScripts/DBS3SetDatasetStatus.py --dataset='${name}' --url=https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter --status=INVALID --recursive=False">>${oFile}

done<"$FILES" 

echo "Making '${oFile}' executable "
chmod u+x ${oFile}

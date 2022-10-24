#!/bin/bash

cmssw_ver=CMSSW_10_6_14
cmssw_base=$HOME
cmssw_location=$cmssw_base/$cmssw_ver
cmssw_src=$cmssw_location/src
logdir=$cmssw_src/logs
report=$logdir/report.log

datasets=$(realpath $1)

if [ -z $datasets ] || [ ! -f $datasets ]; then
  echo "No such file: $datasets";
  exit 1;
fi

nof_datasets=$(cat $datasets | wc -l);
proxy_open=$(( $nof_datasets / 6 + 1 )); # 10 min per dataset

if [ $(which dasgoclient 2>/dev/null | wc -l) = "0" ]; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh;
fi

if [ ! -d $cmssw_location ]; then
  cd $cmssw_base;
  scramv1 project CMSSW $cmssw_ver; # cmsrel $cmssw_ver
fi

cd $cmssw_src;
eval $(scramv1 runtime -sh); # cmsenv

if [ $(which voms-proxy-init 2>/dev/null | wc -l) = "0" ]; then
  source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh;
fi

voms-proxy-init -voms cms --hours $proxy_open;

rm -f lhe_printer.py;
wget -q https://raw.githubusercontent.com/HEP-KBFI/tth-nanoAOD/master/scripts/lhe_printer.py;
chmod +x lhe_printer.py;

rm -f $report;
for sample_name in `cat $datasets`; do
  if [[ "$sample_name" =~ .*/NANOAODSIM$ ]]; then
    miniaod_name=$(dasgoclient -query="parent dataset=$sample_name" | sort | uniq | grep "^/");
  elif [[ "$sample_name" =~ .*/MINIAODSIM$ ]]; then
    miniaod_name=$sample_name;
  else
    continue;
  fi

  dirname=$logdir/$(echo $sample_name | sed 's/^\///g; s/\//__/g');
  if [ ! -d $dirname ]; then
    mkdir -p $dirname;
  fi

  miniaod_files=$(dasgoclient -query="file dataset=$miniaod_name");
  for miniaod_file in $miniaod_files; do
    ./lhe_printer.py -v -d -o $dirname -i $miniaod_file &> $dirname/out.log;
    exit_code=$?;
    if [ "$exit_code" = 0 ]; then
      break;
    else
      rm -f $dirname/*.txt;
    fi
  done

  lhaid=$(grep "Beam 1" $dirname/out.log  | awk '{print $12}');
  if [ -z $lhaid ]; then
    lhaid='N/A';
  fi
  echo -e "$sample_name $lhaid" | tee -a $report;
done

voms-proxy-destroy

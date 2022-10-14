#!/bin/bash

echo "Start"

python crabMonitor.py crab_PhaseII_CMSSW_12_4_9_20220925 resubmit >monitor.log

index=1
while [ ${index} -le 60 ]; do

  sleep 3600;
  echo "`date`";

  python crabMonitor.py crab_PhaseII_CMSSW_12_4_9_20220925 resubmit >monitor.log

  index=$((index+1))
done;

echo "Finish"
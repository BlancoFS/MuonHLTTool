#!/bin/bash

########################
# Phase2HLT-EffFR.sh #
########################
SRC=$1
VER=$2
TAG=$3
DATASET=$4
JobID=$5
doDimuon=$6
doGenMatchForIso=$7

ARGS=("$@")
END=$#
AllArgs=""
for ((i=0;i<END;i++)); do
    AllArgs+="  "
    AllArgs+=${ARGS[$i]}
done

DATE=`date +%Y%m%d`

Name=${VER}_${TAG}_${JobID}

logdir=$SRC/v$DATE
mkdir -p $logdir

subfile=Phase2HLT-EffFR-v$DATE.sub
outfile=Phase2HLT-EffFR_$Name.\$\(Cluster\).\$\(Process\)
jobid=\$\(Cluster\).\$\(Process\)

cat > $subfile << EOF

universe                = vanilla
executable              = ./run_Phase2HLT-EffFR.sh
getenv                  = True
# request_memory          = 4 GB

# should_transfer_files   = YES
# when_to_transfer_output = ON_EXIT_OR_EVICT

initialdir              = $logdir
output                  = $outfile.out
error                   = $outfile.err
log                     = $outfile.log

environment             = CONDOR_ID=Phase2HLT-EffFR_${Name}_$jobid
Arguments               = $AllArgs
Queue
EOF

condor_submit $subfile -batch-name ${TAG}-EffFR
# rm $subfile

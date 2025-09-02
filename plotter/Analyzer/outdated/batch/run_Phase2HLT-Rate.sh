#!/bin/bash

SRC=$1
VER=$2
TAG=$3
DATASET=$4
JobID=$5
doDimuon=$6
doGenMatchForIso=$7

echo "Arguments:"
echo '        '${SRC}
echo '        '${VER}
echo '        '${TAG}
echo '        '${DATASET}
echo '        '${JobID}
echo '        '${doDimuon}
echo '        '${doGenMatchForIso}

# -- set output directory
echo " "

# source /cvmfs/sft.cern.ch/lcg/releases/LCG_95/ROOT/6.16.00/x86_64-centos7-gcc7-opt/ROOT-env.sh
which root

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64/dcap"

export OUTDIR=${SRC}/Outputs_${VER}/Rate
export WORKDIR=/d0/scratch/msoh/Phase2HLT-Rate_${VER}/${TAG}/${JobID}

echo " "
echo "Output directory is : " $OUTDIR
echo " "

mkdir -p ${WORKDIR}

cp ${SRC}/HLTRateAnalyzer.C ${WORKDIR}
# cp ${SRC}/HLTRateAnalyzer_C.* ${WORKDIR}
cp ${SRC}/MuonHLTNtuple.h ${WORKDIR}

cd ${WORKDIR}
ls


# -- run analyzer
# echo root -l -b -q 'HLTRateAnalyzer.C+("'$VER'", "'$TAG'", '$DATASET', "'$JobID'", '$doDimuon', '$doGenMatchForIso')'

echo 'gROOT->LoadMacro("HLTRateAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
root -l -b -q 'HLTRateAnalyzer.C+("'$VER'", "'$TAG'", '$DATASET', "'$JobID'", '$doDimuon', '$doGenMatchForIso')' >tmp.log

ls


# -- Copy output files
if [ ! -d "$OUTDIR" ]; then
mkdir -p ${OUTDIR}
fi
if [ ! -d "$OUTDIR/$TAG" ]; then
mkdir -p ${OUTDIR}/$TAG
fi

echo ""
echo ""
echo "Finished -- ls dir: "
pwd
ls -lrt | tail -n 5
echo "Copying to " ${OUTDIR}/$TAG
cp ${WORKDIR}/*.root ${OUTDIR}/$TAG/

echo ""
echo ""
echo "what is in " ${OUTDIR}/$TAG
ls -lrt ${OUTDIR}/$TAG | tail -n 5

end=`date`
echo ""
echo "Job end `date`"
echo ""

echo "==> " $start $end `uname -n` $processor $rate

# rm -r ${WORKDIR}
exit ${status}

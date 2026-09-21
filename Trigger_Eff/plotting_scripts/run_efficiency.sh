#!/bin/bash

#inputfile for dataset; chamber for ME chamber
#submit_efficiency.py should fill these boys out

INPUT_FILE=$1
CHAMBER=$2
local_dir="output_tmp"

#my hardcoded output path; can change later if not want hardcode
OUTPUT_DIR="root://cmseos.fnal.gov//store/group/lpcmds/ctidmore/all2024"

echo "=========================================="
echo "Job started:  $(date)"
echo "Input file:   $INPUT_FILE"
echo "Chamber:      $CHAMBER"
echo "Output dir:   $OUTPUT_DIR"
echo "=========================================="


source /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc11-opt/setup.sh

#runs script to make roots/efficiencies
mkdir -p "$local_dir"
python3 run_efficiency.py "$INPUT_FILE" "$CHAMBER" "$local_dir"

#takes roots and stores them in my eos directory
DATASET_NAME=$(basename "$INPUT_FILE" .root)
OUTPUT_FILENAME="${DATASET_NAME}_${CHAMBER}.root"

xrdcp -f "$local_dir/${OUTPUT_FILENAME}" "${OUTPUT_DIR}/${OUTPUT_FILENAME}"           

if [ $? -eq 0 ]; then
    echo "Transfer succeeded: ${OUTPUT_DIR}/${OUTPUT_FILENAME}"
else
    echo "ERROR: xrdcp transfer failed!"
    exit 1
fi


echo "Total time ran: $(printf '%02d:%02d:%02d' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)))"
echo "Job finished: $(date)"
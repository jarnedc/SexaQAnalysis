#!/bin/bash

CMSSW_DIR=/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/

source $VO_CMS_SW_DIR/cmsset_default.sh
# shopt -s expand_aliases is needed if you want to use the alias 'cmsenv' created by $VO_CMS_SW_DIR/cmsset_default.sh instead of the less mnemonic eval `scramv1 runtime -sh`
shopt -s expand_aliases


pwd_name=$(echo $PWD)
echo "${pwd_name}"
mkdir Results/"$1"
cd Results/"$1"
mkdir "$2"
cd "$2"
mkdir analyzed_root_files
mkdir error_txt_files
mkdir output_txt_files
cd ../../..

cd $CMSSW_DIR/src
eval `scram runtime -sh`                                         # don't use cmsenv, won't work on batch                                                                                                                                            
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi


for D in `find /pnfs/iihe/cms/store/user/jdeclerc/$1/$2 -type d`
do
   # echo "---------------------------------"
    #echo WILL RUN ON THE FOLLOWING DIRECTORY
   # echo $D

    for filename in $D/*root; do
    #    echo "---------"
	if [[ $filename = *"events_skimmed_2016_trialA"* ]]; then
		
		PART1=$(echo "$filename" | cut -d/ -f8)
		PART2=$(echo "$filename" | cut -d/ -f9)
		PART3=$(echo "$filename" | cut -d/ -f10)
		PART4=$(echo "$filename" | cut -d/ -f11)
		PARTS="${PART2}_${PART3}_${PART4}"

		echo ${PART1} #eg Tau
		echo ${PART2} #eg Tau_Run2016C-07Aug17-v1
		echo ${PART3} #eg 180909_172549
		echo ${PART4} #eg 0000

		bare_filename1=$(basename $filename) 
		bare_filename2="${bare_filename1%%.*}"
		echo $filename
		qsub ${pwd_name}/analyzer_allPDs_2016_step2.sh -v INPUT_ROOTFILE="file://$filename",OUTPUT_ROOTFILE="${pwd_name}/Results/$1/$2/analyzed_root_files/analyzed_${bare_filename2}_${PARTS}.root",OUT_TXT="${pwd_name}/Results/$1/$2/output_txt_files/out_${bare_filename2}_${PARTS}.txt",ERROR_TXT="${pwd_name}/Results/$1/$2/error_txt_files/error_${bare_filename2}_${PARTS}.txt" 

	fi
    done

done

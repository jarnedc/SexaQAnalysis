#!/bin/bash

CMSSW_DIR=/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/
CMSSW_CONFIG_FILE=$CMSSW_DIR/src/V0_angular_correlation_analysis/test/analyzer_AOD_cfg_test.py

source $VO_CMS_SW_DIR/cmsset_default.sh
# shopt -s expand_aliases is needed if you want to use the alias 'cmsenv' created by $VO_CMS_SW_DIR/cmsset_default.sh instead of the less mnemonic eval `scramv1 runtime -sh`
shopt -s expand_aliases


pwd_name=$(echo $PWD)
echo "${pwd_name}"

cd $CMSSW_DIR/src
eval `scram runtime -sh`                                         # don't use cmsenv, won't work on batch                                                                                                                                            
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi

#pwd = $PWD
#echo $pwd
#echo $pwd | awk -F/ '{print 2}'
#echo ./another/example/path | cut -d/ -f3
#echo $PWD | cut -d/ -f3

for D in `find /pnfs/iihe/cms/store/user/jdeclerc/MuonEG/allPDs_2016_multicrab -type d`
do
   # echo "---------------------------------"
    #echo WILL RUN ON THE FOLLOWING DIRECTORY
   # echo $D

    for filename in $D/*root; do
    #    echo "---------"
	if [[ $filename = *"events_skimmed"* ]]; then
		
		PART1=$(echo "$filename" | cut -d/ -f8)
		PART2=$(echo "$filename" | cut -d/ -f10)
		PART3=$(echo "$filename" | cut -d/ -f11)
		PARTS="${PART1}_${PART2}_${PART3}"

		bare_filename1=$(basename $filename) 
		bare_filename2="${bare_filename1%%.*}"
		qsub ${pwd_name}/analyzer_allPDs_2016_step2.sh -v INPUT_ROOTFILE="file://$filename",OUTPUT_ROOTFILE="${pwd_name}/Results/analyzed_${bare_filename2}_${PARTS}.root",OUT_TXT="${pwd_name}/Results/out_${bare_filename2}_${PARTS}.txt",ERROR_TXT="${pwd_name}/Results/error_${bare_filename2}_${PARTS}.txt" 

	fi
    done

done

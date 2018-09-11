#!/bin/bash

CMSSW_DIR=/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/
CMSSW_CONFIG_FILE=$CMSSW_DIR/src/V0_angular_correlation_analysis/test/analyzer_AOD_cfg_test.py

source $VO_CMS_SW_DIR/cmsset_default.sh
# shopt -s expand_aliases is needed if you want to use the alias 'cmsenv' created by $VO_CMS_SW_DIR/cmsset_default.sh instead of the less mnemonic eval `scramv1 runtime -sh`
shopt -s expand_aliases

cd $CMSSW_DIR/src
eval `scram runtime -sh`                                         # don't use cmsenv, won't work on batch                                                                                                                                            
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi


for D in `find /pnfs/iihe/cms/store/user/jdeclerc/MuonEG/allPDs_2016_multicrab -type d`
do
    echo "---------------------------------"
    echo WILL RUN ON THE FOLLOWING DIRECTORY
    echo $D

    for filename in $D/*root; do
        echo "---------"
	if [[ $filename = *"events_skimmed"* ]]; then
		echo will run on the following file: "file://$filename"
		#cmsRun $CMSSW_CONFIG_FILE "file://$D/$filename"  > myout.txt 2>myerr.txt
	fi
    done

done

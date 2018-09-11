###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE

CMSSW_DIR=/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30
CMSSW_CONFIG_FILE=$CMSSW_DIR/src/SexaQAnalysis/V0_angular_correlation_analysis/test/analyzer_AOD_cfg_test.py

source $VO_CMS_SW_DIR/cmsset_default.sh
# shopt -s expand_aliases is needed if you want to use the alias 'cmsenv' created by $VO_CMS_SW_DIR/cmsset_default.sh instead of the less mnemonic eval `scramv1 runtime -sh`
shopt -s expand_aliases 

cd $CMSSW_DIR/src
eval `scramv1 runtime -sh`
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi

#cmsRun $CMSSW_CONFIG_FILE file:///pnfs/iihe/cms/store/user/jdeclerc/MuonEG/allPDs_2016_multicrab/180907_162222/0000/events_skimmed_9.root test_output  > out.txt  2>error.txt
#echo $CMSSW_CONFIG_FILE
#echo ${ROOTFILE}
echo cmsRun
echo $CMSSW_CONFIG_FILE
echo ${ROOTFILE}
#cmsRun /user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/V0_angular_correlation_analysis/test/analyzer_AOD_cfg_test.py file:///pnfs/iihe/cms/store/user/jdeclerc/MuonEG/allPDs_2016_multicrab/180907_162222/0000/events_skimmed_9.root test_output
cmsRun $CMSSW_CONFIG_FILE ${INPUT_ROOTFILE} ${OUTPUT_ROOTFILE}  > ${OUT_TXT}  2>${ERROR_TXT}



#cmsRun /user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/V0_angular_correlation_analysis/test/analyzer_AOD_cfg_test.py file:///pnfs/iihe/cms/store/user/jdeclerc/MuonEG/allPDs_2016_multicrab/180907_162222/0000/events_skimmed_9.root test_output > out.txt  2>error.txt

#echo cmsRun $CMSSW_CONFIG_FILE   > out.txt  2>error.txt
#echo $CMSSW_CONFIG_FILE
echo YOLO
#cmsRun $CMSSW_CONFIG_FILE   > out.txt  2>error.txt
#cmsRun $CMSSW_CONFIG_FILE > myout.txt 2>myerr.txt

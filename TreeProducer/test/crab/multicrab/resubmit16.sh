#!/bin/bash


pwd=$PWD
source $VO_CMS_SW_DIR/cmsset_default.sh                          # make scram available                                                                                                                                                             
cd /user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/                         # your local CMSSW release                                                                                                                                                         
eval `scram runtime -sh`                                         # don't use cmsenv, won't work on batch                                                                                                                                            
cd $pwd

for D in `find ./crab_projects16_trialD -type d` -maxdepth 0
do
    echo WILL RUN ON THE FOLLOWING DIRECTORY
    echo $D
    crab resubmit -d $D
    #crab status -d  $D --verboseErrors
    #crab kill -d $D
done

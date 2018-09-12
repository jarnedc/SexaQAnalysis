#!/bin/bash


pwd=$PWD
source $VO_CMS_SW_DIR/cmsset_default.sh                          # make scram available                                                                                                                                                             
cd /user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/                          # your local CMSSW release                                                                                                                                                         
eval `scram runtime -sh`                                         # don't use cmsenv, won't work on batch                                                                                                                                            
cd $pwd

for D in `find ./crab_projects17 -type d` -maxdepth 0
do
    echo WILL RUN ON THE FOLLOWING DIRECTORY
    echo $D
    crab status -d  $D --verboseErrors
    echo ------------------------------------------------------
done

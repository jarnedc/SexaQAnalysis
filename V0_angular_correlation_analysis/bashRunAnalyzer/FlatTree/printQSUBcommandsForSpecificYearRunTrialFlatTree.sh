#!/bin/bash

#this script will for a certain year, certain run and a certain trial of the running of the treeproducer get all the preFilterInfo_events_skimmed_ files which belong to this year, run and trial and print the output to a .sh file which you can run with bash
#what it does: it hadd all the preFilterInfo_events_skimmed_ files from a specific year, run and trial

pwd_name=$(echo $PWD)
#echo "${pwd_name}"
mkdir Results_FlatTree/"$1"
cd Results_FlatTree/"$1"
mkdir "$2"
cd "$2"
mkdir "$3"
cd "$3"
cd ../../../..

if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi


for D in `find /pnfs/iihe/cms/store/user/jdeclerc/$2/ -type d`
do
   # echo "---------------------------------"
    #echo WILL RUN ON THE FOLLOWING DIRECTORY
   # echo $D
    for filename in $D/*root; do
    #    echo "---------"
	if  [[ $filename = *"preFilterInfo_events_skimmed_"$1"_trial"$3""* ]] && [[ $filename != *"failed"* ]] ; then
		PART1=$(echo "$filename" | cut -d/ -f8)
		PART2=$(echo "$filename" | cut -d/ -f9)
		PART3=$(echo "$filename" | cut -d/ -f10)
		PART4=$(echo "$filename" | cut -d/ -f11)
		PARTS="${PART2}_${PART3}_${PART4}"


		bare_filename1=$(basename $filename) 
		bare_filename2="${bare_filename1%%.*}"

		echo hadd "/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/FlatTree/Results_FlatTree/$1/$2/$3/combined_FlatTree_"${PART2}".root" $D/*"preFilterInfo_events_skimmed_"$1"_trial"$3""* 
		break #once you found a file in a certain directory that is ok, you can go to a next directory now (one directory is a specific year, run and trial)
	fi
    done

done

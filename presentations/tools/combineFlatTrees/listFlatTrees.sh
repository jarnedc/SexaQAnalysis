#!/bin/bash
for D in `find /user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/FlatTree/Results_FlatTree/ -type d`
do
   # echo "---------------------------------"
   # echo WILL RUN ON THE FOLLOWING DIRECTORY
   # echo $D

    for filename in $D/*root; do
    #    echo "---------"
        if [[ $filename = *"combined_FlatTree"* ]] && [[ $filename != *"preFilterInfo"* ]]; 
	then

                PART=$(echo "$filename" | cut -d/ -f14)

		if [[ "$PREV_PART" != "$PART" ]]; 
		then
			echo --NEW_RUN_AND_OR_YEAR--	
		fi

	  	PREV_PART="$PART"

		echo $filename
               # bare_filename1=$(basename $filename)
               # bare_filename2="${bare_filename1%%.*}"
                #echo $filename
               # echo qsub ${pwd_name}/analyzer_allPDs_2016_step2.sh -v INPUT_ROOTFILE="file://$filename",OUTPUT_ROOTFILE="${pwd_name}/Results/$1/$2/$3/analyzed_root_files/analyzed_${bare_filename2}_${PARTS}.root",OUT_TXT="${pwd_name}/Results/$1/$2/$3/output_txt_files/out_${bare_filename2}_${PARTS}.txt",ERROR_TXT="${pwd_name}/Results/$1/$2/$3/error_txt_files/error_${bare_filename2}_${PARTS}.txt" 

        fi
    done

done
echo --NEW_RUN_AND_OR_YEAR--	

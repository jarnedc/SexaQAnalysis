#!/bin/bash

declare -a PDs_array=(
"DisplacedJet"
"SingleElectron"
"MuOnia"
"SinglePhoton"
"BTagMu"
"Tau"
"HTMHT"
"MuonEG"
"Charmonium"
"DoubleMuonLowMass"
"DoubleMuon"
"JetHT"
"SingleMuon"
"MET"
"DoubleEG"
"BTagCSV"
)

for i in "${PDs_array[@]}"
do
  # echo "$i"
   #mkdir Results/"$i"
  # echo "$i"

 #  echo ./printQSUBcommandsForSpecificYearRunTrial.sh $1 "$i" $2 
   echo "FlatTree_to_hadd_$1_"$i"_$2.txt" 
   ./printQSUBcommandsForSpecificYearRunTrialFlatTree.sh $1 "$i" $2 > "FlatTree_to_hadd_$1_"$i"_$2.sh" 
   #./analyzer_allPDs_2016_step1.sh "$i"
   # or do whatever with individual element of the array
done

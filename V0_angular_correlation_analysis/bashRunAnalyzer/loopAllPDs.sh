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
   echo "qsub_to_analyze_$1_"$i"_$2.txt" 
   ./printQSUBcommandsForSpecificYearRunTrial.sh $1 "$i" $2 > "qsub_to_analyze_$1_"$i"_$2.txt" 
done

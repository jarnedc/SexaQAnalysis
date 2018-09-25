#!/bin/bash
#usage example:
# ./haddAll.sh Year Trial 
#eg
# ./haddAll.sh 2016 A
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
   cd Results/$1/"$i"/$2/analyzed_root_files/
   pwd
   rm combined.root
   cd ../../../../..
done

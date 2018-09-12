#!/bin/bash

declare -a PDs_array=(
#"DisplacedJet"
"SingleElectron"
#"MuOnia"
#"SinglePhoton"
#"BTagMu"
#"Tau"
#"HTMHT"
#"MuonEG"
#"Charmonium"
#"DoubleMuonLowMass"
#"DoubleMuon"
#"JetHT"
#"SingleMuon"
#"MET"
#"DoubleEG"
#"BTagCSV"
)

for i in "${PDs_array[@]}"
do
  # echo "$i"
   #mkdir Results/"$i"
   ./analyzer_allPDs_2016_step1.sh "$i"
   # or do whatever with individual element of the array
done

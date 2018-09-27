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
   qsub hadd.sh -v INPUT_PATH="/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/$1/$i/$2/analyzed_root_files/",OUTPUTFILENAME="/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined/combined_"$1"_"$i"_"$2".root"
done


declare -a arr=(
"/BTagCSV/Run2016*/AOD"
"/BTagCSV/Run2017*/AOD"
"/BTagMu/Run2016*/AOD"
"/BTagMu/Run2017*/AOD"
"/Charmonium/Run2016*/AOD"
"/Charmonium/Run2017*/AOD"
"/DisplacedJet/Run2016*/AOD"
"/DisplacedJet/Run2017*/AOD"
"/DoubleEG/Run2016*/AOD"
"/DoubleEG/Run2017*/AOD"
"/DoubleMuon/Run2016*/AOD"
"/DoubleMuon/Run2017*/AOD"
"/DoubleMuonLowMass/Run2016*/AOD"
"/DoubleMuonLowMass/Run2017*/AOD"
"/HTMHT/Run2016*/AOD"
"/HTMHT/Run2017*/AOD"
"/JetHT/Run2016*/AOD"
"/JetHT/Run2017*/AOD"
"/MET/Run2016*/AOD"
"/MET/Run2017*/AOD"
"/MuOnia/Run2016*/AOD"
"/MuOnia/Run2017*/AOD"
"/MuonEG/Run2016*/AOD"
"/MuonEG/Run2017*/AOD"
"/SingleElectron/Run2016*/AOD"
"/SingleElectron/Run2017*/AOD"
"/SingleMuon/Run2016*/AOD"
"/SingleMuon/Run2017*/AOD"
"/SinglePhoton/Run2016*/AOD"
"/SinglePhoton/Run2017*/AOD"
"/Tau/Run2016*/AOD"
"/Tau/Run2017*/AOD")




for i in "${arr[@]}"
do
#   echo "$i"
   dasgoclient -query="$i" -json | cut -d ',' -f 19

   echo "------------------------------------"
   # or do whatever with individual element of the array
done

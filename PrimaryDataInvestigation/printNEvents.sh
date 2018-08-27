#!/bin/bash

dasgoclient -query="/SingleMuon/Run2016B-18Apr2017_ver2-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016C-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016D-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016E-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016F-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016G-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2016H-07Aug17-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017B-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017C-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017D-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017E-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017F-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017G-17Nov2017-v1/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8
dasgoclient -query="/SingleMuon/Run2017H-17Nov2017-v2/AOD" -json | grep "nblocks" | cut -d ',' -f 6,8



declare -a arr=(

)




for i in "${arr[@]}"
do
#   echo "$i"
   dasgoclient -query="$i" -json | grep "nblocks" | cut -d ',' -f 6,8

done

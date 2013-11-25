#!/bin/bash

mkdir -p ${2}
pwd=$PWD
#$ROOTSYS/bin/hadd ${2}/${1}_histo.root ${1}/job*/Filter.Data1.Reco.root
touch ${2}/${1}_total.txt
touch ${2}/${1}_lumi.csv
max=`ls ${1} | grep job | wc -l`
echo "There is ${max} subjobs in task ${1}"
for block in  `seq -s ' ' 1 ${max}`; do
    cd ${1}/job${block}
#    more total.txt
    tail total.txt -n 1 >> ${pwd}/${2}/${1}_total.txt
    cat lumi.csv >> ${pwd}/${2}/${1}_lumi.csv
    cd ${pwd}
done
gsum=0
while read -r line
do
  add=`echo $line`
#  echo $add
  gsum=`expr $gsum + $add`
done < <(more ${2}/${1}_total.txt)
echo $gsum > ${2}/${1}_sumEvents.txt
echo "Total events processed in filtering "${1}" is "${gsum} 
echo "converting csv to json"
csv2json.py ${2}/${1}_lumi.csv --output ${2}/lumi.json
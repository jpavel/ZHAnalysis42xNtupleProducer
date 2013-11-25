#!/bin/sh

time=`date "+%Y%m%d%H%M"`
task=${1}
output=${2}
outputDir=${3}
outputLoc=${4}
inputLoc=${5}

max=`ls ${task} |grep job | wc -l`
echo "There are ${max} jobs in the task ${task}"
source ShowResults_filter_error2.sh ${task} 1> results.out 2> results.err
echo "These jobs need to be re-submitted:"
source ShowResults_filter_good2.sh ${task} 2> resubmit.err | tee resubmit
resubtotal=`more resubmit | wc -l`
successful=`expr ${max} - ${resubtotal}`
echo "Total ${resubtotal}"
mkdir -p ${outputDir}/${output}
mkdir -p ${outputDir}/OLD/${output}_old_${time}
touch ${outputDir}/${output}/file
mv ${outputDir}/${output}/* ${outputDir}/OLD/${output}_old_${time}/
source MergeOutput_filter2.sh ${task} ${outputDir}/${output}/ 
finished=`more ${outputDir}/${output}/*_total.txt | wc -l`
echo "${finished} jobs produced output"
saved=`ls ${outputLoc} | wc -l`
saved_size=`ls -lh ${outputLoc} | head -n 1 | tr -cd [:digit:]`
input_size=`ls -lh ${inputLoc} | head -n 1 | tr -cd [:digit:]`
source removeDuplicate.sh ${outputLoc} 1 ${max}
ratio=`echo "${saved_size} / ${input_size} * 100" | bc -l | awk '{printf("%04.1f",$1);}'`
echo "The size of filtered output is ${saved_size} GB (${ratio} %)"
echo "---------SUMMARY----------"
echo "total: ${max}"
echo "successful: ${successful}"
echo "finished: ${finished}"
echo "saved: ${saved}"

echo "Calculating lumi..."
lumiCalc2.py overview -i ${outputDir}/${output}/lumi.json 1>lumi.out 2>lumi.err
tail lumi.out -n 4
touch Resubmit_now.sh
rm -f Resubmit_now.sh 
touch Resubmit_now_noDelete.sh
rm -f Resubmit_now_noDelete.sh
echo "source Resubmit2.sh resubmit ${task} ${outputLoc} 0" > Resubmit_now.sh
echo "source Resubmit2.sh resubmit ${task} ${outputLoc} 1" > Resubmit_now_noDelete.sh
chmod +x Resubmit_now.sh
chmod +x Resubmit_now_noDelete.sh
echo "To resubmit failed jobs do ./Resubmit_now.sh"
echo "To resubmit failed jobs and do not delete do ./Resubmit_now_noDelete.sh"

#!/bin/bash
    
name=${1}
max=`ls ${name} | grep job | wc -l`
echo "----${name}--------------------"
for block in  `seq -s ' ' 1 ${max}`; do
    echo "job ${block}:"
    more ${name}/job${block}/*stdout* | grep -i '(0)'
    more ${name}/job${block}/*stdout* | grep -i 'ERROR'
    more ${name}/job${block}/*stdout* | grep -i 'FATAL'
    more ${name}_job${block}/*stderr* | grep -i 'error'
    more ${name}_job${block}/*stderr* | grep -i 'segmentation'
    more ${name}_job${block}/*stderr* | grep -i 'fatal'
    ls -ltrh ${name}/job${block} | grep root | grep Filter
    #cd ${name}_job${block}
    more ${name}/job${block}/total.txt
    #    cd ..
done
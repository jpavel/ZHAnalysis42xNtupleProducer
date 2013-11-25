#!/bin/bash
    
name=${1}
max=`ls ${name} | grep job | wc -l`
#echo "----${name}--------------------"
for block in  `seq -s ' ' 1 ${max}`; do
#    echo "job ${block}:"
#    more ${name}_job${block}/*stdout* | grep -i '(0)'
#    more ${name}_job${block}/*stdout* | grep -i 'ERROR'
#    more ${name}_job${block}/*stdout* | grep -i 'FATAL'
    #more ${name}_job${block}/*stderr*
    pass=`ls ${name}/job${block} | grep total | wc -l`
    pass1=`more ${name}/job${block}/*stdout* | grep -i 'ERROR' | wc -l`
    pass2=`more ${name}/job${block}/*stdout* | grep -i 'FATAL' | wc -l`
    pass3=`more ${name}/job${block}/*stderr* | grep -i 'error' | wc -l`
    pass4=`more ${name}/job${block}/*stderr* | grep -i 'segmentation' | wc -l`
    pass1=`expr $pass1 + $pass2`
    pass1=`expr $pass1 + $pass3`
    pass1=`expr $pass1 + $pass4`
#    echo $pass1
#    cd ${name}_job${block}
    if [ $pass -eq 0 ] || [ $pass1 -gt 0 ]; then
      echo "${block}"
    fi
#    cd ..
done
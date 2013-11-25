#!/bin/bash

mydir=$1
filename=output_Ntuples_
before=`ls -1 ${mydir} | wc -l`

echo "Checking duplicates in directory "${mydir}". It has "${before}" files"
 
touch toRemove.list
rm -f toRemove.list
touch toRemove.list

touch toRemove.sh
rm -f toRemove.sh
touch toRemove.sh

touch toKeep.list
rm -f toKeep.list
touch toKeep.list

touch missing.list
rm -f missing.list
touch missing.list

miss=0
for block in  `seq -s ' ' $2 $3`; do
    blocktwo=0
    blocktwo=$(($blocktwo+$block))
    #echo -n ${block}:
    copies=`ls ${mydir} | grep "${filename}${block}_" | wc -l`
#    copies=`ls ${mydir} | grep _${block}_[0-9][_0-9] | grep -v [^${block}]_${block}_....root | wc -l`
    #grep "${filename}${block}_" | wc -l`
    toRemove=`expr $copies - 1`
    if [ $copies -gt 1 ]; then
      echo ${block}:${copies}
      ls -ltrh ${mydir} | grep "${filename}${block}_" | head -n ${toRemove} >> toRemove.list
      echo "--------------------------------------------------------------------------" >> toRemove.list
      ls -1tr ${mydir} | grep "${filename}${block}_" | head -n ${toRemove} >> toRemove.sh
      ls -ltrh ${mydir} | grep "${filename}${block}_" | tail -n 1 >> toKeep.list
    fi
    if [ $copies -lt 1 ]; then
	echo -n ${block}, >> missing.list
	miss=`expr $miss + 1`
    fi
done

touch toRemove.sh.temp
more toRemove.sh |while read line
do
  file=`echo $line`
  echo "srmrm  "$STORAGE${mydir}"/"${file} >> toRemove.sh.temp
done

rm -f toRemove.sh
mv toRemove.sh.temp toRemove.sh
chmod +x toRemove.sh

echo "The files to remove are"
more toRemove.list
echo "The file to keep are"
more toKeep.list
total=`more toRemove.sh | wc -l`
echo "The total number of files to remove is "${total}". To remove old duplicates execute ./toRemove.sh"
echo "There are ${miss} missing results:"
more missing.list
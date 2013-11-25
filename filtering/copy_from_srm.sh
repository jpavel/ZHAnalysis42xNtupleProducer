# script for copying the whole directory from IIHE storage to local directory. Use after establishing a grid proxy. The example command is
#
# source copy_from_srm.sh myDir/mySample
# 
# which would create a directory called myDir/mySample in your homedir and copy there the contents of directory /pnfs/iihe/cms/store/user/$USER/myDir/mySample

input_data=`echo $1`

DIRECTORY=srm://maite.iihe.ac.be//pnfs/iihe/cms/store/user/${USER}/
PNFS=/pnfs/iihe/cms/store/user/${USER}/



total=`ls ${PNFS}${input_data} | grep root | wc -l`
echo "Total number of input root files is" $total
mkdir -p ${input_data}
ls ${PNFS}${input_data} | grep root | while read line
do
 file=`echo $line`
 echo "srmcp ${DIRECTORY}/${input_data}/${file}  file:////user/${USER}/${input_data}/${file}"
 srmcp ${DIRECTORY}/${input_data}/${file}  file:////user/${USER}/${input_data}/${file}
done

more ${1} |while read line
do
  file=`echo $line`
  cd ${2}/job${file}
  echo "Submitting job ${file} ... "
#  sed -i 's/glidein/ZHttNtuples\/53X\/Data/g' ${file}.sh
#  cp *.sh ${file}.sh
  if [ "${4}" != "1" ]; then
    srmrm srm://maite.iihe.ac.be${3}/output_Ntuples_${file}_0_abc.root
  fi
  qsub -q localgrid@cream02 -o script.stdout -e script.stderr ${2}_${file}.sh
  cd ../..
done

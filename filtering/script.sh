#!/bin/bash          

##Some general shell commands
STR="Hello World!"
echo $STR    
echo ">> script.sh is checking where it is"
startDir=$PWD
echo ">> script.sh is checking how much disk space is still available"
df -h
echo ">> script.sh is listing files and directories in the current location"
ls -l
echo ">> script.sh is listing files and directories in userdir on storage element"
ls -l /pnfs/iihe/cms/store/user/$USER

##When accessing files on the storage element it is imporant to execute your code on the /scratch partition of the workernode you are running on. Therefore you need to copy your executable which is accessing/writing root files on
to the /scratch partition and execute it there. This is illustrated below.

echo ">> cd into /scratch partition"
cd /scratch
echo ">> ls of /scratch partition"
ls -l


echo ">> set environment"
pwd=$PWD
source $VO_CMS_SW_DIR/cmsset_default.sh

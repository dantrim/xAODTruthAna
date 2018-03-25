#!/bin/bash

echo "*********** GRID SCRIPT BEGIN **************"

# convert comma separated list into array
list=(`echo $1 | tr ',' ' '`)
nFiles=${#list[@]}
echo $nFiles input files
shift

# write the list of files to the filelist
if [ -f gridFileList.txt ]; then
    rm gridFileList.txt
fi

for file in ${list[@]}; do
    echo $file >> gridFileList.txt
done

# reamining options pass to the job
echo "Passing options to job: $@"

# now run the executable
runWWbb -f gridFileList.txt $@

exitcode=$?
if [[ $exitcode != 0 ]]; then
    echo "*********** GRID SCRIPT FAIL ***********"
    exit 1
fi

echo "*********** GRID SCRIPT END ***********"


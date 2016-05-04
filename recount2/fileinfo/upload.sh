#!/usr/bin/env bash
# Uploads files to Amazon Cloud Drive; copies multiple projects at once
# $1: path to ACD CLI binary
# $2: path to files to upload
# $3: path to destination on Cloud Drive
# $4: number of processes to run simultaneously
ACDCLI=$1
INPUT=$2
OUT=$3
PROCESSES=$4
THREADS=$5

function nrwait() {
    local nrwait_my_arg
    if [[ -z $1 ]] ; then
	nrwait_my_arg=2
    else
	nrwait_my_arg=$1
    fi
    
    while [[ $(jobs -p | wc -l) -ge $nrwait_my_arg ]] ; do
	sleep 0.33;
    done
  }

$ACDCLI sync

cd $INPUT
for project in *
do
	$ACDCLI upload -x $THREADS ${project}/ $OUT/
done
wait

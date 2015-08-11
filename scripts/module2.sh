#!/bin/bash

#liWGS-SV Pipeline: Module 2 (Per-Sample Clustering)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Submit bamstat job per sample
while read ID bam sex; do
  bsub -q normal -o ${OUTDIR}/logs/bamstat.log -e ${OUTDIR}/logs/bamstat.log -sla miket_sc -J ${COHORT_ID}_bamstat "${liWGS_SV}/scripts/bamstat.sh -s 3 -i ${WRKDIR}/${ID}/${ID}.bam -o ${WRKDIR}/${ID}/bamstat/" 
done < ${samples_list}

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_bamstat\|BAMSTAT_GATE" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_bamstat\|BAMSTAT_GATE" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Remove unnecessary pairs files
while read ID bam sex; do
  rm ${WRKDIR}/${ID}/bamstat/*pairs*txt
done < ${samples_list}
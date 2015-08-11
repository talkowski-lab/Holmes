#!/bin/bash

#liWGS-SV Pipeline: Module 3 (Physical Depth Analysis)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Run BinCov on samples
while read ID bam sex; do
  echo -e "${WRKDIR}/${ID}/${ID}.bam\t${ID}" >> ${WRKDIR}/cnMOPS.input
done < ${samples_list}
bsub -q big -R 'rusage[mem=20000]' -sla miket_sc -M 20000 -v 30000 -o ${OUTDIR}/logs/bincov.log -e ${OUTDIR}/logs/bincov.log -J ${COHORT_ID}_binCov "${liWGS_SV}/scripts/binCov.sh ${WRKDIR}/cnMOPS.input ${DICT} physical 1000 ${COHORT_ID} ${WRKDIR}/iCov"

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_binCov" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_binCov" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

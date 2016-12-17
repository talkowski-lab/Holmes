#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Module 2 (Physical Depth Analysis)

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

#Move, bgzip, and tabix individual coverage
while read ID bam sex; do
  mv ${WRKDIR}/iCov/raw_coverages/${ID}.coverage.bed ${WRKDIR}/${ID}/
  bgzip ${WRKDIR}/${ID}/${ID}.coverage.bed
  tabix -p bed ${WRKDIR}/${ID}/${ID}.coverage.bed.gz
done < ${samples_list}
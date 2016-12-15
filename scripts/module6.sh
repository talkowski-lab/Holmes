#!/bin/bash

#liWGS-SV Pipeline: Module 6 (Consensus CNV Merging)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#UPDATE FOR MkIII CONSENSUS: ALWAYS GENOTYPE & IGNORE DNACOPY CALLS (FOR NOW)
# #Set other params
# min_geno=20 #minimum cohort size for incorporating CNV genotype information into consensus script

#Create working directory for consensus CNVs
if ! [ -e ${WRKDIR}/consensusCNV ]; then
  mkdir ${WRKDIR}/consensusCNV
fi

#Merge cnMOPS calls
if [ -e ${WRKDIR}/cnMOPS_del_to_merge.list ]; then
  rm ${WRKDIR}/cnMOPS_del_to_merge.list
fi
if [ -e ${WRKDIR}/cnMOPS_dup_to_merge.list ]; then
  rm ${WRKDIR}/cnMOPS_dup_to_merge.list
fi
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed" >> ${WRKDIR}/cnMOPS_del_to_merge.list
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed" >> ${WRKDIR}/cnMOPS_dup_to_merge.list
done < ${samples_list}
mkdir ${WRKDIR}/consensusCNV/cnMOPS_chrsplit
for contig in $( seq 1 22 ) X Y; do
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/merge_cnMOPS.log -J ${COHORT_ID}_merge "${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/cnMOPS_del_to_merge.list 10000 ${contig} ${COHORT_ID}_cnMOPS_${contig}_dels ${WRKDIR}/consensusCNV/cnMOPS_chrsplit/"
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/merge_cnMOPS.log -J ${COHORT_ID}_merge "${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/cnMOPS_dup_to_merge.list 10000 ${contig} ${COHORT_ID}_cnMOPS_${contig}_dups ${WRKDIR}/consensusCNV/cnMOPS_chrsplit/"
done

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_merge" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_merge" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
    GATEwait=0
  fi
done

if [ -e ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dels.merged.bed ]; then
  rm ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dels.merged.bed
fi
if [ -e ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dups.merged.bed ]; then
  rm ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dups.merged.bed
fi
for contig in $( seq 1 22 ) X Y; do
  cat ${WRKDIR}/consensusCNV/cnMOPS_chrsplit/${COHORT_ID}_cnMOPS_${contig}_dels.merged.${contig}.bed >> ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dels.merged.bed
  cat ${WRKDIR}/consensusCNV/cnMOPS_chrsplit/${COHORT_ID}_cnMOPS_${contig}_dups.merged.${contig}.bed >> ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dups.merged.bed
done

#Launch consensus CNV pipeline
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/consensusCNV.log -J ${COHORT_ID}_CONSENSUS_CNVS "${liWGS_SV}/scripts/consensusCNV.sh ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dels.merged.bed del ${params}"
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/consensusCNV.log -J ${COHORT_ID}_CONSENSUS_CNVS "${liWGS_SV}/scripts/consensusCNV.sh ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/consensusCNV/${COHORT_ID}_cnMOPS_dups.merged.bed dup ${params}"

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_CONSENSUS_CNVS" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_CONSENSUS_CNVS" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
    GATEwait=0
  fi
done

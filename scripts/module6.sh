#!/bin/bash

#liWGS-SV Pipeline: Module 6 (Consensus CNV Merging)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Set other params
min_geno=20 #minimum cohort size for incorporating CNV genotype information into consensus script

#Create working directory for consensus CNVs
mkdir ${WRKDIR}/consensusCNV

#Genotype all CNV intervals if cohort has ≥ ${min_geno} samples
if [ $( cat ${samples_list} | wc - ) -ge ${min_geno} ]; then
  #Create master file of all CNV intervals
  if [ -e ${WRKDIR}/consensusCNV/dels_to_merge.list ]; then
    rm ${WRKDIR}/consensusCNV/dels_to_merge.list
  fi
  if [ -e ${WRKDIR}/consensusCNV/dups_to_merge.list ]; then
    rm ${WRKDIR}/consensusCNV/dups_to_merge.list
  fi
  while read ID bam sex; do
    fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dels.bed
    fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dups.bed
    fgrep del ${WRKDIR}/${ID}/DNAcopy/${ID}.events.tsv | awk -v ID=${ID} -v OFS="\t" '{ if ($3-$2>=10000000 && $1!="X" && $1!="Y") print $1, $2, $3, ID, $4, $7 }' > ${WRKDIR}/${ID}/DNAcopy.dels.bed
    fgrep dup ${WRKDIR}/${ID}/DNAcopy/${ID}.events.tsv | awk -v ID=${ID} -v OFS="\t" '{ if ($3-$2>=10000000 && $1!="X" && $1!="Y") print $1, $2, $3, ID, $4, $7 }' > ${WRKDIR}/${ID}/DNAcopy.dups.bed
    echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dels.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed\n${ID}_DNAcopy\t${WRKDIR}/${ID}/DNAcopy.dels.bed" >> ${WRKDIR}/consensusCNV/dels_to_merge.list
    echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dups.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed\n${ID}_DNAcopy\t${WRKDIR}/${ID}/DNAcopy.dups.bed" >> ${WRKDIR}/consensusCNV/dups_to_merge.list
  done < ${samples_list}
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/mergeCNVinterval.log -e ${OUTDIR}/logs/mergeCNVinterval.log -J ${COHORT_ID}_MERGE "${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensusCNV/dels_to_merge.list 10000 ${COHORT_ID}_del_intervals ${WRKDIR}/consensusCNV/"
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/mergeCNVinterval.log -e ${OUTDIR}/logs/mergeCNVinterval.log -J ${COHORT_ID}_MERGE "${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensusCNV/dups_to_merge.list 10000 ${COHORT_ID}_dup_intervals ${WRKDIR}/consensusCNV/"

  #Gate until complete; 20 sec check; 5 min report
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MERGE" | wc -l )
  GATEwait=0
  until [[ $GATEcount == 0 ]]; do
    sleep 20s
    GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MERGE" | wc -l )
    GATEwait=$[${GATEwait} +1]
    if [[ $GATEwait == 15 ]]; then
      echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
      GATEwait=0
    fi
  done

  #Set rare edge cases where overclustering caused negative size to largest insert in cohort
  scaler=$( fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | cut -f9 | sort -nrk1,1 | head -n1 )
  awk -v scaler=${scaler} -v OFS="\t" '{ if ($3-$2>0) print $0; else print $1, (($3+$2)/2)-(scaler/2), (($3+$2)/2)+(scaler/2), $4, $5, $6, $7, $8, $9 }' ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.bed > ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.bed2
  mv ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.bed2 ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.bed
  awk -v scaler=${scaler} -v OFS="\t" '{ if ($3-$2>0) print $0; else print $1, (($3+$2)/2)-(scaler/2), (($3+$2)/2)+(scaler/2), $4, $5, $6, $7, $8, $9 }' ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.bed > ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.bed2
  mv ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.bed2 ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.bed

  #Genotype all consensus intervals
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log -J ${COHORT_ID}_genotyping "module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV_batch.R ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.bed ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.genotypes.bed"
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log -J ${COHORT_ID}_genotyping "module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV_batch.R ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.bed ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.genotypes.bed"

  #Gate until complete; 20 sec check; 5 min report
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_genotyping" | wc -l )
  GATEwait=0
  until [[ $GATEcount == 0 ]]; do
    sleep 20s
    GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_genotyping" | wc -l )
    GATEwait=$[${GATEwait} +1]
    if [[ $GATEwait == 15 ]]; then
      echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
      GATEwait=0
    fi
  done

  #Split genotype beds per sample
  while read ID bam sex; do
    idx=$( head -n1 ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.genotypes.bed | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' | fgrep -w ${ID} | cut -f1 )
    awk -v idx=${idx} -v OFS="\t" '{ print $2, $3, $4, $(idx), $1 }' ${WRKDIR}/consensusCNV/${COHORT_ID}_del_intervals.merged.genotypes.bed > ${WRKDIR}/${ID}/${ID}.merged_del.genotypes.bed
    awk -v idx=${idx} -v OFS="\t" '{ print $2, $3, $4, $(idx), $1 }' ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_intervals.merged.genotypes.bed > ${WRKDIR}/${ID}/${ID}.merged_dup.genotypes.bed
  done < ${samples_list}

  #Generate consensus CNVs per sample
  while read ID bam sex; do
    ${liWGS_SV}/scripts/consensusCNV_wGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed ${WRKDIR}/${ID}/DNAcopy.dels.bed ${WRKDIR}/${ID}/${ID}.merged_del.genotypes.bed ${ID} del ${params}
    ${liWGS_SV}/scripts/consensusCNV_wGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed ${WRKDIR}/${ID}/DNAcopy.dups.bed ${WRKDIR}/${ID}/${ID}.merged_dup.genotypes.bed ${ID} dup ${params}
  done < ${samples_list}
else
  echo "WARNING [MODULE 6]: Cohort has fewer than ${min_geno} samples; consensus CNVs will not incorporate joint genotyping information" >> ${OUTDIR}/${COHORT_ID}_WARNINGS.txt
  #Generate consensus CNVs per sample
  while read ID bam sex; do
    ${liWGS_SV}/scripts/consensusCNV_noGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed ${WRKDIR}/${ID}/DNAcopy.dels.bed ${ID} del ${params}
    ${liWGS_SV}/scripts/consensusCNV_noGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed ${WRKDIR}/${ID}/DNAcopy.dups.bed ${ID} dup ${params}
  done < ${samples_list}
fi

#Merge all consensus CNVs
if [ -e ${WRKDIR}/consensus_del_to_merge.list ]; then
  rm ${WRKDIR}/consensus_del_to_merge.list
fi
if [ -e ${WRKDIR}/consensus_dup_to_merge.list ]; then
  rm ${WRKDIR}/consensus_dup_to_merge.list
fi
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.consensus.del.bed" >> ${WRKDIR}/consensus_del_to_merge.list
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.consensus.dup.bed" >> ${WRKDIR}/consensus_dup_to_merge.list
done < ${samples_list}
${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensus_del_to_merge.list 10000 ${COHORT_ID}_consensus_dels ${WRKDIR}/consensusCNV/
${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensus_dup_to_merge.list 10000 ${COHORT_ID}_consensus_dups ${WRKDIR}/consensusCNV/

#Genotype merged consensus CNVs
if [ $( cat ${samples_list} | wc - ) -ge ${min_geno} ]; then
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/consensus_genotyping.log -e ${OUTDIR}/logs/consensus_genotyping.log -J ${COHORT_ID}_genotyping "module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV_batch.R ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_dels.merged.bed ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${WRKDIR}/consensusCNV/${COHORT_ID}_del_consensus.merged.genotypes.bed"
  bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/consensus_genotyping.log -e ${OUTDIR}/logs/consensus_genotyping.log -J ${COHORT_ID}_genotyping "module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV_batch.R ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_dups.merged.bed ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${WRKDIR}/consensusCNV/${COHORT_ID}_dup_consensus.merged.genotypes.bed"
done < ${samples_list}

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_genotyping" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_genotyping" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
    GATEwait=0
  fi
done
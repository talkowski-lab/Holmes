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
min_geno=30 #minimum cohort size for incorporating CNV genotype information into consensus script

#Create working directory for consensus CNVs
mkdir ${WRKDIR}/consensusCNV

#Genotype all CNV intervals if cohort has ≥ 30 samples
if [ $( cat ${samples_list} | wc - ) -ge ${min_geno} ]; then
  #Create master file of all CNV intervals
  while read ID bam sex; do
    fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dels.bed
    fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dups.bed
    fgrep del ${WRKDIR}/${ID}/DNAcopy/${ID}.events.tsv | awk -v ID=${ID} -v OFS="\t" '{ if ($3-$2>=10000000 && $1!="X" && $1!="Y") print $1, $2, $3, ID, $4, $7 }' > ${WRKDIR}/${ID}/DNAcopy.dels.bed
    fgrep dup ${WRKDIR}/${ID}/DNAcopy/${ID}.events.tsv | awk -v ID=${ID} -v OFS="\t" '{ if ($3-$2>=10000000 && $1!="X" && $1!="Y") print $1, $2, $3, ID, $4, $7 }' > ${WRKDIR}/${ID}/DNAcopy.dups.bed
    echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dels.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed\n${ID}_DNAcopy\t${WRKDIR}/${ID}/DNAcopy.dels.bed" >> ${WRKDIR}/consensusCNV/dels_to_merge.list
    echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dups.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed\n${ID}_DNAcopy\t${WRKDIR}/${ID}/DNAcopy.dups.bed" >> ${WRKDIR}/consensusCNV/dups_to_merge.list
  done < ${samples_list}
  ${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensusCNV/dels_to_merge.list 10000 ${COHORT_ID}_del_intervals ${WRKDIR}/consensusCNV/
  ${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensusCNV/dups_to_merge.list 10000 ${COHORT_ID}_dup_intervals ${WRKDIR}/consensusCNV/
  ## ****ADD GENOTYPING HERE****
  ## ****ADD REVISED CONSENSUS PIPELINE WITH GENOTYPING INFO HERE****
else
  echo "WARNING [MODULE 6]: Cohort has fewer than ${min_geno} samples; consensus CNVs will not incorporate joint genotyping information" >> ${OUTDIR}/${COHORT_ID}_WARNINGS.txt
  while read ID bam sex; do
    ${liWGS_SV}/scripts/consensusCNV_noGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed ${WRKDIR}/${ID}/DNAcopy.dels.bed ${ID} del ${params}
    ${liWGS_SV}/scripts/consensusCNV_noGeno.sh ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed ${WRKDIR}/${ID}/DNAcopy.dups.bed ${ID} dup ${params}
  done < ${samples_list}
  


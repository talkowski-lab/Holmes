#!/bin/bash

#liWGS-SV Pipeline: Module 6 (Consensus CNV Merging)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Create working directory for consensus CNVs
mkdir ${WRKDIR}/consensusCNV

#Create master file of all CNV intervals
while read ID bam sex; do
  fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dels.bed
  fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe | fgrep -w ${ID} | awk -v ID=${ID} -v OFS="\t" '{ print $1, $3, $5, ID, $7 }' > ${WRKDIR}/${ID}/classifier.dups.bed
  #*******DNAcopy refinement********
  #MUST ADD DNAcopy to list of beds to merge once ready
  # awk -v ID=${ID} -v OFS="\t" '{ if ($3-$2>=1000000 && $1!="X" && $1!="Y") print $1, $2, $3, ID, (CONFIDENCE SCORE HERE) }' (DNAcopy output)
  echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dels.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed" >> ${WRKDIR}/consensusCNV/dels_to_merge.list
  echo -e "${ID}_classifier\t${WRKDIR}/${ID}/classifier.dups.bed\n${ID}_cnMOPS\t${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed" >> ${WRKDIR}/consensusCNV/dups_to_merge.list
done < ${samples_list}
${liWGS_SV}/scripts/mergebeds.sh ${WRKDIR}/consensusCNV/dels_to_merge.list 10000 ${COHORT_ID} ${WRKDIR}/consensusCNV/
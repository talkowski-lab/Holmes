#!/bin/bash

#liWGS-SV Pipeline: Module 5 (Joint Clustering & Classification)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Prepare classifier files
mkdir ${WRKDIR}/classifier
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/bamstat"
done < ${samples_list} > ${WRKDIR}/classifier/classifier.samples.list
while read ID bam sex; do
  echo -e "${ID}\t${WRKDIR}/${ID}/${ID}.coverage.bed.gz"
done < ${samples_list} > ${WRKDIR}/classifier/classifier.icov.list

#Run classifier
cd ${WRKDIR}/classifier
bsub -u nobody -q normal -sla miket_sc -o ${OUTDIR}/logs/classifier.log -e ${OUTDIR}/logs/classifier.log -J ${COHORT_ID}_classifier "${CLASSIFIER_DIR}/run_classify.sh -m ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ${WRKDIR}/classifier/classifier.samples.list ${COHORT_ID} ${WRKDIR}/classifier/classifier.icov.list"

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Apply cluster fix
mkdir ${WRKDIR}/classifier/clusterfix
mkdir ${WRKDIR}/classifier/clusterfix/delsplit
for contig in $( seq 1 22 ) X Y; do
  awk -v OFS="\t" -v contig=${contig} '{ if ($4==contig || $4=="") print $0 }' ${WRKDIR}/classifier/${COHORT_ID}_deletion.clusters | cat -s | sed '1{/^$/d}' > ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.clusters
  bsub -u nobody -q normal -o ${OUTDIR}/logs/classifier_clusterfix.log -e ${OUTDIR}/logs/classifier_clusterfix.log -sla miket_sc -J ${COHORT_ID}_clusterPatch "${CLASSIFIER_DIR}/cleanClusters_patch.sh ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.clusters ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.patched.clusters ${uscore_skip} TRUE"
done
for class in insertion inversion transloc; do
  bsub -u nobody -q normal -o ${OUTDIR}/logs/classifier_clusterfix.log -e ${OUTDIR}/logs/classifier_clusterfix.log -sla miket_sc -J ${COHORT_ID}_clusterPatch "${CLASSIFIER_DIR}/cleanClusters_patch.sh ${WRKDIR}/classifier/${COHORT_ID}_${class}.clusters ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${uscore_skip} FALSE"
done

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_clusterPatch" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_clusterPatch" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Concatenate split deletion patched clusters
for contig in $( seq 1 22 ) X Y; do
  cat ${WRKDIR}/classifier/clusterfix/delsplit/${COHORT_ID}_deletion.chr${contig}.patched.clusters >> ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_deletion.patched.clusters
done

#Complete classifier on fixed clusters
for class in deletion insertion inversion transloc; do
  bsub -u nobody -o ${OUTDIR}/logs/classifier_r2.log -e ${OUTDIR}/logs/classifier_r2.log -q normal -sla miket_sc -J ${COHORT_ID}_classifier_r2 "source /apps/lab/miket/anaconda/4.0.5_py3/bin/deactivate; ${CLASSIFIER_DIR}/rpc_classify.py ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.events.bedpe ${WRKDIR}/classifier/${COHORT_ID}_boot.list ${COHORT_ID}_${class} --cluster-bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.bkpts.bedpe"
done

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier_r2" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_classifier_r2" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Reclassify fixed clusters
mkdir ${WRKDIR}/classifier/clusterfix/newCoords
for class in deletion insertion inversion transloc; do
  bsub -u nobody -q normal -o ${OUTDIR}/logs/reclassify.log -e ${OUTDIR}/logs/reclassify.log -sla miket_sc -J ${COHORT_ID}_reclassify "${CLASSIFIER_DIR}/reclassify_output.sh ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.events.bedpe ${class} ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${class}.patched.clusters ${WRKDIR}/classifier/clusterfix/newCoords FALSE"
done

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_reclassify" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_reclassify" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done
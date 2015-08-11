#!/bin/bash

#liWGS-SV Pipeline: MASTER PIPELINE SCRIPT
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Submit this script to long queue on ERISOne to run entire pipeline

#Read input
samples_list=$1
params=$2

#Source params file
echo -e "STATUS [$(date)]: Loading parameters..."
. ${params}

#Set up output directory tree
echo -e "STATUS [$(date)]: Creating directory trees..."
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
mkdir ${OUTDIR}/QC
mkdir ${OUTDIR}/QC/cohort
mkdir ${OUTDIR}/QC/sample
while read ID bam sex; do
  mkdir ${OUTDIR}/QC/sample/${ID}
done < ${samples_list}
mkdir ${OUTDIR}/data
mkdir ${OUTDIR}/data/seqDepth
mkdir ${OUTDIR}/data/clusters
mkdir ${OUTDIR}/SV_calls
mkdir ${OUTDIR}/logs

#Set up working directory tree
if ! [ -e ${WRKDIR} ]; then
  mkdir ${WRKDIR}
fi
while read ID bam sex; do
  mkdir ${WRKDIR}/${ID}
done < ${samples_list}

#Link & index all bams
echo -e "STATUS [$(date)]: Indexing BAMs..."
while read ID bam sex; do
  ln -s ${bam} ${WRKDIR}/${ID}/${ID}.bam
  bsub -q short -e ${OUTDIR}/logs/bam_index.log -o ${OUTDIR}/logs/bam_index.log -u nobody -sla miket_sc -J ${ID}_index "sambamba index ${WRKDIR}/${ID}/${ID}.bam"
done < ${samples_list}

##STAGE 1: modules 1, 2, and 3
echo -e "STATUS [$(date)]: Beginning PHASE 1..."
#Submit module 1 (QC)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module1.log -e ${OUTDIR}/logs/module1.log -u nobody -J ${COHORT_ID}_MODULE_1 "${liWGS_SV}/scripts/module1.sh ${samples_list} ${params}"
#Submit module 2 (per-sample clustering)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module2.log -e ${OUTDIR}/logs/module2.log -u nobody -J ${COHORT_ID}_MODULE_2 "${liWGS_SV}/scripts/module2.sh ${samples_list} ${params}"
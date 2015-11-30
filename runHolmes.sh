#!/bin/bash

#liWGS-SV Pipeline: MASTER PIPELINE SCRIPT
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Submit this script to long queue on ERISOne to run entire pipeline

#Read input
samples_list=$1
params=$2

#Ensure correct version of gcc for anaconda dependencies
module rm gcc/4.9.0
module rm gcc-4.4
module load gcc/4.9.0

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
mkdir ${OUTDIR}/plots
cp ${params} ${OUTDIR}/${COHORT_ID}.run_parameters.sh

#Set up working directory tree
if ! [ -e ${WRKDIR} ]; then
  mkdir ${WRKDIR}
fi
while read ID bam sex; do
  mkdir ${WRKDIR}/${ID}
done < ${samples_list}
cp ${liWGS_SV}/scripts/SE_largeCNV/* ${WRKDIR}/

#Write start date to temporary file
echo $(date) > ${WRKDIR}/start.tmp

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
#Submit module 2 (physical depth analysis)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module2.log -e ${OUTDIR}/logs/module2.log -u nobody -J ${COHORT_ID}_MODULE_2 "${liWGS_SV}/scripts/module2.sh ${samples_list} ${params}"
#Submit module 3 (per-sample clustering)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module3.log -e ${OUTDIR}/logs/module3.log -u nobody -J ${COHORT_ID}_MODULE_3 "${liWGS_SV}/scripts/module3.sh ${samples_list} ${params}"

#Gate until modules 1 & 2 complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_MODULE_2" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_1\|${COHORT_ID}_MODULE_2" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at PHASE 1a..."
    GATEwait=0
  fi
done

##STAGE 2a: module 4
echo -e "STATUS [$(date)]: PHASE 1a complete; Beginning PHASE 2a..."
#Submit module 4 (physical depth CNV calling)
#If DNAcopy is included in module 4, this job to go to big for creation of DNAcopy reference profiles (requires > 60GB memory for ~300 samples)
#As of Nov 30, 2015, DNAcopy is not run by default because it's too slow, too memory-intensive, and not highly informative. This decision can be revisited later on.
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module4.log -e ${OUTDIR}/logs/module4.log -u nobody -J ${COHORT_ID}_MODULE_4 "${liWGS_SV}/scripts/module4.sh ${samples_list} ${params}"

#Gate until module 3 complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_3" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_3" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at PHASE 1b..."
    GATEwait=0
  fi
done

##STAGE 2b: module 5
echo -e "STATUS [$(date)]: PHASE 1b complete; Beginning PHASE 2b..."
#Submit module 5 (joint clustering)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module5.log -e ${OUTDIR}/logs/module5.log -u nobody -J ${COHORT_ID}_MODULE_5 "${liWGS_SV}/scripts/module5.sh ${samples_list} ${params}"

#Gate until modules 4 & 5 complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_4\|${COHORT_ID}_MODULE_5" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_4\|${COHORT_ID}_MODULE_5" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at PHASE 2..."
    GATEwait=0
  fi
done

##STAGE 3: modules 6 and 7
echo -e "STATUS [$(date)]: Beginning PHASE 3..."
#Submit module 6 (consensus CNV merging)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module6.log -e ${OUTDIR}/logs/module6.log -u nobody -J ${COHORT_ID}_MODULE_6 "${liWGS_SV}/scripts/module6.sh ${samples_list} ${params}"
#Submit module 7 (complex SV categorization)
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module7.log -e ${OUTDIR}/logs/module7.log -u nobody -J ${COHORT_ID}_MODULE_7 "${liWGS_SV}/scripts/module7.sh ${samples_list} ${params}"

#Gate until modules 6 & 7 completes; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_6\|${COHORT_ID}_MODULE_7" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_6\|${COHORT_ID}_MODULE_7" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at PHASE 3..."
    GATEwait=0
  fi
done

##STAGE 4: module 8
echo -e "STATUS [$(date)]: Beginning PHASE 4..."
${liWGS_SV}/scripts/module8.sh ${samples_list} ${params}

##STAGE 5: module 9
echo -e "STATUS [$(date)]: Beginning PHASE 5..."
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/module9.log -e ${OUTDIR}/logs/module9.log -u nobody -J ${COHORT_ID}_MODULE_9 "${liWGS_SV}/scripts/module9.sh ${samples_list} ${params}"

#Gate until module 9 completes; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_9" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_MODULE_9" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at PHASE 5..."
    GATEwait=0
  fi
done

#Move all relevant files to OUTDIR
cp ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${OUTDIR}/data/seqDepth/
cp ${WRKDIR}/final_variants/* ${OUTDIR}/SV_calls
cp ${WRKDIR}/raw_clusters/* ${OUTDIR}/data/clusters
mkdir ${OUTDIR}/SV_calls/annotations
cp ${WRKDIR}/annotations/*_gene_anno.bed ${OUTDIR}/SV_calls/annotations/

#Summarize run outcomes
echo -e "STATUS [$(date)]: Calculating run summary metrics..."
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/finalSummary.log -e ${OUTDIR}/logs/finalSummary.log -u nobody -J ${COHORT_ID}_SUMMARIZE "${liWGS_SV}/scripts/Holmes_summary.sh ${samples_list} ${params}"


#Remove working directory unless specified otherwise
if ! [ ${KEEP_TMP}=="TRUE" ]; then
  rm -rf ${WRKDIR}
fi

#Print final completion notice
echo -e "-----\nSTATUS [$(date)]: Holmes liWGS-SV discovery on ${COHORT_ID} complete\n-----\n\n"

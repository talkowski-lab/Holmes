#!/bin/bash

#liWGS-SV Pipeline: Module 1 (QC)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Submit library metrics QC jobs
while read ID bam sex; do
  bsub -u nobody -o ${OUTDIR}/logs/input_QC.log -e ${OUTDIR}/logs/input_QC.log -sla miket_sc -q normal -J ${COHORT_ID}_QC "sambamba view -h -f bam -F 'not secondary_alignment' ${WRKDIR}/${ID}/${ID}.bam | java -Xmx3g -jar ${PICARD} CollectMultipleMetrics I=/dev/stdin O=${OUTDIR}/QC/sample/${ID}/${ID} AS=true R=${REF} VALIDATION_STRINGENCY=SILENT PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics"
  bsub -u nobody -o ${OUTDIR}/logs/input_QC.log -e ${OUTDIR}/logs/input_QC.log -sla miket_sc -q medium -J ${COHORT_ID}_QC "sambamba view -h -f bam -F 'not secondary_alignment' ${WRKDIR}/${ID}/${ID}.bam | samtools flagstat /dev/stdin > ${OUTDIR}/QC/sample/${ID}/${ID}.flagstat"
  bsub -u nobody -o ${OUTDIR}/logs/input_QC.log -e ${OUTDIR}/logs/input_QC.log -sla miket_sc -q medium -J ${COHORT_ID}_QC "sambamba view -h -f bam -F 'not secondary_alignment' ${WRKDIR}/${ID}/${ID}.bam | bamtools stats -in /dev/stdin -insert > ${OUTDIR}/QC/sample/${ID}/${ID}.stats"
  bsub -u nobody -o ${OUTDIR}/logs/input_QC.log -e ${OUTDIR}/logs/input_QC.log -sla miket_sc -q normal -J ${COHORT_ID}_QC "java -Xmx3g -jar ${PICARD} CollectWgsMetrics I=${WRKDIR}/${ID}/${ID}.bam O=${OUTDIR}/QC/sample/${ID}/${ID}.wgs R=${REF} VALIDATION_STRINGENCY=SILENT"
  bsub -u nobody -o ${OUTDIR}/logs/input_QC.log -e ${OUTDIR}/logs/input_QC.log -sla miket_sc -q big -R 'rusage[mem=25000]' -M 25000 -J ${COHORT_ID}_QC "java -Xmx20g -jar ${PICARD} EstimateLibraryComplexity I=${WRKDIR}/${ID}/${ID}.bam O=${OUTDIR}/QC/sample/${ID}/${ID}.complexity VALIDATION_STRINGENCY=SILENT"
done < ${samples_list}

#Submit dosage QC
while read ID bam sex; do
  bsub -u nobody -o ${OUTDIR}/logs/dosageCheck.log -e ${OUTDIR}/logs/dosageCheck.log -sla miket_sc -q normal -J ${COHORT_ID}_QC "${liWGS_SV}/scripts/WGScheckDosage.sh ${WRKDIR}/${ID}/${ID}.bam ${ID} h37 both ${OUTDIR}/QC/sample/${ID}/"
done < ${samples_list}

#Submit sex check
while read ID bam sex; do
  bsub -u nobody -o ${OUTDIR}/logs/dosageCheck.log -e ${OUTDIR}/logs/dosageCheck.log -sla miket_sc -q short -J ${COHORT_ID}_QC "${liWGS_SV}/scripts/sexCheck.sh ${ID} ${WRKDIR}/${ID}/${ID}.bam ${params}"
done < ${samples_list}
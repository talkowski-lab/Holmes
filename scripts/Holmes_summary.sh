#!/bin/bash

#liWGS-SV Pipeline: End-of-Run Summary Script
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Calculate maximum number of samples to correspond with $polyArt_filter
maxCount=$( echo "(( ${polyArt_filter}*$( cat ${samples_list} | wc -l ) ))" | bc )

#Create report with number of variants per sample
echo -e "##liWGS-SV PIPELINE COHORT SV REPORT\n##RUN DATE: $( echo $(date) | awk '{ print $1, $2, $3, $NF }' )\n\
#ID\tDeletion_wArt\tDuplication_wArt\tInsertion_wArt\tInversion_wArt\tTranslocation_wArt\tComplex_wArt\tUnresolved_wArt\tDeletion_noArt\tDuplication_noArt\tInsertion_noArt\tInversion_noArt\tTranslocation_noArt\tComplex_noArt\tUnresolved_noArt" > ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics
while read ID bam sex; do
  echo ${ID}
  for type in deletion duplication insertion inversion translocation complex unresolved; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.be* | wc -l
  done
  for type in deletion duplication; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bed | awk -v max=${maxCount} '{ if ($7<=max) print $0 }' | wc -l
  done
  for type in insertion inversion translocation; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bedpe | awk -v max=${maxCount} '{ if ($10<=max) print $0 }' | wc -l
  done
  for type in complex unresolved; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bed | awk -v max=${maxCount} '{ if ($6<=max) print $0 }' | wc -l
  done
done < ${samples_list} | paste - - - - - - - - - - - - - - - >> ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics

#Write final report to base level $OUTDIR
echo -e "##################################\n\
#   HOLMES liWGS-SV RUN REPORT   #\n\
##################################\n\n\
+-------------+
| RUN DETAILS |\n\
+-------------+\n\
Cohort: ${COHORT_ID} (n=$( cat ${samples_list} | wc -l))\n\
Started: $( cat ${WRKDIR}/start.tmp | awk '{ print $1, $2, $3, $NF }' )\n\
Finished: $( echo $(date) | awk '{ print $1, $2, $3, $NF }' )\n\
User: ${USER}\n\n\
+-------------------------+\n\
| AVERAGE LIBRARY METRICS |\n\
+-------------------------+" > ${OUTDIR}/${COHORT_ID}.run_summary.txt
Rscript ${liWGS_SV}/scripts/stripQC.R ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ${OUTDIR}/${COHORT_ID}.run_summary.txt
echo -e "\n"

Also include:
clustering (num clusters per sample, outliers)
depth (num depth calls per sample, outliers)
SV calls (num final calls per sample, outliers)
SV sizes (median & IQR per class)
gene annotations (num genes LOF/dup/etc per sample, outliers)

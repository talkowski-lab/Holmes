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
Started: $( cat ${WRKDIR}/start.tmp )\n\
Finished: $( echo $(date) )\n\
User: ${USER}
Output Directory: ${OUTDIR}\n\
Temporary Directory: ${WRKDIR}\n\
Artifact VAF Filter: ${polyArt_filter}\n\
Reference Genome: ${REF}\n\n\
+-------------------------+\n\
| AVERAGE LIBRARY METRICS |\n\
+-------------------------+" > ${OUTDIR}/${COHORT_ID}.run_summary.txt
Rscript ${liWGS_SV}/scripts/stripAndPlotLibraryQC.R ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ${OUTDIR}/${COHORT_ID}.run_summary.txt ${OUTDIR}/plots/${COHORT_ID}.libraryQC.pdf
echo -e "\n\
+----------------+\n\
| RAW SV SIGNALS |\n\
+----------------+" >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo " => Raw Deletion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.events.bedpe | wc -l | paste -d" " <( echo " => Classified Valid Deletion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
while read ID bam sex; do
  cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed | wc -l
done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo " => Average cnMOPS Deletions per Sample:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo " => Raw Insertion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.events.bedpe | wc -l | paste -d" " <( echo " => Classified Valid Insertion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
while read ID bam sex; do
  cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed | wc -l
done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo " => Average cnMOPS Duplications per Sample:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo " => Raw Inversion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.events.bedpe | wc -l | paste -d" " <( echo " => Classified Valid Inversion Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.transloc.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo " => Raw Translocation Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.transloc.events.bedpe | wc -l | paste -d" " <( echo " => Classified Valid Translocation Clusters:" ) - >> ${OUTDIR}/${COHORT_ID}.run_summary.txt


Also include:
clustering (num clusters per sample, outliers)
depth (num depth calls per sample, outliers)
SV calls (num final calls per sample, outliers)
SV sizes (median & IQR per class)
gene annotations (num genes LOF/dup/etc per sample, outliers)

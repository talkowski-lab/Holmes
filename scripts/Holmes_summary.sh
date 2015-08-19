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
#ID\tDeletion_All_wArt\tDeletion_HQ_wArt\tDuplication_All_wArt\tDuplication_HQ_wArt\tInsertion_wArt\tInversion_wArt\tTranslocation_wArt\tComplex_wArt\tUnresolved_wArt\tDeletion_All_noArt\tDeletion_HQ_noArt\tDuplication_All_noArt\tDuplication_HQ_noArt\tInsertion_noArt\tInversion_noArt\tTranslocation_noArt\tComplex_noArt\tUnresolved_noArt" > ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics
while read ID bam sex; do
  echo ${ID}
  for type in deletion duplication insertion inversion translocation complex unresolved; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.be* | wc -l
    if [ ${type} == "deletion" ] || [ ${type} == "duplication" ]; then
      grep -e 'HIGH\|MED' ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.be* | fgrep -w ${ID} | wc -l
    fi
  done
  for type in deletion duplication; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bed | awk -v max=${maxCount} '{ if ($7<=max) print $0 }' | wc -l
    grep -e 'HIGH\|MED' ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.be* | fgrep -w ${ID} | awk -v max=${maxCount} '{ if ($7<=max) print $0 }' | wc -l
  done
  for type in insertion inversion translocation; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bedpe | awk -v max=${maxCount} '{ if ($10<=max) print $0 }' | wc -l
  done
  for type in complex unresolved; do
    fgrep -w ${ID} ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bed | awk -v max=${maxCount} '{ if ($6<=max) print $0 }' | wc -l
  done
done < ${samples_list} | paste - - - - - - - - - - - - - - - - - - - >> ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics

#Get SV Sizes
mkdir ${WRKDIR}/sizes
for type in deletion duplication; do
  awk '{ print $3-$2 }' ${OUTDIR}/SV_calls/${COHORT_ID}.${type}.bed | sort -nk1,1 > ${WRKDIR}/sizes/${type}.size
done
awk '{ print $3-$2 }' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe | sort -nk1,1 > ${WRKDIR}/sizes/insertion_source.size
awk '{ print $6-$5 }' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe | sort -nk1,1 > ${WRKDIR}/sizes/insertion_sink.size

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
+-----------------+\n\
| LIBRARY METRICS |\n\
+-----------------+\n\
Metric\tMin\tQ1\tMed\tMean\tQ3\tMax\tOutliers" > ${OUTDIR}/${COHORT_ID}.run_summary.txt
Rscript ${liWGS_SV}/scripts/stripAndPlotLibraryQC.R ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ${OUTDIR}/${COHORT_ID}.run_summary.txt ${OUTDIR}/plots/${COHORT_ID}.libraryQC.pdf
echo -e "\n\
+----------------+\n\
| RAW SV SIGNALS |\n\
+----------------+" >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
for dummy in 1; do
  echo " => Deletion"
  tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  while read ID bam sex; do
    cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed | wc -l
  done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo "      -Mean cnMOPS per Sample:" ) -
  echo " => Insertion"
  tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  while read ID bam sex; do
    cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed | wc -l
  done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo "      -Mean cnMOPS per Sample:" ) -
  echo " => Inversion"
  tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  echo " => Translocation"
  tail -n2 ${OUTDIR}/data/clusters/${COHORT_ID}.transloc.clusters.txt | head -n1 | cut -f1 | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.transloc.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
done >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
echo -e "\n\
+----------------+\n\
| FINAL SV CALLS |\n\
+----------------+\n\
 => Deletion Variants [All]\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed | wc -l ) ($( awk -v max=${maxCount} '{ if ($7>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f2 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f11 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Deletion Variants [HQ]\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed | grep -e 'HIGH\|MED' | wc -l ) ($( awk -v max=${maxCount} '{ if ($7>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed | grep -e 'HIGH\|MED' | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f3 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f12 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Tandem Duplication Variants [All]\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed | wc -l ) ($( awk -v max=${maxCount} '{ if ($7>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f4 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f13 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Tandem Duplication Variants [HQ]\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed | grep -e 'HIGH\|MED' | wc -l ) ($( awk -v max=${maxCount} '{ if ($7>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed | grep -e 'HIGH\|MED' | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f5 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f14 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Insertion Variants\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe | wc -l ) ($( awk -v max=${maxCount} '{ if ($10>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f6 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f15 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Simple Inversion Variants\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.inversion.bedpe | wc -l ) ($( awk -v max=${maxCount} '{ if ($10>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.inversion.bedpe | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f7 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f16 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Translocation Variants\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.translocation.bedpe | wc -l ) ($( awk -v max=${maxCount} '{ if ($10>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.translocation.bedpe | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f8 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f17 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Complex Variants\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.complex.bed | wc -l ) ($( awk -v max=${maxCount} '{ if ($6>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.complex.bed | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f9 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f18 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
 => Incompletely Resolved Sites\n\
      -Total in cohort: $( sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.unresolved.bed | wc -l ) ($( awk -v max=${maxCount} '{ if ($6>max) print $0 }' ${OUTDIR}/SV_calls/${COHORT_ID}.unresolved.bed | fgrep -v "#" | wc -l ) putative artifacts)\n\
      -Median per sample (w/ artifacts): $( cut -f10 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
      -Median per sample (w/o artifacts): $( cut -f19 ${OUTDIR}/QC/cohort/${COHORT_ID}.SV.metrics | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )\n\
" >> ${OUTDIR}/${COHORT_ID}.run_summary.txt


# Also include:
# SV sizes (median & IQR per class)
# gene annotations (num genes LOF/dup/etc per sample, outliers)

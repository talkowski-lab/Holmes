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
nsamps_cohort=$( cat ${samples_list} | wc -l ) #number of samples in cohort
maxCount=$( echo "(( ${polyArt_filter}*${nsamps_cohort} ))" | bc ) #max for $polyArt_filter

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
awk '{ printf "%.0f\n", (($6+$3)/2)-(($2+$5)/2) }'  ${OUTDIR}/SV_calls/${COHORT_ID}.inversion.bedpe | sort -nk1,1 > ${WRKDIR}/sizes/inversion.size #take average of plus and minus bps for inversions
awk '{ print $3-$2 }' ${OUTDIR}/SV_calls/${COHORT_ID}.complex.bed | sort -nk1,1 > ${WRKDIR}/sizes/complex.size

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
  cut -f1 ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.clusters.txt | sort | uniq | sed '/^$/d' | wc -l | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.deletion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  while read ID bam sex; do
    cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed | wc -l
  done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo "      -Mean cnMOPS per Sample:" ) -
  echo " => Insertion"
  cut -f1 ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.clusters.txt | sort | uniq | sed '/^$/d' | wc -l | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  while read ID bam sex; do
    cat ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed | wc -l
  done < ${samples_list} | awk '{ sum+=$1 } END { print sum/NR }' | cut -f1 -d. | paste -d" " <( echo "      -Mean cnMOPS per Sample:" ) -
  echo " => Inversion"
  cut -f1 ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.clusters.txt | sort | uniq | sed '/^$/d' | wc -l | paste -d" " <( echo "      -Raw Clusters:" ) -
  fgrep Valid ${OUTDIR}/data/clusters/${COHORT_ID}.inversion.events.bedpe | wc -l | paste -d" " <( echo "      -Valid Clusters:" ) -
  echo " => Translocation"
  cut -f1 ${OUTDIR}/data/clusters/${COHORT_ID}.insertion.clusters.txt | sort | uniq | sed '/^$/d' | wc -l | paste -d" " <( echo "      -Raw Clusters:" ) -
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
\n\n\
+----------+\n\
| SV SIZES |\n\
+----------+" >> ${OUTDIR}/${COHORT_ID}.run_summary.txt
#Gather matrix of variants in cohort and per sample, then plot
fgrep -A1 "Deletion Variants [HQ]" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' > ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Tandem Duplication Variants [HQ]" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Insertion Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Simple Inversion Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Translocation Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Complex Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A1 "Incompletely Resolved Sites" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n1 | awk -v OFS="\t" '{ print $4, $5 }' | tr -d "(" | awk -v OFS="\t" '{ print $1-$2, $2 }' >> ${WRKDIR}/varCounts_cohort_toPlot.txt
fgrep -A3 "Deletion Variants [HQ]" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' > ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Tandem Duplication Variants [HQ]" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Insertion Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Simple Inversion Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Translocation Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Complex Variants" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
fgrep -A3 "Incompletely Resolved Sites" ${OUTDIR}/${COHORT_ID}.run_summary.txt | tail -n2 | awk -v OFS="\t" '{ print $NF }' | paste -s | awk -v OFS="\t" '{ print $2, $1-$2 }' >> ${WRKDIR}/varCounts_sample_toPlot.txt
Rscript ${liWGS_SV}/scripts/plotSVCounts.R ${WRKDIR}/varCounts_cohort_toPlot.txt ${WRKDIR}/varCounts_sample_toPlot.txt ${OUTDIR}/plots/${COHORT_ID}.SVcounts.pdf
#Gather sizes & plot
Rscript ${liWGS_SV}/scripts/scrapeSVsizes.R ${WRKDIR}/sizes/ ${OUTDIR}/${COHORT_ID}.run_summary.txt ${OUTDIR}/plots/${COHORT_ID}.SVsizes.pdf
#Gather frequencies & plot
mkdir ${WRKDIR}/frequencies
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed | awk -v nsamps=${nsamps_cohort} '{ print $7/nsamps }' > ${WRKDIR}/frequencies/deletion.list
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed | awk -v nsamps=${nsamps_cohort} '{ print $7/nsamps }' > ${WRKDIR}/frequencies/duplication.list
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.complex.bed | awk -v nsamps=${nsamps_cohort} '{ print $6/nsamps }' > ${WRKDIR}/frequencies/complex.list
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.unresolved.bed | awk -v nsamps=${nsamps_cohort} '{ print $6/nsamps }' > ${WRKDIR}/frequencies/unresolved.list
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe | awk -v nsamps=${nsamps_cohort} '{ print $10/nsamps }' > ${WRKDIR}/frequencies/insertion.list
sed '1d' ${OUTDIR}/SV_calls/${COHORT_ID}.inversion.bedpe | awk -v nsamps=${nsamps_cohort} '{ print $10/nsamps }' > ${WRKDIR}/frequencies/inversion.list
Rscript ${liWGS_SV}/scripts/plotSVFreq.R ${WRKDIR}/frequencies/ ${OUTDIR}/plots/${COHORT_ID}.SVfrequencies.pdf
#continue with report
awk -v max=${maxCount} '{ if ($7<=max) print $4 }' ${OUTDIR}/SV_calls/${COHORT_ID}.deletion.bed > ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($7<=max) print $4 }' ${OUTDIR}/SV_calls/${COHORT_ID}.duplication.bed >> ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($6<=max) print $4 }' ${OUTDIR}/SV_calls/${COHORT_ID}.complex.bed >> ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($6<=max) print $4 }' ${OUTDIR}/SV_calls/${COHORT_ID}.unresolved.bed >> ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($10<=max) print $7 }' ${OUTDIR}/SV_calls/${COHORT_ID}.insertion.bedpe >> ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($10<=max) print $7 }' ${OUTDIR}/SV_calls/${COHORT_ID}.inversion.bedpe >> ${WRKDIR}/nonArtifact.variantIDs.list
awk -v max=${maxCount} '{ if ($10<=max) print $7 }' ${OUTDIR}/SV_calls/${COHORT_ID}.translocation.bedpe >> ${WRKDIR}/nonArtifact.variantIDs.list
rm ${WRKDIR}/annotations/*per_sample.txt
while read ID bam sex; do
  for class in LOF GOF pLOF pD; do
    cat <( cat ${OUTDIR}/SV_calls/${COHORT_ID}*.bedpe | fgrep -w ${ID} | cut -f7 ) <( cat ${OUTDIR}/SV_calls/${COHORT_ID}*.bed | fgrep -w ${ID} | cut -f4 ) | fgrep -wf - ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w ${class} | wc -l >> ${WRKDIR}/annotations/${class}_variants_per_sample.txt
    cat <( cat ${OUTDIR}/SV_calls/${COHORT_ID}*.bedpe | fgrep -w ${ID} | cut -f7 ) <( cat ${OUTDIR}/SV_calls/${COHORT_ID}*.bed | fgrep -w ${ID} | cut -f4 ) | fgrep -wf - ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk -v class=${class} '{ if ($1==class) print $0 }' | wc -l >> ${WRKDIR}/annotations/${class}_genes_per_sample.txt  
    cat <( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/${COHORT_ID}*.bedpe | fgrep -w ${ID} | cut -f7 ) <( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/${COHORT_ID}*.bed | fgrep -w ${ID} | cut -f4 ) | fgrep -wf - ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w ${class} | wc -l >> ${WRKDIR}/annotations/${class}_variants_noArt_per_sample.txt
    cat <( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/${COHORT_ID}*.bedpe | fgrep -w ${ID} | cut -f7 ) <( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/${COHORT_ID}*.bed | fgrep -w ${ID} | cut -f4 ) | fgrep -wf - ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk -v class=${class} '{ if ($1==class) print $0 }' | wc -l >> ${WRKDIR}/annotations/${class}_genes_noArt_per_sample.txt  
  done
done < ${samples_list}
echo -e "\n\
+------------------+\n\
| GENE ANNOTATIONS |\n\
+------------------+\n\
**ALL VARIANTS**
 => Predicted LoF\n\
      -Variants: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w LOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/LOF_variants_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="LOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/LOF_genes_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Predicted Whole-Gene Duplication\n\
      -Variants: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w GOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/GOF_variants_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="GOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/GOF_genes_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Potential LoF\n\
      -Variants: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w pLOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pLOF_variants_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="pLOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pLOF_genes_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Potential Disruption\n\
      -Variants: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w pD | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pD_variants_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( cat ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="pD") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pD_genes_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)
**NON-ARTIFACTUAL VARIANTS (VAF ≤ ${polyArt_filter})**
 => Predicted LoF\n\
      -Variants: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w LOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/LOF_variants_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="LOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/LOF_genes_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Predicted Whole-Gene Duplication\n\
      -Variants: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w GOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/GOF_variants_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="GOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/GOF_genes_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Potential LoF\n\
      -Variants: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w pLOF | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pLOF_variants_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="pLOF") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pLOF_genes_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
 => Potential Disruption\n\
      -Variants: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed -e 's/_/\t/g' -e 's/,/\t/g'  -e '/^$/d' | cut -f$(seq 1 2 10000 | paste -s -d,) | fgrep -w pD | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pD_variants_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)\n\
      -Genes: $( fgrep -wf ${WRKDIR}/nonArtifact.variantIDs.list ${OUTDIR}/SV_calls/annotations/* | cut -f6 | sed 's/,/\n/g' | sed '/^$/d' | sort | uniq | sed 's/_/\t/g' | awk '{ if ($1=="pD") print $0 }' | wc -l ) ($( sort -nk1,1 ${WRKDIR}/annotations/pD_genes_noArt_per_sample.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ) per proband)" >> ${OUTDIR}/${COHORT_ID}.run_summary.txt



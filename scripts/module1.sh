#!/bin/bash

#liWGS-SV Pipeline: Module 1 (QC)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Set additional params
genome_size=2897310462 #n-masked length of grch37

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
  bsub -u nobody -o ${OUTDIR}/logs/sexCheck.log -e ${OUTDIR}/logs/sexCheck.log -sla miket_sc -q short -J ${COHORT_ID}_QC "${liWGS_SV}/scripts/sexCheck.sh ${ID} ${WRKDIR}/${ID}/${ID}.bam ${params}"
done < ${samples_list}

#Submit aneuploidy check
while read ID bam sex; do
  bsub -u nobody -o ${OUTDIR}/logs/aneuploidyCheck.log -e ${OUTDIR}/logs/aneuploidyCheck.log -sla miket_sc -q short -J ${COHORT_ID}_QC "${liWGS_SV}/scripts/chrCopyCount.sh ${ID} ${WRKDIR}/${ID}/${ID}.bam ${params}"
done < ${samples_list}

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep "${COHORT_ID}_QC" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep "${COHORT_ID}_QC" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Parse aneuploidy check
paste <( echo -e "chr\texpected" ) <( cut -f1 ${samples_list} | paste -s ) > ${OUTDIR}/QC/cohort/${COHORT_ID}.aneuploidyCheck.fractions.txt
paste <( echo -e "chr\texpected" ) <( cut -f1 ${samples_list} | paste -s ) > ${OUTDIR}/QC/cohort/${COHORT_ID}.aneuploidyCheck.copies.txt
echo -e "$( seq 1 22 )\nX\nY" > ${TMPDIR}/build.fractions.tmp
echo -e "$( seq 1 22 )\nX\nY" > ${TMPDIR}/build.copies.tmp
ID=$( head -n1 ${samples_list} | cut -f1 )
paste ${TMPDIR}/build.fractions.tmp <( fgrep -v "#" ${OUTDIR}/QC/sample/${ID}/${ID}.aneuploidyCheck | sed '/^$/d' | cut -f5 ) > ${TMPDIR}/build.fractions.tmp2
paste ${TMPDIR}/build.copies.tmp <( perl -E "say \"2\n\" x 24" | sed '/^$/d' ) > ${TMPDIR}/build.copies.tmp2
mv ${TMPDIR}/build.fractions.tmp2 ${TMPDIR}/build.fractions.tmp
mv ${TMPDIR}/build.copies.tmp2 ${TMPDIR}/build.copies.tmp
while read ID bam sex; do
  paste ${TMPDIR}/build.fractions.tmp <( fgrep -v "#" ${OUTDIR}/QC/sample/${ID}/${ID}.aneuploidyCheck | sed '/^$/d' | cut -f4 ) > ${TMPDIR}/build.fractions.tmp2
  paste ${TMPDIR}/build.copies.tmp <( fgrep -v "#" ${OUTDIR}/QC/sample/${ID}/${ID}.aneuploidyCheck | sed '/^$/d' | cut -f6 ) > ${TMPDIR}/build.copies.tmp2
  mv ${TMPDIR}/build.fractions.tmp2 ${TMPDIR}/build.fractions.tmp
  mv ${TMPDIR}/build.copies.tmp2 ${TMPDIR}/build.copies.tmp
done < ${samples_list}
cat ${TMPDIR}/build.fractions.tmp >> ${OUTDIR}/QC/cohort/${COHORT_ID}.aneuploidyCheck.fractions.txt
cat ${TMPDIR}/build.copies.tmp >> ${OUTDIR}/QC/cohort/${COHORT_ID}.aneuploidyCheck.copies.txt
rm ${TMPDIR}/build.fractions.tmp ${TMPDIR}/build.copies.tmp

#Plot aneuploidy check

#Run dosage Z-score
ID=$( head -n1 ${samples_list} | cut -f1 )
cat ${OUTDIR}/QC/sample/${ID}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed > ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed
while read ID bam sex; do
  paste ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed <( cut -f4 ${OUTDIR}/QC/sample/${ID}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed ) > ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed2
  mv ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed2 ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed
done < <( sed '1d' ${samples_list} )
cat <( echo -e "chr\tstart\tend\t$( cut -f1 ${samples_list} | paste -s )" ) ${WRKDIR}/${COHORT_ID}.WGSdosage_ObsVsExp.bed > ${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_ObsVsExp.bed
grep -ve "^[XY]" ${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_ObsVsExp.bed > ${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_ObsVsExp.noXY.bed
Rscript -e "x <- read.table(\"${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_ObsVsExp.noXY.bed\",header=T); x[,4:ncol(x)] <- t(abs(apply(x[,4:ncol(x)],1,scale))); write.table(x,\"${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_absoluteZscores.bed\",row.names=F,col.names=T,sep=\"\\t\",quote=F)"
#get median exp/obs ratio for first 45Mb on chr1
while read ID bam sex; do
  awk '{ if ($1==1 && $3<=45000000) print $4 }' ${OUTDIR}/QC/sample/${ID}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed | awk '{ sum+=$1 }END{ print sum/NR }'
done < ${samples_list} | sort -nk1,1 > ${WRKDIR}/chr1p_45Mb_dos.txt
med_dos=$( cat ${WRKDIR}/chr1p_45Mb_dos.txt | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )

#Collect summary for each sample
echo -e "##liWGS-SV PIPELINE COHORT QC\n##RUN DATE: $( echo $(date) | awk '{ print $1, $2, $3, $NF }' )\n\
#ID\tTotal_Pairs\tRead_Aln_Rate\tPair_Aln_Rate\tProper\tChimera\tRead_Dup\tPair_Dup\tMedian_Insert\tInsert_MAD\tHap_Phys_Cov\tHap_Nuc_Cov\tReported_Sex\tObserved_Sex\tDosage_ZScore" > ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics
while read ID bam sex; do
  total=$( grep '^FIRST_OF_PAIR' ${OUTDIR}/QC/sample/${ID}/${ID}.alignment_summary_metrics | awk '{ print $2 }' ) #total reads in BAM
  rd_aln=$( grep '^PAIR' ${OUTDIR}/QC/sample/${ID}/${ID}.alignment_summary_metrics | awk '{ print $6 }' ) #reads aligned
  rd_aln_rt=$( echo "scale=4;(( ${rd_aln}/(2*${total}) ))" | bc ) #read aln rate
  pr_aln=$( grep '^PAIR' ${OUTDIR}/QC/sample/${ID}/${ID}.alignment_summary_metrics | awk '{ print $17 }' ) #reads aligned in pairs
  pr_aln_rt=$( echo "scale=4;(( ${pr_aln}/(2*${total}) ))" | bc ) #pairwise aln rate
  prop=$( echo "scale=4;(( $( fgrep Proper ${OUTDIR}/QC/sample/${ID}/${ID}.stats | awk '{ print $3 }' | tr -d "()%" )/100 ))" | bc )
  rd_dup=$( echo "scale=4;(( $( grep '^Duplicates' ${OUTDIR}/QC/sample/${ID}/${ID}.stats | awk '{ print $3 }' | tr -d "()%")/100 ))" | bc ) #read dup rate
  pr_dup=$( grep -A1 '^LIBRARY' ${OUTDIR}/QC/sample/${ID}/${ID}.complexity | tail -n1 | awk '{ print $(NF-1) }' ) #pair dup rate
  chim=$( grep '^PAIR' ${OUTDIR}/QC/sample/${ID}/${ID}.alignment_summary_metrics | awk '{ print $21 }' ) #chimera pct
  mis=$( awk '{ if ($8=="FR") print $1 }' ${OUTDIR}/QC/sample/${ID}/${ID}.insert_size_metrics ) #MIS
  ismad=$( awk '{ if ($8=="FR") print $2 }' ${OUTDIR}/QC/sample/${ID}/${ID}.insert_size_metrics ) #ISMAD
  icov=$( echo "scale=2;(( ( ( ${pr_aln}/2 )*${prop}*( 1-${pr_dup} )*${mis} )/${genome_size} ))" | bc ) #calculate insert coverage, uses n-masked size of grch37 defined above
  ncov=$( grep -A1 '^GENOME_TERRITORY' ${OUTDIR}/QC/sample/${ID}/${ID}.wgs | tail -n1 | awk '{ print $2 }' ) #nucleotide cov
  osex=$( fgrep -w "PREDICTED SEX" ${OUTDIR}/QC/sample/${ID}/${ID}.sexCheck | awk '{ print $3 }' ) #predicted sex from sexcheck.sh
  index=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -w ${ID} | cut -f1 )+3" | bc )
  dsign=$( awk '{ if ($1==1 && $3<=45000000) print $4 }' ${OUTDIR}/QC/sample/${ID}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed | awk '{ sum+=$1 }END{ print sum/NR }' | awk -v med_dos=${med_dos} '{ if ($1>=med_dos) print ""; else print "-" }' )
  dosage=$( awk -v idx=${index} '{ print $idx }' ${OUTDIR}/QC/cohort/${COHORT_ID}.WGSdosage_absoluteZscores.bed | sed '1d' | fgrep -v NA | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  echo -e "${ID}\t${total}\t${rd_aln_rt}\t${pr_aln_rt}\t${prop}\t${chim}\t${rd_dup}\t${pr_dup}\t${mis}\t${ismad}\t${icov}\t${ncov}\t${sex}\t${osex}\t${dsign}${dosage}" #print metrics
done < ${samples_list} >> ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics

#Print warnings
while read ID total rd_aln_rt pr_aln_rt prop chim rd_dup pr_dup mis ismad icov ncov sex osex dosage; do
  if [ $( echo "${rd_aln_rt} < 0.8" | bc) -eq 1 ] || [ $( echo "${pr_aln_rt} < 0.8" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} alignment rate < 80%"
  fi
  if [ $( echo "${chim} >= 0.8" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} chimera rate ≥ 30%"
  fi
  if [ $( echo "${ismad} >= 700" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} median insert size is widely distributed (${ismad}bp)"
  fi
  if [ $( echo "${icov} < 30" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} haploid physical coverage < 30X"
  fi
  if [ $( echo "${dosage} > 1" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} WGS dosage Z-score > 1"
  fi
  if [ $( echo "${dosage} < -1" | bc) -eq 1 ]; then
    echo "WARNING [MODULE 1]: ${ID} WGS dosage Z-score < -1"
  fi
  if [ ${sex} != ${osex} ] && [ ${sex} != "U" ]; then
    echo "WARNING [MODULE 1]: ${ID} predicted sex does not match reported sex"
  fi
done < <( fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics ) > ${OUTDIR}/${COHORT_ID}_WARNINGS.txt


#!/bin/bash

#Calculate either nucleotide or physical coverage for a list of input bams

#Reads input
list=$1          #list of bams for analysis.  Note that this should be tab-delimmed with two columns: first col = full path, second col = sample ID
dict=$2          #ref dict
mode=$3          #either "nucleotide" or "physical"
binsize=$4       #bin size for genome segmentation
ID=$5            #group ID
OUTDIR=$6		     #output directory

# checks for appropriate input
if [ $# -eq 6 ]; then

  #loads modules
  module load BEDTools_2.17
  module load samtools-0.1.14

  #creates directory tree
  TMPDIR=`mktemp -d -p /scratch/miket/rlc47temp/tmp.files`
  echo -e "\nWriting temporary files to ${TMPDIR}"
  if [ ! -d "$OUTDIR" ]; then
    mkdir ${OUTDIR}
  fi
  if [ ! -d "${OUTDIR}/raw_coverages" ]; then
    mkdir ${OUTDIR}/raw_coverages
  fi
  if [ ! -d "${OUTDIR}/logfiles" ]; then
    mkdir ${OUTDIR}/logfiles
  fi
  echo -e "Writing output files to ${OUTDIR}"

  #creates contig list
  awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk 'BEGIN{OFS="\t"};{ print $2, $4 }' > ${TMPDIR}/contig.list
  num_contigs=$( cat ${TMPDIR}/contig.list | wc -l )

  #splits proper primary alignment pairs from each input bam
  if [ ${mode} == "physical" ]; then
    echo -e "\nEXTRACTING PROPER PAIRS...\n"
    while read bam sample; do
      mkdir ${TMPDIR}/${sample}
      echo "...submitting ${sample} in parallel..."
      while read contig_line; do
        contig=$( echo ${contig_line} | awk '{ print $1}' )
        bsub -sla miket_sc -u nobody -R 'rusage[mem=6000]' -M 6000 -v 10000 -q normal -e ${OUTDIR}/logfiles/error.txt -o ${OUTDIR}/logfiles/splitPropPairs.out -J ${ID}_step1 "/data/talkowski/tools/bin/sambamba_v0.4.6 view -f bam -F 'paired and proper_pair and not (unmapped or mate_is_unmapped or duplicate or secondary_alignment)' ${bam} ${contig} | /apps/lab/miket/samtools/1.0/bin/samtools sort -T ${TMPDIR}/${sample}.${contig}.screened -O bam -m 3200M /dev/stdin > ${TMPDIR}/${sample}/${sample}.${contig}.screened.bam"
      done < ${TMPDIR}/contig.list
    done < ${list}
  fi

  #creates segmentation bed
  echo -e "\nSEGMENTING REFERENCE...\n"
  awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk 'BEGIN{OFS="\t"};{ print $2, $4 }' > ${TMPDIR}/contig.list
  while read contig length; do
   echo "BINNING ${contig}"
   Rscript -e "options(scipen=1000); write.table(data.frame(as.character(\"${contig}\"),c(seq(0,${length},by=${binsize})),c(seq(${binsize},${length}+${binsize},by=${binsize}))),\"${TMPDIR}/${contig}.bins\",row.names=F,col.names=F,sep=\"\t\",quote=F)"
   cat ${TMPDIR}/${contig}.bins >> ${TMPDIR}/binned_genome.bed
  done < ${TMPDIR}/contig.list

  #get coverage for nucleotide mode
  if [ ${mode} == "nucleotide" ]; then
    echo -e "\nCALCULATING COVERAGE PER CONTIG...\n"
    while read contig length; do
      while read bam sample; do
        bsub -q normal -u nobody -J ${ID}_step1 "/data/talkowski/tools/bin/sambamba_v0.4.6 view -f bam -F 'not secondary_alignment and not duplicate' ${bam} ${contig} | bedtools coverage -counts -abam - -b ${TMPDIR}/${contig}.bins | sort -nk2,2 > ${TMPDIR}/${sample}.${contig}.coverage.bed"
      done < ${list}
    done < ${TMPDIR}/contig.list
  fi

  ##Gate (20 second check; 5 minute report)
  echo -e "\n"
  GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_step1" | wc -l)
  GATEwait=0
  until [[ $GATEcount == 0 ]]; do
   sleep 20s
   GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_step1" | wc -l)
   GATEwait=$[$GATEwait +1]
   if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} step 1 jobs to complete"
    GATEwait=0
   fi
  done

  #merges filtered alignments split by contig (if necessary), otherwise sorts and renames bams
  if [ ${mode} == "physical" ]; then
    if [ $num_contigs -gt 1 ]; then 
     echo -e "\nMERGING & SORTING FILTERED ALIGNMENTS...\n"
     while read bam sample; do
      bsub -sla miket_sc -u nobody -q big -R 'rusage[mem=16000]' -e ${OUTDIR}/logfiles/error.txt -o ${OUTDIR}/logfiles/mergeAlignments.out -M 16000 -v 24000 -J ${ID}_mergefiltered "/data/talkowski/tools/bin//data/talkowski/tools/bin/sambamba_v0.4.6 merge -l 6 ${TMPDIR}/${sample}/${sample}.merged.bam ${TMPDIR}/${sample}/${sample}.*.screened.bam; /data/talkowski/tools/bin//data/talkowski/tools/bin/sambamba_v0.4.6 sort -n --tmpdir=${TMPDIR} -m 12GB -l 6 -o ${TMPDIR}/${sample}.screened.nsort.bam ${TMPDIR}/${sample}/${sample}.merged.bam"
     done < ${list}
    elif [ $num_contigs -eq 1 ]; then
     echo -e "\nSORTING FILTERED ALIGNMENTS...\n"
     contig=$( awk '{ print $1 }' ${TMPDIR}/contig.list )
     while read bam sample; do
      bsub -sla miket_sc -u nobody -q big -R 'rusage[mem=16000]' -e ${OUTDIR}/logfiles/error.txt -o ${OUTDIR}/logfiles/mergeAlignments.out -M 16000 -v 24000 -J ${ID}_mergefiltered "/data/talkowski/tools/bin//data/talkowski/tools/bin/sambamba_v0.4.6 sort -m 4GB -n -l 6 -o ${TMPDIR}/${sample}.screened.nsort.bam ${TMPDIR}/${sample}/${sample}.${contig}.screened.bam"
     done < ${list}
    fi
    ##Gate (20 second check; 5 minute report)
    echo -e "\n"
    GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_mergefiltered" | wc -l)
    GATEwait=0
    until [[ $GATEcount == 0 ]]; do
     sleep 20s
     GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_mergefiltered" | wc -l)
     GATEwait=$[$GATEwait +1]
     if [[ $GATEwait == 15 ]]; then
      echo "$(date): INCOMPLETE"
      echo "$(date): Waiting on ${GATEcount} alignment mergers to complete"
      GATEwait=0
     fi
    done
  fi

  #creates fragments
  if [ ${mode} == "physical" ]; then
    echo -e "\nSYNTHESIZING FRAGMENTS...\n"
    while read bam sample; do
      bsub -sla miket_sc -e ${OUTDIR}/logfiles/error.txt -o ${OUTDIR}/logfiles/fragmentSynthesis.out -u nobody -q normal -J ${ID}_synthesizefragments "bedtools bamtobed -bedpe -i ${TMPDIR}/${sample}.screened.nsort.bam | awk -v OFS=\"\t\" '{ if (\$2>=0 && \$6>\$2) print \$1, \$2, \$6 }' > ${TMPDIR}/${sample}.fragments.bed"
    done < ${list}
    ##Gate (20 second check; 5 minute report)
    echo -e "\n"
    GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_synthesizefragments" | wc -l)
    GATEwait=0
    until [[ $GATEcount == 0 ]]; do
     sleep 20s
     GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_synthesizefragments" | wc -l)
     GATEwait=$[$GATEwait +1]
     if [[ $GATEwait == 15 ]]; then
      echo "$(date): INCOMPLETE"
      echo "$(date): Waiting on ${GATEcount} fragment synthesis jobs to complete"
      GATEwait=0
     fi
    done
  fi

  #calculates coverage
  if [ ${mode} == "physical" ]; then
    echo -e "\nCALCULATING COVERAGE...\n"
    while read bam sample; do
     bsub -sla miket_sc -u nobody -e ${OUTDIR}/logfiles/error.txt -o ${OUTDIR}/logfiles/calculateCoverage.out -q big  -R 'rusage[mem=40000]' -M 40000 -v 60000 -J ${ID}_getcoverage "bedtools coverage -counts -a ${TMPDIR}/${sample}.fragments.bed -b ${TMPDIR}/binned_genome.bed > ${TMPDIR}/${sample}.unsorted.coverage.bed; sort -V -k1,1 -k2,2 ${TMPDIR}/${sample}.unsorted.coverage.bed > ${TMPDIR}/${sample}.coverage.bed"
    done < ${list}

    ##Gate (20 second check; 5 minute report)
    echo -e "\n"
    GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_getcoverage" | wc -l)
    GATEwait=0
    until [[ $GATEcount == 0 ]]; do
     sleep 20s
     GATEcount=$(bjobs -w -u ${USER} | awk '{ print $7 }' | grep "${ID}_getcoverage" | wc -l)
     GATEwait=$[$GATEwait +1]
     if [[ $GATEwait == 15 ]]; then
      echo "$(date): INCOMPLETE"
      echo "$(date): Waiting on ${GATEcount} coverage calculations to complete"
      GATEwait=0
     fi
    done
  fi

  #Merges output
  echo -e "\nCONSOLIDATING OUTPUT...\n"
  if [ ${mode} == "nucleotide" ]; then
    while read bam sample; do
      while read contig length; do
        cat ${TMPDIR}/${sample}.${contig}.coverage.bed >> ${TMPDIR}/${sample}.coverage.bed
      done < ${TMPDIR}/contig.list
    done < ${list}
  fi
  sed 1d ${list} > ${TMPDIR}/appendedsamples.tmp
  first=$( head -n1 ${list} | awk '{ print $2 }' )
  echo -e "Chr\nStart\nStop\n${first}" > ${TMPDIR}/${ID}.results.header.txt
  cp ${TMPDIR}/${first}.coverage.bed ${TMPDIR}/${ID}.${mode}.cov_matrix.bed
  while read bam sample; do
    echo "...appending ${sample}..."
    awk '{ print $4 }' ${TMPDIR}/${sample}.coverage.bed > ${TMPDIR}/${sample}.cov_vals.txt
    paste ${TMPDIR}/${ID}.${mode}.cov_matrix.bed ${TMPDIR}/${sample}.cov_vals.txt > ${TMPDIR}/matrix.build.tmp
    mv ${TMPDIR}/matrix.build.tmp ${TMPDIR}/${ID}.${mode}.cov_matrix.bed
    echo "${sample}" >> ${TMPDIR}/${ID}.results.header.txt
  done < ${TMPDIR}/appendedsamples.tmp

  #Cleans and exits
  echo -e "\nCLEANING AND EXITING...\n"
  paste -s ${TMPDIR}/${ID}.results.header.txt > ${TMPDIR}/${ID}.results.header.reformatted.txt
  cat ${TMPDIR}/${ID}.results.header.reformatted.txt ${TMPDIR}/${ID}.${mode}.cov_matrix.bed > ${OUTDIR}/${ID}.${mode}.cov_matrix.bed
  while read bam sample; do
    mv ${TMPDIR}/${sample}.coverage.bed ${OUTDIR}/raw_coverages/${sample}.coverage.bed
  done < ${list}
  # rm -rf ${TMPDIR}

  echo -e "COMPLETE\nFinal output written to:"
  echo "${OUTDIR}/${ID}.${mode}.cov_matrix.bed"

else
  echo -e "\n\nBinwise Coverage Calculator\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
  echo "Usage:"
  echo "  binCov.sh [samples.list] [ref.fa.dict] [nucleotide/physical] [binsize] [group_ID] [OUTDIR]"
  echo ""
fi



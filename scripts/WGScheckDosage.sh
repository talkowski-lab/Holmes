#!/bin/bash

#liWGS-SV Pipeline: WGS dosage assessment script
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read Input
bam=$1     #full path to indexed bam file
ID=$2      #string, ID for library
REF=$3     #string, specify either "h37" or "hg19"
check=$4   #string, either "summary", "genome", or "both"
OUTDIR=$5  #full path to output directory

# checks for appropriate input
if [ $# -eq 5 ]; then

#Hard-code variable paths
if [ ${REF} == "h37" ]; then
  intervals=/data/talkowski/Samples/SFARI/miscFiles/dosage_intervals.bed
  dict=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.mainContigs.dict
else
  intervals=/data/talkowski/Samples/SFARI/miscFiles/dosage_intervals.hg19.bed
  dict=/data/talkowski/tools/ref/hg19_bwa07/hg19.lex.mainContigs.dict
fi

#Set params
h37=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.fa
binsize=500000 #500kb bins by default
TMPDIR=/scratch/miket/rlc47temp/tmp.files
module load bedtools/2.22.1
if [ ${check} == "genome" ] || [ ${check} == "both" ]; then
  bins=`mktemp`
fi

#Create output directory
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
if ! [ -e ${OUTDIR}/${ID}_WGSdosageCheck ]; then
  mkdir ${OUTDIR}/${ID}_WGSdosageCheck
fi

#Get library size
total=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} )

#different output for different specifications of "check"
if [ ${check} == "summary" ] || [ ${check} == "both" ]; then
  #Get counts
  echo -e "NA\tNA\tNA\tTOTAL_LIBRARY\t${total}" > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed
  sambamba view -f bam -F 'not secondary_alignment and not duplicate' ${bam} | bedtools coverage -counts -abam - -b <( cut -f1-4 ${intervals} ) | sort -nk1,1 -k2,2 >> ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed

  #Calculate fractions
  sed '1d' ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed | awk -v OFS="\t" -v total=${total} '{ print $1, $2, $3, $4, $5/total }' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.fractions.bed

  #Calculate summary fractions
  echo -e "IntervalClass\tFractionOfLibrary" > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.WGS_dosageSummary.txt
  for color in red blue yellow; do
    fgrep ${color} ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.counts.bed | awk -v OFS="\t" -v total=${total} -v color=${color} '{ sum+=$5 } END { print color, sum/total }' >> ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.WGS_dosageSummary.txt
  done
fi
if [ ${check} == "genome" ] || [ ${check} == "both" ]; then
  #Bin genome
  contig_dir=`mktemp -d`
  while read contig length; do
    Rscript -e "options(scipen=1000); write.table(data.frame(as.character(\"${contig}\"),c(seq(0,${length},by=${binsize})),c(seq(${binsize},${length}+${binsize},by=${binsize}))),\"${contig_dir}/${contig}.bins\",row.names=F,col.names=F,sep=\"\t\",quote=F)"
    cat ${contig_dir}/${contig}.bins >> ${bins}
  done < <( awk '{ print $2, $3 }' ${dict} | sed 's/\:/\t/g' | sed '1d' | awk 'BEGIN{OFS="\t"};{ print $2, $4 }' )

  #Get counts
  exp_per_bin=$( echo "${total}/$( cat ${bins} | wc -l )" | bc )
  if [ ${REF} == "h37" ]; then
    sambamba view -f bam -F 'not secondary_alignment and not duplicate' ${bam} | bedtools coverage -counts -abam - -b ${bins} | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed
  else
    sambamba view -f bam -F 'not secondary_alignment and not duplicate' ${bam} | bedtools coverage -counts -abam - -b ${bins} | sed 's/^chr//g' | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' | sed 's/^/chr/g' > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed
  fi

  #Calculate fractions
  awk -v OFS="\t" -v exp_per_bin=${exp_per_bin} '{ print $1, $2, $3, $4/exp_per_bin }' ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.counts.bed > ${OUTDIR}/${ID}_WGSdosageCheck/${ID}.genome.ObsVsExp.bed

  #Clean up
  rm ${bins}
  rm -r ${contig_dir}
fi

#If input inappropriate, displays usage and exits
else
 echo -e "\n\nChromatin Dosage Anomaly Diagnostic Test\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
 echo "Usage:"
 echo "  WGScheckDosage.sh [path/to/bam] [ID] [h37/hg19] [summary/genome] [OUTDIR]"
 echo ""
fi

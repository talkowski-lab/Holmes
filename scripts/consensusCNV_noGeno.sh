#!/bin/bash

#################################
#     Consensus CNV Pipeline    #
#          No Genotyping        #
#          Talkowski Lab        #
#################################

# Consensus Groups without Genotyping:
#     A [HIGH]: Valid cluster, cnMOPS support, <30% blacklist
#     B [MED]: cnMOPS call, ≥50kb, <30% blacklist, no clustering overlap
#     C [MED]: valid cluster, cnMOPS support, ≥30% blacklist
#     D [LOW]: cnMOPS call, ≥50kb, ≥30% blacklist
#     E [LOW]: valid cluster, <25kb, no cnMOPS support

#Input
events=$1               #full path to classifier output (events.bedpe format required)
cnMOPS=$2               #full path to cnMOPS file (bed format required)
DNAcopy_raw=$3          #full path to DNAcopy file (bed format required)
ID=$4                   #sample ID.  Must match sample ID used in classifier and cnMOPS input files
cnvtype=$5              #del or dup
params=$6               #params file from liWGS-SV

#Source params file
. ${params}

#Set params
TMPDIR=/scratch/miket/rlc47temp/tmp.files
module load bedtools/2.22.1
BL=`mktemp`
DNAcopy=`mktemp`
premerge=`mktemp`
fgrep -v GL ${CNV_BLACKLIST} > ${BL}
cnMOPS_m=`mktemp`
cnMOPS_cutoff=50000 #size below which cnMOPS calls without clustering support become low-qual
clustering_cutoff=25000 #size below which clusters without depth support become low-qual

#Pushes classifier output to bed
##NOTE: NO LONGER USES format_calls.py BECAUSE OF COORDINATE CHANGES IN CLUSTERING PATCH##
case ${cnvtype} in
  del)
    fgrep Valid ${events} | awk -v OFS="\t" -v ID=${ID} '$20 ~ ID { print $1, $3, $5, $7 }' | awk -v OFS="\t" '{ if ($3-$2<1) print $1, $3-1, $3, $4, $5; else print $0 }' | sed 's/[\t]+/\t/g' > ${TMPDIR}/${ID}_classifier.${cnvtype}.bed
    ;;
  dup)
    fgrep Valid ${events} | awk -v OFS="\t" -v ID=${ID} '$20 ~ ID { print $1, $2, $6, $7 }' | awk -v OFS="\t" '{ if ($3-$2<1) print $1, $3-1, $3, $4, $5; else print $0 }' | sed 's/[\t]+/\t/g' > ${TMPDIR}/${ID}_classifier.${cnvtype}.bed
    ;;
esac

#Removes DNAcopy intervals with >30% coverage by N-masked reference regions, also excludes X & Y
bedtools coverage -a ${NMASK} -b ${DNAcopy_raw} | awk -v OFS="\t" '{ if ($NF<0.3 && $1!="X" && $1!="Y") print $1, $2, $3, $4, $5, $6 }' > ${DNAcopy}

#Merges cnMOPS calls & DNAcopy calls, spans across assembly gaps but only if gap is <50% size of existant call
cat <( cut --complement -f6 ${cnMOPS} ) <( awk -v OFS="\t" '{ print $1, $2, $3, "DNAcopy_"$5, $6, "DNAcopy" }' ${DNAcopy} ) > ${premerge}
while read chr start end A B C; do
  bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b ${NMASK} | awk -v OFS="\t" -v end=${end} -v start=${start} '{ if (($6-$5)<((end-start)/2)) print $4, $5, $6, "GAP_EXTENSION", "GAP_EXTENSION", "GAP_EXTENSION" }' | cat <( echo -e "${chr}\t${start}\t${end}\t${A}\t${B}\t${C}" ) - | sort -nk1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o collapse >> ${cnMOPS_m}
done < ${premerge}

#GROUP A, HIGH - Valid cluster w/ cnMOPS support, <30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | bedtools intersect -v -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupA_noGeno", "HIGH" }' > ${TMPDIR}/${ID}.${cnvtype}.A.bed

#GROUP B, MEDIUM - cnMOPS calls larger than max size required for clustering overlap, <30% blacklist, no clustering overlap
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} -v OFS="\t" '{ if ($3-$2>=min) print $1, $2, $3, $6, "GroupB_noGeno", "MED" }' ${cnMOPS_m} ) -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -v -f 0.3 -a - -b ${BL} > ${TMPDIR}/${ID}.${cnvtype}.B.bed

#GROUP C, MEDIUM - Valid cluster w/ cnMOPS support, ≥30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | bedtools intersect -u -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupC_noGeno", "MED" }' > ${TMPDIR}/${ID}.${cnvtype}.C.bed

#GROUP D, LOW - cnMOPS calls larger than max size required for clustering overlap, ≥30% blacklist, no clustering overlap
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} -v OFS="\t" '{ if ($3-$2>=min) print $1, $2, $3, $6, "GroupD_noGeno", "LOW" }' ${cnMOPS_m} ) -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -u -f 0.3 -a - -b ${BL} > ${TMPDIR}/${ID}.${cnvtype}.D.bed

#GROUP E, LOW - Valid cluster w/ no cnMOPS support
bedtools intersect -v -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupE_noGeno", "LOW" }' > ${TMPDIR}/${ID}.${cnvtype}.E.bed

##Sort & write out
cat ${TMPDIR}/${ID}.${cnvtype}.A.bed ${TMPDIR}/${ID}.${cnvtype}.B.bed ${TMPDIR}/${ID}.${cnvtype}.C.bed ${TMPDIR}/${ID}.${cnvtype}.D.bed ${TMPDIR}/${ID}.${cnvtype}.E.bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${WRKDIR}/${ID}/${ID}.consensus.${cnvtype}.bed

##CLEAN UP
rm ${TMPDIR}/${ID}.${cnvtype}.*.bed ${BL} ${DNAcopy} ${cnMOPS_m} ${premerge}
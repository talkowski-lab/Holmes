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
DNAcopy=$3              #full path to DNAcopy file (bed format required)
ID=$4                   #sample ID.  Must match sample ID used in classifier and cnMOPS input files
cnvtype=$5              #del or dup
params=$6               #params file from liWGS-SV

#Source params file
. ${params}

#Set params
TMPDIR=/scratch/miket/rlc47temp/tmp.files
module load bedtools/2.22.1
BL=`mktemp`
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

#Merges cnMOPS calls & DNAcopy calls, spans across assembly gaps
cat <( cut --complement -f6 ${cnMOPS} ) <( awk -v OFS="\t" '{ print $1, $2, $3, "DNAcopy_"$5, $6, "DNAcopy" }' ${DNAcopy} ) <( awk -v OFS="\t" '{ print $0, "GAP_EXTENSION", "GAP_EXTENSION" }' /data/talkowski/rlc47/src/GRCh37_assemblyGaps.bed ) | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | bedtools merge -c 4,5,6 -o collapse -i - | sed -e 's/^23/X/g' -e 's/^24/Y/g'  | awk -v OFS="\t" -v ID=${ID} '$4 ~ /sampleName|DNAcopy/ { print $1, $2, $3+1, $4, $5, $6 }' > ${cnMOPS_m}

#GROUP A, HIGH - Valid cluster w/ cnMOPS support, <30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | bedtools intersect -v -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupA_NoGeno", "HIGH" }' > ${TMPDIR}/${ID}.${cnvtype}.A.bed

#GROUP B, MEDIUM - cnMOPS calls larger than max size required for clustering overlap, <30% blacklist, no clustering overlap
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} -v OFS="\t" '{ if ($3-$2>=min) print $1, $2, $3, $6, "GroupB_NoGeno", "MED" }' ${cnMOPS_m} ) -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -v -f 0.3 -a - -b ${BL} > ${TMPDIR}/${ID}.${cnvtype}.B.bed

#GROUP C, MEDIUM - Valid cluster w/ cnMOPS support, ≥30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | bedtools intersect -u -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupC_NoGeno", "MED" }' > ${TMPDIR}/${ID}.${cnvtype}.C.bed

#GROUP D, LOW - cnMOPS calls larger than max size required for clustering overlap, ≥30% blacklist, no clustering overlap
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} -v OFS="\t" '{ if ($3-$2>=min) print $1, $2, $3, $6, "GroupD_NoGeno", "LOW" }' ${cnMOPS_m} ) -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -u -f 0.3 -a - -b ${BL} > ${TMPDIR}/${ID}.${cnvtype}.D.bed

#GROUP E, LOW - Valid cluster w/ no cnMOPS support
bedtools intersect -v -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b ${cnMOPS_m} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupE_NoGeno", "LOW" }' > ${TMPDIR}/${ID}.${cnvtype}.E.bed

##Sort & write out
cat ${TMPDIR}/${ID}.${cnvtype}.A.bed ${TMPDIR}/${ID}.${cnvtype}.B.bed ${TMPDIR}/${ID}.${cnvtype}.C.bed ${TMPDIR}/${ID}.${cnvtype}.D.bed$ {TMPDIR}/${ID}.${cnvtype}.E.bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${WRKDIR}/${ID}/${ID}.consensus.${cnvtype}.bed

##CLEAN UP
rm ${TMPDIR}/${ID}.${cnvtype}.*.bed*
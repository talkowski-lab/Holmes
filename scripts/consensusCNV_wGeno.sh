#!/bin/bash

#################################
#     Consensus CNV Pipeline    #
#          w/ Genotyping        #
#          Talkowski Lab        #
#################################

# Consensus Groups without Genotyping:
#   A [HIGH]: Valid cluster, cnMOPS or genotyping support, <30% blacklist
#   B [HIGH]: cnMOPS call, ≥50kb, <30% blacklist, genotyping pass, no clustering overlap
#   C [MED]: cnMOPS call, <50kb, genotyping pass, <30% blacklist
#   D [MED]: valid cluster, genotyping or cnMOPS support, ≥30% blacklist
#   E [MED]: cnMOPS call, ≥50kb, genotyping pass, ≥30% blacklist
#   F [LOW]: cnMOPS call, ≥50kb, no clustering support, no genotyping support
#   G [LOW]: cnMOPS call, <50kb, genotyping pass, ≥30% blacklist
#   H [LOW]: valid cluster, <25kb, no cnMOPS or genotyping support

#Input
events=$1               #full path to classifier output (events.bedpe format required)
cnMOPS=$2               #full path to cnMOPS file (bed format required)
DNAcopy_raw=$3          #full path to DNAcopy file (bed format required)
genotypes=$4            #full path to sample genotype bed (cols 1-3 bed, col 4 = predicted genotype, col 5 = interval ID)
ID=$5                   #sample ID.  Must match sample ID used in classifier and cnMOPS input files
cnvtype=$6              #del or dup
params=$7               #params file from liWGS-SV

#Source params file
. ${params}

#Set params
TMPDIR=/scratch/miket/rlc47temp/tmp.files
module load bedtools/2.22.1
BL=`mktemp`
gcnv=`mktemp`
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

#Gets intervals predicted to be CN-variable in right direction by genotyping
if [ ${cnvtype}=="del" ]; then
  awk -v OFS="\t" '{ if ($4<2) print $0 }' ${genotypes} > ${gcnv}
else
  awk -v OFS="\t" '{ if ($4>2) print $0 }' ${genotypes} > ${gcnv}
fi

#Removes genotyped intervals with >30% coverage by N-masked reference regions, also excludes X & Y
bedtools coverage -a ${NMASK} -b ${gcnv} | awk -v OFS="\t" '{ if ($NF<0.3 && $1!="X" && $1!="Y") print $1, $2, $3, $4, $5 }' > ${gcnv}2
mv ${gcnv}2 ${gcnv}

#Removes DNAcopy intervals with >30% coverage by N-masked reference regions, also excludes X & Y
bedtools coverage -a ${NMASK} -b ${DNAcopy_raw} | awk -v OFS="\t" '{ if ($NF<0.3 && $1!="X" && $1!="Y") print $1, $2, $3, $4, $5, $6 }' > ${DNAcopy}

#Merges cnMOPS calls & DNAcopy calls, spans across assembly gaps but only if gap is <50% size of existant call
cat <( cut --complement -f6 ${cnMOPS} ) <( awk -v OFS="\t" '{ print $1, $2, $3, "DNAcopy_"$5, $6, "DNAcopy" }' ${DNAcopy} ) > ${premerge}
while read chr start end A B C; do
  bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b ${NMASK} | awk -v OFS="\t" -v end=${end} -v start=${start} '{ if (($6-$5)<((end-start)/2)) print $4, $5, $6, "GAP_EXTENSION", "GAP_EXTENSION", "GAP_EXTENSION" }' | cat <( echo -e "${chr}\t${start}\t${end}\t${A}\t${B}\t${C}" ) - | sort -nk1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o collapse >> ${cnMOPS_m}
done < ${premerge}

#GROUP A, HIGH - Valid cluster w/ cnMOPS OR genotyping support, <30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b <( cat ${cnMOPS_m} ${gcnv} | cut -f1-5 ) | bedtools intersect -v -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupA_wGeno", "HIGH" }' > ${TMPDIR}/${ID}.${cnvtype}.A.bed

#GROUP B, HIGH - cnMOPS ≥ 50kb, genotyping support, <30% blacklist, no clustering overlap
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} '{ if ($3-$2>=min) print $0 }' ${cnMOPS_m} ) -b ${gcnv} | bedtools intersect -v -f 0.51 -a - -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -v -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupB_wGeno", "HIGH" }' > ${TMPDIR}/${ID}.${cnvtype}.B.bed

#GROUP C, MED - cnMOPS < 50kb, genotyping support, <30% blacklist, no clustering overlap
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} '{ if ($3-$2<min) print $0 }' ${cnMOPS_m} ) -b ${gcnv} | bedtools intersect -v -f 0.51 -a - -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -v -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupC_wGeno", "MED" }' > ${TMPDIR}/${ID}.${cnvtype}.C.bed

#GROUP D, MED - Valid cluster w/ cnMOPS OR genotyping support, ≥30% blacklist
bedtools intersect -u -r -f 0.51 -a ${TMPDIR}/${ID}_classifier.${cnvtype}.bed -b <( cat ${cnMOPS_m} ${gcnv} | cut -f1-5 ) | bedtools intersect -u -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupD_wGeno", "MED" }' > ${TMPDIR}/${ID}.${cnvtype}.D.bed

#GROUP E, MED - cnMOPS ≥ 50kb, genotyping support, ≥30% blacklist, no clustering overlap
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} '{ if ($3-$2>=min) print $0 }' ${cnMOPS_m} ) -b ${gcnv} | bedtools intersect -v -f 0.51 -a - -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -u -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupE_wGeno", "MED" }' > ${TMPDIR}/${ID}.${cnvtype}.E.bed

#GROUP F, LOW - cnMOPS ≥ 50kb, no genotyping support, no clustering overlap
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} '{ if ($3-$2>=min) print $0 }' ${cnMOPS_m} ) -b ${gcnv} | bedtools intersect -v -f 0.51 -a - -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupF_wGeno", "LOW" }' > ${TMPDIR}/${ID}.${cnvtype}.F.bed

#GROUP G, LOW - cnMOPS < 50kb, genotyping support, ≥30% blacklist, no clustering overlap
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${cnMOPS_cutoff} '{ if ($3-$2<min) print $0 }' ${cnMOPS_m} ) -b ${gcnv} | bedtools intersect -v -f 0.51 -a - -b ${TMPDIR}/${ID}_classifier.${cnvtype}.bed | bedtools intersect -u -f 0.3 -a - -b ${BL} | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupG_wGeno", "LOW" }' > ${TMPDIR}/${ID}.${cnvtype}.G.bed

#GROUP H, LOW - Valid cluster, no cnMOPS support, no genotyping support, <25kb
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${clustering_cutoff} '{ if ($3-$2<min) print $0 }' ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ) -b <( cat ${cnMOPS_m} ${gcnv} | cut -f1-5 ) | awk -v OFS="\t" '{ print $1, $2, $3, $4, "GroupH_wGeno", "LOW" }' > ${TMPDIR}/${ID}.${cnvtype}.H.bed

##Sort & write out
cat ${TMPDIR}/${ID}.${cnvtype}.A.bed ${TMPDIR}/${ID}.${cnvtype}.B.bed ${TMPDIR}/${ID}.${cnvtype}.C.bed ${TMPDIR}/${ID}.${cnvtype}.D.bed ${TMPDIR}/${ID}.${cnvtype}.E.bed ${TMPDIR}/${ID}.${cnvtype}.F.bed ${TMPDIR}/${ID}.${cnvtype}.G.bed ${TMPDIR}/${ID}.${cnvtype}.H.bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${WRKDIR}/${ID}/${ID}.consensus.${cnvtype}.bed

##CLEAN UP
rm ${BL} ${gcnv} ${cnMOPS_m} ${DNAcopy} ${premerge} ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ${TMPDIR}/${ID}.${cnvtype}.*.bed
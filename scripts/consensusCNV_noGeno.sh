#!/bin/bash

#################################
#     Consensus CNV Pipeline    #
#          Talkowski Lab        #
#################################

#Input
events=$1               #full path to classifier output (events.bedpe format required)
cnMOPS=$2               #full path to cnMOPS file (bed-style required)
OUTDIR=$3               #output = ${OURDIR}/${ID}.consensus.${cnvtype}.bed
ID=$4                   #sample ID.  Must match sample ID used in classifier and cnMOPS input files
cnvtype=$5              #del or dup
minsize_cnMOPS=$6       #minimum size for reliable cnMOPS overlap
minsize_global=$7       #minimum size for any confident CNV call
maxsize_clustering=$8   #maximum size to use cnMOPS/clustering overlap before just allowing cnMOPS calls
minsize_NAHR=$9         #minimum size for a cnMOPS call to be considered for NAHR
allosomes=$10           #allow cnMOPS on allosomes? (TRUE or FALSE)

#Set params
TMPDIR=/scratch/miket/rlc47temp/tmp.files
module load bedtools/2.22.1
cnMOPS_m=`mktemp`

#Pushes classifier output to bed
##NOTE: NO LONGER USES format_calls.py BECAUSE OF COORDINATE CHANGES IN CLUSTERING PATCH##
case ${cnvtype} in
  del)
    fgrep Valid ${events} | awk -v OFS="\t" -v ID=${ID} '$20 ~ ID { print $0 }' > ${TMPDIR}/${ID}.${cnvtype}.class.bedpe
    awk -v OFS="\t" '{ print $1, $3, $5, $7 }' ${TMPDIR}/${ID}.${cnvtype}.class.bedpe > ${TMPDIR}/${ID}_classifier.${cnvtype}.bed
    #python /data/talkowski/rlc47/code/dna-scripts/format_calls.py -s ${cnvtype} -o ${TMPDIR} classifier ${TMPDIR}/${ID}.${cnvtype}.class.bedpe ${ID}
    ;;
  dup)
    fgrep Valid ${events} | awk -v OFS="\t" -v ID=${ID} '$20 ~ ID { print $0 }' > ${TMPDIR}/${ID}.${cnvtype}.class.bedpe
    awk -v OFS="\t" '{ print $1, $2, $6, $7 }' ${TMPDIR}/${ID}.${cnvtype}.class.bedpe > ${TMPDIR}/${ID}_classifier.${cnvtype}.bed
    #python /data/talkowski/rlc47/code/dna-scripts/format_calls.py -s ${cnvtype} -o ${TMPDIR} classifier ${TMPDIR}/${ID}.${cnvtype}.class.bedpe ${ID}
    ;;
esac

#Merges cnMOPS calls across assembly gaps & removes allosomes
if [ ${allosomes} == "FALSE" ]; then
    cat ${cnMOPS} <( awk -v OFS="\t" '{ print $0, "GAP_EXTENSION", "GAP_EXTENSION", "GAP_EXTENSION" }' /data/talkowski/rlc47/src/GRCh37_assemblyGaps.bed ) | awk -v OFS="\t" '{ if ($1!="X" && $1!="Y") print $0 }' | sort -nk1,1 -k2,2n | bedtools merge -c 4,5,6,7 -o collapse -i - | awk -v OFS="\t" -v ID=${ID} '$4 ~ ID { print $1, $2, $3+1, $4, $5, $6, $7 }' > ${cnMOPS_m}
else
    cat ${cnMOPS} <( awk -v OFS="\t" '{ print $0, "GAP_EXTENSION", "GAP_EXTENSION", "GAP_EXTENSION" }' /data/talkowski/rlc47/src/GRCh37_assemblyGaps.bed ) | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | bedtools merge -c 4,5,6,7 -o collapse -i - | sed -e 's/^23/X/g' -e 's/^24/Y/g'  | awk -v OFS="\t" -v ID=${ID} '$4 ~ ID { print $1, $2, $3+1, $4, $5, $6, $7 }' > ${cnMOPS_m}
fi

#GROUP A, HIGH - Valid cluster w/ depth support above global min
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${minsize_global} -v OFS="\t" '{ if ($3-$2>=min) print $0, "GroupA", "HIGH" }' ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ) -b ${cnMOPS_m} > ${TMPDIR}/${ID}.${cnvtype}.A.bed

#GROUP B, MEDIUM - cnMOPS calls larger than max size required for clustering overlap
awk -v min=${maxsize_clustering} -v OFS="\t" '{ if ($3-$2>=min) print $1, $2, $3, $7, "GroupB", "MED" }' ${cnMOPS_m} > ${TMPDIR}/${ID}.${cnvtype}.B.bed

#GROUP C, MEDIUM - Valid cluster w/ depth support below global min
bedtools intersect -u -r -f 0.51 -a <( awk -v min=${minsize_global} -v OFS="\t" '{ if ($3-$2<min && $3-$2>0) print $0, "GroupC", "MED" }' ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ) -b ${cnMOPS_m} > ${TMPDIR}/${ID}.${cnvtype}.C.bed

#GROUP D, LOW - Valid cluster above min cnMOPS threshold w/o depth support
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${minsize_cnMOPS} -v OFS="\t" '{ if ($3-$2>=min) print $0, "GroupD", "LOW" }' ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ) -b ${cnMOPS_m} > ${TMPDIR}/${ID}.${cnvtype}.D.bed

#GROUP E, LOW - Valid cluster below min cnMOPS threshold w/o depth support
bedtools intersect -v -r -f 0.51 -a <( awk -v min=${minsize_cnMOPS} -v OFS="\t" '{ if ($3-$2<min && $3-$2>0) print $0, "GroupE", "LOW" }' ${TMPDIR}/${ID}_classifier.${cnvtype}.bed ) -b ${cnMOPS_m} > ${TMPDIR}/${ID}.${cnvtype}.E.bed

#GROUP F, LOW - cnMOPS calls greater than minimum size for NAHR calling, but also below max size required for clustering overlap, with both ends of the interval within 15kb of a seg dup
bedtools pairtopair -type both -is -a <( awk -v min=${minsize_NAHR} -v max=${maxsize_clustering} -v OFS="\t" '{ if ($3-$2>=min && $3-$2<max) print $1, $2-7500, $2+7500, $1, $3-7500, $3+7500, $7 }' ${cnMOPS_m} ) -b <( awk -v OFS="\t" '{ print $1, $2, $3, $7, $8, $9 }' /data/talkowski/rlc47/src/genomicSuperDups.txt | sed 's/chr//g' | awk -v OFS="\t" '{ if ($1==$4) print $0 }' ) | awk '{ print $1"+"$2+7500"+"$5+7500"+"$7"+GroupF+LOW" }' | sort | uniq | sed 's/\+/\t/g' > ${TMPDIR}/${ID}.${cnvtype}.F.bed

#MERGE GROUPS THAT USED DEPTH SUPPORT (A, B, C, F)
bedtools merge -c 4,5,6 -o collapse -i <( cat ${TMPDIR}/${ID}.${cnvtype}.A.bed ${TMPDIR}/${ID}.${cnvtype}.B.bed ${TMPDIR}/${ID}.${cnvtype}.C.bed ${TMPDIR}/${ID}.${cnvtype}.F.bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) > ${TMPDIR}/${ID}.consensus.${cnvtype}.bed

#ADD GROUPS THAT DID NOT USE DEPTH SUPPORT (D, E)
cat ${TMPDIR}/${ID}.consensus.${cnvtype}.bed ${TMPDIR}/${ID}.${cnvtype}.D.bed ${TMPDIR}/${ID}.${cnvtype}.E.bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${TMPDIR}/${ID}.consensus.${cnvtype}.bed2
mv ${TMPDIR}/${ID}.consensus.${cnvtype}.bed2 ${TMPDIR}/${ID}.consensus.${cnvtype}.bed 

##SCREEN OUT ABPARTS OVERLAP > 50% FOR CALLS SMALLER THAN 2MB
cat <( bedtools coverage -a /data/talkowski/rlc47/src/abParts.bed -b <( awk -v OFS="\t" '{ if ($3-$2<=2000000) print $0 }' ${TMPDIR}/${ID}.consensus.${cnvtype}.bed ) | awk -v OFS="\t" '{ if ($NF<=0.5) print $1, $2, $3, $4, $5, $6 }' ) <( awk -v OFS="\t" '{ if ($3-$2>2000000) print $0 }' ${TMPDIR}/${ID}.consensus.${cnvtype}.bed ) | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -k2,2n | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${OUTDIR}/${ID}.consensus.${cnvtype}.bed

##CLEAN UP
rm ${TMPDIR}/${ID}*.bed*
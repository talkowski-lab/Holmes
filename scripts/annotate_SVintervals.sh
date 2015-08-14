#!/bin/bash

#Read input
input=$1 #input bed file, must be 5-column format: chr, start, end, type, event_ID
type=$2 #type: choose from DEL, DUP, INV, INS_SOURCE, and INS_SINK
OUTFILE=$3 #full path to output file
params=$4 #liWGS-SV params

#Loads liWGS-SV params
. ${params}

#Set other parameters
TMPDIR=`mktemp -d`
module load bedtools/2.22.1

#Process variants
case ${type} in
DEL)
  while read chr start end class eID; do
    if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
      rm ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
    fi
    ${liWGS_SV}/scripts/build_exons.sh hg19 "chr${chr}" ${start} ${end} ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
    if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
    genes=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( sed 's/^chr//g' ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ) | awk '{ print "LOF_"$NF }' | sort | uniq | awk -v ORS="," '{ print $0 }' )
   else
    genes=""
   fi
   echo -e "${chr}\t${start}\t${end}\t${class}\t${eID}\t${genes}"
  done < <( awk -v OFS="\t" '{ if ($4=="DEL") print $0 }' ${input} ) | fgrep -v "WARNING:" > ${OUTFILE}
  ;;

DUP)
  while read chr start end class eID; do
    GOF=$( awk -v OFS="\t" -v chr=${chr} -v start=${start} -v end=${end} '{ if ($1==chr && $3<=end && $2>=start) print "GOF_"$4 }' ${refFlat} | sort | uniq | awk -v ORS="," '{ print $0 }' )
    bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t$( echo "${start}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b1.tmp
    bedtools intersect -wb -a <( echo -e "${chr}\t${end}\t$( echo "${end}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b2.tmp
    fgrep -wf ${TMPDIR}/b1.tmp ${TMPDIR}/b2.tmp > ${TMPDIR}/ovr.tmp
    if [ $( cat ${TMPDIR}/ovr.tmp | wc -l ) -gt 0 ]; then
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        rm ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      fi
      ${liWGS_SV}/scripts/build_exons.sh hg19 "chr${chr}" ${start} ${end} ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        LOF=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( sed 's/^chr//g' ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ) | awk '{ print $NF }' | sort | uniq | fgrep -wf ${TMPDIR}/ovr.tmp | awk -v ORS="," '{ print "LOF_"$1 }' )
      else
        LOF=""
      fi
    else
      LOF=""
    fi
    genes=$( echo -e "${GOF}${LOF}" )
   echo -e "${chr}\t${start}\t${end}\t${class}\t${eID}\t${genes}"
  done < <( awk -v OFS="\t" '{ if ($4=="DUP") print $0 }' ${input} ) | fgrep -v "WARNING:" >> ${OUTFILE}
  ;;

INS_SINK)
  while read chr start end class eID; do
    genes=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b ${refFlat} | cut -f7 | sort | uniq | awk -v ORS="," '{ print "pD_"$0 }' )
    echo -e "${chr}\t${start}\t${end}\t${class}\t${eID}\t${genes}"
  done < <( awk -v OFS="\t" '{ if ($4=="INS_SINK") print $0 }' ${input} ) >> ${OUTFILE}
  ;;

INS_SOURCE)
  while read chr start end class eID; do
    bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t$( echo "${start}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b1.tmp
    bedtools intersect -wb -a <( echo -e "${chr}\t${end}\t$( echo "${end}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b2.tmp
    LOFa=$( echo "$( fgrep -wvf ${TMPDIR}/b1.tmp ${TMPDIR}/b2.tmp | awk -v ORS="," '{ print "pLOF_"$1 }' )$( fgrep -wvf ${TMPDIR}/b2.tmp ${TMPDIR}/b1.tmp | awk -v ORS="," '{ print "pLOF_"$1 }' )" )
    fgrep -wf ${TMPDIR}/b1.tmp ${TMPDIR}/b2.tmp > ${TMPDIR}/ovr.tmp
    if [ $( cat ${TMPDIR}/ovr.tmp | wc -l ) -gt 0 ]; then
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        rm ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      fi
      ${liWGS_SV}/scripts/build_exons.sh hg19 "chr${chr}" ${start} ${end} ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        LOFb=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( sed 's/^chr//g' ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ) | awk '{ print $NF }' | sort | uniq | fgrep -wf ${TMPDIR}/ovr.tmp | awk -v ORS="," '{ print "pLOF_"$1 }' )
      else
        LOFb=""
      fi
    else
      LOFb=""
    fi
    genes=$( echo -e "${LOFa}${LOFb}" )
   echo -e "${chr}\t${start}\t${end}\t${class}\t${eID}\t${genes}"
  done < <( awk -v OFS="\t" '{ if ($4=="INS_SOURCE") print $0 }' ${input} ) | fgrep -v "WARNING:" >> ${OUTFILE}
  ;;

INV)
  while read chr start end class eID; do
    bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t$( echo "${start}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b1.tmp
    bedtools intersect -wb -a <( echo -e "${chr}\t${end}\t$( echo "${end}+1" | bc )" ) -b ${refFlat} | cut -f7 | sort | uniq > ${TMPDIR}/b2.tmp
    LOFa=$( echo "$( fgrep -wvf ${TMPDIR}/b1.tmp ${TMPDIR}/b2.tmp | awk -v ORS="," '{ print "LOF_"$1 }' )$( fgrep -wvf ${TMPDIR}/b2.tmp ${TMPDIR}/b1.tmp | awk -v ORS="," '{ print "LOF_"$1 }' )" )
    fgrep -wf ${TMPDIR}/b1.tmp ${TMPDIR}/b2.tmp > ${TMPDIR}/ovr.tmp
    if [ $( cat ${TMPDIR}/ovr.tmp | wc -l ) -gt 0 ]; then
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        rm ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      fi
      ${liWGS_SV}/scripts/build_exons.sh hg19 "chr${chr}" ${start} ${end} ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp
      if [ -e ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ]; then
        LOFb=$( bedtools intersect -wb -a <( echo -e "${chr}\t${start}\t${end}" ) -b <( sed 's/^chr//g' ${TMPDIR}/${ID}.${chr}_${start}_${end}.buildExons.tmp ) | awk '{ print $NF }' | sort | uniq | fgrep -wf ${TMPDIR}/ovr.tmp | awk -v ORS="," '{ print "LOF_"$1 }' )
      else
        LOFb=""
      fi
    else
      LOFb=""
    fi
    genes=$( echo -e "${LOFa}${LOFb}" )
   echo -e "${chr}\t${start}\t${end}\t${class}\t${eID}\t${genes}"
  done < <( awk -v OFS="\t" '{ if ($4=="INV") print $0 }' ${input} ) | fgrep -v "WARNING:" >> ${OUTFILE}
  ;;

esac

#Cleanup
rm -rf ${TMPDIR}


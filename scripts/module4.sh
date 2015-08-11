#!/bin/bash

#liWGS-SV Pipeline: Module 4 (Physical Depth CNV Calling)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}
cd ${WRKDIR}

#Load necessary modules for SE large CNV caller
module load intel_parallel_xe/xe
module load R/3.1.0
Rscript -e "if(\"DNAcopy\" %in% rownames(installed.packages()) == FALSE){source(\"http://bioconductor.org/biocLite.R\"); biocLite(\"DNAcopy\")}; suppressPackageStartupMessages(library(DNAcopy))"

#Submit cnMOPS - autosomes
mkdir ${WRKDIR}/cnMOPS
cov=${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed
for contig in $( seq 1 22 ); do
  mkdir ${WRKDIR}/cnMOPS/${contig}
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed
  for binsize in 1 3 10 30; do
    bsub -q big -M 30000 -sla miket_sc -u nobody -o ${OUTDIR}/logs/cnMOPS.log -e ${OUTDIR}/logs/cnMOPS.log -R 'rusage[mem=30000]' -v 40000 -J ${COHORT_ID}_cnMOPS "Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID} ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed ${WRKDIR}/cnMOPS/${contig}/"
  done
done

#Submit cnMOPS - allosomes
if [ ${other_assign} == "MALE" ]; then
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="M" || $14=="O") print $1 }' > ${WRKDIR}/males.list
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="F") print $1 }' > ${WRKDIR}/females.list
else
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="M") print $1 }' > ${WRKDIR}/males.list
  fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | awk '{ if ($14=="F" || $14=="O") print $1 }' > ${WRKDIR}/females.list
fi
Midx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/males.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )
Fidx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/females.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )
for contig in X Y; do
  mkdir ${WRKDIR}/cnMOPS/${contig}
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) | cut -f${Midx} > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.M.bed
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) | cut -f${Fidx} > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.F.bed
  for binsize in 1 3 10 30; do
    bsub -u nobody -q big -M 30000 -o ${OUTDIR}/logs/cnMOPS.log -e ${OUTDIR}/logs/cnMOPS.log -sla miket_sc -u rlc47 -R 'rusage[mem=30000]' -v 40000 -J ${COHORT_ID}_cnMOPS "Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID}_M ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.M.bed ${WRKDIR}/cnMOPS/${contig}/"
    bsub -u nobody q big -M 30000 -o ${OUTDIR}/logs/cnMOPS.log -e ${OUTDIR}/logs/cnMOPS.log -sla miket_sc -u rlc47 -R 'rusage[mem=30000]' -v 40000 -J ${COHORT_ID}_cnMOPS "Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID}_F ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.F.bed ${WRKDIR}/cnMOPS/${contig}/"
  done
done  

#Submit SE large CNV caller
Rscript ${liWGS_SV}/scripts/SE_largeCNV/getwithinlibrarynorm_query.R ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ${COHORT_ID} ${WRKDIR}
while read ID bam sex; do
  bsub -q normal -M 6000 -R 'rusage[mem=6000]' -v 10000 -sla miket_sc -u nobody -o ${OUTDIR}/logs/DNAcopy.log -e ${OUTDIR}/logs/DNAcopy.log -J ${COHORT_ID}_DNAcopy "${liWGS_SV}/scripts/SE_largeCNV/largescaleCNVpipeline.sh ${samples_list} ${params} ${ID}"
done < ${samples_list}

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cnMOPS" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_cnMOPS" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Merge cnMOPS calls
for contig in $( seq 1 22 ); do
  echo ${contig}
  for binsize in 1 3 10 30; do
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
  done
done
for contig in X Y; do
  echo ${contig}
  for binsize in 1 3 10 30; do
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_M.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
    fgrep -v "#" ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}_F.${binsize}000bpBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff
  done
done

#Split cnMOPS Calls by sample
mkdir ${WRKDIR}/cnMOPS/cnMOPS_calls
while read bam ID; do
  echo ${ID}
  mkdir ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}
  fgrep -w "${ID}" ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff | grep 'CN[0-1]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed
  fgrep -w "${ID}" ${WRKDIR}/cnMOPS/${COHORT_ID}.cnMOPS_master.cnMOPS.gff | grep 'CN[3-9]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  ${WRKDIR}/${ID}/${ID}.cnMOPS.dels.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' ${WRKDIR}/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  ${WRKDIR}/${ID}/${ID}.cnMOPS.dups.bed
done < ${WRKDIR}/cnMOPS.input

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_DNAcopy" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_DNAcopy" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Move remaining DNAcopy samples to sample dirs
while read ID bam sex; do
  mv ${WRKDIR}/${ID}.SD.profile ${WRKDIR}/${ID}/DNAcopy/
done < ${samples_list}




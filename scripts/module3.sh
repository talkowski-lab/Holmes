#!/bin/bash

#liWGS-SV Pipeline: Module 3 (Physical Depth Analysis)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Run BinCov on samples
while read ID bam sex; do
  echo -e "${WRKDIR}/${ID}/${ID}.bam\t${ID}" >> ${WRKDIR}/cnMOPS.input
done < ${samples_list}
bsub -q big -R 'rusage[mem=20000]' -sla miket_sc -M 20000 -v 30000 -o ${OUTDIR}/logs/bincov.log -e ${OUTDIR}/logs/bincov.log -J ${COHORT_ID}_binCov "${liWGS_SV}/scripts/binCov.sh ${WRKDIR}/cnMOPS.input ${DICT} physical 1000 ${COHORT_ID} ${WRKDIR}/iCov"

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_binCov" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_binCov" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo "$(date): INCOMPLETE"
    echo "$(date): Waiting on ${GATEcount} jobs to complete"
    GATEwait=0
  fi
done

#Submit cnMOPS - autosomes
mkdir ${WRKDIR}/cnMOPS
cov=${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed
for contig in $( seq 1 22 ); do
  mkdir ${WRKDIR}/cnMOPS/${contig}
  cat <( head -n1 ${cov} ) <( awk -v chr=${contig} '{ if ($1==chr) print }' ${cov} ) > ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed
  for binsize in 1 3 10 30; do
    bsub -q big -M 30000 -sla miket_sc -u rlc47 -R 'rusage[mem=30000]' -v 40000 -J ${COHORT_ID}_cnMOPS "Rscript ${liWGS_SV}/scripts/cnMOPS_postcoverage.R -m insert -r ${binsize} -b ${binsize}000 -I ${COHORT_ID} ${WRKDIR}/cnMOPS/${contig}/${COHORT_ID}.rawCov.chr${contig}.bed ${WRKDIR}/cnMOPS/${contig}/"
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
for contig in X Y; do
  
done

for contig in $( seq 1 22 ); do
  echo ${contig}
  for binsize in 1 3 10 30; do
    fgrep -v "#" /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/${contig}/PE250.${binsize}00kbBins.cnMOPS.gff | awk -v OFS="\t" '{ print $1, $4, $5, $9, $10, $11, $12 }' | sed 's/^chr//g' >> /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/PE250.cnMOPS_master.cnMOPS.gff
  done
done
mkdir /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls
while read bam ID; do
  echo ${ID}
  mkdir /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}
  fgrep -w "${ID}" /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/PE250.cnMOPS_master.cnMOPS.gff | grep 'CN[0-1]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed
  fgrep -w "${ID}" /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/PE250.cnMOPS_master.cnMOPS.gff | grep 'CN[3-9]' | sed 's/median\=//g' | sed 's/mean\=//g' | sed 's/CN\=//g' | sed 's/\;//g' > /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dels.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.dels.bed
  bedtools merge -d 1 -c 4,5,6,7 -o distinct,mean,mean,distinct -i <( sed -e 's/^X/23/g' -e 's/^Y/24/g' /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.preMerge.dups.bed | sort -nk1,1 -k2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' ) >  /data/talkowski/Collaboration/JaffeAssembly/cnMOPS/cnMOPS_calls/${ID}/${ID}.cnMOPS.dups.bed
done < /data/talkowski/Collaboration/JaffeAssembly/bams.list2

#Move, bgzip, and tabix individual coverage
while read ID bam sex; do
  mv ${WRKDIR}/iCov/raw_coverages/${ID}.coverage.bed ${WRKDIR}/${ID}/
  bgzip ${WRKDIR}/${ID}/${ID}.coverage.bed
  tabix -p bed ${WRKDIR}/${ID}/${ID}.coverage.bed.gz
done < ${samples_list}






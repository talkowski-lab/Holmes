#!/bin/bash

#liWGS-SV Pipeline: Module 8 (Call reformatting & variant consolidation)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Make temp files
pool=`mktemp`
clusters=`mktemp`
agg=`mktemp`

#Make output directories
mkdir ${WRKDIR}/final_variants
mkdir ${WRKDIR}/raw_clusters

#Move clusters, events, and breakpoints to clusters file
for type in deletion insertion inversion transloc; do
  cp ${WRKDIR}/classifier/clusterfix/newCoords/${type}.events.reclassified.bedpe ${WRKDIR}/raw_clusters/${COHORT_ID}.${type}.events.bedpe
  cp ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${type}.patched.clusters ${WRKDIR}/raw_clusters/${COHORT_ID}.${type}.clusters.txt
  cp ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_${type}.patched.bkpts.bedpe ${WRKDIR}/raw_clusters/${COHORT_ID}.${type}.bkpts.bedpe
done

#Write out consensus CNVs
for type in deletion duplication; do
  echo -e "#chr\tstart\tend\tID\tmergedist_5\tmergedist_3\tobservations\tsamples\tinfo" > ${WRKDIR}/final_variants/${COHORT_ID}.${type}.bed
done
cat ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_dels.merged.bed > ${WRKDIR}/final_variants/${COHORT_ID}.deletion.bed
cat ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_dups.merged.bed > ${WRKDIR}/final_variants/${COHORT_ID}.duplication.bed

#Write out simple inversions
echo -e "#chr\tminplus\tmaxplus\tchr\tminminus\tmaxminus\tID\tstA\tstB\tobservations\tsamples\tclusters" > ${WRKDIR}/final_variants/${COHORT_ID}.inversion.bedpe
fgrep -w SIMPLE ${WRKDIR}/classifier/clusterfix/newCoords/inversion.classifications.list | fgrep -w Valid > ${pool}
while read clusterlist chr minp maxp minm maxm trash; do
  if [ $( echo ${clusterlist} | sed 's/,/\n/g' | wc -l ) -gt 1 ]; then
    echo ${clusterlist} | sed 's/,/\n/g' > ${clusters}
    samples=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | paste -s -d, )
    observations=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | wc -l )
  else
    samples=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f20 | tr -d "[]" )
    observations=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f19 )
  fi
  echo -e "${chr}\t${minp}\t${maxp}\t${chr}\t${minm}\t${maxm}\t${COHORT_ID}_simpleinversion\t+\t-\t${observations}\t[${samples}]\t[${clusterlist}]"
done < ${pool} | awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, $6, $7"_"NR, $8, $9, $10, $11, $12 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.inversion.bedpe

#Write out chromosomal translocations
echo -e "#chr\tminplus\tmaxplus\tchr\tminminus\tmaxminus\tID\tarmA\tarmB\tobservations\tsamples\tclusters" > ${WRKDIR}/final_variants/${COHORT_ID}.translocation.bedpe
fgrep RecipCTX ${WRKDIR}/classifier/clusterfix/newCoords/translocation.classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cluster type; do
  observations=$( fgrep -w ${cluster} ${WRKDIR}/final_variants/${COHORT_ID}.transloc.bedpe | cut -f19 )
  samples=$( fgrep -w ${cluster} ${WRKDIR}/final_variants/${COHORT_ID}.transloc.bedpe | cut -f20 | tr -d "[]" )
  armA=${type:0:1}
  armB=${type:1:1}
echo -e "${chrA}\t${minA}\t${maxA}\t${chrB}\t${minB}\t${maxB}\t${COHORT_ID}_CTX\t${armA}\t${armB}\t${observations}\t[${samples}]\t[${cluster}]"
done < ${pool} | awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, $6, $7"_"NR, $8, $9, $10, $11, $12 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.translocation.bedpe

#Write out insertions
echo -e "#source_chr\tsource_start\tsource_end\tsink_chr\tsink_min\tsink_max\tID\tsource_orientation\tsink_orientation\tobservations\tsamples\tclusters" > ${WRKDIR}/final_variants/${COHORT_ID}.insertion.bedpe
fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/inversion.classifications.list | fgrep INS5 > ${pool}
while read cluster chr minp maxp minm maxm trash; do
  samples=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f20 | tr -d "[]" )
  observations=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f19 )
  echo -e "${chr}\t${maxm}\t${maxp}\t${chr}\t${minp}\t${minm}\t${COHORT_ID}_insertion\t-\t+\t${observations}\t[${samples}]\t[${cluster}]"
done < ${pool} > ${agg}
fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/inversion.classifications.list | fgrep INS3 > ${pool}
while read cluster chr minp maxp minm maxm trash; do
  samples=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f20 | tr -d "[]" )
  observations=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.inversion.events.bedpe | cut -f19 )
  echo -e "${chr}\t${minm}\t${minp}\t${chr}\t${maxp}\t${maxm}\t${COHORT_ID}_insertion\t-\t+\t${observations}\t[${samples}]\t[${cluster}]"
done < ${pool} >> ${agg} 
##ADD MORE SOURCES OF INSERTION HERE (TRANSLOC & COMPLEX LINKED)
sort -nk1,1 -k2,2n ${agg} | awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, $6, $7"_"NR, $8, $9, $10, $11, $12 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.insertion.bedpe







#Clean up
rm ${pool}
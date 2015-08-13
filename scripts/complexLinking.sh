#!/bin/bash

#liWGS-SV Pipeline: Complex Linking Script
#August 2015
#Contact: hbrand@mgh.harvard.edu or rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Load modules
module load anaconda/2.0.1
module load bedtools/2.22.1

#Set other parameters
min_sample_concordance=0.2 #minimum overlap between clusters to consider potential complex site

#Determine clustering distance
while read ID bam sex; do
  basename ${WRKDIR}/${ID}/bamstat/deletion_clusters*txt | cut -f3 -d_ | tr -d "d"
done < ${samples_list} | sort -nk1,1 > ${WRKDIR}/clustering_distances.list
clustdist=$( tail -n1 ${WRKDIR}/clustering_distances.list )

#Create output directory
mkdir ${WRKDIR}/classifier/clusterfix/newCoords/Complex

#Get valid IDs & breakpoints
cat ${WRKDIR}/classifier/clusterfix/newCoords/*.reclassified.bedpe | fgrep Valid | awk '{print $7}' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/valid.ids.txt
sed 's/\.0[1-4]//g' ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_*.bkpts.bedpe | fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/valid.ids.txt | awk -v OFS="\t" '{print $1, $2, $3, $7, $20"\n"$4, $5, $6, $7, $20}' | sort -k1,1 -k2,2n > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/All.SV.valid.bed

#Cluster valid breakpoints
bedtools cluster -d ${clustdist} -i ${WRKDIR}/classifier/clusterfix/newCoords/Complex/All.SV.valid.bed > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/All.SV.valid.cluster.txt

#Pull out only clusters with multiple events
uniq -D -f 5 ${WRKDIR}/classifier/clusterfix/newCoords/Complex/All.SV.valid.cluster.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.valid.clusters.txt
cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.valid.clusters.txt | awk '{print $6}' | sort -nu > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.ids.txt

#Iterate and do stuff (?)
i=0
if [ -e ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int.list.txt ]; then
  rm ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int.list.txt
fi
while read p; do
  clusterid=$( echo ${p} )
  awk -v var=${clusterid} '{if ($NF==var) print $0}' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.valid.clusters.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt
  countint=1
  while [  ${countint} -gt 0 ]; do
    count1=1
    count=0
    awk 'NR==1 { print $4 }' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt 
    while [ ${count} -lt ${count1} ]; do
      fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | awk '{ print $5 }' | tr -d "[]" | tr ',' '\n' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.ids.txt  
      fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.ids.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | awk '{ print $4 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt
      fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | awk '{ print $5 }' | tr -d "[]" | tr ',' '\n' | fgrep -wf - ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | awk '{ print $4 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt
      count=$( sort -u ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt | wc -l )
      count1=$( sort -u ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt | wc -l )
    done
    sort -u ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt | tr '\n' '\t' | awk '{ print $0 }' >> ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int.list.txt 
    fgrep -wvf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int1.txt
    cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int1.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt
    countint=$( cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | wc -l )
  done
done < ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.ids.txt

#Step 2
awk '{ if ($2!="") print $0 }' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int.list.txt | awk -v cohort=${COHORT_ID} '{ print cohort"_complex_"NR "\t" $0 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt
cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt
rm ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int2.list.txt
countint=1
while [ ${countint} -gt 0 ]; do
  awk 'NR==1 { print $1 }' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt
  count1=1
  count=0
  while [  ${count} -lt ${count1} ]; do
    fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt | cut -f2- | tr '\t' '\n' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.ids.txt  
    fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.ids.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt | awk '{ print $1 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt
    fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt | cut -f2- | tr '\t' '\n' | fgrep -wf - ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt | awk '{ print $1 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt
    count=$( sort -u ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt | wc -l )
    count1=$( sort -u ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt | wc -l )
  done
  fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist1.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int1.list.txt | cut -f2- | tr '\t' '\n' | sort -u | tr '\n' '\t' | awk '{ print $0 }' >> ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int2.list.txt
  fgrep -wvf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterlist.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int1.txt
  cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int1.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt
  countint=$( cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/int.txt | wc -l )
done
awk -v cohort=${COHORT_ID} '{print cohort"_complex_"NR $0}' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complex.int2.list.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ComplexEvents.txt 

#Percent for complex events
if [ -e ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complexEventsSampleOverlap.txt ]; then
  rm ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complexEventsSampleOverlap.txt
fi
while read line; do 
  id=$( echo ${line} | awk '{ print $1 }' )
  echo -e "${line}" | cut -f2- | tr '\t' '\n' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt
  #Look for highest shared cluster percent for each individual cluster
  while read cluster; do 
    fgrep -w ${cluster} ${WRKDIR}/classifier/clusterfix/newCoords/*.events.reclassified.bedpe | awk '{ print $20 }' | tr -d "[]" | tr ',' '\n' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid.txt
    while read cluster1; do
      if [ "${cluster}" != "${cluster1}" ]; then
        fgrep -w ${cluster1} ${WRKDIR}/classifier/clusterfix/newCoords/*.events.reclassified.bedpe | awk '{ print $20 }' | tr -d "[]" | tr ',' '\n' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid1.txt
        clustercount=$( cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid1.txt | sort | uniq -c | awk '{ if ($1>1) print }' | wc -l )
        clustertotal=$( cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid.txt ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterid1.txt | sort | uniq -c | wc -l )
        echo "scale=3; ${clustercount}/${clustertotal}" | bc -l >> ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterpercent.txt
      fi
    done < ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt
  done < ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt
  clustermax=$( sort -nrk1,1 ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterpercent.txt | head -n1 )
  rm ${WRKDIR}/classifier/clusterfix/newCoords/Complex/clusterpercent.txt
  #Look for samples seen in at least two of the clusters
  multiple=$( fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt ${WRKDIR}/classifier/clusterfix/newCoords/*.events.reclassified.bedpe | awk '{ print $20 }' | tr -d "[]" | tr ',' '\n' | sort | uniq -c | awk '{ if ($1>1) print }' | wc -l )
  total=$( fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt ${WRKDIR}/classifier/clusterfix/newCoords/*.events.reclassified.bedpe | awk '{ print $20 }' | tr -d "[]" | tr ',' '\n' | sort | uniq -c | wc -l )
  percent=$( echo "scale=3; ${multiple}/${total}" | bc -l )
  #Overlap between two clusters
  overlap=$( fgrep -wf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/listcheck.txt ${WRKDIR}/classifier/clusterfix/newCoords/*.events.reclassified.bedpe | awk '{ print $20 }' | sort | uniq -c | awk '{ if ($1>1) print "yes"; else print "no" }' | sort -r | head -n1 )
  echo -e "${id}\t${percent}\t${overlap}\t${clustermax}" >> ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complexEventsSampleOverlap.txt
done < ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ComplexEvents.txt 

#CNV within CNV
cat ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe | fgrep Valid | cut -f1-7,20 > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.txt
awk -v OFS="\t" '{print $1, $2, $6, $0}' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.txt > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.bed

#Bedtools intersect while removing abparts
cat ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.bed | bedtools intersect -wa -wb -a - -b ${abParts} | awk '{ print $10 }' > ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ids.abparts.txt
# bedtools intersect -wa -wb -a ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.bed -b ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ALL.CNV.valid.bed | fgrep -vwf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ids.abparts.txt | awk '{ if ($10!=$21 && $6<$16 && $9>$20 && $11==$22) print }'

#Write final output & clean up
awk -v min=${min_sample_concordance} '{ if ($2>=min) print $1 }' ${WRKDIR}/classifier/clusterfix/newCoords/Complex/complexEventsSampleOverlap.txt | fgrep -wf - ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ComplexEvents.txt | fgrep -wvf ${WRKDIR}/classifier/clusterfix/newCoords/Complex/ids.abparts.txt > ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list
rm -rf ${WRKDIR}/classifier/clusterfix/newCoords/Complex

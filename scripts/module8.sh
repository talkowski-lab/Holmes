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
fgrep RecipCTX ${WRKDIR}/classifier/clusterfix/newCoords/transloc.classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cluster type; do
  observations=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f19 )
  samples=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f20 | tr -d "[]" )
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
fgrep Insertion ${WRKDIR}/classifier/clusterfix/newCoords/transloc.classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cluster class; do
  observations=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f19 )
  samples=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f20 | tr -d "[]" )
  if [ ${class} == "Insertion" ]; then
    echo -e "${chrA}\t${minA}\t${maxA}\t${chrB}\t${minB}\t${maxB}\t${COHORT_ID}_insertion\t+\t+\t${observations}\t[${samples}]\t[${cluster}]"
  else
    echo -e "${chrA}\t${minA}\t${maxA}\t${chrB}\t${minB}\t${maxB}\t${COHORT_ID}_insertion\t-\t+\t${observations}\t[${samples}]\t[${cluster}]"
  fi
done < ${pool} >> ${agg}
fgrep Insertion ${WRKDIR}/classifier/clusterfix/newCoords/complex_classifications.list | fgrep -v DupFlanked > ${pool}
while read chrA minA maxA chrB minB maxB cxcluster class flag; do
  clusters=$( fgrep -w ${cxcluster} ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list | sed 's/\t/,/g' | sed 's/,$//g' )
  samples=$( while read cluster; do fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g'; done < <( echo ${clusters} | sed 's/,/\n/g' ) | sort | uniq | paste -s -d, )
  observations=$( echo ${samples} | sed 's/,/\n/g' | wc -l )
  echo -e "${chrA}\t${minA}\t${maxA}\t${chrB}\t${minB}\t${maxB}\t${COHORT_ID}_insertion\t+\t+\t${observations}\t[${samples}]\t[${clusters}]"
done < ${pool} >> ${agg}
sort -nk1,1 -k4,4n -k2,2n -k5,5n ${agg} | awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, $6, $7"_"NR, $8, $9, $10, $11, $12 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.insertion.bedpe

#Write out complex (resolved) events as bed; bed intervals correspond to min and max coordinates. If inter-chromosomal, chromosome="Multiple" and coordinates will be ".".
echo -e "#chr\tmin\tmax\tID\tpredicted_event_class\tobservations\tsamples\tclusters" > ${WRKDIR}/final_variants/${COHORT_ID}.complex.bed
fgrep Valid ${WRKDIR}/classifier/clusterfix/newCoords/inversion.classifications.list | fgrep -v SE | fgrep -v SIMPLE | fgrep -v INS | fgrep -v FLAG > ${pool}
while read clusterlist chr minp maxp minm maxm class_old trash; do
  case ${class_old} in
    DEL3)
      class="INVdel"
      ;;
    DEL5)
      class="delINV"
      ;;
    DEL35|LINKED_del3INVdel5)
      class="delINVdel"
      ;;
    DUP3_DEL5|LINKED_dup3INVdel5)
      class="delINVdup"
      ;;
    DUP3)
      class="INVdup"
      ;;
    DUP5_DEL3|LINKED_del3INVdup5)
      class="dupINVdel"
      ;;
    DUP5)
      class="dupINV"
      ;;
    dupINVdup|LINKED_dup3INVdup5)
      class="dupINVdup"
      ;;
    LINKED_INVdel5_OTHER)
      class="delINV_other"
      ;;
    LINKED_INVdup5_OTHER)
      class="dupINV_other"
      ;;    
    LINKED_del3INV_OTHER)
      class="INVdel_other"
      ;;    
    LINKED_dup3INV_OTHER)
      class="INVdup_other"
      ;;  
    TANDEM_DUP_INV)
      class="IR"
      ;;  
    LINKED_OTHER|*)
      class="INV_other"
      ;;  
  esac
  min=$( echo -e "${minp}\n${maxp}\n${minm}\n${maxm}" | fgrep -v "NA" | sort -nk1,1 | head -n1 )
  max=$( echo -e "${minp}\n${maxp}\n${minm}\n${maxm}" | fgrep -v "NA" | sort -nk1,1 | tail -n1 )
  if [ $( echo ${clusterlist} | sed 's/,/\n/g' | wc -l ) -gt 1 ]; then
    echo ${clusterlist} | sed 's/,/\n/g' > ${clusters}
    samples=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | paste -s -d, )
    observations=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | wc -l )
  else
    samples=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" )
    observations=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f19 )
  fi
  echo -e "${chr}\t${min}\t${max}\t${COHORT_ID}_complexSV\t${class}\t${observations}\t[${samples}]\t[${clusters}]"
done < ${pool} > ${agg}
grep -e 'DupFlankedInsertion\|Inversion_FlankingDel3\|Inversion_FlankingDel5' ${WRKDIR}/classifier/clusterfix/newCoords/complex_classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cxcluster class_old flag; do
  clusters=$( fgrep -w ${cxcluster} ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list | sed 's/\t/,/g' | sed 's/,$//g' )
  samples=$( while read cluster; do fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g'; done < <( echo ${clusters} | sed 's/,/\n/g' ) | sort | uniq | paste -s -d, )
  observations=$( echo ${samples} | sed 's/,/\n/g' | wc -l )
  if [ ${chrA} == ${chrB} ]; then
    min=$( echo -e "${minA}\n${maxA}\n${minB}\n${maxB}" | fgrep -v "NA" | sort -nk1,1 | head -n1 )
    max=$( echo -e "${minA}\n${maxA}\n${minB}\n${maxB}" | fgrep -v "NA" | sort -nk1,1 | tail -n1 )
    case ${class_old} in
      Inversion_FlankingDel3)
        class="INVdel"
        ;;
      Inversion_FlankingDel5)
        class="delINV"
        ;;
      DupFlankedInsertion)
        class="complexDUPins"
        ;;
    esac
    echo -e "${chrA}\t${min}\t${max}\t${COHORT_ID}_complexSV\t${class}\t${observations}\t[${samples}]\t[${clusters}]"
  else
    case ${class_old} in
      Inversion_FlankingDel3)
        class="INVdel"
        ;;
      Inversion_FlankingDel5)
        class="delINV"
        ;;
      DupFlankedInsertion)
        class="complexDUPins"
        ;;
    esac
    echo -e "Multiple\t.\t.\t${COHORT_ID}_complexSV\t${class}\t${observations}\t[${samples}]\t[${clusters}]"
  fi
done < ${pool} >> ${agg}
sort -nk1,1 -k4,4n -k2,2n -k5,5n ${agg} | awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, $7, $8 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.complex.bed

#Write out all unresolved events as bed; bed intervals correspond to min and max coordinates. If inter-chromosomal, chromosome="Multiple" and coordinates will be ".".
echo -e "#chr\tmin\tmax\tID\tobservations\toutput_class\tsamples\tclusters" > ${WRKDIR}/final_variants/${COHORT_ID}.unresolved.bed
grep -e 'Unresolved\|SKIP_MissingBPs' ${WRKDIR}/classifier/clusterfix/newCoords/complex_classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cxcluster class_old flag; do
  clusters=$( fgrep -w ${cxcluster} ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list | sed 's/\t/,/g' | sed 's/,$//g' )
  samples=$( while read cluster; do fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g'; done < <( echo ${clusters} | sed 's/,/\n/g' ) | sort | uniq | paste -s -d, )
  observations=$( echo ${samples} | sed 's/,/\n/g' | wc -l )
  nchrs=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | sed 's/\:/\t/g' | awk -v OFS="\n" '{ print $2, $5 }' | sort | uniq | wc -l )
  if [ ${nchrs} -eq 1 ]; then
    chr=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | sed 's/\:/\t/g' | awk -v OFS="\n" '{ print $2, $5 }' | sort | uniq )
    min=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | awk -v OFS="\n" '{ print $2, $3, $5, $6 }' | sort -nk1,1 | head -n1 )
    max=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | awk -v OFS="\n" '{ print $2, $3, $5, $6 }' | sort -nk1,1 | tail -n1 )
    echo -e "${chr}\t${min}\t${max}\t${COHORT_ID}_unresolved\tcomplex\t${observations}\t[${samples}]\t[${clusters}]"
  else
    echo -e "Multiple\t.\t.\t${COHORT_ID}_unresolved\tcomplex\t${observations}\t[${samples}]\t[${clusters}]"
  fi
done < ${pool} > ${agg}
cut -f7 ${WRKDIR}/classifier/clusterfix/newCoords/complex_classifications.list | fgrep -wvf - ${WRKDIR}/classifier/clusterfix/newCoords/${COHORT_ID}.putative_complex_sites.list | sed 's/\t/,/g' | sed 's/$,//g' > ${pool}
while read clusters; do
  samples=$( while read cluster; do fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g'; done < <( echo ${clusters} | sed 's/,/\n/g' ) | sort | uniq | paste -s -d, )
  observations=$( echo ${samples} | sed 's/,/\n/g' | wc -l )
  nchrs=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | sed 's/\:/\t/g' | awk -v OFS="\n" '{ print $2, $5 }' | sort | uniq | wc -l )
  if [ ${nchrs} -eq 1 ]; then
    chr=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | sed 's/\:/\t/g' | awk -v OFS="\n" '{ print $2, $5 }' | sort | uniq )
    min=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | awk -v OFS="\n" '{ print $2, $3, $5, $6 }' | sort -nk1,1 | head -n1 )
    max=$( echo -e "${clusters}" | sed 's/,/\n/g' | fgrep -wf - ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | awk -v OFS="\n" '{ print $2, $3, $5, $6 }' | sort -nk1,1 | tail -n1 )
    echo -e "${chr}\t${min}\t${max}\t${COHORT_ID}_unresolved\tcomplex\t${observations}\t[${samples}]\t[${clusters}]"
  else
    echo -e "Multiple\t.\t.\t${COHORT_ID}_unresolved\tcomplex\t${observations}\t[${samples}]\t[${clusters}]"
  fi
done < ${pool} >> ${agg}
fgrep Unresolved ${WRKDIR}/classifier/clusterfix/newCoords/transloc.classifications.list > ${pool}
while read chrA minA maxA chrB minB maxB cluster type; do
  observations=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f19 )
  samples=$( fgrep -w ${cluster} ${WRKDIR}/raw_clusters/${COHORT_ID}.transloc.events.bedpe | cut -f20 | tr -d "[]" )
  echo -e "Multiple\t.\t.\t${COHORT_ID}_unresolved\ttransloc\t${observations}\t[${samples}]\t[${cluster}]"
done < ${pool} >> ${agg}
grep -e 'SE\|FLAG' ${WRKDIR}/classifier/clusterfix/newCoords/inversion.classifications.list | fgrep Valid > ${pool}
while read clusterlist chr minp maxp minm maxm class_old trash; do
  min=$( echo -e "${minp}\n${maxp}\n${minm}\n${maxm}" | fgrep -v "NA" | sort -nk1,1 | head -n1 )
  max=$( echo -e "${minp}\n${maxp}\n${minm}\n${maxm}" | fgrep -v "NA" | sort -nk1,1 | tail -n1 )
  if [ $( echo ${clusterlist} | sed 's/,/\n/g' | wc -l ) -gt 1 ]; then
    echo ${clusterlist} | sed 's/,/\n/g' > ${clusters}
    samples=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | paste -s -d, )
    observations=$( fgrep -wf ${clusters} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" | sed 's/,/\n/g' | sort | uniq | wc -l )
  else
    samples=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f20 | tr -d "[]" )
    observations=$( fgrep -w ${clusterlist} ${WRKDIR}/raw_clusters/${COHORT_ID}.*.events.bedpe | cut -f19 )
  fi
  echo -e "${chr}\t${min}\t${max}\t${COHORT_ID}_unresolved\tinversion\t${observations}\t[${samples}]\t[${clusterlist}]"
done < ${pool} >> ${agg}
sort -nk1,1 -k2,2n ${agg} | awk -v OFS="\t" '{ print $1, $2, $3, $4"_"NR, $5, $6, $7, $8 }' >> ${WRKDIR}/final_variants/${COHORT_ID}.unresolved.bed

#Clean up
rm ${pool} ${clusters} ${agg}
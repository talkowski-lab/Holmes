#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Script to screen classifier events for retrotransposons

#Rules:
#1) only considers valid events
#2) requires at least 30% of all reads to be in an exon to report event\
#3) only considers primary contigs (no MT either)

#Read input
events=$1              #two-column tab-delimmed list.  First column: deletion/insertion/inversion/transloc. Second column: path to events.bedpe. Preferrably reclassified with reclassify_output.sh
clusters=$2            #two-column tab-delimmed list.  First column: deletion/insertion/inversion/transloc. Second column: path to clusters.bedpe. Must be filtered using cleanClusters_patch.sh
OUTDIR=$3              #output directory

#Sets params
refFlat=/data/talkowski/tools/bin/TGDB/BACKUP/hg19_refFlat_3_13_2014.bed
min_pct=0.3
TMPDIR=`mktemp -d`
TMPPOS=`mktemp`
module load bedtools/2.22.1
module load samtools/1.0
del_events=$( fgrep deletion ${events} | cut -f2 )
ins_events=$( fgrep insertion ${events} | cut -f2 )
inv_events=$( fgrep inversion ${events} | cut -f2 )
tloc_events=$( fgrep transloc ${events} | cut -f2 )
del_clusters=$( fgrep deletion ${clusters} | cut -f2 )
ins_clusters=$( fgrep insertion ${clusters} | cut -f2 )
inv_clusters=$( fgrep inversion ${clusters} | cut -f2 )
tloc_clusters=$( fgrep transloc ${clusters} | cut -f2 )

#Converts refFlat to list of exons
echo -e "\n\n***BUILDING EXON LIST***\n\n"
for contig in $( seq 1 22 ) X Y; do
  /data/talkowski/rlc47/code/dna-scripts/build_exons.sh hg19 chr${contig} 1 300000000 ${TMPDIR}/${contig}_exons.bed
done
cat ${TMPDIR}/*_exons.bed | sed 's/^chr//g' | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -nk2,2 | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${TMPDIR}/all_exons.bed

#Prepares outfile
echo -e "eID\tTotal_Mappings\tMappings_in_Exons\tPct_Exonic\tGene(s)" > ${OUTDIR}/candidate_retrotransposons.list

#Deletion
echo -e "\n\n***CHECKING DELETIONS***\n\n"
while read eID; do
  eidx=$( echo ${eID} | sed 's/deletion_/\t/g' | cut -f2 )
  awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $4, $6, $6+1"\n"$5, $7, $7+1 }' ${del_clusters} > ${TMPPOS} #Gets mapping coordinates of all reads in cluster
  num_rds=$( cat ${TMPPOS} | wc -l ) #counts number of mappings in cluster
  min_rds=$( echo "(( ${num_rds} * ${min_pct} ))" | bc | cut -d. -f1 ) #sets 15% read cutoff
  num_ovr=$( bedtools intersect -u -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | wc -l )
  if [[ ${num_ovr} -ge ${min_rds} ]]; then
    gene_ovr=$( bedtools intersect -wao -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | cut -f7 | sort | uniq | grep -ve '^$' | paste -s | sed 's/\t/,/g' )
    pct_exonic=$( echo "scale=2;(( 100*${num_ovr}/${num_rds} ))" | bc )
    echo -e "${eID}\t${num_rds}\t${num_ovr}\t${pct_exonic}%\t${gene_ovr}"
  fi
done < <( fgrep Valid ${del_events} | cut -f7 ) >> ${OUTDIR}/candidate_retrotransposons.list

#Insertion
echo -e "\n\n***CHECKING INSERTIONS***\n\n"
while read eID; do
  eidx=$( echo ${eID} | sed 's/insertion_/\t/g' | cut -f2 )
  awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $4, $6, $6+1"\n"$5, $7, $7+1 }' ${ins_clusters} > ${TMPPOS} #Gets mapping coordinates of all reads in cluster
  num_rds=$( cat ${TMPPOS} | wc -l ) #counts number of mappings in cluster
  min_rds=$( echo "(( ${num_rds} * ${min_pct} ))" | bc | cut -d. -f1 ) #sets 15% read cutoff
  num_ovr=$( bedtools intersect -u -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | wc -l )
  if [[ ${num_ovr} -ge ${min_rds} ]]; then
    gene_ovr=$( bedtools intersect -wao -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | cut -f7 | sort | uniq | grep -ve '^$' | paste -s | sed 's/\t/,/g' )
    pct_exonic=$( echo "scale=2;(( 100*${num_ovr}/${num_rds} ))" | bc )
    echo -e "${eID}\t${num_rds}\t${num_ovr}\t${pct_exonic}%\t${gene_ovr}"
  fi
done < <( fgrep Valid ${ins_events} | cut -f7 ) >> ${OUTDIR}/candidate_retrotransposons.list

#Inversion
echo -e "\n\n***CHECKING INVERSIONS***\n\n"
while read eID; do
  eidx=$( echo ${eID} | sed 's/inversion_/\t/g' | cut -f2 )
  awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $4, $6, $6+1"\n"$5, $7, $7+1 }' ${inv_clusters} > ${TMPPOS} #Gets mapping coordinates of all reads in cluster
  num_rds=$( cat ${TMPPOS} | wc -l ) #counts number of mappings in cluster
  min_rds=$( echo "(( ${num_rds} * ${min_pct} ))" | bc | cut -d. -f1 ) #sets 15% read cutoff
  num_ovr=$( bedtools intersect -u -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | wc -l )
  if [[ ${num_ovr} -ge ${min_rds} ]]; then
    gene_ovr=$( bedtools intersect -wao -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | cut -f7 | sort | uniq | grep -ve '^$' | paste -s | sed 's/\t/,/g' )
    pct_exonic=$( echo "scale=2;(( 100*${num_ovr}/${num_rds} ))" | bc )
    echo -e "${eID}\t${num_rds}\t${num_ovr}\t${pct_exonic}%\t${gene_ovr}"
  fi
done < <( fgrep Valid ${inv_events} | cut -f7 ) >> ${OUTDIR}/candidate_retrotransposons.list

#Translocation
echo -e "\n\n***CHECKING TRANSLOCATIONS***\n\n"
while read eID; do
  eidx=$( echo ${eID} | sed 's/transloc_/\t/g' | cut -f2 )
  awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $4, $6, $6+1"\n"$5, $7, $7+1 }' ${tloc_clusters} > ${TMPPOS} #Gets mapping coordinates of all reads in cluster
  num_rds=$( cat ${TMPPOS} | wc -l ) #counts number of mappings in cluster
  min_rds=$( echo "(( ${num_rds} * ${min_pct} ))" | bc | cut -d. -f1 ) #sets 15% read cutoff
  num_ovr=$( bedtools intersect -u -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | wc -l )
  if [[ ${num_ovr} -ge ${min_rds} ]]; then
    gene_ovr=$( bedtools intersect -wao -a ${TMPPOS} -b ${TMPDIR}/all_exons.bed | cut -f7 | sort | uniq | grep -ve '^$' | paste -s | sed 's/\t/,/g' )
    pct_exonic=$( echo "scale=2;(( 100*${num_ovr}/${num_rds} ))" | bc )
    echo -e "${eID}\t${num_rds}\t${num_ovr}\t${pct_exonic}%\t${gene_ovr}"
  fi
done < <( fgrep Valid ${tloc_events} | cut -f7 ) >> ${OUTDIR}/candidate_retrotransposons.list

#Cleanup
rm -rf ${TMPDIR}


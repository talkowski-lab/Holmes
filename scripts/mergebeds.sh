#!/bin/bash

######################################
# Talkowski Lab
# BED-format CNV Variant Merger
######################################

#Read input
in_list=$1             #input list. Format: sample ID <TAB> full/path/to/bedfile/to/be/merged
cluster_dist=$2        #maximum distance between breakpoints for concordance, in bases. 10kb for Eco with cnMOPS
contig=$3              #Contig to merge
name=$4                #output name. Will be appended to output file
OUTDIR=$5              #output directory. Will write to ${OUTDIR}/${name}.merged.bed

#9-col bed output format (tab-delimmed):
#1) chr
#2) mean bp A
#3) mean bp B
#4) consensus ID
#5) Merge span A
#6) Merge span B
#7) observations (count)
#8) samples (bracket-enclosed, comma-delimmed)

#Set params
module load bedtools/2.22.1
TMPDIR=`mktemp -d /tmp/tmp.XXXXXXXX`
samp_list=`mktemp`
clust=`mktemp`
remove=`mktemp`
prefinal=`mktemp`
add=`mktemp`

#Create $OUTDIR if necessary
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Write all samples to list
cut -f1 ${in_list} > ${samp_list}

#Clean sample beds to same format
while read ID bed; do
  paste <( awk -v OFS="\t" -v contig=${contig} -v ID=${ID} '{ if ($1==contig) print $1, $2, $3, ID, ID"_"NR }' ${bed} ) <( awk -v contig=${contig} '{ if ($1==contig) print $0 }' ${bed} | cut --complement -f1-3  | sed 's/\t/\//g' ) > ${TMPDIR}/${ID}.bed
done < ${in_list}

#cat & sort master bed
cat ${TMPDIR}/*bed | sed -e 's/^X/23/g' -e 's/^Y/24/g' | sort -nk1,1 -nk2,2 -nk3,3 > ${TMPDIR}/all.bed

#read through master bed & merge
while read chr start end ID eID info; do
  #skip if event has already been clustered
  if [ $( fgrep -w ${eID} ${remove} | wc -l ) -eq 0 ]; then
    awk -v chr=${chr} -v start=${start} -v end=${end} -v dist=${cluster_dist} '{ if ($1==chr && $2>=(start-dist) && $2<=(start+dist) && $3>=(end-dist) && $3<=(end+dist)) print $0 }' ${TMPDIR}/all.bed | fgrep -wvf ${remove} > ${clust}
    if [ $( cat ${clust} | wc -l ) -gt 1 ]; then
      #Calculates average breaks
      midA=$( cut -f2 ${clust} | awk '{ sum+=$1 } END { printf "%.0f\n", sum/NR }' )
      midB=$( cut -f3 ${clust} | awk '{ sum+=$1 } END { printf "%.0f\n", sum/NR }' )
      rep=1
      #searches cluster_dist around average breaks to avoid arbitrary cluster nucleation issue. Does not allow hits from samples already in cluster
      until [ ${rep} -eq 0 ]; do
        awk -v chr=${chr} -v midA=${midA} -v midB=${midB} -v dist=${cluster_dist} '{ if ($1==chr && $2>=(midA-dist) && $2<=(midA+dist) && $3>=(midB-dist) && $3<=(midB+dist)) print $0 }' ${TMPDIR}/all.bed | fgrep -wvf <( cut -f5 ${clust} ) | fgrep -wvf ${remove} | fgrep -wvf <( cut -f4 ${clust} ) > ${add}
        if [ $( cat ${add} | wc -l ) -gt 0 ]; then
          cat ${add} >> ${clust}
          rep=1
        else
          rep=0
        fi
      done
      echo ${chr}
      cut -f2 ${clust} | awk '{ sum+=$1 } END { printf "%.0f\n", sum/NR }'
      cut -f3 ${clust} | awk '{ sum+=$1 } END { printf "%.0f\n", sum/NR }'
      echo "blank"
      echo "$( cut -f2 ${clust} | sort -nk1,1 | tail -n1 )-$( cut -f2 ${clust} | sort -nk1,1 | head -n1 )" | bc
      echo "$( cut -f3 ${clust} | sort -nk1,1 | tail -n1 )-$( cut -f3 ${clust} | sort -nk1,1 | head -n1 )" | bc
      cat ${clust} | wc -l #prints number of observations
      cut -f4 ${clust} | sort | uniq | paste -s -d, | sed -e 's/^/[/g' -e 's/$/]/g' #prints samples
      cut -f6 ${clust} | paste -s -d\| | sed -e 's/^/[/g' -e 's/$/]/g' #prints info
      cut -f5 ${clust} >> ${remove}
    else
      awk -v OFS="\n" '{ print $1, $2, $3, "blank", "0", "0", "1", "["$4"]", "["$6"]" }' ${clust}
      cut -f5 ${clust} >> ${remove}
    fi
  fi
done < ${TMPDIR}/all.bed | paste - - - - - - - - - > ${prefinal}

#Print final output
awk -v OFS="\t" -v contig=${contig} -v name=${name} '{ print $1, $2, $3, name"_"contig"_"NR, $5, $6, $7, $8, $9 }' ${prefinal} | sort -nk1,1 -nk2,2 -nk3,3 | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${OUTDIR}/${name}.merged.${contig}.bed

#Clean up
rm -rf ${TMPDIR}
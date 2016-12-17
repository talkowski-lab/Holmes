#!/bin/bash

#Cluster filter patch
#Ryan Collins
#Talkowski Lab
#March 2015

##NOTE: RUN THIS AFTER JOINT CLUSTERING BUT BEFORE CLASSIFICATION

#Read Input
clusters=$1
outfile=$2
underscore_skip=$3 #number of underscores to include as part of the sample ID; "STUDY_ID" format would be underscore_skip=1, etc.
del=$4 #TRUE if you want to split cluster where individual sample don't overlap; only necessary for deletions.  ASSUMES DELETION CLUSTERS (+/-)

#Set params
module load bedtools/2.22.1
TMPCLUST=`mktemp`
exclusion_list=`mktemp`
samples=`mktemp`
work=`mktemp`
groups=`mktemp`
abParts=/data/talkowski/rlc47/src/abParts.bed
fields=$( echo "1+${underscore_skip}" | bc )

#Get last cluster idx
last=$( tail -n2 ${clusters} | cut -f1 | head -n1 )

#Remove pairs where read and mate map to identical position
awk -v OFS="\t" '{ if ($1=="" || ($4==$5 && $6!=$7) || $4!=$5) print $0 }' ${clusters} | grep -A1 . | grep -v "^--$" > ${TMPCLUST}

#Remove pairs that map inside abParts
while read chr start end; do
  awk -v OFS="\t" -v chr=${chr} -v start=${start} -v end=${end} '{ if ( ($4==chr && $6<=end && $6>=start) || ($5==chr && $7<=end && $7>=start) ) print $3 }' ${TMPCLUST} | sort -nk1,1 | uniq >> ${exclusion_list}
done < ${abParts}
fgrep -wvf ${exclusion_list} ${TMPCLUST} > ${TMPCLUST}2
mv ${TMPCLUST}2 ${TMPCLUST}

#Split del clusters
if [ ${del} == "TRUE" ]; then
  newclust=1
  for cidx in $( cut -f1 ${TMPCLUST} | sort -nk1,1 | uniq ); do
    awk -v OFS="\t" -v cidx=${cidx} '{ if ($1==cidx) print $0 }' ${TMPCLUST} > ${work}
    cut -f3 ${work} | cut -d_ -f1,${fields} | sort | uniq > ${samples}
    while read ID; do
      awk '{ print $4 }' ${work} | head -n1
      fgrep ${ID} ${work} | cut -f6 | sort -nk1,1 | head -n1
      fgrep ${ID} ${work} | cut -f6 | sort -nk1,1 | tail -n1
      echo ${ID}
    done < ${samples} | paste - - - - | sort -nk1,1 -k2,2 | bedtools merge -d -1 -c 4 -o distinct -i - > ${groups}
    fgrep -f <( head -n1 ${groups} | cut -f4 | sed 's/\,/\n/g' ) ${work} >> ${TMPCLUST}_del
    echo -e "" >> ${TMPCLUST}_del
    if [ $( cat ${groups} | wc -l ) -gt 1 ]; then
      while read chr start stop group; do
        fgrep -f <( echo ${group} | sed 's/\,/\n/g' ) ${work} | cut --complement -f1 | awk -v OFS="\t" -v idx=${newclust} -v chr=${chr} '{ print "R"chr"_"idx, $0 }' >> ${TMPCLUST}_del
        echo -e "" >> ${TMPCLUST}_del
        newclust=$( echo "${newclust}+1" | bc )
      done < <( sed '1d' ${groups} )
    fi
  done
  mv ${TMPCLUST}_del ${TMPCLUST}
fi

#Write to outfile
mv ${TMPCLUST} ${outfile}

#Clean up
rm ${exclusion_list}

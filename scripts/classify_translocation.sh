#!/bin/bash

###################
# TRANSLOCATION CLASSIFICATION SCRIPT
# RYAN COLLINS
# TALKOWSKI LAB
# MARCH 2015
###################

# checks for appropriate input
if [ $# -eq 3 ]; then

  #Read input
  events=$1              #events.bedpe. Preferrably reclassified with reclassify_output.sh
  clusters=$2            #clusters.bedpe. Must be filtered using cleanClusters_patch.sh
  OUTDIR=$3              #output directory.  Will write to ${OUTDIR}/transloc.classifications.list 

  #Set params
  module load bedtools/2.22.1
  TMPDIR=/scratch/miket/rlc47temp/tmp.files
  TMPCLST=`mktemp`
  VALID=`mktemp`
  ORIENTATIONS=`mktemp`
  SE_IDs=`mktemp`
  DE_IDs=`mktemp`
  MULTIALLELIC_TMP=`mktemp`
  hetchrom=/data/talkowski/rlc47/src/GRCh37_heterochromatin.bed
  SEPP=`mktemp`
  SEMM=`mktemp`
  SEPM=`mktemp`
  SEMP=`mktemp`
  match_IDs=`mktemp`
  samps=`mktemp`
  matches=`mktemp`
  usedSE=`mktemp`

  #Write all valid events to VALID
  fgrep Valid ${events} > ${VALID}

  #Write header to outfile
  echo -e "#chrSource\tstartSource\tendSource\tchrSink\tminSink\tmaxSink\teID\tcategorization" > ${OUTDIR}/transloc.classifications.list

  #Iterate over clusters to clean unusual cases and assign processing groups
  echo "SORTING CLUSTERS INTO PRE-CATEGORIZATION CLASSES..."
  while read eID; do
    #get cluster and write to TMPCLST
    eidx=$( echo "${eID}" | sed 's/transloc_/\t/g' | cut -f2 )
    awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $0 }' ${clusters} > ${TMPCLST}
    #write proportions of pairwise orientations to ORIENTATIONS
    cut -f12-13 $TMPCLST | sed -e 's/\(113\|177\|145\|81\)/M/g' -e 's/\(129\|65\|97\|161\)/P/g' -e 's/\t//g' | sort | uniq -c | awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1 > ${ORIENTATIONS}
    types=$( cat ${ORIENTATIONS} | wc -l )
    #Write to either SE, DE, or MULTIALLELIC files based on breakdown of orientations
    case ${types} in
      #Write SE to SE file
      1)
        echo ${eID} >> ${SE_IDs}
        ;;
      #write DE to DE file if their orientations match; if not, write to SE file
      2)
        PM=$( fgrep "PM" ${ORIENTATIONS} | wc -l )
        MP=$( fgrep "MP" ${ORIENTATIONS} | wc -l )
        PP=$( fgrep "PP" ${ORIENTATIONS} | wc -l )
        MM=$( fgrep "MM" ${ORIENTATIONS} | wc -l )
        if [ ${PM} -gt 0 ] && [ ${MP} -gt 0 ]; then
          echo ${eID} >> ${DE_IDs}
        elif [ ${PP} -gt 0 ] && [ ${MM} -gt 0 ]; then
          echo ${eID} >> ${DE_IDs}
        else
          echo ${eID} >> ${SE_IDs}
        fi
        ;;
      #Write 3-enders to DE file
      3)
        echo ${eID} >> ${DE_IDs}
        ;;
      #If most dominant pairing ≥ 80% of reads in 4-ender, write to DE file. Else, write to MULTIALLELIC file
      4)
        PM=$( fgrep "PM" ${ORIENTATIONS} | cut -f1 )
        MP=$( fgrep "MP" ${ORIENTATIONS} | cut -f1 )
        PP=$( fgrep "PP" ${ORIENTATIONS} | cut -f1 )
        MM=$( fgrep "MM" ${ORIENTATIONS} | cut -f1 )
        ins=$( echo "${PM}+${MP}" | bc )
        inv=$( echo "${PP}+${MM}" | bc )
        total=$( echo "${ins}+${inv}" | bc )
        insp=$( echo "$( echo "scale=3;(( ${ins}/${total} ))" | bc ) * 1000" | bc | cut -d. -f1 )
        invp=$( echo "$( echo "scale=3;(( ${inv}/${total} ))" | bc ) * 1000" | bc | cut -d. -f1 )
        if [ ${insp} -ge 800 ] || [ ${invp} -ge 800 ]; then
          echo ${eID} >> ${DE_IDs}
        else
          echo ${eID} >> ${MULTIALLELIC_TMP}
        fi
        ;;
    esac
  done < <( cut -f7 ${VALID} )

  #CATEGORIZE CLUSTERED DOUBLE ENDERS
  echo "CATEGORIZING CLUSTERED DOUBLE-ENDERS..."
  while read eID; do
    #get cluster and write to TMPCLST
    eidx=$( echo "${eID}" | sed 's/transloc_/\t/g' | cut -f2 )
    awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $0 }' ${clusters} > ${TMPCLST}
    #Get dominant pairing
    cut -f12-13 $TMPCLST | sed -e 's/\(113\|177\|145\|81\)/M/g' -e 's/\(129\|65\|97\|161\)/P/g' -e 's/\t//g' | sort | uniq -c | awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1 > ${ORIENTATIONS}
    if [ $( fgrep "PM" ${ORIENTATIONS} | wc -l ) -gt 0 ]; then
      PM=$( fgrep "PM" ${ORIENTATIONS} | cut -f1 )
    else
      PM=0
    fi
    if [ $( fgrep "MP" ${ORIENTATIONS} | wc -l ) -gt 0 ]; then
      MP=$( fgrep "MP" ${ORIENTATIONS} | cut -f1 )
    else
      MP=0
    fi
    if [ $( fgrep "PP" ${ORIENTATIONS} | wc -l ) -gt 0 ]; then
      PP=$( fgrep "PP" ${ORIENTATIONS} | cut -f1  )
    else
      PP=0
    fi
    if [ $( fgrep "MM" ${ORIENTATIONS} | wc -l ) -gt 0 ]; then
      MM=$( fgrep "MM" ${ORIENTATIONS} | cut -f1 )
    else
      MM=0
    fi
    ins=$( echo "${PM}+${MP}" | bc )
    inv=$( echo "${PP}+${MM}" | bc )
    #Classify inverted insertion
    if [ ${inv} -gt ${ins} ]; then
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get +/+ coordinates
      plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==65 || $12==129) ) print $6; else if ($5==chrA && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==65 || $12==129) ) print $6; else if ($5==chrB && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      #Get -/- coordinates
      minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==113 || $12==177) ) print $6; else if ($5==chrA && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==113 || $12==177) ) print $6; else if ($5==chrB && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #ins B -> A if pA < mA and pB > mB
      if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tInverted_Insertion"
      #ins B -> A if pA> mA and pB < mB
      elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tInverted_Insertion"
      #potential reciprocal if pA < mA and pB < mB
      elif [ ${plusA} -le ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        #recheck with actual "closest coordinates"
        plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==65 || $12==129) ) print $6; else if ($5==chrA && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | tail -n1 )
        plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==65 || $12==129) ) print $6; else if ($5==chrB && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | tail -n1 )
        minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==113 || $12==177) ) print $6; else if ($5==chrA && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | head -n1 )
        minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==113 || $12==177) ) print $6; else if ($5==chrB && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | head -n1 )
         #ins B -> A if pA < mA and pB > mB
        if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
          ##WRITE OUTPUT -- FORMAT TBD
          echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tInverted_Insertion"
        #ins B -> A if pA> mA and pB < mB
        elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
          ##WRITE OUTPUT -- FORMAT TBD
          echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tInverted_Insertion"
        #otherwise, real recip TLOC
        else
          plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==65 || $12==129) ) print $6; else if ($5==chrA && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
          plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==65 || $12==129) ) print $6; else if ($5==chrB && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
          minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==113 || $12==177) ) print $6; else if ($5==chrA && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
          minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==113 || $12==177) ) print $6; else if ($5==chrB && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
          #Check P-Q assignments
          if [ $( awk -v chrA=${chrA} -v posA=${plusA} '{ if ($1==chrA && $2>=posA) print $0 }' <( fgrep centromere ${hetchrom} ) | wc -l ) -gt 0 ]; then
            echo -e "${chrA}\t${plusA}\t${minusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tPQ_RecipCTX"
          else
            echo -e "${chrB}\t${plusB}\t${minusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tPQ_RecipCTX"
          fi
        fi
      else
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${minusB}\t${plusB}\t${eID}\tComplex_Unresolved"
      fi
    #Else treat as same-strand insertion
    else
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get A+/- coordinates
      plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==97 || $12==161) ) print $6; else if ($5==chrA && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==81 || $12==145) ) print $6; else if ($5==chrA && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #Get B+/- coordinates
      plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==97 || $12==161) ) print $6; else if ($5==chrB && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==81 || $12==145) ) print $6; else if ($5==chrB && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #ins B -> A if pA < mA and pB > mB
      if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tInsertion"
      #ins B -> A if pA> mA and pB < mB
      elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tInsertion"
      # potential reciprocal if pA < mA and pB < mB
      elif [ ${plusA} -le ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        #recheck with actual "closest coordinates"
        plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==97 || $12==161) ) print $6; else if ($5==chrA && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | tail -n1 )
        minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==81 || $12==145) ) print $6; else if ($5==chrA && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | head -n1 )
        plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==97 || $12==161) ) print $6; else if ($5==chrB && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | tail -n1 )
        minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==81 || $12==145) ) print $6; else if ($5==chrB && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | head -n1 )
        #ins B -> A if pA < mA and pB > mB
        if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
          ##WRITE OUTPUT -- FORMAT TBD
          echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tInsertion"
        #ins B -> A if pA> mA and pB < mB
        elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
          ##WRITE OUTPUT -- FORMAT TBD
          echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tInsertion"
        #otherwise, real recip TLOC
        else
          plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==97 || $12==161) ) print $6; else if ($5==chrA && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
          minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==81 || $12==145) ) print $6; else if ($5==chrA && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
          plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==97 || $12==161) ) print $6; else if ($5==chrB && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
          minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==81 || $12==145) ) print $6; else if ($5==chrB && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
          #Check P/Q assignments
          if [ $( awk -v chrA=${chrA} -v posA=${plusA} '{ if ($1==chrA && $2>=posA) print $0 }' <( fgrep centromere ${hetchrom} ) | wc -l ) -gt 0 ]; then
            echo -e "${chrA}\t${plusA}\t${minusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tPP_RecipCTX"
          else
            echo -e "${chrA}\t${plusA}\t${minusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tQQ_RecipCTX"
          fi
        fi
      else
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${minusB}\t${plusB}\t${eID}\tComplex_Unresolved"
      fi
    fi
  done < ${DE_IDs} >> ${OUTDIR}/transloc.classifications.list

  #SPLIT SINGLE-ENDERS BY ORIENTATION
  while read eID; do
    #get cluster and write to TMPCLST
    eidx=$( echo "${eID}" | sed 's/transloc_/\t/g' | cut -f2 )
    awk -v OFS="\t" -v eidx=${eidx} '{ if ($1==eidx) print $0 }' ${clusters} > ${TMPCLST}
    #Get dominant pairing
    cut -f12-13 $TMPCLST | sed -e 's/\(113\|177\|145\|81\)/M/g' -e 's/\(129\|65\|97\|161\)/P/g' -e 's/\t//g' | sort | uniq -c | awk -v OFS="\t" '{ print $1, $2 }' | sort -nrk1,1 > ${ORIENTATIONS}
    dom=$( head -n1 ${ORIENTATIONS} | cut -f2 )
    #Write PP to SEPP
    if [ ${dom} == "PP" ]; then
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get nearest coordinates
      nearA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==65 || $12==129) ) print $6; else if ($5==chrA && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      nearB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==65 || $12==129) ) print $6; else if ($5==chrB && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      #write out
      echo -e "${chrA}\t${nearA}\t${chrB}\t${nearB}\t${eID}" >> ${SEPP}
    fi
    #Write MM to SEMM
    if [ ${dom} == "MM" ]; then
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get nearest coordinates
      nearA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==113 || $12==177) ) print $6; else if ($5==chrA && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      nearB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==113 || $12==177) ) print $6; else if ($5==chrB && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #write out
      echo -e "${chrA}\t${nearA}\t${chrB}\t${nearB}\t${eID}" >> ${SEMM}
    fi
    #Write PM to SEMP
    if [ ${dom} == "PM" ]; then
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get nearest coordinates
      nearA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==97 || $12==161) ) print $6; else if ($5==chrA && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      nearB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==81 || $12==145) ) print $6; else if ($5==chrB && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #write out
      echo -e "${chrA}\t${nearA}\t${chrB}\t${nearB}\t${eID}" >> ${SEPM}
    fi
    #Write MP to SEMP
    if [ ${dom} == "MP" ]; then
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get nearest coordinates
      nearA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==81 || $12==145) ) print $6; else if ($5==chrA && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      nearB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==97 || $12==161) ) print $6; else if ($5==chrB && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      #write out
      echo -e "${chrA}\t${nearA}\t${chrB}\t${nearB}\t${eID}" >> ${SEMP}
    fi
  done < ${SE_IDs}

  #LINK INVERTED SE THAT CAN DRAW OUT AS INV INS
  while read pchrA plusA pchrB plusB peID; do
    #Write cluster samples to samps
    fgrep -w ${peID} ${VALID} | cut -f20 | tr -d "[]" | sed 's/,/\n/g' > ${samps}
    #Get all corresponding -/- clusters on same two contigs
    awk -v pchrA=${pchrA} -v pchrB=${pchrB} '{ if ( ($1==pchrA || $1==pchrB) && ($3==pchrA || $3==pchrB) ) print $0 }' ${SEMM} | cut -f5 > ${match_IDs}
    #Check for at least one sample overlap between clusters
    #If multiple candidates, take one with most samples overlapping
    while read meID; do
      echo ${meID}
      fgrep -w ${meID} ${VALID} | awk '{ print $20 }' | tr -d "[]" | sed 's/,/\n/g' | fgrep -wf - ${samps} | wc -l
    done < ${match_IDs} | paste - - | sort -nrk2,2 > ${matches}
    if [ $( awk '{ if ($2>0) print $1 }' ${matches} | wc -l ) -gt 0 ]; then
      meID=$( head -n1 ${matches} | cut -f1 )
      echo ${meID} >> ${usedSE}
      #if match, test coordinates
      pidx=$( echo ${peID} | sed 's/transloc_/\t/g' | cut -f2 )
      midx=$( echo ${meID} | sed 's/transloc_/\t/g' | cut -f2 )
      awk -v pidx=${pidx} -v midx=${midx} -v OFS="\t" '{ if ($1==pidx || $1==midx) print $0 }' ${clusters} > ${TMPCLST}
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get +/+ coordinates
      plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==65 || $12==129) ) print $6; else if ($5==chrA && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==65 || $12==129) ) print $6; else if ($5==chrB && ($13==65 || $13==129) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      #Get -/- coordinates
      minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==113 || $12==177) ) print $6; else if ($5==chrA && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==113 || $12==177) ) print $6; else if ($5==chrB && ($13==113 || $13==177) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #ins B -> A if pA < mA and pB > mB
      if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${eID}\tInverted_Insertion"
      #ins B -> A if pA> mA and pB < mB
      elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${eID}\tInverted_Insertion"
      else
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${minusB}\t${plusB}\t${peID},${meID}\tComplex_Unresolved"
      fi
    else
      fgrep -w ${peID} ${VALID} | cut -f1-7 | paste - <( echo "SE_Unresolved" )
    fi
  done < ${SEPP} >> ${OUTDIR}/transloc.classifications.list
  #All remaining MM SEs are SE unresolved
  while read a b c d eID; do
    fgrep -w ${eID} ${VALID} | cut -f1-7 | paste - <( echo "SE_Unresolved" )
  done < <( fgrep -wvf ${usedSE} ${SEMM} ) >> ${OUTDIR}/transloc.classifications.list

  #LINK INVERTED SE THAT CAN DRAW OUT AS REG INS
  while read achrA plusA achrB minusB aeID; do
    #Write cluster samples to samps
    fgrep -w ${aeID} ${VALID} | cut -f20 | tr -d "[]" | sed 's/,/\n/g' > ${samps}
    #Get all corresponding -/- clusters on same two contigs
    awk -v achrA=${achrA} -v achrB=${achrB} '{ if ( ($1==achrA || $1==achrB) && ($3==achrA || $3==achrB) ) print $0 }' ${SEMP} | cut -f5 > ${match_IDs}
    #Check for at least one sample overlap between clusters
    #If multiple candidates, take one with most samples overlapping
    while read beID; do
      echo ${beID}
      fgrep -w ${beID} ${VALID} | awk '{ print $20 }' | tr -d "[]" | sed 's/,/\n/g' | fgrep -wf - ${samps} | wc -l
    done < ${match_IDs} | paste - - | sort -nrk2,2 > ${matches}
    if [ $( awk '{ if ($2>0) print $1 }' ${matches} | wc -l ) -gt 0 ]; then
      beID=$( head -n1 ${matches} | cut -f1 )
      echo ${beID} >> ${usedSE}
      #if match, test coordinates
      aidx=$( echo ${aeID} | sed 's/transloc_/\t/g' | cut -f2 )
      bidx=$( echo ${beID} | sed 's/transloc_/\t/g' | cut -f2 )
      awk -v aidx=${aidx} -v bidx=${bidx} -v OFS="\t" '{ if ($1==aidx || $1==bidx) print $0 }' ${clusters} > ${TMPCLST}
      #Get chromosomes
      chrA=$( head -n1 ${TMPCLST} | cut -f4 )
      chrB=$( head -n1 ${TMPCLST} | cut -f5 )
      #Get A+/- coordinates
      plusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==97 || $12==161) ) print $6; else if ($5==chrA && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      minusA=$( awk -v chrA=${chrA} '{ if ($4==chrA && ($12==81 || $12==145) ) print $6; else if ($5==chrA && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #Get B+/- coordinates
      plusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==97 || $12==161) ) print $6; else if ($5==chrB && ($13==97 || $13==161) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      minusB=$( awk -v chrB=${chrB} '{ if ($4==chrB && ($12==81 || $12==145) ) print $6; else if ($5==chrB && ($13==81 || $13==145) ) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #ins B -> A if pA < mA and pB > mB
      if [ ${plusA} -le ${minusA} ] && [ ${plusB} -gt ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrB}\t${minusB}\t${plusB}\t${chrA}\t${plusA}\t${minusA}\t${aeID},${beID}\tInsertion"
      #ins B -> A if pA> mA and pB < mB
      elif [ ${plusA} -gt ${minusA} ] && [ ${plusB} -le ${minusB} ]; then
        ##WRITE OUTPUT -- FORMAT TBD
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${plusB}\t${minusB}\t${aeID},${beID}\tInsertion"
      else
        echo -e "${chrA}\t${minusA}\t${plusA}\t${chrB}\t${minusB}\t${plusB}\t${aeID},${beID}\tComplex_Unresolved"
      fi
    else
      fgrep -w ${aeID} ${VALID} | cut -f1-7 | paste - <( echo "SE_Unresolved" )
    fi
  done < ${SEPM} >> ${OUTDIR}/transloc.classifications.list
  #All remaining MM SEs are SE unresolved
  while read a b c d eID; do
    fgrep -w ${eID} ${VALID} | cut -f1-7 | paste - <( echo "SE_Unresolved" )
  done < <( fgrep -wvf ${usedSE} ${SEMP} ) >> ${OUTDIR}/transloc.classifications.list

  #Add Multiallelic Complex to outfile
  while read eID; do
    fgrep -w ${eID} ${VALID} | cut -f1-7 | paste - <( echo "Multiallelic_Unresolved" )
  done < ${MULTIALLELIC_TMP} >> ${OUTDIR}/transloc.classifications.list

  #Clean up
  rm ${TMPCLST} ${VALID} ${ORIENTATIONS} ${SE_IDs} ${DE_IDs} ${MULTIALLELIC_TMP} ${SEPP} ${SEMM} ${SEPM} ${SEMP} ${match_IDs} ${samps} ${matches} ${usedSE}

# if inappropriate input; displays useage
else
  echo -e "\n\nTranslocation classification script\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
  echo "Usage:"
  echo "  classify_translocation.sh [tloc.events.bedpe] [tloc.clusters.bedpe] [OUTDIR]"
  echo ""
fi






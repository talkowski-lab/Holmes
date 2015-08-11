#!/bin/bash

###################
# INVERSION CLASSIFICATION SCRIPT
# RYAN COLLINS
# TALKOWSKI LAB
# MARCH 2015
###################

# checks for appropriate input
if [ $# -eq 8 ]; then

  #Read input
  events=$1              #events.bedpe. Preferrably reclassified with reclassify_output.sh
  clusters=$2            #clusters.bedpe. Must be filtered using cleanClusters_patch.sh
  complex_size=$3        #minimum size (bp) of breakpoint complexity to report. Should correspond to max clustering size
  deletions=$4           #deletion.events.bedpe run on same cohort.  Sample IDs must match.  Preferrably reclassified with reclassify_output.sh
  insertions=$5          #insertion.events.bedpe run on same cohort.  Sample IDs must match.  Preferrably reclassified with reclassify_output.sh
  cnMOPS_dels=$6         #list of cnMOPS del calls per sample.  Format: ID TAB path_to_cnMOPS.del.bed
  cnMOPS_dups=$7         #list of cnMOPS dup calls per sample.  Format: ID TAB path_to_cnMOPS.dup.bed
  OUTDIR=$8              #output directory.  Will write to ${OUTDIR}/inversion.classifications.list

  #Set params
  module load bedtools/2.22.1
  TMPDIR=/scratch/miket/rlc47temp/tmp.files
  TMPCLST=`mktemp`
  SAMPS=`mktemp`
  cnDELS=`mktemp`
  cnDUPS=`mktemp`
  peDELS=`mktemp`
  peDUPS=`mktemp`
  preOVR=`mktemp`
  cnDEL_OVR=`mktemp`
  cnDUP_OVR=`mktemp`
  peDEL_OVR=`mktemp`
  peDUP_OVR=`mktemp`
  OVR_table=`mktemp`
  RES=`mktemp`
  SEP=`mktemp`
  SEM=`mktemp`
  SELINK=`mktemp`
  LINKED=`mktemp`

  #Write header to output file
  echo -e "#EVENT_ID\tCHR\tMIN_PLUS\tMAX_PLUS\tMIN_MINUS\tMAX_MINUS\tCLASS\tINNER_SIZE\tOUTER_SIZE\tCALL\tDELS\tPRED_LOSS\tDUPS\tPRED_GAIN" > ${OUTDIR}/inversion.classifications.list

  #Iterate over all IDs in events file
  while read eID; do
    #store call & chr
    call=$( awk -v eID=${eID} '{ if ($7==eID) print $21 }' ${events} )
    chr=$( awk -v eID=${eID} '{ if ($7==eID) print $1 }' ${events} )
    #write original cluster to TMPCLST
    cidx=$( echo "${eID}" | awk -v FS="_" '{ print $NF }' )
    awk -v OFS="\t" -v cidx=${cidx} '{ if ($1==cidx) print $0 }' ${clusters} > ${TMPCLST}
    #write samples to SAMPS
    awk -v eID=${eID} '{ if ($7==eID) print $20 }' ${events} | tr -d "[]" | sed 's/,/\n/g' > ${SAMPS}
    #Gather dels & dups
    rm ${cnDELS} ${cnDUPS} ${peDELS} ${peDUPS}
    while read ID; do
      cut -f1-3 $( awk -v ID=${ID} '{ if ($1==ID) print $2 }' ${cnMOPS_dels} ) >> ${cnDELS}
      cut -f1-3 $( awk -v ID=${ID} '{ if ($1==ID) print $2 }' ${cnMOPS_dups} ) >> ${cnDUPS}
      fgrep Valid ${deletions} | awk -v ID=${ID} -v OFS="\t" '$20 ~ ID { print $1, $3, $3+1, $7"\n"$1, $5-1, $5, $7 }' >> ${peDELS}
      fgrep Valid ${deletions} | awk -v ID=${ID} -v OFS="\t" '$20 ~ ID { print $1, $2-1, $2, $7"\n"$1, $6, $6+1, $7 }'  >> ${peDUPS}
    done < ${SAMPS}
    
    #categorize +/+ single-enders
    if [ $( awk -v OFS="\t" -v eID=${eID} '{ if ($7==eID && $9=="+") print $0 }' ${events} | wc -l ) -gt 0 ]; then
      #Overlap with dosage alterations
      awk -v OFS="\t" -v eID=${eID} -v buffer=${complex_size} '{ if ($7==eID && $9=="+") print $1, $3-buffer, $3+buffer, $1, $6-buffer, $6+buffer }' ${events} > ${preOVR}
      bedtools pairtobed -type either -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
      bedtools pairtobed -type either -a ${preOVR} -b ${cnDUPS} > ${cnDUP_OVR}
      bedtools pairtobed -type either -a ${preOVR} -b ${peDELS} > ${peDEL_OVR}
      bedtools pairtobed -type either -a ${preOVR} -b ${peDUPS} > ${peDUP_OVR}
      nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
      if [ ${nDel} -gt 0 ]; then
        predLoss=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | awk '{ sum+=$1} END {print sum}' )
      else
        predLoss=0
      fi
      nDup=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
      if [ ${nDup} -gt 0 ]; then
        predGain=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | awk '{ sum+=$1} END {print sum}' )
      else
        predGain=0
      fi
      for file in ${cnDEL_OVR} ${cnDUP_OVR} ${peDEL_OVR} ${peDUP_OVR}; do
        wc -l ${file}
      done | sed 's/\s/\t/g' > ${OVR_table}
      #Write output depending on overlaps
      if [ $( awk '{ if ($1>0) print $2 }' ${OVR_table} | wc -l ) -eq 0 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_UNRESOLVED", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      elif [ $( cat ${cnDEL_OVR} | wc -l ) -gt 0 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_DEL", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      elif [ $( cat ${cnDUP_OVR} | wc -l ) -gt 0 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_DUP", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      elif [ $( cat ${cnDUP_OVR} | wc -l ) -eq 1 ] && [ $( cat ${cnDEL_OVR} | wc -l ) -eq 1 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_DEL_DUP", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      elif [ $( cat ${peDEL_OVR} | wc -l ) -gt 0 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_EXCISION", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      elif [ $( cat ${peDUP_OVR} | wc -l ) -gt 0 ]; then
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_INSERTION", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      else
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, $3, $6, "NA", "NA", "SE5_COMPLEX_OTHER", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      fi
    
    #categorize -/- single-enders
    elif [ $( awk -v OFS="\t" -v eID=${eID} '{ if ($7==eID && $9=="-") print $0 }' ${events} | wc -l ) -gt 0 ]; then
      #Check for inv tandem dup
      mmin=$( awk -v eID=${eID} '{ if ($7==eID) print $2 }' ${events} )
      mmax=$( awk -v eID=${eID} '{ if ($7==eID) print $5 }' ${events} )
      awk -v OFS="\t" -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmin-buffer, mmax+buffer}' ${events} > ${preOVR}
      bedtools intersect -r -f 0.51 -a ${preOVR} -b ${cnDUPS} > ${cnDUP_OVR}
      if [ $( echo "${mmax}-${mmin}" | bc ) -ge ${complex_size} ] && [ $( cat ${cnDUP_OVR} | wc -l ) -gt 0 ]; then
        nDup=$( cat ${cnDUP_OVR} | wc -l )
        predGain=$( echo "${mmax}-${mmin}" | bc )
        awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "TANDEM_DUP_INV", "NA", "NA", $21, "0", "0", nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
      else
        #Overlap with dosage alterations
        awk -v OFS="\t" -v eID=${eID} -v buffer=${complex_size} '{ if ($7==eID && $9=="-") print $1, $2-buffer, $2+buffer, $1, $5-buffer, $5+buffer }' ${events} > ${preOVR}
        bedtools pairtobed -type either -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
        bedtools pairtobed -type either -a ${preOVR} -b ${cnDUPS} > ${cnDUP_OVR}
        bedtools pairtobed -type either -a ${preOVR} -b ${peDELS} > ${peDEL_OVR}
        bedtools pairtobed -type either -a ${preOVR} -b ${peDUPS} > ${peDUP_OVR}
        nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
        if [ ${nDel} -gt 0 ]; then
          predLoss=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | awk '{ sum+=$1} END {print sum}' )
        else
          predLoss=0
        fi
        nDup=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
        if [ ${nDup} -gt 0 ]; then
          predGain=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | awk '{ sum+=$1} END {print sum}' )
        else
          predGain=0
        fi
        for file in ${cnDEL_OVR} ${cnDUP_OVR} ${peDEL_OVR} ${peDUP_OVR}; do
          wc -l ${file}
        done | sed 's/\s/\t/g' > ${OVR_table}
        #Write output depending on overlaps
        if [ $( awk '{ if ($1>0) print $2 }' ${OVR_table} | wc -l ) -eq 0 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_UNRESOLVED", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        elif [ $( cat ${cnDEL_OVR} | wc -l ) -gt 0 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_DEL", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        elif [ $( cat ${cnDUP_OVR} | wc -l ) -gt 0 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_DUP", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        elif [ $( cat ${cnDUP_OVR} | wc -l ) -eq 1 ] && [ $( cat ${cnDEL_OVR} | wc -l ) -eq 1 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="+") print eID, chr, "NA", "NA", $2, $5, "SE3_DEL_DUP", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        elif [ $( cat ${peDEL_OVR} | wc -l ) -gt 0 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_EXCISION", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        elif [ $( cat ${peDUP_OVR} | wc -l ) -gt 0 ]; then
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_INSERTION", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        else
          awk -v OFS="\t" -v eID=${eID} -v chr=${chr} -v nDel=${nDel} -v predLoss=${predLoss} -v nDup=${nDup} -v predGain=${predGain} '{ if ($7==eID && $9=="-") print eID, chr, "NA", "NA", $2, $5, "SE3_COMPLEX_OTHER", "NA", "NA", $21, nDel, predLoss, nDup, predGain }' ${events} >> ${OUTDIR}/inversion.classifications.list
        fi
      fi
    
    #categorize double-enders
    else      
      #get closest coordinates of +/+ break (3rd quartile)
      pA=$( awk '{ if ($12==129 || $12==65) print $6 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      pB=$( awk '{ if ($13==129 || $13==65) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
      #get closest coordinates of -/- break (1st quartile)
      mA=$( awk '{ if ($12==113 || $12==177) print $6 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      mB=$( awk '{ if ($13==113 || $13==177) print $7 }' ${TMPCLST} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
      #determine min and max coordinates of +/+ and -/- breaks
      pmin=$( echo -e "${pA}\n${pB}" | sort -n | head -n1 )
      pmax=$( echo -e "${pA}\n${pB}" | sort -n | tail -n1 )
      mmin=$( echo -e "${mA}\n${mB}" | sort -n | head -n1 )
      mmax=$( echo -e "${mA}\n${mB}" | sort -n | tail -n1 )

      #PMIN<MMIN<PMAX<MMAX
      if [ ${pmin} -lt ${mmin} ] && [ ${mmin} -lt ${pmax} ] && [ ${pmax} -lt ${mmax} ]; then
        isize=$( echo -e "${pmax}-${mmin}" | bc )
        osize=$( echo -e "${mmax}-${pmin}" | bc )
        #Simple if MMIN-PMIN<Size & MMAX-PMAX<Size
        if [ $( echo "${mmin}-${pmin}" | bc ) -lt ${complex_size} ] && [ $( echo "${mmax}-${pmax}" | bc ) -lt ${complex_size} ]; then
          echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tSIMPLE\t${isize}\t${osize}\t${call}\t0\t0\t0\t0" >> ${OUTDIR}/inversion.classifications.list
        #3' DEL if MMIN-PMIN>=Size & MMAX-PMAX<Size
        elif [ $( echo "${mmin}-${pmin}" | bc ) -ge ${complex_size} ] && [ $( echo "${mmax}-${pmax}" | bc ) -lt ${complex_size} ]; then
          #Overlap with dosage alterations
          awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmin-buffer, pmin+buffer, $1, mmin-buffer, mmin+buffer}' ${events} > ${preOVR}
          bedtools pairtobed -type either -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
          nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
          predLoss=$( echo "${mmin}-${pmin}" | bc )
          echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDEL3\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
        #5' DEL if MMIN-PMIN<Size & MMAX-PMAX>=Size
        elif [ $( echo "${mmin}-${pmin}" | bc ) -lt ${complex_size} ] && [ $( echo "${mmax}-${pmax}" | bc ) -ge ${complex_size} ]; then
          #Overlap with dosage alterations
          awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmax-buffer, pmax+buffer, $1, mmax-buffer, mmax+buffer}' ${events} > ${preOVR}
          bedtools pairtobed -type either -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
          nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
          predLoss=$( echo "${mmax}-${pmax}" | bc )
          echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDEL5\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
        #5' AND 3' DEL if MMIN-PMIN>=Size & MMAX-PMAX>=Size
        elif [ $( echo "${mmin}-${pmin}" | bc ) -ge ${complex_size} ] && [ $( echo "${mmax}-${pmax}" | bc ) -ge ${complex_size} ]; then
          #Overlap with dosage alterations
          awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmin-buffer, pmin+buffer, $1, mmin-buffer, mmin+buffer"\n"$1, pmax-buffer, pmax+buffer, $1, mmax-buffer, mmax+buffer }' ${events} > ${preOVR}
          bedtools pairtobed -type either -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
          nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
          predLoss=$( echo "(${mmax}-${pmax})+(${mmin}-${pmin})" | bc )
          echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDEL35\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
        fi

      #MMIN<PMIN<PMAX<MMAX
      elif [ ${mmin} -lt ${pmin} ] && [ ${pmin} -lt ${pmax} ] && [ ${pmax} -lt ${mmax} ]; then
        #Overlap with dosage alterations
        awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmin-buffer, pmin+buffer}' ${events} > ${preOVR}
        bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDUPS} > ${cnDUP_OVR}
        nDup=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
        predGain=$( echo "${pmin}-${mmin}" | bc )
        #Check for excision
        awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmin-buffer, pmin+buffer}' ${events} > ${preOVR}
        excision=$( bedtools intersect -f 0.1 -r -u -a ${preOVR} -b ${peDELS} | wc -l )
        #if dup match and no excision match, 3' DUP
        if [ ${nDup} -gt 0 ] && [ ${excision} -eq 0 ]; then
          isize=$( echo "${pmax}-${mmin}" | bc )
          osize=$( echo "${mmax}-${pmin}" | bc )
          #check for 5' DEL
          if [ $( echo -e "${mmax}-${pmax}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmax-buffer, mmax+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmax}-${pmax}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDUP3_DEL5\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t${nDup}\t${predGain}" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDUP3\t${isize}\t${osize}\t${call}\t0\t0\t${nDup}\t${predGain}" >> ${OUTDIR}/inversion.classifications.list
          fi
        #if excision match and no dup match, 3' EXCISION-INS
        elif [ ${nDup} -eq 0 ] && [ ${excision} -gt 0 ]; then
          isize=$( echo "${pmin}-${mmin}" | bc )
          osize=$( echo "${mmax}-${pmax}" | bc )
          #check for insertion site DEL (iDEL)
          if [ $( echo -e "${mmax}-${pmax}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmax-buffer, mmax+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmax}-${pmax}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tEXINS3_iDEL\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tEXINS3\t${isize}\t${osize}\t${call}\t0\t0\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          fi
        # no excision match nor dup match, then assume 3' INS with unknown origin, could be duped or excised
        else
          isize=$( echo "${pmin}-${mmin}" | bc )
          osize=$( echo "${mmax}-${pmax}" | bc )
          #check for insertion site DEL (iDEL)
          if [ $( echo -e "${mmax}-${pmax}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmax-buffer, mmax+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmax}-${pmax}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tINS3_UNK_iDEL\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tINS3_UNK\t${isize}\t${osize}\t${call}\t0\t0\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          fi
        fi

      #PMIN<MMIN<MMAX<PMAX
      elif [ ${pmin} -lt ${mmin} ] && [ ${mmin} -lt ${mmax} ] && [ ${mmax} -lt ${pmax} ]; then
        #Overlap with dosage alterations
        awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmax-buffer, pmax+buffer}' ${events} > ${preOVR}
        bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDUPS} > ${cnDUP_OVR}
        nDup=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
        predGain=$( echo "${pmax}-${mmax}" | bc )
        #Check for excision
        awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmax-buffer, pmax+buffer}' ${events} > ${preOVR}
        excision=$( bedtools intersect -f 0.1 -r -u -a ${preOVR} -b ${peDELS} | wc -l )
        #if dup match and no excision match, 5' DUP
        if [ ${nDup} -gt 0 ] && [ ${excision} -eq 0 ]; then
          isize=$( echo "${pmax}-${mmax}" | bc )
          osize=$( echo "${mmin}-${pmin}" | bc )
          #check for 3' DEL
          if [ $( echo -e "${mmin}-${pmin}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmin-buffer, mmin+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmin}-${pmin}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDUP5_DEL3\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t${nDup}\t${predGain}" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tDUP5\t${isize}\t${osize}\t${call}\t0\t0\t${nDup}\t${predGain}" >> ${OUTDIR}/inversion.classifications.list
          fi
        #if excision match and no dup match, 3' EXCISION-INS
        elif [ ${nDup} -eq 0 ] && [ ${excision} -gt 0 ]; then
          isize=$( echo "${pmax}-${mmax}" | bc )
          osize=$( echo "${mmin}-${pmin}" | bc )
          #check for insertion site DEL (iDEL)
          if [ $( echo -e "${mmin}-${pmin}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmin-buffer, pmin+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmin}-${pmin}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tEXINS5_iDEL\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tEXINS5\t${isize}\t${osize}\t${call}\t0\t0\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          fi
        # no excision match nor dup match, then assume 3' INS with unknown origin, could be duped or excised
        else
          isize=$( echo "${pmax}-${mmax}" | bc )
          osize=$( echo "${mmin}-${pmin}" | bc )
          #check for insertion site DEL (iDEL)
          if [ $( echo -e "${mmin}-${pmin}" | bc ) -gt ${complex_size} ]; then
            #overlap with del dosage
            awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, pmin-buffer, mmin+buffer}' ${events} > ${preOVR}
            bedtools intersect -r -f 0.1 -a ${preOVR} -b ${cnDELS} > ${cnDEL_OVR}
            nDel=$( cut -f7-9 ${cnDEL_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
            predLoss=$( echo "${mmin}-${pmin}" | bc )
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tINS5_UNK_iDEL\t${isize}\t${osize}\t${call}\t${nDel}\t${predLoss}\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          else
            echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tINS5_UNK\t${isize}\t${osize}\t${call}\t0\t0\t0\t0" >> ${OUTDIR}/inversion.classifications.list
          fi
        fi

      #MMIN<PMIN<MMAX<PMAX
      elif [ ${mmin} -lt ${pmin} ] && [ ${pmin} -lt ${mmax} ] && [ ${mmax} -lt ${pmax} ]; then
        isize=$( echo "${pmax}-${mmax}" | bc )
        osize=$( echo "${pmin}-${mmin}" | bc )
        #Overlap with dosage alterations
        awk -v OFS="\t" -v pmin=${pmin} -v pmax=${pmax} -v mmin=${mmin} -v mmax=${mmax} -v buffer=${complex_size} -v eID=${eID} '{ if ($7==eID) print $1, mmin-buffer, pmin+buffer"\n"$1, mmax-buffer, pmax+buffer}' ${events} > ${preOVR}
        bedtools intersect -wao -r -f 0.1 -a ${preOVR} -b ${cnDUPS} | cut -f4-6 | sort -nk2,2 -k3,3 | uniq > ${cnDUP_OVR}
        nDup=$( cut -f7-9 ${cnDUP_OVR} | sort -nk2,2 | bedtools merge -i - | awk '{ print $3-$2 }' | wc -l )
        predGain=$( echo "${pmax}-${mmax}" | bc )
        echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tdupINVdup\t${isize}\t${osize}\t${call}\t0\t0\t${nDup}\t${predGain}" >> ${OUTDIR}/inversion.classifications.list

      #FLAG - for DEs that don't fit above categories, inspect by hand
      else
        echo -e "${eID}\t${chr}\t${pmin}\t${pmax}\t${mmin}\t${mmax}\tFLAG\tNA\tNA\t${call}\tNA\tNA\tNA\tNA" >> ${OUTDIR}/inversion.classifications.list
      fi
    fi
  done < <( cut -f7 ${events} )

  #Split remaining SE & attempt linking OF VALID ONLY; WILL NOT LINK BREAKPOINTS WITH > 5MB SEPARATING PMIN/MMIN OR PMAX/MMAX
  fgrep SE5 ${OUTDIR}/inversion.classifications.list | awk -v OFS="\t" '{ print $2, $3, $4, $1, ".", "+", $7, $10 }' | sort -nk1,2 > ${SEP}
  fgrep SE3 ${OUTDIR}/inversion.classifications.list | awk -v OFS="\t" '{ print $2, $5, $6, $1, ".", "-", $7, $10 }' | sort -nk1,2 > ${SEM}
  fgrep -v SE5 ${OUTDIR}/inversion.classifications.list | fgrep -v SE3 > ${RES}
  bedtools intersect -S -wao -r -f 0.05 -a <( fgrep Valid ${SEP} ) -b <( fgrep Valid ${SEM} ) | awk '{ if ($9!="." && sqrt(($10-$2)*($10-$2))<=5000000 && sqrt(($11-$3)*($11-$3))<=5000000) print $0 }' > ${SELINK}
  awk -v OFS="\n" '{ print $4, $12 }' ${SELINK} | fgrep -wvf - <( cat ${SEP} ${SEM} ) | cut -f4 | fgrep -wf - ${OUTDIR}/inversion.classifications.list >> ${RES}

  #Categorize linked SE pairs
  while read chrP minP maxP eIDP x stP catP callP chrM minM maxM eIDM x stM catM callM ovr; do
    #Recategorize
    if [ ${catP} == "SE5_DEL" ] && [ ${catM} == "SE3_DEL" ]; then
      catN="LINKED_del3INVdel5"
    elif [ ${catP} == "SE5_DEL" ] && [ ${catM} == "SE3_DUP" ]; then
      catN="LINKED_dup3INVdel5"
    elif [ ${catP} == "SE5_DUP" ] && [ ${catM} == "SE3_DEL" ]; then
      catN="LINKED_del3INVdup5"
    elif [ ${catP} == "SE5_DUP" ] && [ ${catM} == "SE3_DUP" ]; then
      catN="LINKED_dup3INVdup5"
    elif [ ${catP} == "SE5_DEL" ] && [ ${catM} != "SE3_DEL" ] && [ ${catM} != "SE3_DUP" ]; then
      catN="LINKED_INVdel5_OTHER"
    elif [ ${catP} == "SE5_DUP" ] && [ ${catM} != "SE3_DEL" ] && [ ${catM} != "SE3_DUP" ]; then
      catN="LINKED_INVdup5_OTHER"
    elif [ ${catM} == "SE3_DEL" ] && [ ${catP} != "SE5_DEL" ] && [ ${catP} != "SE5_DUP" ]; then
      catN="LINKED_del3INV_OTHER"
    elif [ ${catM} == "SE3_DUP" ] && [ ${catP} != "SE5_DEL" ] && [ ${catP} != "SE5_DUP" ]; then
      catN="LINKED_dup3INV_OTHER"    
    else
      catN="LINKED_OTHER"
    fi
    #Write original SE lines to tmp file
    fgrep -w ${eIDP} ${OUTDIR}/inversion.classifications.list > ${LINKED}
    fgrep -w ${eIDM} ${OUTDIR}/inversion.classifications.list >> ${LINKED}
    #Gather necessary additional info
    nDelP=$( fgrep -w ${eIDP} ${LINKED} | cut -f11 )
    predLossP=$( fgrep -w ${eIDP} ${LINKED} | cut -f12 )
    nDupP=$( fgrep -w ${eIDP} ${LINKED} | cut -f13 )
    predGainP=$( fgrep -w ${eIDP} ${LINKED} | cut -f14 )
    nDelM=$( fgrep -w ${eIDM} ${LINKED} | cut -f11 )
    predLossM=$( fgrep -w ${eIDM} ${LINKED} | cut -f12 )
    nDupM=$( fgrep -w ${eIDM} ${LINKED} | cut -f13 )
    predGainM=$( fgrep -w ${eIDM} ${LINKED} | cut -f14 )
    #Write new output
    echo -e "${eIDP},${eIDM}\t${chrP}\t${minP}\t${maxP}\t${minM}\t${maxM}\t${catN}\tNA\tNA\t${callP},${callM}\t${nDelP}/${nDelM}\t${predLossP}/${predLossM}\t${nDupP}/${nDupM}\t${predGainP}/${predGainM}" 
  done < ${SELINK} >> ${RES}

  #Clean up
  sort -nk2,2 -nk3,3 ${RES} > ${OUTDIR}/inversion.classifications.list
  rm -rf ${TMPCLST} ${SAMPS} ${cnDELS} ${cnDUPS} ${peDELS} ${peDUPS} ${preOVR} ${cnDEL_OVR} ${cnDUP_OVR} ${peDEL_OVR} ${peDUP_OVR} ${OVR_table} ${SEP} ${SEM} ${SELINK} ${LINKED} ${RES}

# if inappropriate input; displays useage
else
  echo -e "\n\nInversion classification script\n\nContact: Ryan Collins (rcollins@chgr.mgh.harvard.edu)\n\n"
  echo "Usage:"
  echo "  classify_inversion.sh [inv.events.bedpe] [inv.clusters.bedpe] [complex_size] [del.events.bedpe] [ins.events.bedpe] [cnMOPS_dels.bed] [cnMOPS_dups.bed] [OUTDIR]"
  echo ""
fi
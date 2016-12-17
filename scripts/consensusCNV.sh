#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Consensus linking with genotyping

#Consensus group descriptions:
#  A1: valid cluster, 100% cnMOPS concordance, genotyping confirm
#  A2: valid cluster, 100% cnMOPS concordance, genotyping reject -- MANUALLY CHECK ALL SITES
#  A3: valid cluster, imperfect (0%<x<100%) cnMOPS concordance, genotyping confirm
#  A4: valid cluster, imperfect (0%<x<100%) cnMOPS concordance, genotyping reject -- MANUALLY CHECK ALL SITES
#  B1: valid cluster, no cnMOPS concordance, genotyping confirm -- MANUALLY CHECK ALL SITES
#  B2: valid cluster, no cnMOPS concordance, genotyping reject
#  C1: cnMOPS call, no cluster, genotyping confirm
#  C2: cnMOPS call, no cluster, genotyping reject

#Group designation: consensus groups will have a "B" appended to the group if there is ≥30% overlap with CNV blacklisted sites (e.g. A1 vs A1B, etc)

#Genotyping confirmation: genotyping is considered a "confirmation" if at least half of samples in the prior call (either clustering or cnMOPS call) return a concordant genotype.
#  -Note that sites with >80% power to be resolved by a two-sample T-test of priors will by design have 100% confirmation or 0% confirmation. Underpowered sites are probably worth checking by hand.

#Note: males and females are genotyped independently on X, and females are not genotyped on Y (all sites automatically reported "-9", or "no data" encoding)

#Consensus group confidence levels:
# HIGH: A1, A2, A3, C1
# MED: A1B, A2B, A3B, A4, B1, C1B
# LOW: A4B, B1B, B2, B2B, C2, C2B

#Output format in tab-delimited columns, as follows:
# [1] chr
# [2] start
# [3] end
# [4] consensus ID
# [5] consensus group
# [6] confidence
# [7] number of samples
# [8] list of samples
# [9] cluster IDs
# [10] count of samples with cluster
# [11] list of samples with cluster
# [12] cnMOPS IDs
# [13] count of samples with cnMOPS
# [14] list of samples with cnMOPS
# [15] cnMOPS concordance vs. union prior (fraction)
# [16] genotyping result (PASS or FAIL)
# [17] genotyping concordance vs. union prior (fraction)
# [18] genotyping p-value (mean of per-sample p-values if two-sample t-test < 0.8 power)
# [19] genotyping power (will be set to 0 if cluster has single sample)
# [20] predicted homozygoes (deletions only, dups will be "." here)
# [21] flag for manual check, if necessary (values: CHECK or "." )

#Read input
events=$1               #full path to classifier output (events.bedpe format required)
cnMOPS_raw=$2           #full path to merged cnMOPS file (mergebeds.sh format required)
cnvtype=$3              #del or dup
params=$4               #params file from liWGS-SV

#Source params file
. ${params}

#Set other params
val_clust=`mktemp`
cnMOPS=`mktemp`
overlap=`mktemp`
cnMOPS_to_remove=`mktemp`
samp_union=`mktemp`
samp_cnMOPS=`mktemp`
samp_clust=`mktemp`
gA_tmp=${WRKDIR}/consensusCNV/${cnvtype}.gA.prelim.bed
gB_tmp=${WRKDIR}/consensusCNV/${cnvtype}.gB.prelim.bed
gC_tmp=${WRKDIR}/consensusCNV/${cnvtype}.gC.prelim.bed
preOut=`mktemp`

#Convert classifier output to extended bed file
case ${cnvtype} in
  del)
    fgrep Valid ${events} | awk -v OFS="\t" '{ print $1, $3, $5, $7, $20 }' | awk -v OFS="\t" '{ if ($3-$2<1) print $1, $3-1, $3, $4, $5, $6; else print $0 }' | sed 's/[\t]+/\t/g' > ${val_clust}
    ;;
  dup)
    fgrep Valid ${events} | awk -v OFS="\t" '{ print $1, $2, $6, $7, $20 }' | awk -v OFS="\t" '{ if ($3-$2<1) print $1, $3-1, $3, $4, $5, $6; else print $0 }' | sed 's/[\t]+/\t/g' > ${val_clust}
    ;;
esac

#Convert merged cnMOPS to extended bed file
cut -f1-4,8 ${cnMOPS_raw} | sed 's/[\t]+/\t/g' > ${cnMOPS}

#Get max insert size in cohort; will automatically force cluster-based coordinates to be at least this distance apart if smaller
max_insert=$( fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | cut -f9 | sort -nrk1,1 | head -n1 )

#Get number of samples in cohort
nsamp_all=$( cat ${samples_list} | wc -l )

#Get sex indexes; used for genotyping X and Y
Midx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/males.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )
Fidx=$( echo "$( awk -v OFS="\t" '{ print NR, $1 }' ${samples_list} | fgrep -wf ${WRKDIR}/females.list | awk '{ print ($1)+3 }' )" | cat <( echo -e "1\n2\n3" ) - | paste -s -d, )

#Get counts of males and females; used for genotyping X and Y
nmale=$( cat ${WRKDIR}/males.list | wc -l )
nfem=$( cat ${WRKDIR}/females.list | wc -l )

#Screen for clusters with ≥ 50% cnMOPS reciprocal overlap: groups A & B
while read chr start end cID samples; do
  #correct coordinates if necessary
  if [ $( echo "scale=0;(( ${end}-${start} ))" | bc ) -lt ${max_insert} ]; then
    mid=$( echo "scale=0;(( ( ${end}+${start} )/2 ))" | bc | cut -f1 -d. )
    if [ ${mid} -gt $( echo "scale=0;(( ${max_insert}/2 ))" | bc ) ]; then
      start=$( echo "scale=0;(( ${mid} - (${max_insert}/2) ))" | bc )
      end=$( echo "scale=0;(( ${mid} + (${max_insert}/2) ))" | bc )
    else
      start=0
      end=${max_insert}
    fi
  fi
  echo ${samples} | tr -d "[]" | sed 's/,/\n/g' | sort | uniq > ${samp_clust}
  bedtools intersect -wa -wb -r -f 0.5 -a <( echo -e "${chr}\t${start}\t${end}" ) -b ${cnMOPS} | fgrep -wf ${samp_clust} > ${overlap}
  #Group A (>0% cnMOPS concordance)
  if [ $( cat ${overlap} | wc -l ) -gt 0 ]; then
    #add cnMOPS IDs to list of calls to exclude from group C
    cut -f7 ${overlap} >> ${cnMOPS_to_remove}
    #generate ancillary info/metrics
    cut -f8 ${overlap} | tr -d "[]" | sed 's/,/\n/g' > ${samp_cnMOPS}
    cat ${samp_clust} ${samp_cnMOPS} | sort | uniq > ${samp_union}
    cnMOPS_concordance=$( echo "scale=2;(( $( fgrep -wf ${samp_cnMOPS} ${samp_clust} | wc -l ) / $( cat ${samp_clust} | wc -l ) ))" | bc )
    #report event with clustering coordinates but union of all samples between cnMOPS/clustering, and also report cnMOPS concordance, how many samples had clustering, and how many samples had cnMOPS
    #temp file format: chr, start, stop, placeholder for ID, count of union samples, union samples, cluster ID, count of cluster samples, cluster samples, cnMOPS ID(s), count of cnMOPS samples, cnMOPS samples (union), cnMOPS concordance
    echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_union} | wc -l )\t[$( paste -s -d, ${samp_union} )]\t${cID}\t$( cat ${samp_clust} | wc -l )\t[$( paste -s -d, ${samp_clust} )]\t$( cut -f7 ${overlap} | paste -s -d, )\t$( cat ${samp_cnMOPS} | wc -l )\t[$( paste -s -d, ${samp_cnMOPS} )]\t${cnMOPS_concordance}" >> ${gA_tmp}
  #Group B (0% cnMOPS concordance)
  else
    #report event with clustering coordinates and which samples had clustering
    #temp file format: chr, start, stop, placeholder for ID, count of union samples, union samples, cluster ID, count of cluster samples, cluster samples, cnMOPS ID(s), count of cnMOPS samples, cnMOPS samples (union), cnMOPS concordance
    echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_clust} | wc -l )\t${samples}\t${cID}\t$( cat ${samp_clust} | wc -l )\t${samples}\t.\t0\t.\t0.00" >> ${gB_tmp}
  fi
done < ${val_clust}

#Remove intersected cnMOPS IDs & print all cnMOPS calls with no clustering to gC_tmp. Also restricts to ≥ 10kb to speed up genotyping
while read chr start end cID samples; do
  if [ $( echo "${end}-${start}" | bc ) -ge 10000 ]; then
    #report event with cnMOPS coordinates and which samples had cnMOPS
    #temp file format: chr, start, stop, placeholder for ID, count of union samples, union samples, cluster ID, count of cluster samples, cluster samples, cnMOPS ID(s), count of cnMOPS samples, cnMOPS samples (union), cnMOPS concordance
    echo ${samples} | tr -d "[]" | sed 's/,/\n/g' | sort | uniq > ${samp_cnMOPS}
    echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t.\t0\t.\t${cID}\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t1.00" >> ${gC_tmp}
  fi
done < <( fgrep -wvf ${cnMOPS_to_remove} ${cnMOPS} )

#Assign temporary group-based CNV IDs to all intervals pre-genotyping
awk -v OFS="\t" -v class=${cnvtype} '{ $4="GroupA_"class"_"NR; print }' ${gA_tmp} > ${gA_tmp}2
mv ${gA_tmp}2 ${gA_tmp}
awk -v OFS="\t" -v class=${cnvtype} '{ $4="GroupB_"class"_"NR; print }' ${gB_tmp} > ${gB_tmp}2
mv ${gB_tmp}2 ${gB_tmp}
awk -v OFS="\t" -v class=${cnvtype} '{ $4="GroupC_"class"_"NR; print }' ${gC_tmp} > ${gC_tmp}2
mv ${gC_tmp}2 ${gC_tmp}

#Merge all groups, filter coverage matrix on N-masked reference regions, and genotype all CNV intervals
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots
#Genotype autosomes
for contig in $( seq 1 22 ); do
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  cat ${gA_tmp} ${gB_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed
  cut -f1-4,6 ${gC_tmp} | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed
  head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  awk -v OFS="\t" -v contig=${contig} '{ if ($1==contig) print $0 }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  if [ ${cnvtype} == "del" ]; then
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
  else
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE FALSE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE FALSE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
  fi
done
#Genotype chrX by sex in both M and F
for contig in X; do
  if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig} ]; then
    rm -r ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  fi
  if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig} ]; then
    rm -r ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  fi
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  for sex in M F; do
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}
  done
  awk -v chr=${contig} '{ if ($1==chr) print }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} | cat <( head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ) - > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Midx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Fidx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cat ${gA_tmp} ${gB_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed
  cut -f1-4,6 ${gC_tmp} | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed
  #split genotyping sites by sex
  for grouping in C notC; do
    while read chr start end ID samples; do
      #if > 1 sample is male, add to male list and remove female samples from priors
      if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/males.list - | wc -l ) -gt 0 ]; then
        newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/females.list - | paste -d, -s )
        echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed
      fi
      #if > 1 sample is female, add to female list and remove male samples from priors
      if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/females.list - | wc -l ) -gt 0 ]; then
        newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/males.list - | paste -d, -s )
        echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed
      fi
    done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed
  done
  # Submit genotyping
  for sex in M F; do
    if [ ${cnvtype} == "del" ]; then
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
    else
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE FALSE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE FALSE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
    fi
  done
done
#Genotype Y only in males
for contig in Y; do
  if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig} ]; then
    rm -r ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  fi
  if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig} ]; then
    rm -r ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  fi
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  for sex in M; do
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}
  done
  awk -v chr=${contig} '{ if ($1==chr) print }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} | cat <( head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ) - > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Midx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cat ${gA_tmp} ${gB_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed
  cut -f1-4,6 ${gC_tmp} | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed
  #split genotyping sites by sex
  for grouping in C notC; do
    while read chr start end ID samples; do
      #if > 1 sample is male, add to male list and remove female samples from priors
      if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/males.list - | wc -l ) -gt 0 ]; then
        newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/females.list - | paste -d, -s )
        echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed
      fi
    done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed
  done
  #Submit genotyping
  for sex in M; do
    if [ ${cnvtype} == "del" ]; then
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
    else
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.notC.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${contig}.bed TRUE Z TRUE FALSE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.notC.err"
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.C.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.C.${contig}.bed TRUE Z TRUE FALSE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.C.err"
    fi
  done
  #Create fake genotyping matrix for Y in females (all -9, "no data" encoding)
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F
  for grouping in notC C; do
    echo -e "CNV_ID\tchr\tstart\tend\t$( paste -s ${WRKDIR}/females.list )" > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed
    while read chr start end ID samples; do
      echo -e "${ID}\t${chr}\t${start}\t${end}"
      perl -E "say \"-9\t\" x ${nfem}"
    done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed | paste - - >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed
  done
done

#Gate until complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_${cnvtype}_genotypeCNV" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${COHORT_ID}_${cnvtype}_genotypeCNV" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Waiting on ${GATEcount} jobs..."
    GATEwait=0
  fi
done

#Merge & complete allosome genotyping
for contig in X Y; do
  for grouping in notC C; do
    if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.${grouping}.log ]; then
      rm ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.${grouping}.log
    fi
    paste <( head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed ) <( head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed | cut -f5- ) > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed
    while read chr start end ID samples; do
      echo -e "${ID}\t${chr}\t${start}\t${end}"
      #Set all males to "-9" genotype if site not genotyped in males
      if [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed | wc -l ) -eq 0 ]; then
        perl -E "say \"-9\t\" x ${nmale}"
      else
        fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed | cut -f5-
      fi
      #Set all females to "-9" genotype if site not genotyped in females
      if [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed | wc -l ) -eq 0 ]; then
        perl -E "say \"-9\t\" x ${nfem}"
      else
        fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed | cut -f5-
      fi   
    done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${grouping}.${contig}.bed | paste - - - | sed 's/[\t]+/\t/g' >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${grouping}.${contig}.bed
    #cat both M and F logfiles
    cat ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}.genotyping.${grouping}.log ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}.genotyping.${grouping}.log > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.${grouping}.log
  done
done

#Overlap group A with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
    gConcordance=$( echo "scale=2;(( $( cat ${overlap} | wc -l ) / $( cat ${samp_union} | wc -l ) ))" | bc )
  else
    gConcordance="0"
  fi
  #Determine genotyping outcome
  if [ $( echo "${gConcordance} >= 0.5" | bc ) -eq 1 ]; then
    gRes="PASS"
  else
    gRes="FAIL"
  fi
  #Get predicted homozygous deletions, print "CHECK" if any found at duplication allele
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep "Data" | wc -l ) -gt 0 ]; then
    #no homozygotes called with insufficient data
    homo="."
  else
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep "[" | awk '{ print $NF }' | tr -d "[]" | sed 's/,/\n/g' | sed '/^$/d' | paste -s -d, | awk '{ print "["$1"]" }' )
    if [ "${homo}" == "[]" ]; then
      homo="."
    fi
    if [ ${cnvtype} == "dup" ] && ! [ ${homo} == "." ]; then
      homo="CHECK"
    fi
  fi
  #reset genotyping metadata
  unset pval
  unset power
  if [ $( cat ${samp_union} | wc -l ) -ge ${nsamp_all} ] || [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w "No Data" | wc -l ) -gt 0 ]; then
    #genotyping is impossible if call appears in all samples or if call has "no data"
    pval="."
    power="."
  elif [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | wc -l ) -gt 0 ] && [ $( cat ${samp_union} | wc -l ) -gt 1 ]; then
    #Get permuted p-value if power >= 0.8
    power=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | sed -e 's/median\ power\:/\n/g' -e 's/median\ perm\:/GREPTHIS\n/g' | fgrep -w GREPTHIS | awk '{ print $1 }' | sort -nrk1,1 | head -n1 )
    if [ $( echo "${power} >= 0.8" | bc ) -eq 1 ]; then
      pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | sed -e 's/median\ perm\:/\n/g' -e 's/separation=/GREPTHIS\n/g' | fgrep -w GREPTHIS | awk '{ print $1 }' | sort -nk1,1 | head -n1 )
    #Otherwise get median pval of all samples
    else
      pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
    fi
  else
    #parse genotyping output with 1 sample
    power="0"
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w Cannot | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  fi
  #print pvalue and power errors if not previously set
  if [ -z ${pval} ]; then
    pval="ERR"
  fi
  if [ -z ${power} ]; then
    power="ERR"
  fi
  #Groups A1 & A2 have 100% cnMOPS concordance
  if [ $( echo "${cnMOPS_consensus} == 1" | bc ) -eq 1 ]; then
    #Group A1 has genotype confirmation
    if [ ${gRes} == "PASS" ]; then
      #Group A1B has ≥ 30% blacklist overlap
      if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
        group="A1B"
        conf="MED"
        gCheck="."
      #Group A1 has < 30% blacklist overlap
      else
        group="A1"
        conf="HIGH"
        gCheck="."
        #mark to check by hand if not all samples are genotyped the same
        if [ $( echo "${gConcordance} < 1" | bc ) -eq 1 ]; then
          gCheck="CHECK"
        else
          gCheck="."
        fi
      fi
    #Group A2 has genotype failure
    else
      #Group A2B has ≥ 30% blacklist overlap
      if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
        group="A2B"
        conf="MED"
        gCheck="."
      #Group A2 has < 30% blacklist overlap
      else
        group="A2"
        conf="HIGH"
        gCheck="CHECK"
      fi
    fi
  #Groups A3 & A4 have >0% cnMOPS concordance but also < 100% cnMOPS concordance
  else
    #Group A3 has genotype confirmation
    if [ ${gRes} == "PASS" ]; then
      #Group A3B has ≥ 30% blacklist overlap
      if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
        group="A3B"
        conf="MED"
        gCheck="."
      #Group A3 has < 30% blacklist overlap
      else
        group="A3"
        conf="HIGH"
        gCheck="."
        #mark to check by hand if not all samples are genotyped the same
        if [ $( echo "${gConcordance} < 1" | bc ) -eq 1 ]; then
          gCheck="CHECK"
        else
          gCheck="."
        fi
      fi
    #Group A4 has genotype failure
    else
      #Group A4B has ≥ 30% blacklist overlap
      if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
        group="A4B"
        conf="LOW"
        gCheck="."
      #Group A4 has < 30% blacklist overlap
      else
        group="A4"
        conf="MED"
        gCheck="CHECK"
      fi
    fi
  fi
  echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
done < ${gA_tmp} > ${preOut}A

#Overlap group B with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.notC.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
    gConcordance=$( echo "scale=2;(( $( cat ${overlap} | wc -l ) / $( cat ${samp_union} | wc -l ) ))" | bc )
  else
    gConcordance="0"
  fi
  #Determine genotyping outcome
  if [ $( echo "${gConcordance} >= 0.5" | bc ) -eq 1 ]; then
    gRes="PASS"
  else
    gRes="FAIL"
  fi
  #Get predicted homozygous deletions, print "CHECK" if any found at duplication allele
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep "Data" | wc -l ) -gt 0 ]; then
    #no homozygotes called with insufficient data
    homo="."
  else
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep "[" | awk '{ print $NF }' | tr -d "[]" | sed 's/,/\n/g' | sed '/^$/d' | paste -s -d, | awk '{ print "["$1"]" }' )
    if [ "${homo}" == "[]" ]; then
      homo="."
    fi
    if [ ${cnvtype} == "dup" ] && ! [ ${homo} == "." ]; then
      homo="CHECK"
    fi
  fi
#reset genotyping metadata
  unset pval
  unset power
  if [ $( cat ${samp_union} | wc -l ) -ge ${nsamp_all} ] || [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w "No Data" | wc -l ) -gt 0 ]; then
    #genotyping is impossible if call appears in all samples or if call has "no data"
    pval="."
    power="."
  elif [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | wc -l ) -gt 0 ] && [ $( cat ${samp_union} | wc -l ) -gt 1 ]; then
    #Get permuted p-value if power >= 0.8
    power=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | sed -e 's/median\ power\:/\n/g' -e 's/median\ perm\:/GREPTHIS\n/g' | fgrep -w GREPTHIS | awk '{ print $1 }' | sort -nrk1,1 | head -n1 )
    if [ $( echo "${power} >= 0.8" | bc ) -eq 1 ]; then
      pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | sed -e 's/median\ perm\:/\n/g' -e 's/separation=/GREPTHIS\n/g' | fgrep -w GREPTHIS | awk '{ print $1 }' | sort -nk1,1 | head -n1 )
    #Otherwise get median pval of all samples
    else
      pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w ncol | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
    fi
  else
    #parse genotyping output with 1 sample
    power="0"
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.notC.log | fgrep -w Cannot | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  fi
  #print pvalue and power errors if not previously set
  if [ -z ${pval} ]; then
    pval="ERR"
  fi
  if [ -z ${power} ]; then
    power="ERR"
  fi
  #group B1 is genotyping pass
  if [ ${gRes} == "PASS" ]; then
    #Group B1B has ≥ 30% blacklist overlap
    if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
      group="B1B"
      conf="LOW"
      gCheck="."
    #Group B1 has < 30% blacklist overlap
    else
      group="B1"
      conf="MED"
      gCheck="CHECK"
    fi
    echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
  #group B2 is genotyping fail
  else
     #Group B2B has ≥ 30% blacklist overlap
    if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
      group="B2B"
      conf="LOW"
      gCheck="."
    #Group B2 has < 30% blacklist overlap
    else
      group="B2"
      conf="LOW"
      gCheck="."
    fi
    echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
  fi
done < ${gB_tmp} > ${preOut}B

#Overlap group C with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.C.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.C.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}_Genotypes.C.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
    gConcordance=$( echo "scale=2;(( $( cat ${overlap} | wc -l ) / $( cat ${samp_union} | wc -l ) ))" | bc )
  else
    gConcordance="0"
  fi
  #Determine genotyping outcome -- different for group C than other groups, as only samples with concordant genotyping are kept in calls in C1
  if [ $( echo "${gConcordance} > 0" | bc ) -eq 1 ]; then
    gRes="PASS"
  else
    gRes="FAIL"
  fi
  #Get predicted homozygous deletions, print "CHECK" if any found at duplication allele
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.C.log | fgrep "Data" | wc -l ) -gt 0 ]; then
    #no homozygotes called with insufficient data
    homo="."
  else
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.C.log | fgrep "[" | awk '{ print $NF }' | tr -d "[]" | sed 's/,/\n/g' | sed '/^$/d' | paste -s -d, | awk '{ print "["$1"]" }' )
    if [ "${homo}" == "[]" ]; then
      homo="."
    fi
    if [ ${cnvtype} == "dup" ] && ! [ ${homo} == "." ]; then
      homo="CHECK"
    fi
  fi
  #parse genotyping output -- report median p among all samples with CN!=2 && CN!=-9, unless genotyping is CN=2 for all samples, in which case report the median p among all samples
  unset pval
  power="." #power is always NA for groupC
  if [ ${gRes} == "PASS" ]; then
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.C.log | grep -e 'ncol\|Cannot' | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3 && $1!=2 && $2!=-9) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  else
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/${COHORT_ID}.${cnvtype}.genotyping.C.log | grep -e 'ncol\|Cannot' | awk '{ print $NF }' | sed -e 's/\:/\t/g' -e 's/\#/\t/g' -e 's/\//\n/g' | sed '/^$/d' | awk '{ if (NF==3) print $3 }' | sort -nk1,1 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' )
  fi
  #print error if pval is unset
  if [ -z ${pval} ]; then
    pval="ERR"
  fi
  #group C1 is genotyping pass, but only samples with pass genotyping are retained in call (due to individual z-testing)
  if [ ${gRes} == "PASS" ]; then
    #Only report samples with GT pass
    samples_union=$( cat ${overlap} | paste -s -d, | awk '{ print "["$1"]" }' )
    count_union=$( cat ${overlap} | wc -l )
    #Group C1B has ≥ 30% blacklist overlap
    if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
      group="C1B"
      conf="MED"
      gCheck="."
    #Group C1 has < 30% blacklist overlap
    else
      group="C1"
      conf="HIGH"
      #mark to check by hand if not all samples are genotyped the same
      if [ $( echo "${gConcordance} < 1" | bc ) -eq 1 ]; then
        gCheck="CHECK"
      else
        gCheck="."
      fi
    fi
    if [ ${count_union} -gt 0 ]; then
      echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
    fi
  #group C2 is genotyping fail
  else
     #Group C2B has ≥ 30% blacklist overlap
    if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
      group="C2B"
      conf="LOW"
      gCheck="."
    #Group C2 has < 30% blacklist overlap
    else
      group="C2"
      conf="LOW"
      gCheck="CHECK"
    fi
    if [ ${count_union} -gt 0 ]; then
      echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t1.00\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
    fi
  fi
done < ${gC_tmp} > ${preOut}C

#Merge groups
cat ${preOut}A ${preOut}B ${preOut}C > ${preOut}

#Size Filters
#All A groups have no size filtering
awk -v OFS="\t" '$5 ~ /A/ { print $0 }' ${preOut} > ${preOut}2
#Group B1 has no size filtering
awk -v OFS="\t" '$5 ~ /B1/ { print $0 }' ${preOut} >> ${preOut}2
#Group B2 < 25kb
awk -v OFS="\t" '$5 ~ /B2/ { print $0 }' ${preOut} | awk '{ if ($3-$2<25000) print $0 }' >> ${preOut}2
#Group C1 ≥ 10kb
awk -v OFS="\t" '$5 ~ /C1/ { print $0 }' ${preOut} | awk '{ if ($3-$2>=10000) print $0 }' >> ${preOut}2
#Group C2 ≥ 50kb
awk -v OFS="\t" '$5 ~ /C2/ { print $0 }' ${preOut} | awk '{ if ($3-$2>=50000) print $0 }' >> ${preOut}2

#Write header to output file
echo -e "#Chr\tStart\tEnd\tID\tGroup\tConfidence\tObservations\tSamples\tCluster_IDs\tCluster_Observations\tCluster_Samples\tcnMOPS_IDs\tcnMOPS_Observations\tcnMOPS_Samples\tcnMOPS_Concordance\tGenotyping\tGenotyping_Concordance\tGenotyping_pVal\tGenotyping_Power\tHomozygotes\tManual_Check" > ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_${cnvtype}s.bed

#Sort pre-output file
sort -Vk1,1 -k2,2n -k3,3n ${preOut}2 >> ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_${cnvtype}s.bed

#Move prefilter file to WRKDIR
mv ${preOut} ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_${cnvtype}s.preFilter.bed

#Clean up
rm ${val_clust} ${cnMOPS} ${overlap} ${cnMOPS_to_remove} ${samp_union} ${samp_cnMOPS} ${samp_clust}
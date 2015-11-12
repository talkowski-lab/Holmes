#!/bin/bash

#################################
#     Consensus CNV Pipeline    #
#          w/ Genotyping        #
#             Mk III            #
#          Talkowski Lab        #
#################################

#Consensus group descriptions:
#  A1: valid cluster, 100% cnMOPS concordance, genotyping confirm
#  A2: valid cluster, 100% cnMOPS concordance, genotyping reject -- MANUALLY CHECK ALL SITES
#  A3: valid cluster, imperfect (0%<x<100%) cnMOPS concordance, genotyping confirm
#  A4: valid cluster, imperfect (0%<x<100%) cnMOPS concordance, genotyping reject -- MANUALLY CHECK ALL SITES
#  B1: valid cluster, no cnMOPS concordance, genotyping confirm -- MANUALLY CHECK ALL SITES
#  B2: valid cluster, no cnMOPS concordance, genotyping reject
#  C1: cnMOPS call, no cluster, genotyping confirm

#Group designation: consensus groups will have a "B" appended to the group if there is ≥30% overlap with CNV blacklisted sites (e.g. A1 vs A1B, etc)

#Genotyping confirmation: genotyping is considered a "confirmation" if at least half of samples in the prior call (either clustering or cnMOPS call) return a concordant genotype.
#  -Note that sites with >80% power to be resolved by a two-sample T-test of priors will by design have 100% confirmation or 0% confirmation. Underpowered sites are probably worth checking by hand.

#Note: males and females are genotyped independently on X, and females are not genotyped on Y (all sites automatically reported "-9", or "no data" encoding)

#Consensus group confidence levels:
# HIGH: A1, A2, A3, C1
# MED: A1B, A2B, A3B, A4, B1, C1B, C2
# LOW: A4B, B1B, B2, B2B, C2B

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
    echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t.\t0\t.\t${cID}\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t." >> ${gC_tmp}
  fi
done < <( fgrep -wvf ${cnMOPS_to_remove} ${cnMOPS} )

#Assign temporary group-based CNV IDs to all intervals pre-genotyping
awk -v OFS="\t" '{ $4="GroupA_"NR; print }' ${gA_tmp} > ${gA_tmp}2
mv ${gA_tmp}2 ${gA_tmp}
awk -v OFS="\t" '{ $4="GroupB_"NR; print }' ${gB_tmp} > ${gB_tmp}2
mv ${gB_tmp}2 ${gB_tmp}
awk -v OFS="\t" '{ $4="GroupC_"NR; print }' ${gC_tmp} > ${gC_tmp}2
mv ${gC_tmp}2 ${gC_tmp}

#Merge all groups, filter coverage matrix on N-masked reference regions, and genotype all CNV intervals
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit
mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots
#Genotype autosomes
for contig in $( seq 1 22 ); do
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  cat ${gA_tmp} ${gB_tmp} ${gC_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  awk -v OFS="\t" -v contig=${contig} '{ if ($1==contig) print $0 }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  if [ ${cnvtype} == "del" ]; then
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.err"
  else
    bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.err"
  fi
done
#Genotype chrX by sex in both M and F
for contig in X; do
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  for sex in M F; do
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}
  done
  awk -v chr=${contig} '{ if ($1==chr) print }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} | cat <( head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ) - > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Midx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Fidx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cat ${gA_tmp} ${gB_tmp} ${gC_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  #split genotyping sites by sex
  while read chr start end ID samples; do
    #if > 1 sample is male, add to male list and remove female samples from priors
    if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/males.list - | wc -l ) -gt 0 ]; then
      newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/females.list - | paste -d, -s )
      echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
    fi
    #if > 1 sample is female, add to female list and remove male samples from priors
    if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/females.list - | wc -l ) -gt 0 ]; then
      newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/males.list - | paste -d, -s )
      echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
    fi
  done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  # Submit genotyping
  for sex in M F; do
    if [ ${cnvtype} == "del" ]; then
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.err"
    else
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.err"
    fi
  done
done
#Genotype Y only in males
for contig in Y; do
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  for sex in M; do
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}
    mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}
  done
  awk -v chr=${contig} '{ if ($1==chr) print }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} | cat <( head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed ) - > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cut -f ${Midx} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  cat ${gA_tmp} ${gB_tmp} ${gC_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  #split genotyping sites by sex
  while read chr start end ID samples; do
    #if > 1 sample is male, add to male list and remove female samples from priors
    if [ $( echo "${samples}" | sed 's/,/\n/g' | fgrep -wf ${WRKDIR}/males.list - | wc -l ) -gt 0 ]; then
      newsamp=$( echo "${samples}" | sed 's/,/\n/g' | fgrep -wvf ${WRKDIR}/females.list - | paste -d, -s )
      echo -e "${chr}\t${start}\t${end}\t${ID}\t${newsamp}" >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
    fi
  done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  #Submit genotyping
  for sex in M; do
    if [ ${cnvtype} == "del" ]; then
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.err"
    else
      bsub -q big -R 'rusage[mem=20000]' -M 20000 -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}/${sex}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${sex}/${COHORT_ID}.${cnvtype}.genotyping.err"
    fi
  done
  #Create fake genotyping matrix for Y in females (all -9, "no data" encoding)
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F
  echo -e "CNV_ID\tchr\tstart\tend\t$( paste -s ${WRKDIR}/females.list )" > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed
  while read chr start end ID samples; do
    echo -e "${ID}\t${chr}\t${start}\t${end}"
    perl -E "say \"-9\t\" x ${nfem}"
  done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed | paste - - >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed
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
  if [ -e ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log ]; then
    rm ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log
  fi
  paste <( head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed ) <( head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed | cut -f5- ) > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed
  while read chr start end ID samples; do
    echo -e "${ID}\t${chr}\t${start}\t${end}"
    #Set all males to "-9" genotype if site not genotyped in males
    if [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed | wc -l ) -eq 0 ]; then
      perl -E "say \"-9\t\" x ${nmale}"
    else
      fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed | cut -f5-
    fi
    #Set all females to "-9" genotype if site not genotyped in females
    if [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed | wc -l ) -eq 0 ]; then
      perl -E "say \"-9\t\" x ${nfem}"
    else
      fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed | cut -f5-
    fi
    #Decide which logfile to use (for pval & power)
    if [ ${contig} != "Y" ] && [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}.genotyping.log | wc -l ) -gt 0 ]; then
      if [ $( fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}.genotyping.log | wc -l ) -gt 0 ]; then
        if [ ${nfem} -ge ${nmale} ]; then
          fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}.genotyping.log >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log
        else
          fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}.genotyping.log >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log
        fi
      else
        fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/F/${COHORT_ID}.${cnvtype}.genotyping.log >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log
      fi
    else
      fgrep -w ${ID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/M/${COHORT_ID}.${cnvtype}.genotyping.log >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log
    fi      
  done < ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed | paste - - - | sed 's/[\t]+/\t/g' >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed
done

#Overlap group A with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
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
  #Get predicted homozygoes (deletions only)
  if [ ${cnvtype} == "del" ] && [ ${gRes} == "PASS" ]; then
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep "[" | awk '{ print $NF }' | sed 's/,\]/\]/g' )
    if [ ${homo} == "[]" ]; then
      homo="."
    fi
  else
    homo="."
  fi
  #parse genotyping output with >2 samples
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | wc -l ) -gt 0 ]; then
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ p\:/\n/g' | tail -n1 | awk '{ print $1 }' )
    power=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ power\:/\n/g' | tail -n1 | awk '{ print $1 }' )
  #parse genotyping output with 1 sample
  elif [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w "No Data" | wc -l ) -gt 0 ]; then
    pval="."
    power="."
  else
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w Cannot | awk '{ print $NF }' | cut -d\# -f2 | tr -d "/" )
    if [ -z ${pval} ]; then
      pval="."
    fi
    power=0
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
done < ${gA_tmp} > ${preOut}

#Overlap group B with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
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
  #Get predicted homozygoes (deletions only)
  if [ ${cnvtype} == "del" ] && [ ${gRes} == "PASS" ]; then
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep "[" | awk '{ print $NF }' | sed 's/,\]/\]/g' )
    if [ ${homo} == "[]" ]; then
      homo="."
    fi
  else
    homo="."
  fi
  #parse genotyping output with >2 samples
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | wc -l ) -gt 0 ]; then
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ p\:/\n/g' | tail -n1 | awk '{ print $1 }' )
    power=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ power\:/\n/g' | tail -n1 | awk '{ print $1 }' )
  #parse genotyping output with 1 sample
  elif [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w "No Data" | wc -l ) -gt 0 ]; then
    pval="."
    power="."
  else
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w Cannot | awk '{ print $NF }' | cut -d\# -f2 | tr -d "/" )
    if [ -z ${pval} ]; then
      pval="."
    fi
    power=0
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
  #group B2 is genotyping fail and limited to < 25kb
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
    if [ $( echo "${end}-${start}" | bc ) -le 25000 ]; then
      echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
    fi
  fi
done < ${gB_tmp} >> ${preOut}

#Overlap group C with genotyping
while read chr start end eID count_union samples_union cluster_ID count_cluster samples_cluster cnMOPS_ID count_cnMOPS samples_cnMOPS cnMOPS_consensus; do
  echo "${samples_union}" | tr -d "[]" | sed 's/,/\n/g' > ${samp_union}
  #calculate genotyping concordance
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | wc -l ) -gt 0 ]; then
    head -n1 ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | cut -f$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}_Genotypes.${chr}.bed | sed 's/\t/\n/g' | sed '1,4d' | awk '{ if ($1!=2 && $1!=-9) print NR+4 }' | paste -s -d, ) | sed 's/\t/\n/g' | fgrep -wf - ${samp_union} > ${overlap}
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
  #Get predicted homozygoes (deletions only)
  if [ ${cnvtype} == "del" ] && [ ${gRes} == "PASS" ]; then
    homo=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep "[" | awk '{ print $NF }' | sed 's/,\]/\]/g' )
    if [ ${homo} == "[]" ]; then
      homo="."
    fi
  else
    homo="."
  fi
  #parse genotyping output with >2 samples
  if [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | wc -l ) -gt 0 ]; then
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ p\:/\n/g' | tail -n1 | awk '{ print $1 }' )
    power=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w ncol | sed 's/median\ power\:/\n/g' | tail -n1 | awk '{ print $1 }' )
  #parse genotyping output with 1 sample
  elif [ $( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w "No Data" | wc -l ) -gt 0 ]; then
    pval="."
    power="."
  else
    pval=$( fgrep -w ${eID} ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${chr}/SFARIEcoV5.${cnvtype}.genotyping.log | fgrep -w Cannot | awk '{ print $NF }' | cut -d\# -f2 | tr -d "/" )
    if [ -z ${pval} ]; then
      pval="."
    fi
    power=0
  fi
  #group C1 is genotyping pass and restricted to ≥ 10kb
  if [ ${gRes} == "PASS" ]; then
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
    if [ $( echo "${end}-${start}" | bc ) -ge 10000 ]; then
      echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
    fi
  #group C2 is genotyping fail and restricted to ≥ 50kb
  else
     #Group C2B has ≥ 30% blacklist overlap
    if [ $( echo "$( bedtools coverage -a ${CNV_BLACKLIST} -b <( echo -e "${chr}\t${start}\t${end}" ) | awk '{ print $NF }' ) >= 0.3" | bc ) -eq 1 ]; then
      group="C2B"
      conf="LOW"
      gCheck="."
    #Group C2 has < 30% blacklist overlap
    else
      group="C2"
      conf="MED"
      gCheck="CHECK"
    fi
    if [ $( echo "${end}-${start}" | bc ) -ge 50000 ]; then
      echo -e "${chr}\t${start}\t${end}\t${eID}\t${group}\t${conf}\t${count_union}\t${samples_union}\t${cluster_ID}\t${count_cluster}\t${samples_cluster}\t${cnMOPS_ID}\t${count_cnMOPS}\t${samples_cnMOPS}\t${cnMOPS_consensus}\t${gRes}\t${gConcordance}\t${pval}\t${power}\t${homo}\t${gCheck}"
    fi
  fi
done < ${gC_tmp} >> ${preOut}

#Write header to output file
echo -e "#Chr\tStart\tEnd\tID\tGroup\tConfidence\tObservations\tSamples\tCluster_IDs\tCluster_Observations\tCluster_Samples\tcnMOPS_IDs\tcnMOPS_Observations\tcnMOPS_Samples\tcnMOPS_Concordance\tGenotyping\tGenotyping_Concordance\tGenotyping_pVal\tGenotyping_Power\tHomozygotes\tManual_Check" > ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_${cnvtype}s.bed

#Sort pre-output file
sort -Vk1,1 -k2,2n -k3,3n ${preOut} >> ${WRKDIR}/consensusCNV/${COHORT_ID}_consensus_${cnvtype}s.bed

#Clean up
rm ${val_clust} ${cnMOPS} ${overlap} ${cnMOPS_to_remove} ${samp_union} ${samp_cnMOPS} ${samp_clust} ${preOut}
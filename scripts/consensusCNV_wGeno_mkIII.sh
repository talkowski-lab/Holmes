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

#Consensus group confidence levels:
# HIGH: A1, A2, A3, C1
# MED: A1B, A2B, A3B, A4, B1, C1B, C2
# LOW: A4B, B1B, B2, B2B, C2B

#Output format: 18 tab-delimited columns
# [1] chr
# [2] start
# [3] end
# [4] consensus ID
# [5] consensus group
# [6] number of samples
# [7] list of samples
# [8] cluster IDs
# [9] count of samples with cluster
# [10] list of samples with cluster
# [11] cnMOPS IDs
# [12] count of samples with cnMOPS
# [13] list of samples with cnMOPS
# [14] cnMOPS concordance vs. clustering (fraction)
# [15] genotyping result (PASS or FAIL)
# [16] genotyping concordance vs. union prior (fraction)
# [17] genotyping p-value (mean of per-sample p-values if two-sample t-test < 0.8 power)
# [18] genotyping power (will be set to 0 if cluster has single sample)

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
gA_tmp=`mktemp`
gB_tmp=`mktemp`
gC_tmp=`mktemp`

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

#Screen for clusters with ≥ 50% cnMOPS reciprocal overlap: groups A & B
while read chr start end cID samples; do
  #correct coordinates if necessary
  if [ $( echo "scale=0;(( ${end}-${start} ))" | bc ) -lt ${max_insert} ]; then
    mid=$( echo "${end}-${start}" | bc )
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

#Remove intersected cnMOPS IDs & print all cnMOPS calls with no clustering to gC_tmp: group C
while read chr start end cID samples; do
  #report event with cnMOPS coordinates and which samples had cnMOPS
  #temp file format: chr, start, stop, placeholder for ID, count of union samples, union samples, cluster ID, count of cluster samples, cluster samples, cnMOPS ID(s), count of cnMOPS samples, cnMOPS samples (union), cnMOPS concordance
  echo ${samples} | tr -d "[]" | sed 's/,/\n/g' | sort | uniq > ${samp_cnMOPS}
  echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t.\t0\t.\t${cID}\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t." >> ${gC_tmp}
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
for contig in $( seq 1 20 ) X Y; do
  # mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}
  mkdir ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}
  # cat ${gA_tmp} ${gB_tmp} ${gC_tmp} | sort -k1,1V -k2,2n -k3,3n | cut -f1-4,6 | tr -d "[]" | awk -v contig=${contig} '{ if ($1==contig) print $0 }' > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed
  # head -n1 ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  # awk -v OFS="\t" -v contig=${contig} '{ if ($1==contig) print $0 }' ${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed | bedtools intersect -v -a - -b ${NMASK} >> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed
  if [ ${cnvtype} == "del" ]; then
    bsub -q normal -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0; Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE TRUE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.err"
  else
    bsub -q normal -sla miket_sc -J ${COHORT_ID}_${cnvtype}_genotypeCNV -o ${OUTDIR}/logs/genotyping.log -e ${OUTDIR}/logs/genotyping.log "cd ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/plots/${contig}; module load R/3.1.0 Rscript ${liWGS_SV}/scripts/genotypeCNV.R ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_toGenotype.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.physical.cov_matrix.${contig}.bed ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}_Genotypes.${contig}.bed TRUE Z TRUE FALSE > ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.log 2> ${WRKDIR}/consensusCNV/${cnvtype}_genotyping/chrsplit/${contig}/${COHORT_ID}.${cnvtype}.genotyping.err"
  fi
done



#Clean up
rm ${val_clust} ${cnMOPS} ${overlap} ${cnMOPS_to_remove} ${samp_union} ${samp_cnMOPS} ${samp_clust} ${gA_tmp} ${gB_tmp} ${gC_tmp}


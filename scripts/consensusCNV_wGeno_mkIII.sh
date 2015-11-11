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

#Input
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

#Screen for clusters with ≥ 50% cnMOPS reciprocal overlap: groups A & B
while read chr start end cID samples; do
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
  echo -e "${chr}\t${start}\t${end}\tblank\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t.\t0\t.\t.\t${cID}\t$( cat ${samp_cnMOPS} | wc -l )\t${samples}\t." >> ${gC_tmp}
done < <( fgrep -wvf ${cnMOPS_to_remove} ${cnMOPS} )

#Merge all groups & genotype intervals





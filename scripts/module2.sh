#!/bin/bash

#liWGS-SV Pipeline: Module 2 (Per-Sample Clustering)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Submit bamstat job per sample
bsub -q normal -o ${OUTDIR}/logs/bamstat.log -e ${OUTDIR}/logs/bamstat.log -sla miket_sc -J ${COHORT_ID}_bamstat "/data/talkowski/rlc47/code/dna-scripts/bamstat.sh -s 3 -i ${WRKDIR}/quad/SFARIQuad${lib}/SFARIQuad${lib}.merged.mkdup.nsort.bam -o ${WRKDIR}/quad/SFARIQuad${lib}/SV/bamstat_3cluster" 

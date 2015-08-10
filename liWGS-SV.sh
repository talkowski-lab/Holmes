#!/bin/bash

##MASTER liWGS-SV PIPELINE SCRIPT
#Contact: rcollins@chgr.mgh.harvard.edu

#Submit this script to long queue on ERISOne to run entire pipeline

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Set up output directory tree
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
mkdir ${OUTDIR}/QC
mkdir ${OUTDIR}/cohort
while read ID bam sex; do
  mkdir ${OUTDIR}/QC/${ID}
done < ${samples_list}
mkdir ${OUTDIR}/data
mkdir ${OUTDIR}/data/seqDepth
mkdir ${OUTDIR}/data/clusters
mkdir ${OUTDIR}/SV_calls

#Set up working directory tree
if ! [ -e ${WRKDIR} ]; then
  mkdir ${WRKDIR}
fi
while read ID bam sex; do
  mkdir ${WRKDIR}/${ID}
done < ${samples_list}
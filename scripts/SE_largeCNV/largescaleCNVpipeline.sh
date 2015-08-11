#!/bin/bash

#liWGS-SV Pipeline: SE Large CNV Pipeline
#August 2015
#Contact: Serkan Erdin (serdin@cmgh.harvard.edu)

#Read input
samples_list=$1
params=$2
ID=$3

#Source params file
. ${params}
cd ${WRKDIR}

#Load modules
module load intel_parallel_xe/xe
module load R/3.1.0
Rscript -e "if(\"DNAcopy\" %in% rownames(installed.packages()) == FALSE){source(\"http://bioconductor.org/biocLite.R\"); biocLite(\"DNAcopy\")}; suppressPackageStartupMessages(library(DNAcopy))"

#Runs
mkdir ${WRKDIR}/${ID}/DNAcopy
perl normalizeJLcoverage.pl ${WRKDIR}/${COHORT_ID}.query.bindata.txt ${DNAcopy_ref} ${ID} > ${WRKDIR}/${ID}/DNAcopy/${ID}.log2.txt
Rscript --no-environ JumpingLibrary.CNV.R ${WRKDIR}/${ID}/DNAcopy/${ID}.log2.txt ${ID} ${WRKDIR}/${ID}/DNAcopy
Rscript --no-environ JumpingLibrary.plotchromosomecoverage.R ${WRKDIR}/${ID}/DNAcopy/${ID}.log2.txt ${ID} ${WRKDIR}/${ID}/DNAcopy
perl ${liWGS_SV}/scripts/SE_largeCNV/mergeSegments.pl ${ID}.log2.segments.p_value --ref-arm-sizes ${liWGS_SV}/data/h37.chr.arm.bed --output-basename ${ID} --amp-threshold 0.3 --del-threshold -0.3
mv ${WRKDIR}/${ID}.SD.profile ${WRKDIR}/${ID}.log2.segments.p_value ${WRKDIR}/${ID}.events.tsv ${WRKDIR}/${ID}.summary.txt ${WRKDIR}/${ID}/DNAcopy/

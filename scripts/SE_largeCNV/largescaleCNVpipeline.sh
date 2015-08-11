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

#Sets variables
queryinsertcov=${WRKDIR}/iCov/${COHORT_ID}.physical.cov_matrix.bed
outdir="/data/talkowski/Samples/NF/LargeScaleCNV_JumpingLibrary"
refarmsize="chr.arm.hg19.txt"

mkdir ${WRKDIR}/${ID}/DNAcopy
perl normalizeJLcoverage.pl ${WRKDIR}/${COHORT_ID}.query.bindata.txt ${DNAcopy_ref} ${ID} > ${WRKDIR}/${ID}/DNAcopy/${sampleid}.log2.txt
Rscript --no-environ JumpingLibrary.CNV.R ${sampleid}.log2.txt ${sampleid} ${outdir}
perl ~/PerlScripts/mergeSegments.pl ${sampleid}.log2.segments.p_value --ref-arm-sizes ${refarmsize} --output-basename ${sampleid} --amp-threshold 0.3 --del-threshold -0.3
mkdir ${sampleid}
mv ${sampleid}.* ${sampleid}/
fi

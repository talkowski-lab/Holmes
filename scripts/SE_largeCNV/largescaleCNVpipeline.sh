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
queryinsertcov="/data/nf/MTa425_150316/Merged/NF2_cnMobs_human_filtered/NF2_human_filtered.insert.cov_matrix.bed"
referencefileid="Sfari-all"
queryfileid="Sfari-all-query"
outdir="/data/talkowski/Samples/NF/LargeScaleCNV_JumpingLibrary"
refarmsize="chr.arm.hg19.txt"

Rscript getwithinlibrarynorm_query.R ${queryinsertcov} ${queryfileid}  ${outdir} keep 

perl normalizeJLcoverage.pl ${queryfileid}.query.bindata.txt ${liWGS_SV}/data/DNAcopy_reference.bindata.bed ${sampleid} > ${sampleid}.log2.txt
Rscript --no-environ JumpingLibrary.CNV.R ${sampleid}.log2.txt ${sampleid} ${outdir}
Rscript --no-environ JumpingLibrary.plotchromosomecoverage.R ${sampleid}.log2.txt ${sampleid} ${outdir}
perl ~/PerlScripts/mergeSegments.pl ${sampleid}.log2.segments.p_value --ref-arm-sizes ${refarmsize} --output-basename ${sampleid} --amp-threshold 0.3 --del-threshold -0.3
mkdir ${sampleid}
mv ${sampleid}.* ${sampleid}/
fi

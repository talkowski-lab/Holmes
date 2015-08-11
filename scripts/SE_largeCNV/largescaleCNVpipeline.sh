#!/bin/bash

###############################################
#
# running read-depth based CNV pipeline
# Serkan Erdin (serdin@chgr.mgh.harvard.edu)
# Talkowski Lab
#
###############################################


### loading necessary modules
module load intel_parallel_xe/xe
module load R/3.1.0

#Checks for appropriate input
if [ $# -eq 1 ]; then


### packages needed
### DNAcopy


referenceinsertcov="/data/talkowski/Samples/SFARI/batch1/cnMOPS/SFARI.batch1.covMatrix.bed"
queryinsertcov="/data/nf/MTa425_150316/Merged/NF2_cnMobs_human_filtered/NF2_human_filtered.insert.cov_matrix.bed"
referencefileid="Sfari-all"
queryfileid="Sfari-all-query"
outdir="/data/talkowski/Samples/NF/LargeScaleCNV_JumpingLibrary"
refarmsize="chr.arm.hg19.txt"



sampleid=$1

#Rscript getwithinlibrarynorm_reference.R ${referenceinsertcov} ${referencefileid}  ${outdir} keep 
#Rscript getwithinlibrarynorm_query.R ${queryinsertcov} ${queryfileid}  ${outdir} keep 
perl normalizeJLcoverage.pl ${queryfileid}.query.bindata.txt ${referencefileid}.reference.bindata.txt ${sampleid} > ${sampleid}.log2.txt
Rscript --no-environ JumpingLibrary.CNV.R ${sampleid}.log2.txt ${sampleid} ${outdir}
Rscript --no-environ JumpingLibrary.plotchromosomecoverage.R ${sampleid}.log2.txt ${sampleid} ${outdir}
perl ~/PerlScripts/mergeSegments.pl ${sampleid}.log2.segments.p_value --ref-arm-sizes ${refarmsize} --output-basename ${sampleid} --amp-threshold 0.3 --del-threshold -0.3
mkdir ${sampleid}
mv ${sampleid}.* ${sampleid}/
fi

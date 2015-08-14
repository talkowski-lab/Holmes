#!/bin/bash

#liWGS-SV Pipeline: Example parameters file
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Copy this format exactly to create your own parameters file

#Run preferences
export COHORT_ID=test #cohort ID (string)
export OUTDIR=/scratch/miket/rlc47temp/liWGS-SV_${COHORT_ID} #output directory - only final outputs will be written here
export WRKDIR=/scratch/miket/rlc47temp/tmp.files/${COHORT_ID}_liWGS-SV #working directory - all temp files will be written here. May temporarily need lots of storage, so best to avoid writing to /tmp
export other_assign=MALE #sex group to assign "other" sex for depth-based calling. must be either "MALE" or "FEMALE"
export uscore_skip=1 #underscore skip; determines how many underscores in clustering read ID are assumed to be a part of sample ID; e.g. format like "STUDY_SAMPLE" should have uscore_skip=1, etc
export KEEP_TMP=TRUE #Boolean (all caps) to indicate whether working directory should be kept once pipeline is complete. Individual modules may be rerun with identical parameters file if working directory is saved.
export pre_bamstat=FALSE #boolean (all caps) to indicate if bamstat on all samples has already been run. If set as TRUE, bamstat_paths MUST BE CORRECTLY SPECIFIED otherwise the whole pipeline will fail
export GENOTYPE_OVERRIDE=FALSE #boolean (all caps) to indicate if you want to manually disable CNV genotyping. Any value other than "TRUE" will revert to default operation (genotyping used if cohort size >= ${min_geno}, set in mod. 6)

#Full paths for reference files & executables 
export liWGS_SV=/data/talkowski/rlc47/code/liWGS-SV #path to liWGS-SV git repo
export REF=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.fa #reference fasta
export DICT=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71_reord_bwa07/Ensembl_hgGRCh37_71_ERCC_reord.mainContigs.dict #reference dictionary, restricted to chromosomes where calls should be made
export CLASSIFIER_DIR=/data/talkowski/rlc47/code/classifier #classifier git repo
export PYCLUSTER_DIR=/data/talkowski/rlc47/code/pycluster #pycluster git repo
export PICARD=/data/talkowski/tools/bin/picard-tools-1.137/picard.jar #picard.jar executable
export sambamba=/data/talkowski/tools/bin/sambamba_v0.4.6 #sambamba executable
export DNAcopy_ref=/data/talkowski/rlc47/src/DNAcopy_reference.bindata.bed #DNAcopy reference profile
export CNV_BLACKLIST=/data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed #blacklist file, restrict CNV calling on >30% coverage of features in list 
export bamstat_paths=NA #file containing paths to pre-run bamstat directories. First column: ID, second column: full path to bamstat directory; tab-delimited. Will be ignored unless pre_bamstat="TRUE"
export abParts=/data/talkowski/rlc47/src/abParts.bed #path to antibody parts annotation file (available from UCSC), used for exclusion of putative complex site
export NMASK=/data/talkowski/rlc47/src/GRCh37_Nmask.bed #bed file corresponding to N-masked regions of reference genome
export refFlat=/data/talkowski/tools/bin/TGDB/BACKUP/hg19_refFlat_3_13_2014.bed
export GTF=/data/talkowski/tools/ref/Ensembl_hgGRCh37_71/Homo_sapiens_GRCh37_71_ERCC.gtf

#Update user paths
export PATH=${PATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}
export PYTHOPATH=${PYTHONPATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}
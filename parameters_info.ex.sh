#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Example parameters file
#Copy this format exactly for your own parameters file

###Run preferences###
#cohort ID (string)
export COHORT_ID=your_cohort_id_here 
#output directory - only final outputs will be written here
export OUTDIR=your_output_directory_here
#working directory - all temp files will be written here. May temporarily need 
#lots of storage, so best to avoid writing to /tmp on restricted systems
export WRKDIR=${OUTDIR}/TMPFILES
#sex group to assign "other" sex for depth-based calling. Must be either "MALE" 
#or "FEMALE"
export other_assign=FEMALE
#underscore skip; determines how many underscores in clustering read ID are 
#assumed to be a part of sample ID; e.g. format like "STUDY_SAMPLE" 
#should have uscore_skip=1, etc
export uscore_skip=1
#Boolean (all caps) to indicate whether working directory should be kept once 
#pipeline is complete. Individual modules may be rerun with identical parameters
#file if working directory is saved.
export KEEP_TMP=TRUE
#boolean (all caps) to indicate if bamstat on all samples has already been run. 
#If set as TRUE, bamstat_paths MUST BE CORRECTLY SPECIFIED otherwise the 
#whole pipeline will fail
export pre_bamstat=FALSE 
#boolean (all caps) to indicate if you want to manually disable CNV genotyping. 
#Any value other than "TRUE" will revert to default operation (genotyping used 
#if cohort size >= ${min_geno}, set in mod. 6)
export GENOTYPE_OVERRIDE=FALSE
#maximum MAF/VAF before variant is considered artifactual and/or reference variant
export polyArt_filter=0.5

###Full paths for reference files & executables###
#path to Holmes git repo
export liWGS_SV=/path/to/Holmes/
#reference fasta
export REF=/path/to/reference.fa
#reference dictionary, restricted to chromosomes where calls should be made
export DICT=/path/to/reference.fa.dict
#classifier git repo
export CLASSIFIER_DIR=${liWGS_SV}/classifier
#pycluster git repo
export PYCLUSTER_DIR=${liWGS_SV}/pycluster
#picard.jar executable
export PICARD=/path/to/picard.jar
#sambamba executable
export sambamba=/path/to/sambamba_v0.4.6
#blacklist file, restrict CNV calling on >30% coverage of features in list 
export CNV_BLACKLIST=/path/to/CNV_blacklist.bed
#file containing paths to pre-run bamstat directories. First column: ID, 
#second column: full path to bamstat directory; tab-delimited. Will be 
#ignored unless pre_bamstat="TRUE"
export bamstat_paths=${OUTDIR}/bamstat.list
#path to antibody parts annotation file (available from UCSC), 
#used for exclusion of putative complex site
export abParts=/data/talkowski/rlc47/src/abParts.bed
#bed file corresponding to N-masked regions of reference genome
export NMASK=/path/to/GRCh37_Nmask.bed
#UCSC refFlat (genes)
export refFlat=/path/to/refFlat.bed
#Gencode GTF
export GTF=/path/to/gencode.gtf

#Update user paths
export PATH=${PATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}
export PYTHONPATH=${PYTHONPATH}:${CLASSIFIER_DIR}:${PYCLUSTER_DIR}
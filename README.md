[![DOI](https://zenodo.org/badge/40485508.svg)](https://zenodo.org/badge/latestdoi/40485508)
# Holmes
Pipeline for structural variation detection from liWGS libraries  

Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski  
Massachusetts General Hospital, The Broad Institute, and Harvard Medical School
Contact: Ryan L. Collins <rlcollins@g.harvard.edu>  

If you use this software, please cite:
Collins RL, et al. *Defining the spectrum of large inversions, complex structural variation, and chromothripsis in the morbid genome.* In Press (2016)  

Code development credits: Ryan L. Collins, Harrison Brand, Matthew R. Stone, Vamsee Pillalamarri, Joseph T. Glessner, Claire Redin, Colby Chiang, Ian Blumenthal, Adrian Heilbut

### Table of Contents  
1. Overview  
2. Modules
  Module 1: Cohort QC
  Module 2: Physical Depth Analysis  
  Module 3: Per-Sample Read Pair Clustering  
  Module 4: Physical Depth-Based CNV Discovery  
  Module 5: Cohort-Wide Joint Read Pair Clustering & Classification
  Module 6: Consensus CNV Categorization  
  Module 7: Balanced & Complex SV Resolution  
  Module 8: Call Consolidation & Reformatting  
  Module 9: SV Annotation
3. Miscellanea & FAQs  

## Overview  
Holmes is a pipeline for SV discovery and annotation in cohorts (n>20-30) of human liWGS libraries. This tool can detect, resolve, and classify all known forms of SV, from canonical copy number variants (CNVs) to balanced rearrangements (e.g. inversions, reciprocal translocations) to more complex chromosomal rearrangements. Currently, this tool only supports human samples and the GRCh37 reference assembly.  

The pipeline is run in nine independent modules split into five sequential stages, as shown below:  


**Note: the codebase for this tool has hard-coded variable paths and other dependencies particular to the Partners Healthcare computing cluster. These dependencies are in the process of being resolved so the tool can be more easily deployed on other clusters. For now, please contact <rlcollins@g.harvard.edu> with any questions or issues.**  

Preliminary description:  

Execution command:  
runHolmes.sh samples.list parameters_info.sh  

##INPUT##
samples.list: three columns, tab delimited, col 1) sample ID, col 2) full path to sample bam, col 3) expected sex (M=XY, F=XX, O=other, U=unknown (will defer to predicted sex by sexcheck))  
parameters_info.sh: shell script to export all parameters for pipeline run  

##PRE-MODULE STEPS##
Symlinks & indexes all bams  
Creates working and output directory trees  
Loads necessary modules  

##MODULE 1: QC##
Runs the following:  
	Picard EstimateLibraryComplexity  
	Picard CollectAlignmentSummaryMetrics  
	Picard CollectInsertSizeMetrics  
	Picard CollectWgsMetrics  
	Samtools flagstat  
	Bamtools stats  
	Sex Check  
	WGS Dosage Bias Check  
Checks for nominal QC values, reports errors to ${OUTDIR}/${COHORT_ID}_WARNINGS.txt  
Writes master QC table to ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics  

##MODULE 2: PHYSICAL DEPTH ANALYSES##
Runs binCov to generate 1kb binned physical depth for each library  
BGZips & tabix indexes each coverage file (for classifier)  

##MODULE 3: PER-SAMPLE CLUSTERING##
**Rate-limiting step of entire pipeline**
If ${pre_bamstat} isn't set as "TRUE", bamstat is run at min cluster size = 3 for each sample  
If ${pre_bamstat}="TRUE", bamstat clusters and stats.file are copied from preexisting paths to ${WRKDIR}  
Removes *pairs.txt and *pairs.sorted.txt to save space  

##MODULE 4: PHYSICAL DEPTH CNV CALLING##
Runs cnMOPS on autosomes on all samples  
Runs cnMOPS on allosomes on samples split by M/F. "Other" sex samples pooled with either M or F depending on ${other_assign}  
Merges cnMOPS calls per sample  

##MODULE 5: JOINT RECLUSTERING & CLASSIFICATION##
Runs classifier  
Patches clusters  
Reclassifies patched clusters  
Applies final classification labels & sets coordinate reporting to be 1st or 3rd quartile of reads, respectively (to avoid overclustering/negative sizes)  

##MODULE 6: CONSENSUS CNV CALLING##
Runs in one of two modes: with or without genotyping information  
Mode chosen by parameter ${min_geno}, set in module6.sh, which corresponds to the minimum number of samples in the cohort to use genotyping  
Consensus Groups with Genotyping:  
	A [HIGH]: Valid cluster, cnMOPS or genotyping support, <30% blacklist  
	B [HIGH]: cnMOPS call, ≥50kb, <30% blacklist, genotyping pass, no clustering overlap  
	C [MED]: cnMOPS call, <50kb, genotyping pass, <30% blacklist  
	D [MED]: valid cluster, genotyping or cnMOPS support, ≥30% blacklist  
	E [MED]: cnMOPS call, ≥50kb, genotyping pass, ≥30% blacklist  
	F [LOW]: cnMOPS call, ≥50kb, no clustering support, no genotyping support  
	G [LOW]: cnMOPS call, <50kb, genotyping pass, ≥30% blacklist  
	H [LOW]: valid cluster, <25kb, no cnMOPS or genotyping support  
Consensus Groups without Genotyping:  
    A [HIGH]: Valid cluster, cnMOPS support, <30% blacklist  
    B [MED]: cnMOPS call, ≥50kb, <30% blacklist, no clustering overlap  
    C [MED]: valid cluster, cnMOPS support, ≥30% blacklist  
    D [LOW]: cnMOPS call, ≥50kb, ≥30% blacklist  
    E [LOW]: valid cluster, <25kb, no cnMOPS support  
Returns single merged file each for consensus dels and consensus dups  

##MODULE 7: COMPLEX SV CATEGORIZATION##
Runs inversion classification script  
Runs translocation classification script  
Runs complex linking script  
Runs complex parsing script  

##MODULE 8: VARIANT CONSOLIDATION & REFORMATTING
Outputs the following seven variant files:  
-Deletion  
-Duplication  
-Inversion  
-Insertion  
-Translocation  
-Complex  
-Unresolved  

#!/bin/bash

#liWGS-SV Pipeline: Module 7 (Complex SV Categorization)
#August 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
samples_list=$1
params=$2

#Source params file
. ${params}

#Dynamically determine min complex size
minCXsize=$( fgrep -v "#" ${OUTDIR}/QC/cohort/${COHORT_ID}.QC.metrics | cut -f9 | sort -nrk1,1 | head -n1 )

#Make cnMOPS list
while read ID bam sex; do
	echo -e "${ID}\t${WRKDIR}\t${ID}\t${ID}.cnMOPS.dels.bed" > ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list
	echo -e "${ID}\t${WRKDIR}\t${ID}\t${ID}.cnMOPS.dups.bed" > ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list
done < ${samples_list}

#Submit complex linking script
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/complex_linking.log -e ${OUTDIR}/logs/complex_linking.log -J ${COHORT_ID}_cxLINK "${liWGS_SV}/scripts/complexLinking.sh ${samples_list} ${params}"

#Submit inversion categorization
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/inversion_categorization.log -e ${OUTDIR}/logs/inversion_categorization.log -J ${COHORT_ID}_classification "${liWGS_SV}/scripts/classify_inversion.sh ${WRKDIR}/classifier/clusterfix/newCoords/inversion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_inversion.patched.clusters ${minCXsize} ${WRKDIR}/classifier/clusterfix/newCoords/deletion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/newCoords/insertion.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dels.list ${WRKDIR}/classifier/clusterfix/newCoords/cnMOPS_dups.list ${WRKDIR}/classifier/clusterfix/newCoords/"

#Submit translocation categorization
bsub -q normal -sla miket_sc -o ${OUTDIR}/logs/transloc_categorization.log -e ${OUTDIR}/logs/transloc_categorization.log -J ${COHORT_ID}_classification "${liWGS_SV}/scripts/classify_translocation.sh ${WRKDIR}/classifier/clusterfix/newCoords/transloc.events.reclassified.bedpe ${WRKDIR}/classifier/clusterfix/${COHORT_ID}_transloc.patched.clusters ${WRKDIR}/classifier/clusterfix/newCoords/"


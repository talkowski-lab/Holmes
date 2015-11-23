#!/bin/bash

#liWGS-SV Pipeline: Full-Chromosome Aneuploidy Check
#November 2015
#Contact: rcollins@chgr.mgh.harvard.edu

#Read input
ID=$1
bam=$2
params=$3

#Source params file
. ${params}

#Set additional params
genome_size=2897310462 #n-masked length of grch37

#Write header to outfile
echo -e "#HOLMES ANEUPLOIDY CHECK\n#${ID}\n\n#contig\tlibrary_count\tcontig_count\tobserved_fraction\texpected_fraction\tpredicted_copies" > ${OUTDIR}/QC/sample/${ID}/${ID}.aneuploidyCheck

#Get primary alignment count in library
total=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} )

#Iterate over primary contigs
for contig in $( seq 1 22 ) X Y; do
  count=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} ${contig} )
  frac_ex=$( echo -e "scale=6; (( $( fgrep -w "SN:${contig}" ${DICT} | cut -f3 | cut -d\: -f2 ) / ${genome_size} ))" | bc )
  frac_obs=$( echo -e "scale=6; (( ${count} / ${total} ))" | bc )
  copies=$( echo -e "scale=2; (( 2 * ${frac_obs} / ${frac_ex} ))" | bc | awk '{printf "%.1f\n",$1}' )
  echo -e "${contig}\t${total}\t${count}\t${frac_obs}\t${frac_ex}\t${copies}"
done >> ${OUTDIR}/QC/sample/${ID}/${ID}.aneuploidyCheck
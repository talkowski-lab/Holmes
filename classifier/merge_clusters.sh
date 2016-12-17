#!/bin/bash

# set -e 

# sample_list=$1
# prefix=$2
# svtype=$3
outdir=$PWD

usage(){
cat <<EOF
Usage: merge_clusters.sh SAMPLE_LIST PREFIX SVTYPE

SAMPLE_LIST List of samples to merge and classify.
	    Column format: sample_ID bamstat_directory cluster_distance
            cluster_size_cutoff gcov icov_file
PREFIX	    Prefix to use when writing merged output. Generally cohort name 
-t SVTYPE	    SV type to filter and merge {deletion,insertion,inversion,transloc}

-o OUTDIR   Output directory [CWD] 

-h		Print this message
EOF
}

# Check input
if [[ $# -lt 3 ]]; then
	usage
	exit 1
fi

# Help flag
while getopts i:p:t:o:h opt; do
	case $opt in
        i)
            sample_list=$OPTARG
            ;;
        p)
            prefix=$OPTARG
            ;;
        t)
            svtype=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
		h)
			usage
			exit 0
			;;
	esac
done

# Create output directory if it doesn't exist
if [ ! -d "$outdir" ]; then
    mkdir $outdir
fi

while read sample bstat_dir dist cutoff null; do

    # Old naming convention
    fin=${bstat_dir}/*${svtype}_clusters_d*_q-1.txt

    # New naming convention includes cluster size
    if [ ! -f "$fin" ]; then
        fin=${bstat_dir}/*${svtype}_clusters_d*_q-1_s*.txt
    fi

    joboutput=$(bsub -q normal -o ${outdir}/${svtype}_filter.out -sla miket_sc "
filter_sample_clusters.py ${fin} ${outdir}/${sample}_${svtype}.reads ${sample} -c 3 \
    --exclude /data/talkowski/rlc47/src/b37.lumpy.exclude.4-13.bed;
sort -k2,2 -k5,5 -k3n,3n ${outdir}/${sample}_${svtype}.reads > ${outdir}/${sample}_${svtype}.sorted.reads")

	expr "$joboutput" : 'Job <\([0-9]\+\)>' >> ${outdir}/filter_${svtype}.joblist
done < <( sed '/^#/d' ${sample_list} )

while [[ $(fgrep -f ${outdir}/filter_${svtype}.joblist <(bjobs | awk '{print $1}') | wc -l) -ne 0 ]]
do
        sleep 20s
done

joboutput=$(bsub -q medium -o ${outdir}/merge.out -sla miket_sc "
sort --merge -k2,2 -k5,5 -k3n,3n ${outdir}/*_${svtype}.sorted.reads > ${outdir}/${prefix}_${svtype}.reads")

jobid=$(expr "$joboutput" : 'Job <\([0-9]\+\)>')

while [[ $(fgrep ${jobid} <(bjobs | awk '{print $1}') | wc -l) -ne 0 ]]; do
    sleep 20s
done

# clean up
rm ${outdir}/filter_${svtype}.joblist

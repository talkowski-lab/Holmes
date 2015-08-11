#!/bin/bash

# BAMFILE=$1
# STATSDIR=$2

set -e

usage(){
cat <<EOF
Usage: bamstat.sh [-s SIZE] [-n] -i bamfile -o output_dir

-i bamfile      Name-sorted bamfile to process.
-o output_dir   Directory to save discordant read pairs to.
-s SIZE         (optional) Minimum cluster size used by readPairCluster.
                Defaults to 7% of estimated haploid coverage in the bamfile.
-n              (optional) Name sort bamfile before running bamstat.

-h              Print this message
EOF
}

# Check input
if [[ $# -lt 2 ]]; then
	usage
	exit 1
fi

if [[ -d "$STATSDIR" ]]
then
	echo "Output directory $STATSDIR already exists"
	exit 1
fi

# Set size if specified by the user
SIZE=0
SORT=false
while getopts ns:i:o:h opt; do
	case $opt in
		s)
			SIZE=$OPTARG
			;;
		n)
			SORT=true
			;;
		i)
			BAMFILE=$OPTARG
			;;
		o)
			STATSDIR=$OPTARG
			;;
		h)
			usage
			exit 0
			;;
	esac
done

# Run Bamstat
mkdir $STATSDIR
cd $STATSDIR

if $SORT; then
	# TODO: Bamstat currently closes input file after sampling, 
	# preventing use of pipes. Fix this.
	# sambamba sort -n -o /dev/stdout $BAMFILE | bamstat -i - -d 7 -b 
	/data/talkowski/tools/bin/sambamba_v0.4.6 view -f bam -F 'not (secondary_alignment or duplicate)' -o ${BAMFILE}.nosecondary.tmp_bamstat.bam ${BAMFILE}
	/data/talkowski/tools/bin/sambamba_v0.4.6 sort -n -m 4GB -p -o ${BAMFILE}.tmp_nsort.bam ${BAMFILE}.nosecondary.tmp_bamstat.bam
	bamstat -i ${BAMFILE}.tmp_nsort.bam -d 7 -b
	rm ${BAMFILE}.tmp_nsort.bam ${BAMFILE}.nosecondary.tmp_bamstat.bam
else
	sambamba view -f bam -F 'not (secondary_alignment or duplicate)' -o ${BAMFILE}.nosecondary.tmp_bamstat.bam ${BAMFILE}
	bamstat -i ${BAMFILE}.nosecondary.tmp_bamstat.bam -d 7 -b
	rm ${BAMFILE}.nosecondary.tmp_bamstat.bam
fi

# Compute readpaircluster parameters
# 1) Minimum mapq for clustering = -1
MAPQ=-1

# 2) Compute clustering distance as median insert size + 7*MAD
dist_cutoff=7
insert=$(grep -P "Actual FR median insert size:" stats.file | grep -oP "\d+")
mad=$(grep -P "Actual FR median absolute deviation:" stats.file | grep -oP "\d+")
DISTANCE=$(echo "$insert + $dist_cutoff * $mad" | bc)

# 3) Compute minimum cluster size as 7% of estimated haploid coverage
if [[ $SIZE -eq 0 ]]; then
	size_cutoff=0.07
	hg19_usable_size=2897310462
	# proper_pairs=$(samtools view -c -f 2 -F 1024 ${BAMFILE})
    proper_pairs=$(fgrep "Proper pairs" stats.file | grep -oP "(\d+)")
    coverage=$(python -c "print $insert * $proper_pairs / ${hg19_usable_size}.0")
    # coverage=$(scrape_bamstat.py 
	SIZE=$(python -c "print int(round($size_cutoff * $coverage))")
fi

echo ReadPairCluster parameters:
echo Minimum mapq $MAPQ
echo Max insert $DISTANCE
echo Minimum cluster size $SIZE

# Rename for consistency
mv translocation_pairs.txt transloc_pairs.txt

# Coordinate sort each set of reads then pass to readpaircluster
for svtype in "deletion" "insertion" "inversion" "transloc"
do
	bsub -sla miket_sc -q normal -o rpc.out -J BAMSTAT_GATE "
sort -k2,2 -k5,5 -k3n,3n ${svtype}_pairs.txt > ${svtype}_pairs.sorted.txt;
readPairCluster -d $DISTANCE -s $SIZE -q $MAPQ -r ${svtype}_pairs.sorted.txt > ${svtype}_clusters_d${DISTANCE}_q${MAPQ}_s${SIZE}.txt"
done


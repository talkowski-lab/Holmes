#!/bin/bash
#
# rpc.sh
# Copyright (C) 2015 Matthew Stone <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.
#


set -e

usage(){
cat <<EOF
Usage: rpc_single.sh [-b] [-n N_MAD] (-d DIST | -m METRICS) SAMPLE BAM

SAMPLE, BAM, and one of DIST and METRICS are required arguments. 

If DIST is not provided, median insert size and MAD will be parsed from METRICS
to compute clustering distance. Distance is computed as median + N_MAD * MAD.

SAMPLE      Sample ID
BAM         Bam file (coordinate sorted)
-d DIST     Clustering distance (generally median insert + 7 MAD)
-m METRICS  Metrics file output by Picard's InsertSizeMetrics
-n N_MAD    Number of MAD added to median insert to compute minimum
            clustering distance [7]

-h          Print this message
EOF
}
# -b          Convert rpc cluster files to bedpe format [FALSE]
            # (Requires compress_clusters.py script from pycluster repo)

dist=""
metrics=""
bedpe=false
n_mad=7
while getopts ":d:m:bn:h" opt; do
  case "$opt" in
    d)
      dist=$OPTARG
      ;;
    m)
      metrics=$OPTARG
      ;;
    # b)
      # bedpe=true
      # ;;
    n)
      n_mad=$OPTARG
      ;;
    h)
      usage
      exit 0
      ;;
  esac
done
shift $(( OPTIND - 1))

sample=$1
bam=$2

# If distance wasn't specified, calculate it from metrics file
if [[ -z $dist ]]; then
  if [[ -z $metrics ]]; then
    echo "ERROR: Must specify either distance or metrics file"
    usage
    exit 1
  fi

  read median mad <<<$(head -n8 $metrics | tail -n1 | cut -f1-2)
  dist=$(($median + $n_mad * $mad))
fi

# Extract discordant reads from bam
# Equivalent to `samtools view -F 3342`
sambamba view -f bam -F "((paired) and
                          (not (proper_pair or unmapped or mate_is_unmapped or 
                                secondary_alignment or supplementary or 
                                duplicate)))" $bam > ${sample}.disc.bam

# Convert discordant reads to bamstat format
# (Filters pairs where mapq=0 on both sides)
disc_to_rpc.py ${sample}.disc.bam ${sample}

for sv in del dup inv tloc; do
  bsub \
    -q big \
    -M 96000 \
    -R 'rusage[mem=96000]' \
    -v 96000 \
    -sla miket_sc \
    -J rpc_${sample}_${sv} \
    -o rpc.out "
  rpc.py \
    -d $dist \
    -x /data/talkowski/rlc47/src/b37.lumpy.exclude.4-13.bed.gz \
    ${sample}.${sv}.pairs.txt \
    ${sample}.${sv}.clusters.txt;
  compress_clusters.py \
      --no-split \
      ${sample}.${sv}.clusters.txt \
      ${sample}.${sv}.rpc.bedpe \
      ${sample}_${sv}
  "
done

# if [[ $bedpe ]]; then
  # for sv in del dup inv tloc; do
    # compress_clusters.py \
      # --no-split \
      # ${sample}.${sv}.clusters.txt \
      # ${sample}.${sv}.rpc.bedpe \
      # ${sample}_${sv}
  # done
# fi

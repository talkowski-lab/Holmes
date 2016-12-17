#!/bin/bash

set -e

usage(){
cat <<EOF
Usage: run_classify.sh SAMPLE_LIST PREFIX ICOV_LIST

SAMPLE_LIST List of samples to merge and classify.
  Column format: sample_ID bamstat_directory
PREFIX  Prefix to use when assigning cluster names and writing output.
  Generally cohort name (e.g. LOGIC, SFARI)
ICOV_LIST (optional) List of sample insert coverage files.
  Insert coverage files must be bgzipped and tabix indexed
  bed files.
  Column format: sample_ID insert_coverage_filepath

-m  Take Holmes output QC metrics file for insert sizes

-h  Print this message
EOF
}

# Help flag
bootstrap="bootstrap_classifier.py"
while getopts ":m:h" opt; do
 case ${opt} in
  m)
   bootstrap="bootstrap_classifier.py -m ${OPTARG}"
   ;;
  h)
   usage
   exit 0
   ;;
 esac
done
shift $(( OPTIND - 1))

#Read positional args
sample_list=$1
prefix=$2
icov_file=$3

# Check input
if [[ $# -lt 3 ]]; then
 usage
 exit 1
fi
#unset positional arguments so conda source deactivate will work
set --

#Ensure default python env (python 2)
source /apps/lab/miket/anaconda/4.0.5/envs/collins_py3/bin/deactivate
source /apps/lab/miket/anaconda/4.0.5/bin/activate collins_py2

# Ensure correct gcc version for anaconda dependencies
module rm gcc/4.9.0
module rm gcc-4.4
module load gcc/4.9.0

# If no insert coverage file specified, create it
if [[ $icov_file == "" ]]; then
 pull_icov.sh ${sample_list} ${prefix}
 icov_file=${prefix}_icov.list
fi

# Generate all the data needed by the classifier
sed '/^#/d' ${sample_list} | ${bootstrap} | paste - <(awk '{print $2}' $icov_file) > ${prefix}_boot.list
sample_list=${prefix}_boot.list

# Compute max clustering distance
cutoff=$(awk 'BEGIN {max = 0} {if ($3 > max) max=$3} END {print max}' <(sed '/^#/d' ${sample_list}))

# Require cluster size to be at least 3
cluster_size=3

# 1) For each sample in $SAMPLE_LIST, merge reads from all samples together for reclustering
# 2) Cluster reads from all samples
# 3) Calculate cluster mapq, uniq, global coverage, and local coverage, then classify clusters
for svtype in inversion insertion deletion transloc; do
 bsub -J ${prefix}_MERGE -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -o ${svtype}_classify.out "merge_clusters.sh -i ${sample_list} -p ${prefix} -t ${svtype} -o filtered"
  # Otherwise making filtered directory might fail
  sleep 3s
done

#Gate until merge jobs complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_MERGE" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_MERGE" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at read merging..."
    GATEwait=0
  fi
done

#Split merged reads by chromosome and submit RPC
mkdir chrsplit_clustering
for svtype in inversion deletion; do
 for chr in $( seq 1 22 ) X Y; do
  awk -v chr=${chr} '{ if ($2==chr) print $0 }' filtered/${prefix}_${svtype}.reads | sed -e "s/b'\([ACTGN]\{25\}\)'/\1/g" > chrsplit_clustering/${prefix}_${svtype}.${chr}.reads
  #Old C++ RPC implementation
  # bsub -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -J ${prefix}_cluster -o ${svtype}_classify.out "readPairCluster -d ${cutoff} -q -1 -s 3 -r chrsplit_clustering/${prefix}_${svtype}.${chr}.reads > chrsplit_clustering/${prefix}_${svtype}.${chr}.clusters"
  #Matt's python implementation
  bsub -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -J ${prefix}_cluster -o ${svtype}_classify.out "rpc.py -d ${cutoff} -q -1 -s ${cluster_size} chrsplit_clustering/${prefix}_${svtype}.${chr}.reads chrsplit_clustering/${prefix}_${svtype}.${chr}.clusters"
 done
done
#exclude all insertion pairs not separated by > 1kb to protect against small sequenced linear fragments
for chr in $( seq 1 22 ) X Y; do
  awk -v chr=${chr} '{ if ($2==chr && $6-$3>1000) print $0 }' filtered/${prefix}_insertion.reads | sed -e "s/b'\([ACTGN]\{25\}\)'/\1/g" > chrsplit_clustering/${prefix}_insertion.${chr}.reads
  #Old C++ RPC implementation
  # bsub -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -J ${prefix}_cluster -o ${svtype}_classify.out "readPairCluster -d ${cutoff} -q -1 -s 3 -r chrsplit_clustering/${prefix}_insertion.${chr}.reads > chrsplit_clustering/${prefix}_insertion.${chr}.clusters"
  #Matt's python implementation
  bsub -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -J ${prefix}_cluster -o insertion_classify.out "rpc.py -d ${cutoff} -q -1 -s ${cluster_size} chrsplit_clustering/${prefix}_insertion.${chr}.reads chrsplit_clustering/${prefix}_insertion.${chr}.clusters"
done
for chrA in $( seq 1 22 ) X Y; do
 for chrB in $( seq 1 22 ) X Y; do
  awk -v chrA=${chrA} -v chrB=${chrB} '{ if ($2==chrA && $5==chrB) print $0 }' filtered/${prefix}_transloc.reads | sed -e "s/b'\([ACTGN]\{25\}\)'/\1/g" > chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.reads
  if [ -e chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.reads ] && [ -s chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.reads ]; then
   #Old C++ RPC implementation
   # bsub -q big -R 'rusage[mem=20000]' -M 20000 -v 28000 -sla miket_sc -J ${prefix}_cluster -o ${svtype}_classify.out "readPairCluster -d ${cutoff} -q -1 -s 3 -r chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.reads > chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.clusters"
   #Matt's python implementation
   bsub -q normal -R 'rusage[mem=6000]' -M 6000 -v 10000 -sla miket_sc -J ${prefix}_cluster -o transloc_classify.out "rpc.py -d ${cutoff} -q -1 -s ${cluster_size} chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.reads chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.clusters"
  fi
 done
done

#Gate until clustering jobs complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_cluster" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_cluster" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at joint clustering..."
    GATEwait=0
  fi
done

#Merge clusters and classify
for svtype in inversion insertion deletion; do
 for chr in $( seq 1 22 ) X Y; do
  awk -v chr=${chr} '{ if ($1!="") print chr"_"$0; else print $0 }' chrsplit_clustering/${prefix}_${svtype}.${chr}.clusters >> ${prefix}_${svtype}.clusters
 done
 bsub -q normal -sla miket_sc -J ${prefix}_classify -o ${svtype}_classify.out "rpc_classify.py ${prefix}_${svtype}.clusters ${prefix}_${svtype}.events.bedpe ${sample_list} ${prefix}_${svtype} --cluster-bedpe ${prefix}_${svtype}.bkpts.bedpe"
done
for chrA in $( seq 1 22 ) X Y; do
 for chrB in $( seq 1 22 ) X Y; do
  if [ -e chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.clusters ]; then
    awk -v chrA=${chrA} -v chrB=${chrB} '{ if ($1!="") print chrA"_"chrB"_"$0; else print $0 }' chrsplit_clustering/${prefix}_transloc.${chrA}_${chrB}.clusters >> ${prefix}_transloc.clusters
  fi
 done
done
bsub -q normal -sla miket_sc -J ${prefix}_classify -o transloc_classify.out "rpc_classify.py ${prefix}_transloc.clusters ${prefix}_transloc.events.bedpe ${sample_list} ${prefix}_transloc --cluster-bedpe ${prefix}_transloc.bkpts.bedpe"

#Gate until classification jobs complete; 20 sec check; 5 min report
GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_classify" | wc -l )
GATEwait=0
until [[ $GATEcount == 0 ]]; do
  sleep 20s
  GATEcount=$( bjobs -w | awk '{ print $7 }' | grep -e "${prefix}_classify" | wc -l )
  GATEwait=$[${GATEwait} +1]
  if [[ $GATEwait == 15 ]]; then
    echo -e "STATUS [$(date)]: Gated at cluster classification..."
    GATEwait=0
  fi
done

# clean up
if [ ${KEEP_TMP} != "TRUE" ]; then
  rm -r chrsplit_clustering
fi

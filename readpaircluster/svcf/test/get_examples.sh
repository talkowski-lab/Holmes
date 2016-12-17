#!/bin/bash
#
# get_examples.sh
#
# 
#
# Copyright (C) 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.


member=$1
family=SFARI_d11194
sample=${family}${member}

svdir=/data/talkowski/Samples/SFARI/samples/${sample}/SV

TEST_COUNT=3
example_lumpy=example.lumpy.${member}.vcf
example_delly=example.delly.${member}.vcf
example_rpc=example.rpc.${member}.bedpe

function print_header {
  fin=$1
  grep -P -e "^#" $fin
}

function sample_vcf {
  vcf=$1
  sv=$2

  sed '/^#/d' $vcf \
    | fgrep $sv \
    | head -n $TEST_COUNT \
    | sort -k1,1V -k2,2n
}

# Get test LUMPY VCF
# Keep header and first $test_count of each svtype
lumpy_vcf=${svdir}/lumpy/${sample}.vcf
print_header $lumpy_vcf > $example_lumpy

for sv in DEL DUP INV BND; do
  sample_vcf $lumpy_vcf $sv >> $example_lumpy
done
gzip -f $example_lumpy


# Get test DELLY VCF
delly_vcf=${svdir}/delly/${sample}.DEL.vcf
print_header $delly_vcf > $example_delly

for sv in DEL DUP INV TRA; do
  delly_vcf=${svdir}/delly/${sample}.${sv}.vcf
  sample_vcf $delly_vcf $sv >> $example_delly
done

# Get test RPC bedpe
if [[ -f "$example_rpc" ]]; then
  rm $example_rpc
fi

for sv in del dup inv tloc; do
  rpc_file=${svdir}/rpc/mapq_filter/${sample}.del.clusters.txt
  bedpe=$(mktemp -p $SCRATCH_TMP)
  compress_clusters.py $rpc_file $bedpe ${sample}_${sv} --no-split
  head -n $TEST_COUNT $bedpe >> $example_rpc
  rm $bedpe
done
sort -k1,1V -k4,4V -k2,2n $example_rpc > ${example_rpc}.sorted
mv ${example_rpc}.sorted $example_rpc



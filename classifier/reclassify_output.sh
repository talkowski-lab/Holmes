#!/bin/bash

############################
#     Reclassify Output    #
#       Talkowski Lab      #
############################

#Input
events=$1        #events.bedpe
svtype=$2        #deletion, insertion, inversion, transloc
clusters=$3      #corresponding cluster file
OUTDIR=$4        #output = ${outdir}/${svtype}.events.reclassified.bedpe
JS=$5            #if TRUE, will REMOVE all clusters on span <=50 (not just mark them as invalid)

#Set params
INVALID=`mktemp`
VALID=`mktemp`
POOL=`mktemp`
precoords=`mktemp`
postcoords=`mktemp`
cluster=`mktemp`

#Type-specific filtering
case ${svtype} in
  deletion)
    #Global filtering (5 total filters)
    if [ ${JS} == "TRUE" ]; then
        fgrep deletion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3 && $3-$2>50 && $6-$5>50) print $0 }' > ${POOL}
    else
        fgrep deletion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3) print $0 }' > ${POOL}
    fi
    awk -v OFS="\t" '{ if ( ($15>=$14 && $15!=0 && ($15-$14)/$15>0.304) || ($14>$15 && ($14-$15)/$14>0.304) || ($12+$13<27) || ($17+$18<0.039) || $3-$2<52 || $6-$5<52 || (($3-$2)+($6-$5)>61591) ) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Invalid"; print }' > ${INVALID}
    cut -f7 ${INVALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Deletions ≥ 25kb are valid pending cnMOPS 51% reciprocal overlap 
    awk -v OFS="\t" '{ if ($5-$3>=25000) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid_pending_cnMOPS"; print }' > ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Deletions < 25kb must have min(mapq)>10 and min(uniq)≥85%
    awk -v OFS="\t" '{ if ($5-$3<25000 && $12>10 && $13>10 && $14>=0.85 && $15>=0.85) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    awk -v OFS="\t" '{ $21="Invalid"; print }' ${POOL} >> ${INVALID}
    ;;

  insertion)
    #Global filtering (5 total filters)
    if [ ${JS} == "TRUE" ]; then
        fgrep insertion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3 && $3-$2>50 && $6-$5>50) print $0 }' > ${POOL}
    else
        fgrep insertion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3) print $0 }' > ${POOL}
    fi
    awk -v OFS="\t" '{ if ( ($15>=$14 && $15!=0 && ($15-$14)/$15>0.304) || ($14>$15 && ($14-$15)/$14>0.304) || ($12+$13<27) || ($17+$18<0.039) || $3-$2<52 || $6-$5<52 || (($3-$2)+($6-$5)>61591) ) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Invalid"; print }' > ${INVALID}
    cut -f7 ${INVALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Insertions ≥ 40kb with sum(local) ≥ 40% and global ≥ 35% are valid, with or without cnMOPS support
    awk -v OFS="\t" '{ if ($6-$2>=40000 && $17+$18>=0.4 && $16>=0.35) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid_pending_cnMOPS"; print }' > ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Insertions ≥ 40kb with sum(local) < 40% are valid pending cnMOPS 51% reciprocal overlap 
    awk -v OFS="\t" '{ if ($6-$2>=40000) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid_pending_cnMOPS"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Insertions < 40kb must have min(mapq)>10 and min(uniq)≥85% and max(span)≤16kb
    awk -v OFS="\t" '{ if ($6-$2<40000 && $12>10 && $13>10 && $14>=0.85 && $15>=0.85 && ($3-$2)<=16000 && ($6-$5)<=16000 ) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    awk -v OFS="\t" '{ $21="Invalid"; print }' ${POOL} >> ${INVALID}
    ;;

  inversion)
    #Global filtering (5 total filters)
    if [ ${JS} == "TRUE" ]; then
        fgrep inversion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3 && $3-$2>50 && $6-$5>50) print $0 }' > ${POOL}
    else
        fgrep inversion ${events} | awk -v OFS="\t" '{ if ($11/$19>=3) print $0 }' > ${POOL}
    fi
    awk -v OFS="\t" '{ if ( ($15>=$14 && $15!=0 && ($15-$14)/$15>0.304) || ($14>$15 && ($14-$15)/$14>0.304) || ($12+$13<27) || ($17+$18<0.039) || $3-$2<52 || $6-$5<52 || (($3-$2)+($6-$5)>61591) ) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Invalid"; print }' > ${INVALID}
    cut -f7 ${INVALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Single-enders are only valid if they have min(mapq) ≥ 2.7 and sum(local) ≥ 16.4
    awk -v OFS="\t" '{ if ($9!="." && $12>=2.7 && $13>=2.7 && $17+$18>=0.164) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' > ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    awk -v OFS="\t" '{ if ($9!=".") print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Invalid"; print }' >> ${INVALID}
    cut -f7 ${INVALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #High-confidence double-enders have global ≥ 5.6%, min(mapq) ≥ 6.6, min(uniq) ≥ 93.9%, min(span) ≥ 3.1kb, min(local) ≥ 2.7%, and sum(local)>14%
    awk '{ if ($9=="." && $10=="." && $16>=0.056 && $17>=0.027 && $18>=0.027 && $12>=6.6 && $13>=6.6 && $14>=0.939 && $15>=0.939 && $3-$2>=3100 && $6-$5>=3100 && $17+$18>.14) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Low-confidence double-enders have global ≥ 5.6%, min(mapq) ≥ 6.6, min(uniq) ≥ 93.9%, min(span) ≥ 3.1kb, min(local) ≥ 2.7%, max(uniq)=1, and max(span)≤25kb
    awk '{ if ($9=="." && $10=="." && $16>=0.056 && $17>=0.027 && $18>=0.027 && $12>=6.6 && $13>=6.6 && $14>=0.939 && $15>=0.939 && $3-$2>=3100 && $6-$5>=3100 && ($14==1 || $15==1) && $3-$2<=25000 && $6-$5<=25000) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    awk -v OFS="\t" '{ $21="Invalid"; print }' ${POOL} >> ${INVALID}
    ;;

  transloc)
    #Global filtering (5 total filters)
    if [ ${JS} == "TRUE" ]; then
        fgrep transloc ${events} | awk -v OFS="\t" '{ if ($11/$19>=3 && $3-$2>50 && $6-$5>50) print $0 }' > ${POOL}
    else
        fgrep transloc ${events} | awk -v OFS="\t" '{ if ($11/$19>=3) print $0 }' > ${POOL}
    fi
    awk -v OFS="\t" '{ if ( ($15>=$14 && $15!=0 && ($15-$14)/$15>0.304) || ($14>$15 && ($14-$15)/$14>0.304) || ($12+$13<27) || ($17+$18<0.039) || $3-$2<52 || $6-$5<52 || (($3-$2)+($6-$5)>61591) ) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Invalid"; print }' > ${INVALID}
    cut -f7 ${INVALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Valid double-enders have global ≥ 7.1% and sum(mapq) ≥ 45.83
    awk '{ if ($9=="." && $10=="." && $16>=0.071 && $12+$13>=45.83) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' > ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    #Valid single-enders have sum(mapq) ≥ 53.72
    awk '{ if ($9!="." && $10!="." && $12+$13>=53.72) print $0 }' ${POOL} | awk -v OFS="\t" '{ $21="Valid"; print }' >> ${VALID}
    cut -f7 ${VALID} | fgrep -wvf - ${POOL} > ${POOL}.tmp
    mv ${POOL}.tmp ${POOL}
    awk -v OFS="\t" '{ $21="Invalid"; print }' ${POOL} >> ${INVALID}
    ;;
esac

#Change coordinate reporting in events.bedpe output
case ${svtype} in
    deletion)
        cat ${VALID} ${INVALID} > ${precoords}
        while read idx; do
            echo ${idx}
            awk -v idx=${idx} '{ if ($1==idx) print $0 }' ${clusters} > ${cluster}
            if [ $( egrep "_deletion_${idx}\b" ${events} | cut -f9 ) == "+" ]; then
                as=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                as=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            fi
            if [ $( egrep "_deletion_${idx}\b" ${events} | cut -f10 ) == "+" ]; then
                bs=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                bs=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            fi
            egrep "_deletion_${idx}\b" ${precoords} | awk -v OFS="\t" -v as=${as} -v ae=${ae} -v bs=${bs} -v be=${be} '{ print $1, as, ae, $4, bs, be, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' >> ${postcoords}
        done < <( cut -f7 ${precoords} | sed 's/deletion_/\t/g' | cut -f2 )
        ;;

    insertion)
        cat ${VALID} ${INVALID} > ${precoords}
        while read idx; do
            echo ${idx}
            awk -v idx=${idx} '{ if ($1==idx) print $0 }' ${clusters} > ${cluster}
            if [ $( egrep "_insertion_${idx}\b" ${events} | cut -f9 ) == "+" ]; then
                as=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                as=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            fi
            if [ $( egrep "_insertion_${idx}\b" ${events} | cut -f10 ) == "+" ]; then
                bs=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                bs=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            fi
            egrep "_insertion_${idx}\b" ${precoords} | awk -v OFS="\t" -v as=${as} -v ae=${ae} -v bs=${bs} -v be=${be} '{ print $1, as, ae, $4, bs, be, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' >> ${postcoords}
        done < <( cut -f7 ${precoords} | sed 's/insertion_/\t/g' | cut -f2 )
        ;;

    inversion)
        cat ${VALID} ${INVALID} > ${precoords}
        while read idx; do
            echo ${idx}
            awk -v idx=${idx} '{ if ($1==idx) print $0 }' ${clusters} > ${cluster}
            if [ $( egrep "_inversion_${idx}\b" ${events} | cut -f9 ) == "+" ]; then
                as=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            elif [ $( egrep "_inversion_${idx}\b" ${events} | cut -f9 ) == "-" ]; then
                as=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                ae=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                #cols 2 & 3 in double-ender inversion are 3rd quartile reads nearest +/+ breakpoint on both sides
                #cols 5 & 6 in double-ender inversion are 1rd quartile reads nearest -/- breakpoint on both sides
                as=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
                bs=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
            fi
            if [ $( egrep "_inversion_${idx}\b" ${events} | cut -f10 ) == "+" ]; then
                bs=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            elif [ $( egrep "_inversion_${idx}\b" ${events} | cut -f10 ) == "-" ]; then
                bs=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
            else
                ae=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
                be=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
            fi
            egrep "_inversion_${idx}\b" ${precoords} | awk -v OFS="\t" -v as=${as} -v ae=${ae} -v bs=${bs} -v be=${be} '{ print $1, as, ae, $4, bs, be, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' >> ${postcoords}
        done < <( cut -f7 ${precoords} | sed 's/inversion_/\t/g' | cut -f2 )
        ;;

    #Remove transloc coordinate switching because it's arbitrary
    transloc)
        cat ${VALID} ${INVALID} > ${postcoords}
        # cat ${VALID} ${INVALID} > ${precoords}
        # while read idx; do
        #     echo ${idx}
        #     awk -v idx=${idx} '{ if ($1==idx) print $0 }' ${clusters} > ${cluster}
        #     if [ $( egrep "_transloc_${idx}\b" ${events} | cut -f9 ) == "+" ]; then
        #         as=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         ae=$( awk '{ if ($12==97 || $12==161 || $12==129 || $12==65) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     elif [ $( egrep "_transloc_${idx}\b" ${events} | cut -f9 ) == "-" ]; then
        #         as=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         ae=$( awk '{ if ($12==81 || $12==145 || $12==113 || $12==177) print $6 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     else
        #         #cols 2 & 3 correspond to Q1 & Q3 positions on chrA for double-ender translocations
        #         as=$( cut -f6 ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         ae=$( cut -f6 ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     fi
        #     if [ $( egrep "_transloc_${idx}\b" ${events} | cut -f10 ) == "+" ]; then
        #         bs=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         be=$( awk '{ if ($13==97 || $13==161 || $13==129 || $13==65) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     elif [ $( egrep "_transloc_${idx}\b" ${events} | cut -f10 ) == "-" ]; then
        #         bs=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         be=$( awk '{ if ($13==81 || $13==145 || $13==113 || $13==177) print $7 }' ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     else
        #         #cols 5 & 6 correspond to Q1 & Q3 positions on chrB for double-ender translocations
        #         bs=$( cut -f7 ${cluster} | sort -n | perl -e '$d=.25;@l=<>;print $l[int($d*$#l)]' )
        #         be=$( cut -f7 ${cluster} | sort -n | perl -e '$d=.75;@l=<>;print $l[int($d*$#l)]' )
        #     fi
        #     egrep "_transloc_${idx}\b" ${precoords} | awk -v OFS="\t" -v as=${as} -v ae=${ae} -v bs=${bs} -v be=${be} '{ print $1, as, ae, $4, bs, be, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' >> ${postcoords}
        # done < <( cut -f7 ${precoords} | sed 's/transloc_/\t/g' | cut -f2 )
        ;;
esac

#Cleanup
sed -e 's/^X/23/g' -e 's/^Y/24/g' ${postcoords} | sort -nk1,1 -nk3,3 -nk2,2 -nk4,4 | sed -e 's/^23/X/g' -e 's/^24/Y/g' > ${OUTDIR}/${svtype}.events.reclassified.bedpe
rm -rf ${INVALID} ${VALID} ${POOL} ${precoords} ${postcoords} ${cluster}


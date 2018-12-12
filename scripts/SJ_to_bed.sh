#!/bin/bash
##run from results directory
##Output BED file is sorted for indexing and loading into IGV
##awk for converting SJ.out.tab to bed12 format
##based on code originally published by frymor at http://seqanswers.com/forums/showthread.php?t=62896
for sj in */*SJ.out.tab; do
    echo ${sj}
    echo "Converting..."
    awk \
        {'if($4=="2") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$8"\t-\t"$2-$9-1"\t"$3+$9"\t255,0,0\t2\t"$9","$9"\t","0,"$3-$2+$9+1; \
        else \
        if($4=="1") print ""$1"\t"$2-$9-1"\t"$3+$9"\tJUNC000"NR"\t"$8"\t+\t"$2-$9-1"\t"$3+$9"\t0,0,255\t2\t"$9","$9"\t","0,"$3-$2+$9+1'} \
    ${sj} > ${sj%.*}.bed12
    echo "Sorting..."
    sort -V -o ${sj%.*}.sort.bed ${sj%.*}.bed12
done
echo "Complete"

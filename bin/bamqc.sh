#!/bin/bash
#trimgalore.sh

usage() {
    echo "-h Help documentation for bamqc.sh"
    echo "-r  --Reference Genome: GRCh38 or GRCm38"
    echo "-b  --BAM File"
    echo "-p  --Prefix for output file name"
    echo "Example: bash bamqc.sh -p prefix -r ref.tar.gz -b SRR1551047.bam"
    exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:c:p:h opt
do
    case $opt in
    b) sbam=$OPTARG;;
    p) prefix=$OPTARG;;
    h) usage;;
    esac
done
shift $(($OPTIND -1))

if [[ -z $threads ]]
then
    threads=`nproc`
fi

issorted=$(samtools stats -@ ${threads} ${sbam} | grep "is sorted: 1")
if [[ -z $issorted ]]
then
    samtools sort -T ./ -@ ${threads} -o ${prefix}.sort.bam ${sbam}
    sbam="${prefix}.sort.bam"
fi

samtools view -@ ${threads} -F 1024 -1 -o ${prefix}.filt.bam ${sbam}
sbam="${prefix}.filt.bam"
samtools index -@ ${threads} ${sbam}

fastqc -f bam ${sbam} 
samtools flagstat -@ $threads ${sbam} > ${prefix}.flagstat.txt
samtools idxstats -@ $threads ${sbam} > ${prefix}.idxstat.txt
samtools stats -@ $threads ${sbam} > ${prefix}.stat.txt
samtools coverage -o ${prefix}.coverage.txt ${sbam}

qualimap bamqc -bam ${sbam} --java-mem-size=18G -outdir ${prefix}_qualimap > ${prefix}_bamqc_log.txt

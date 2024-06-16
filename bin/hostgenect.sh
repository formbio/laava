#!/bin/bash
#bamcount.sh

usage() {
  echo "-h Help documentation for bamcount.sh"
  echo "-p --Prefix for output file name"
  echo "-t --threads (number of CPUs)"
  echo "Example: bash bamcount.sh -p prefix SRR1551047.bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :p:t:h opt
do
  case $opt in
    p) prefix=$OPTARG;;
    t) threads=$OPTARG;;
    h) usage;;
  esac
done

shift $(($OPTIND -1))

if [[ -z $threads ]]
then
  threads=`nproc`
fi
if [[ -z $prefix ]]
then
  prefix=$(echo $sbam | cut -d. -f1)
fi
sbam=$@
bedtools merge -i rnaref/exons.bed -o distinct -c 4 > exons.merged.bed
bedtools coverage -b ${sbam} -a exons.merged.bed > ${prefix}.bedout.txt
merge_bedcov.pl $prefix ${prefix}.bedout.txt

#!/bin/bash
module load samtools

> full.list
for bamfile in `ls *.trimmed.realigned.indelqual.readsrem.bam`
do
    prefix=`cut -d"." -f1 <(echo "${bamfile}")`
    echo $prefix
    samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f ../NC_045512.2.fasta ${bamfile} | cut -f1-4 > ${prefix}.depth
    cp ${prefix}_trimmed_union_snpEff_filtered.vcf ${prefix}.vcf
    echo ${prefix} >> full.list
done


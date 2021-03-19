#!/usr/bin/env bash

if [ $# != 7 ]
    then
        echo "Input1 : fastq 1 file"
        echo "Input2 : fastq 2 file"
        echo "Input3 : output directory ( bitseq output will be made in the directory )"
        echo "Input4 : bowtie2 index directory"
        echo "Input5 : fasta file"
        echo "Input6 : threads"
        echo "Input7 : name"
        exit 1
fi

#fastq1=$(ls $1*1.fastq)
#fastq2=$(ls $1*2.fastq)

# bitseq step 1. bowtie2 
#bowtie2 -q -k 5 --threads $6 --no-mixed --no-discordant -x $4 -1 $1 -2 $2 -S $3$7.sam
/home/seokju/New/build/bowtie2-2.3.5.1-linux-x86_64/bowtie2 -q -k 100 --threads $6 --no-mixed --no-discordant -x $4 -1 $1 -2 $2 | samtools view -hSb - -o $3$7.bam

# bitseq step 2. parseAlignment
/home/seokju/New/build/BitSeq/parseAlignment $3$7.bam -o $3$7.prob --trSeqFile $5 --trInfoFile $3$7.tr --uniform --verbose

# bitseq step 3. estimateVBExpression
/home/seokju/New/build/BitSeq/estimateVBExpression -o $3$7 -O RPKM -P $6 -t $3$7.tr $3$7.prob


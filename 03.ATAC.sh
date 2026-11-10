#!/bin/bash
#SBATCH --job-name=ATAC
#SBATCH --partition=low
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --error=%j.err
#SBATCH --output=%j.out

cd work_pahh
rawdata=cleandata_path
bowtie2 -x /ref/GRCm38/GRCm38 -1 $rawdata/sample_1.clean.fq.gz -2 $rawdata/sample/sample_2.clean.fq.gz | samtools sort -@ 5  -O bam -o sample.bam -
samtools index sample.bam
samtools flagstat sample.bam >sample.raw.stat
sambamba markdup -r sample.bam sample.rmdup.bam
samtools index sample.rmdup.bam
samtools flagstat sample.rmdup.bam >sample.rmdup.stat
samtools view -h -f 2 -q 30 sample.rmdup.bam | grep -v chrM | samtools sort -O bam -@ 5 -o - >sample.last.bam
samtools index sample.last.bam
samtools flagstat sample.last.bam >sample.last.stat
bedtools bamtobed -i sample.last.bam >sample.bed
macs2 callpeak -t sample.bed -g mm --nomodel --shift -100 --extsize 200 -n sample --outdir /home/fanxinhao/result/gut_muscle/ATAC/peak

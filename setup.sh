#!/bin/bash

# Step 1: Install dependencies
pip install pysam
pip install pybwa

# Step 2: Index the reference genome using BWA and Samtools
bwa index data/reference.fasta
samtools faidx data/reference.fasta

# Step 3 Align the reads to the reference genome using BWA and generate SAM files
bwa mem data/reference.fasta data/variant1.fasta > data/aligned1_reads.sam
bwa mem data/reference.fasta data/variant2.fasta > data/aligned2_reads.sam
bwa mem data/reference.fasta data/variant3.fasta > data/aligned3_reads.sam

# Step 4: Convert the SAM files to BAM files
samtools view -S -b data/aligned1_reads.sam > data/aligned1_reads.bam
samtools view -S -b data/aligned2_reads.sam > data/aligned2_reads.bam
samtools view -S -b data/aligned3_reads.sam > data/aligned3_reads.bam

# Step 5: Sort the BAM files
samtools sort data/aligned1_reads.bam -o data/sorted_aligned1_reads.bam
samtools sort data/aligned2_reads.bam -o data/sorted_aligned2_reads.bam
samtools sort data/aligned3_reads.bam -o data/sorted_aligned3_reads.bam

# Step 6: Index the BAM files
samtools index data/sorted_aligned1_reads.bam
samtools index data/sorted_aligned2_reads.bam
samtools index data/sorted_aligned3_reads.bam

# We have now generated the BAM files with aligned sequenced reads that are sorted and indexed.
# We will use these BAM files to identify single nucleotide substitutions using our tool, SNVC.

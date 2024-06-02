#!/bin/bash

# Step 1: Create and activate a virtual environment (optional but recommended)
python3 -m venv venv
source venv/bin/activate

# Step 2: Upgrade pip to the latest version
python3 -m pip install --upgrade pip

# Step 3: Install dependencies
python3 -m pip install pysam

# Step 4: Index the reference genome using BWA and Samtools
bwa index data/reference.fasta
samtools faidx data/reference.fasta

# Step 5: Align the reads to the reference genome using BWA and generate SAM files
bwa mem data/reference.fasta data/variant1.fasta > data/aligned1_reads.sam
bwa mem data/reference.fasta data/variant2.fasta > data/aligned2_reads.sam
bwa mem data/reference.fasta data/variant3.fasta > data/aligned3_reads.sam

# Step 6: Convert the SAM files to BAM files
samtools view -S -b data/aligned1_reads.sam > data/aligned1_reads.bam
samtools view -S -b data/aligned2_reads.sam > data/aligned2_reads.bam
samtools view -S -b data/aligned3_reads.sam > data/aligned3_reads.bam

# Step 7: Sort the BAM files
samtools sort data/aligned1_reads.bam -o data/sorted_aligned1_reads.bam
samtools sort data/aligned2_reads.bam -o data/sorted_aligned2_reads.bam
samtools sort data/aligned3_reads.bam -o data/sorted_aligned3_reads.bam

# Step 8: Index the BAM files
samtools index data/sorted_aligned1_reads.bam
samtools index data/sorted_aligned2_reads.bam
samtools index data/sorted_aligned3_reads.bam

# We have now generated the BAM files with aligned sequenced reads that are sorted and indexed.
# We will use these BAM files to identify single nucleotide substitutions using our tool, SNVC.
"setup.sh" 38L, 1582C                                                                                      38,95         All

# Single Nucleotide Variant Caller (SNVC) 

## Installation instructions
### Note: We assume you already have python installed since these installations require python. Please install python on your laptop if you haven't already. 

### Step 1: Clone our github repository and cd into the directory
``` git clone https://github.com/archxta/cse185_project ```
``` cd cse185_project ```

### Step 2: Install BWA and Samtools
  ``` pip install pysam ``` to install samtools

  ``` pip install pybwa ``` to install bwa

### Step 3: Index the reference genome using BWA and Samtools
``` bwa index reference.fa ```

``` samtools faidx reference.fasta ```

### Step 4: Align the reads to reference genome using BWA and generate the SAM file
``` bwa mem reference.fasta variant1.fasta > aligned1_reads.sam ```

``` bwa mem reference.fasta variant2.fasta > aligned2_reads.sam ```

``` bwa mem reference.fasta variant3.fasta > aligned3_reads.sam ```

### Step 5: Convert this sam file to the BAM file
``` samtools view -S -b aligned1_reads.sam > aligned1_reads.bam ```

``` samtools view -S -b aligned2_reads.sam > aligned2_reads.bam ```

``` samtools view -S -b aligned3_reads.sam > aligned3_reads.bam ```

### Step 6: Sort the BAM file
``` samtools sort aligned1_reads.bam -o sorted_aligned1_reads.bam ```

``` samtools sort aligned2_reads.bam -o sorted_aligned2_reads.bam ```

``` samtools sort aligned3_reads.bam -o sorted_aligned3_reads.bam ```

### Step 7: Index the BAM file
``` samtools index sorted_aligned1_reads.bam ```

``` samtools index sorted_aligned2_reads.bam ```

``` samtools index sorted_aligned3_reads.bam ```

We have now generated the BAM files with aligned sequenced reads that are sorted and indexed. We will use these BAM files to identify single nucleotide substitutions using our tool, SNVC.  

 

















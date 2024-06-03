# Single Nucleotide Variant Caller (SNVC) 

## Description


## Installation instructions
### Note: We assume you already have python3 installed since these installations require python. You can check if you have python3 installed by doing ```python3 --version```. Please install python3 on your laptop if you haven't already. Another way to run this tool  would be on Juypter Notebook terminal. 

### For Mac users: if you don't already have brew installed please install it using ``` /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" ``` and follow the instructions in the script. You can check if brew is installed by running this command ``` brew --version ```. 



Step 1: Clone our github repository and cd into the directory

``` git clone https://github.com/archxta/cse185_project ```

``` cd cse185_project ```

Step 2: Install the necessary dependents

``` brew install samtools ```

``` brew install bwa ```

``` pip3 install pysam ```

```python3 -m pip install matplotlib ```

Step 3: Run our setup script in your terminal to generate the BAM files (Note: ignore any error generated by this command, it should still create the BAM files successfully)

```bash setup.sh```

Step 4: We have now generated the BAM files with aligned sequenced reads that are sorted and indexed. We will use these BAM files to identify single nucleotide substitutions using our tool, SNVC.  To run our variant caller on our example datasets, run the following commands

```python3 code.py data/reference.fasta data/variant1.fasta data/sorted_aligned1_reads.bam results/output1.vcf```

```python3 code.py data/reference.fasta data/variant2.fasta data/sorted_aligned2_reads.bam results/output2.vcf```

```python3 code.py data/reference.fasta data/variant3.fasta data/sorted_aligned3_reads.bam results/output3.vcf```

This should have generated the final VCF file which contains information about the genetic variants detected in our DNA sequencing data that we obtained from a sickle cell mutation on the HBB gene. 

You can view the VCF file contents using the ```cat``` command followed by the filepath and filename (ex. ```cat results/output1.vcf```









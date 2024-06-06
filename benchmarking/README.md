# Benchmarking Single Nucleotide Variant Caller (SNVC) against VarScan

For convinience and datasize limitation purposes, we are using the .mpileup file that we already generated dunring Lab1 to perform benchmarking. We will be using the NA12878_child.mpileup - this can be downloaded from the Juypter Notebook folder corresponding to Lab1. We have also included it in this GitHub Repository for easier access. 

## Step 1: Running our tool on data from Lab1

python3 mpileup.py NA12878_child.mpileup snvoutput.vcf 

## Step 2: Run VarScan on the same input file by following the steps below:
Run the folloring commands:

curl -L https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download > VarScan.jar

java -jar VarScan.jar mpileup2snp NA12878_child.mpileup --min-var-frequency 0.2 --min-freq-for-hom 0.8 --p-value 0.01 --output-vcf 1 --variants-only > varscanoutput.vcf


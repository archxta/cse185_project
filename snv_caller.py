import argparse
import os
from src.index import index_reference_genome
from src.align import align_reads_to_reference
from src.call_variants import index_bam_file, analyze_reads
from src.utils import write_vcf

def main():
    parser = argparse.ArgumentParser(description='Single Nucleotide Variant Caller')
    parser.add_argument('-r', '--reference', required=True, help='Path to the reference genome in FASTA format')
    parser.add_argument('-f', '--fastq', required=True, help='Path to the reads in FASTQ format')
    parser.add_argument('-o', '--output', required=True, help='Path to the output VCF file')
    args = parser.parse_args()

    reference_fasta = args.reference
    reads_fastq = args.fastq
    output_vcf = args.output

    print("Indexing the reference genome...")
    bwa_index = index_reference_genome(reference_fasta)

    print("Aligning reads to the reference genome...")
    alignment = align_reads_to_reference(bwa_index, reads_fastq)

    bam_file = alignment.to_bam()
    bam_file_name = "aligned_reads.bam"
    bam_file.save(bam_file_name)

    print("Indexing the BAM file...")
    index_bam_file(bam_file_name)

    print("Analyzing reads for variants...")
    bam = pysam.AlignmentFile(bam_file_name, "rb")
    reference_name = os.path.basename(reference_fasta).split('.')[0]
    variants = analyze_reads(bam, reference_name)

    print("Writing variants to VCF file...")
    write_vcf(variants, output_vcf)

    print("Variant calling complete. Results saved to:", output_vcf)

if __name__ == "__main__":
    main()

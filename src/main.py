import sys
from src.indexing import index_reference_genome
from src.alignment import align_reads_to_reference
from src.bam_indexing import index_bam_file
from src.variant_analysis import analyze_reads
from src.write_vcf import write_vcf

def main():
    if len(sys.argv) != 5:
        print("Usage: python src/main.py <reference_fasta> <reads_fastq> <output_bam> <output_vcf>")
        sys.exit(1)

    reference_fasta, reads_fastq, output_bam, output_vcf = sys.argv[1:]

    # Index reference genome
    index_reference_genome(reference_fasta)

    # Align reads to reference genome
    output_sam = output_bam.replace('.bam', '.sam')
    align_reads_to_reference(reference_fasta, reads_fastq, output_sam)

    # Index BAM file
    index_bam_file(output_bam)

    # Analyze reads and call variants
    variants = analyze_reads(output_bam, reference_fasta)

    # Write VCF file
    write_vcf(variants, reference_fasta, output_vcf)
    print(f"Variants written to {output_vcf}")

if __name__ == "__main__":
    main()

import pysam
import os
import sys
from pysam import FastaFile

def index_reference_genome(reference_fasta):
    os.system(f'bwa index {reference_fasta}')
    os.system(f'samtools faidx {reference_fasta}')
    print(f"Indexed {reference_fasta}")

def align_reads_to_reference(reference_fasta, reads_fastq, output_sam):
    os.system(f'bwa mem {reference_fasta} {reads_fastq} > {output_sam}')
    print(f"Aligned {reads_fastq} to {reference_fasta}")

def index_bam_file(bam_file):
    os.system(f'samtools index {bam_file}')
    print(f"Indexed {bam_file}")

def analyze_reads(bam_file, reference_fasta):
    bam = pysam.AlignmentFile(bam_file, "rb")
    reference = FastaFile(reference_fasta)
    variants = []

    for pileupcolumn in bam.pileup():
        ref_pos = pileupcolumn.reference_pos
        ref_name = bam.get_reference_name(pileupcolumn.tid)
        ref_base = reference.fetch(ref_name, ref_pos, ref_pos+1)
        base_counts = {}
        total_reads = 0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base_counts[base] = base_counts.get(base, 0) + 1
                total_reads += 1

        af = 0.0
        dp = 0
        if total_reads > 0:
            af = sum(base_counts.values()) / total_reads
            dp = total_reads

        for base, count in base_counts.items():
            if base != ref_base and count / total_reads > 0.2:
                variants.append((ref_name, ref_pos+1, ref_base, base))

    return variants

def write_vcf(variants, reference_fasta, output_vcf):
    vcf_header = f"""##fileformat=VCFv4.2
##reference={reference_fasta}
#CHROM\tPOS\tREF\tALT"""
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    with open(output_vcf, 'w') as vcf:
        vcf.write(vcf_header + '\n')
        for variant in variants:
            chrom, pos, ref, alt = variant
            vcf.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python scripts/code.py <reference_fasta> <reads_fastq> <output_bam> <output_vcf>")
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
